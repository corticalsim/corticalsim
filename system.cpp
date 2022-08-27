#include "corticalSimReal.h"
#include "meshImport.h"

System::System(char* parFile) :
    EventDescriptorID(0),
    eventID(0),
    memUsage(0),
    totalLength(0),
    systemTime(0),
    systemTimeOffset(0),
    stopSignal(false),
    timeQueue(this, &System::identity, &System::identity),
    vPlusQueue(this, &System::vPlusToTime, &System::timeToVPlus),
    totalSEventCount(0),
    totalInvalidDEventCount(0),
    totalValidDEventCount(0),
    boundaryCrossingCount(0),
    totalZipperCount(0),
    totalCrossoverCount(0),
    totalInducedCatastropheCount(0),
    totalLengthSeveringCount(0),
    totalIntersectionSeveringCount(0),
    countSegments(0),
    countTrajectories(0),
    countIntersections(0),
    dataSaved(false),
    currentTimeTag(0),
    wallClockStartTime(time(0)),
    nextSnapshotEventTime(0),
    nextStatusEventTime(0),
    p(this),geometry(NULL)
{
     	/// read parameters
     	p.initialize(parFile);

        // initialize number of growing tips in different faces
        for(int dTag =0;dTag<=p.faceNumber;dTag++)
        growingTipsReg.push_back(0);

	// load geometry
	geometry = new TriMeshGeometry(this);

     	/// initialize the random seed with the desired one
     	randomGen.seed(static_cast<int>(p.seed));

	/// used in multiPcross & multiPzip events
        makeBinomialTable();

        /// perform further initialization
        initializeOutput();
        p.writeToFile();

	/// start creating MT events
        determineStochasticEvent();
        timeQueue.pushGlobal(p.measurementInterval, status);
        nextStatusEventTime = p.measurementInterval;

        /// if snapShot events are requested
        if (p.movieEnabled)
        {
            timeQueue.pushGlobal(p.movieFrameInterval, snapshot);
            nextSnapshotEventTime = p.movieFrameInterval;
        }

        /// if parameter file is updated on run
        if (p.PPB)
        {
            timeQueue.pushGlobal(p.PPBformingTime, signalPPB);
            nextPPBEventTime = p.PPBformingTime;
        }

	/// if parameter file is updated on run
        if (p.newParameterReadInterval > ZERO_CUTOFF)
        {
            timeQueue.pushGlobal(p.newParameterReadInterval, parameter_change);
            nextParameterEventTime = p.newParameterReadInterval;
        }

    return;
}

System::~System()
{
    // destroy all MTs
    growing_mts.removeAll();
    shrinking_mts.removeAll();

    return;
}

void System::run(double simTime,string& cS, string & mS)
{
    int wallTimePollingCounter(0);
    int memoryPollingCounter(0);
    timeQueue.pushGlobal(simTime - systemTimeOffset, stop);

    if(p.showOutput == 1)
    {
    cout <<"==============================================="<<endl;
    cout << " Time stamped simulation output " <<endl;
    cout <<"==============================================="<<endl;
    }

    while (!stopSignal)
    {
        if(++wallTimePollingCounter == CLOCK_POLLING_INTERVAL)
        {
            if (difftime(time(0),wallClockStartTime) >= p.wallClockLimit)
            {
                //time's up; send a stop signal
                cerr << "Wall clock limit (" << p.wallClockLimit << " seconds) reached. Stopping.\n";
                timeQueue.pushGlobal(systemTime, stop);

                // write final snapshot when running out of time.
                if (p.movieEnabled)
                geometry->outputSnapshot(movieFile);

		// emergency exit the program
		emergencyBreak();
		exit(-1);
            }
            wallTimePollingCounter = 0;
        }

        if (++memoryPollingCounter == MEMORY_POLLING_INTERVAL)
        {
            if (estimateMemoryFootprint() > 1024*1024*p.memoryLimit)
            {
                cerr << "Memory limit (approx. " << p.memoryLimit << " MB) reached. Stopping.\n";
                timeQueue.pushGlobal(systemTime, stop);

		// emergency exit the program
		emergencyBreak();
		exit(-1);
            }
            memoryPollingCounter = 0;
        }

        if (eventID > QUEUE_FLUSH_INTERVAL)
        {
            #ifdef DBG_QMGMT
            cout << "Flushing queues to reduce numerical drift.\n";
            #endif

	    //moved forward (used to be after flush and reload)
            if(!integrityCheck())
                exit(-1);

            flushAndReload();
            if(!integrityCheck())
                exit(-1);

            #ifdef DBG_QMGMT
            cout << "Queue reloaded.\n";
            #endif
        }

	// continue to next event
        nextEvent();
    }

    // stop signal received : write down necessary output(s)
    writeMeasurementsToFile(0);
    closeFiles();

    return;
}

void System::emergencyBreak()
{
    // write current output(s) and close all opened files
    writeMeasurementsToFile(0);
    closeFiles();

    return;
}

int System::estimateMemoryFootprint()
{
    // calculate memory in use
    int total(0);
    total += (growing_mts.size() + shrinking_mts.size())*(sizeof(Microtubule) + 2*sizeof(MTTip*));
    total += EventDescriptorMap.size() * (sizeof(pair<EventDescriptorIndex, EventDescriptor*>) + 3*sizeof(void*));
    total += countSegments * (sizeof(Segment) + sizeof(Segment*));
    total += countTrajectories * (sizeof(Trajectory));
    total += countIntersections * (2*(sizeof(pair<double, Intersection>) + 3*sizeof(void*)));
    total += (timeQueue.queue.size() + vPlusQueue.queue.size()) * (sizeof(DeterministicEvent) + 2*sizeof(void*));
    total += measurementHistory.size() * sizeof(Measurement);

    return total;
}

bool System::integrityCheck()
{
    bool valid(true);

    // integrity check of the gemoetry embeded components
    if (!geometry->integrityCheck())
        valid = false;

    // integrity check of each growing MTs
    Microtubule* mt = growing_mts.first();
    while (mt != NULL)
    {
        if (!mt->integrityCheck())
            valid = false;

        mt = mt->next();
    }

    // integrity check of each shrinking MTs
    mt = shrinking_mts.first();
    while (mt != NULL)
    {
        if (!mt->integrityCheck())
            valid = false;

        mt = mt->next();
    }

    return valid;
}

void System::flushAndReload(bool refreshParameterEvent)
{
    double timeOffset(systemTime);

    // update all microtubule lengths
    updateAll();

    // recalculate total lengths to avoid drift (in DBG mode, also check for large deviations)
    double regionLength(0.);
    double systemLength(0);
    int ridx(0);
    Trajectory* tptr(NULL);

    for (ridx = 0; ridx < geometry->regions.size(); ridx++)
    {
        regionLength = 0;
        tptr = geometry->regions[ridx]->trajectories.first();

        while (tptr != NULL)
        {
            regionLength += tptr->segmentLength();
            tptr = tptr->next();
        }

        if (abs(geometry->regions[ridx]->totalLength - regionLength) > ZERO_CUTOFF*max(100.,regionLength))
        {
            cerr << "Unacceptable drift in region length (drift = " << abs(geometry->regions[ridx]->totalLength - regionLength) << "). Exiting. [try lowering flush interval]\n";
            emergencyBreak();
            exit(-1);
        }
        geometry->regions[ridx]->totalLength = regionLength;
        systemLength += regionLength;
    }

    if(abs(totalLength - systemLength) > ZERO_CUTOFF*max(1.,totalLength))
    {
        cerr << "Unacceptable drift in system length (drift = " << abs(totalLength - systemLength) << "). Exiting. [try lowering flush interval]\n";
        emergencyBreak();
        exit(-1);
    }

    totalLength = systemLength;

    // reset master timer
    systemTimeOffset += timeOffset;

    systemTime = 0;
    timeQueue.flush();
    vPlusQueue.flush();
    EventDescriptorMap.clear();
    eventID = 0;
    EventDescriptorID = 0;

    // now, reload all events (this is a PROBLEM for global events that have not been saved. They need to be re-inserted manually)
    determineStochasticEvent();
    timeQueue.pushGlobal(p.stopTime - systemTimeOffset, stop);
    nextStatusEventTime -= timeOffset;
    timeQueue.pushGlobal(nextStatusEventTime, status);

    if (refreshParameterEvent)
    {
        nextParameterEventTime -= timeOffset;
        timeQueue.pushGlobal(nextParameterEventTime, parameter_change);
    }

    if (p.movieEnabled)
    {
        nextSnapshotEventTime -= timeOffset;
        timeQueue.pushGlobal(nextSnapshotEventTime, snapshot);
    }

    if(p.PPB)
    {
        nextPPBEventTime -= timeOffset;
	timeQueue.pushGlobal(nextPPBEventTime,signalPPB);
    }

    Microtubule* mt(NULL);
    Segment* seg(NULL);
    mt = growing_mts.first();

    while (mt != NULL)
    {
        // this can go into a change time function
        mt->nucleationTime -= timeOffset;
        seg = mt->segments.first();

        while (seg  != NULL)
        {
            seg->nucleationTime -= timeOffset;
            seg = seg->next();
        }

        mt->plus.event.index = registerEventDescriptor(&(mt->plus.event));
        mt->minus.event.index = registerEventDescriptor(&(mt->minus.event));
        mt->disappearEvent.index = registerEventDescriptor(&(mt->disappearEvent));
        mt->plus.determineEvent();
        mt->minus.determineEvent();
        mt->setDisappearEvent();
        mt = mt->next();
    }

    mt = shrinking_mts.first();

    while (mt != NULL)
    {
        mt->nucleationTime -= timeOffset;
        seg = mt->segments.first();

        {
            seg->nucleationTime -= timeOffset;
            seg = seg->next();
        }

        mt->plus.event.index = registerEventDescriptor(&(mt->plus.event));
        mt->minus.event.index = registerEventDescriptor(&(mt->minus.event));
        mt->disappearEvent.index = registerEventDescriptor(&(mt->disappearEvent));
        mt->plus.determineEvent();
        mt->minus.determineEvent();
        mt->setDisappearEvent();
        mt = mt->next();
    }

    return;
}

#ifndef NO_INLINE
inline
#endif

double LambertW(double x)
{
    double out(0.0);
    const double expMin1 = exp(-1.0);
    const double zeroCutoff = pow(2.0, -DBL_MANT_DIG/6.0);
    const double etaCutoff = 34.0 * pow(2.0, 2.0*(1.0-DBL_MANT_DIG)/7.0);

    // LambertW function (Barry et al., ACM Trans. on Math. Software, Vol. 21, No. 2, June 1995, pp. 161--171)
    if (x < -expMin1) return -VERY_LARGE;

    if (abs(x) < zeroCutoff)
    {
        out = x*(60 + 114*x + 17*x*x)/(60 + 174*x + 101*x*x);
    }
    else if (x > 20.0)
    {
        double h = exp(- 1.124491989777808 /(0.4225028202459761 + log(x)));
        out = log(x /log(x/pow(log(x),h)));

        double z_n = log(x/out) - out;
        double temp = 2.0*(1.0 + out)*(1.0 + out + z_n*2.0/3.0);
        double e_n = (z_n/(1.0 + out))*(temp - z_n)/(temp - 2.0*z_n);
        out = out*(1.0 + e_n);
    }
    else
    {
        double eta = 2.0 + 2.0*exp(1.0)*x;
        double sqrt_eta = sqrt(eta);
        double D(0.0);
        if (sqrt_eta >= etaCutoff)
        {
            double n2 = 4.612634277343749*sqrt(sqrt(sqrt_eta + 1.09556884765625));
            double n1 = (4.0 - 3.0*sqrt(2.0) + (2.0*sqrt(2.0) - 3.0)*n2)/(sqrt(2.0) - 2.0);
            D = n1*sqrt_eta/(n2 + sqrt_eta);

            out = -1.0 + sqrt_eta/(1.0 + sqrt_eta/(3.0 + D));

            double z_n = log(x/out) - out;
            double temp = 2.0*(1.0 + out)*(1.0 + out + z_n*2.0/3.0);
            double e_n = (z_n/(1.0 + out))*(temp - z_n)/(temp - 2.0*z_n);
            out = out*(1.0 + e_n);
        }
        else
        {
            D = sqrt_eta/(8.0/3.0 + sqrt_eta/(135.0/83.0 + sqrt_eta/(166.0/39.0 + sqrt_eta*3167.0/3549.0)));
            out = -1.0 + sqrt_eta/(1.0 + sqrt_eta/(3.0 + D));
        }
    }

    #ifdef DBG_ACID_TEST
    double test = (x - out*exp(out))/((1.0 + out)*exp(out));
    if (abs(test) > ZERO_CUTOFF*abs(x))
        cerr << "LambertW() accuracy [" << test << "] below threshold. Output may contain inaccuracies.\n";
    #endif

    return out;
}

double System::vPlusToTime(double dist)
{
    // for infinite pool (implimented as identity function)
    if (p.restrictedPool == 0)
        return dist;

    // for finite tubulin pool
    double Lmax = (p.poolDensity*geometry->area);
    double invLMax = 1.0/(p.poolDensity*geometry->area);

    double alpha = (p.vPlus - p.vTM)*growing_mts.size() + (p.vMin - p.vTM)*shrinking_mts.size();

    if  (growing_mts.size() == 0)
    {
        if (shrinking_mts.size() == 0)
            return dist/(1.0 - totalLength*invLMax);
        else
            return ((Lmax - totalLength) - sqrt(pow(Lmax - totalLength,2) - 2.*alpha*Lmax*dist))/alpha;
    }

    double beta = p.vPlus*growing_mts.size()*invLMax;

    if ((shrinking_mts.size() == 0) && ((p.vTM/p.vPlus) < ZERO_CUTOFF))
        return -log(1.0 - growing_mts.size()*p.vPlus*dist/(Lmax-totalLength))/beta;

    double pv = (alpha/beta - totalLength)/(Lmax - alpha/beta);
    double qv = growing_mts.size()*p.vPlus*dist/(Lmax - alpha/beta);

    // verify whether this is an optimal trade-off
    if (pv - qv < 250)
        return (qv-pv + LambertW(pv*exp(pv-qv)))/beta;

    // if not, do asymptotic approximation
    double L1 = pv-qv + log(pv);
    double L2 = log(L1);
    return (1.0/beta)*(qv - pv + L1 - L2 + L2/L1 + L2*(L2-2)/(2*pow(L1,2))+ L2*(6-9*L2+2*pow(L2,2))/(6*pow(L1,3))
            + L2*(-12+36*L2-22*pow(L2,2)+3*pow(L2,3))/(12*pow(L1,4)));

}

double System::timeToVPlus(double time)
{
    // for infinite pool (implimented as identity function)
    if (p.restrictedPool == 0)
        return time;

    // for finite tubulin pool
    double invLMax = 1.0/(p.poolDensity*geometry->area);
    double alpha = (p.vPlus - p.vTM)*growing_mts.size() + (p.vMin - p.vTM)*shrinking_mts.size();

    if  (growing_mts.size() == 0)
        return ((1.0 - totalLength*invLMax) - 0.5*alpha*invLMax*time)*time;

    double beta = p.vPlus*growing_mts.size()*invLMax;
    double lambda_inf = invLMax*alpha/beta;
    double lambda_0 = totalLength*invLMax;

    return (1.0 - lambda_inf)*time - (lambda_0 - lambda_inf)*(1.0 - exp(-beta*time))/beta;
}

void System::advanceTime(double newTime)
{
    #ifdef DBG_SYSTEM
    cout << "DBG/SYSTEM: advancing time by " << newTime - systemTime << " seconds. New system time is " << newTime << ".\n";
    #endif

    // make a jump in time step (taking advantage of using event-driven algorithm)
    timeQueue.advanceTime(newTime);
    vPlusQueue.advanceTime(newTime);
    systemTime = newTime;
    totalLength += timeQueue.progression(currentTimeTag)*(p.vMin*shrinking_mts.size() - p.vTM*(growing_mts.size() + shrinking_mts.size()) )
                   + vPlusQueue.progression(currentTimeTag)*(p.vPlus*growing_mts.size());

    if (++currentTimeTag == POSITION_CACHE_SIZE)
    {
        currentTimeTag = 0;
        updateAll(true);
    }

    timeQueue.storeTime(currentTimeTag);
    vPlusQueue.storeTime(currentTimeTag);

    return;
}

void System::updateAll(bool forceUpdate)
{
    // always nice to be updated in lengths of all regions and microtubules (convenient for snapshots, etc.)
    for (int ridx = 0; ridx < geometry->regions.size(); ridx++)
        geometry->regions[ridx]->updateRegionLength(forceUpdate);

    Microtubule* mt;
    mt = growing_mts.first();

    while (mt!= NULL)
    {
        mt->updateLength(forceUpdate);

        #ifdef DBG_ACID_TEST
        if (mt->segments.size() == 1)
        {
            if ((mt->segments.first()->end - mt->segments.first()->start)*mt->segments.first()->dir < -ZERO_CUTOFF)
            {
                cerr << "ERROR: negative segment length after updating.\n";
                exit(-1);
            }
        }
        #endif

        mt = mt->next();
    }

    mt = shrinking_mts.first();
    while (mt!= NULL)
    {
        mt->updateLength(forceUpdate);

        #ifdef DBG_ACID_TEST
        if (mt->segments.size() == 1)
        {
            if ((mt->segments.first()->end - mt->segments.first()->start)*mt->segments.first()->dir < -ZERO_CUTOFF)
            {
                cerr << "ERROR: negative segment length after updating.\n";
                exit(-1);
            }
        }
        #endif
        mt = mt->next();
    }

    return;
}

void System::nextEvent()
{
    double tqValue = numeric_limits<double>::max();
    double dqValue = numeric_limits<double>::max();

    // maybe clean up a little
    if (!timeQueue.empty())
        tqValue = timeQueue.firstEventTime();

    if (!vPlusQueue.empty())
        dqValue = vPlusQueue.firstEventTime();

    if ((nextStochasticEventTime < tqValue) && (nextStochasticEventTime < dqValue))
    {
        // execute stochastic event
        totalSEventCount++;
        #ifdef DBG_EVENT
        cout << "DBG/Event: executing stochastic event [time=" << nextStochasticEventTime << ", type=" << nextStochasticEventType << "]\n";
        #endif

        advanceTime(nextStochasticEventTime);

        // execute next stochastic event
        switch(nextStochasticEventType)
        {
        case catastrophe:
            handleCatastropheEvent();
            break;
        case rescue:
            handleRescueEvent();
            break;
        case katanin:
            handleSeveringEvent();
            break;
        #ifdef CROSS_SEV
        case severingAtCross:
            handleSeveringAtCrossEvent();
            break;
        #endif
        case nucleation:
            handleNucleationEvent();
            break;
        }
        determineStochasticEvent();
    }

    else
    {
        double nextEventTime(0.0);

        // execute deterministic event
        DeterministicEvent event;
        if (tqValue < dqValue)
        {
            event = timeQueue.pop();
            nextEventTime = tqValue;
        }
        else
        {
            event = vPlusQueue.pop();
            nextEventTime = dqValue;
        }

        if (event.infoIdx != -1)
        {
            map<EventDescriptorIndex, EventDescriptor*>::iterator eventItr = EventDescriptorMap.find(event.infoIdx);

            // whether the tip still exists and the tag still matches (valid event)
            if ((eventItr != EventDescriptorMap.end())
                    && (eventItr->second->tag == event.tag))
            {
                #ifdef DBG_EVENT
                cout << "DBG/Event: executing deterministic event [time=" << event.eventTimeDist << ", type=" << eventItr->second->type << ", tag=" << event.tag << "]\n";
                #endif

                totalValidDEventCount++;
                advanceTime(nextEventTime);
                eventItr->second->mt->handleEvent(eventItr->second);
                determineStochasticEvent();
            }

            // invalid event
            else
            {
                totalInvalidDEventCount++;
                #ifdef DBG_EVENT
                cout << "DBG/Event: Event is no longer valid. [tag= " << event.tag << "]\n";
                #endif
            }
        }

        // execute global event
        else
        {
            #ifdef DBG_EVENT
            cout << "DBG/Event: executing global event [time=" << event.eventTimeDist << ", type=" << event.global_type << "]\n";
            #endif

            advanceTime(event.eventTimeDist);
            handleGlobalEvent(event);
        }
    }

    return;
}

#ifndef NO_INLINE
inline
#endif

void System::handleCatastropheEvent()
{

    int tipNumber(0), totalRelevantTips(0),ridx(0),domainType(0);
    RegionMTTipTag tipTag;

    double norm(0.0);
    for(int dTag = 0; dTag<= p.faceNumber;dTag++)
    norm+= growingTipsReg[dTag]*p.RegionKcatMultiplier[dTag];

    double randGen = randomGen.rand();

    double normFraction(0.0);
    for(int dTag = 0; dTag<= p.faceNumber;dTag++)
    {
	normFraction+= growingTipsReg[dTag]*p.RegionKcatMultiplier[dTag];
        if(randGen < (double)(normFraction)/norm)
        {
		domainType = dTag;
		totalRelevantTips = growingTipsReg[dTag];
		break;
	}
    }

    tipNumber = randomGen.randInt(totalRelevantTips - 1);

    // determine the region this tip is located in
    // depending on whether it's a high or low number, start scanning from the beginning or end
    if (tipNumber < totalRelevantTips/2)
    {
        while (true)
        {
            if (geometry->regions[ridx]->faceTag == domainType)
            {
                if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size())
                    break;

                tipNumber -= geometry->regions[ridx]->growingPlusTipList.size();
            }
            ridx++;
        }
    }

    else
    {
        // flip the direction of indexing: back to front
        tipNumber = totalRelevantTips- 1 - tipNumber;
        ridx = geometry->regions.size()-1;
        while (true)
        {
            if (geometry->regions[ridx]->faceTag == domainType)
            {
                if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size())
                    break;
                tipNumber -= geometry->regions[ridx]->growingPlusTipList.size();
            }
            ridx--;
        }
        // flip the direction of the result for the next step
        tipNumber = geometry->regions[ridx]->growingPlusTipList.size() - 1 - tipNumber;
    }

    // Now, select the tip with # 'tipNumber' from the current region.
    // Again, if the number is high within the number of tips in the region, start from the end
    if (tipNumber < geometry->regions[ridx]->growingPlusTipList.size()/2)
    {
        tipTag = geometry->regions[ridx]->growingPlusTipList.begin();
        while (tipNumber >0)
        {
            tipTag++;
            tipNumber--;
        }
    }
    else
    {
        tipTag = --(geometry->regions[ridx]->growingPlusTipList.end());
        while (tipNumber < geometry->regions[ridx]->growingPlusTipList.size()-1)
        {
            tipTag--;
            tipNumber++;
        }
    }

    (*tipTag)->mt->catastrophe();

    return;
}

#ifndef NO_INLINE
inline
#endif

void System::handleRescueEvent()
{
    int ridx(0);
    int tipNumber(0);
    RegionMTTipTag tipTag;

    // select the number of the tip to be rescued
    tipNumber = randomGen.randInt(shrinking_mts.size() - 1);

    // determine the region this tip is located : binary search algorithm
    if (tipNumber < shrinking_mts.size()/2)
    {
        while (tipNumber >= geometry->regions[ridx]->shrinkingPlusTipList.size())
        {
            tipNumber -= geometry->regions[ridx]->shrinkingPlusTipList.size();
            ridx++;
        }
    }

    else
    {
        // flip the direction of indexing: back to front
        tipNumber = shrinking_mts.size()- 1 - tipNumber;
        ridx = geometry->regions.size()-1;
        while (tipNumber >= geometry->regions[ridx]->shrinkingPlusTipList.size())
        {
            tipNumber -= geometry->regions[ridx]->shrinkingPlusTipList.size();
            ridx--;
        }
        // flip the direction of the result for the next step
        tipNumber = geometry->regions[ridx]->shrinkingPlusTipList.size() - 1 - tipNumber;
    }

    // determine which tip in the selected region : binary search algorithm
    if (tipNumber < geometry->regions[ridx]->shrinkingPlusTipList.size()/2)
    {
        tipTag = geometry->regions[ridx]->shrinkingPlusTipList.begin();
        while (tipNumber > 0)
        {
            tipTag++;
            tipNumber--;
        }
    }
    else
    {
        tipTag = --(geometry->regions[ridx]->shrinkingPlusTipList.end());
        while (tipNumber < geometry->regions[ridx]->shrinkingPlusTipList.size()-1)
        {
            tipTag--;
            tipNumber++;
        }
    }

    // rescue the selected tip
    (*tipTag)->mt->rescue();

    return;
}

#ifndef NO_INLINE
inline
#endif

void System::handleSeveringEvent()
{
    double cutLength(0.0);
    Segment* cutSeg(NULL);
    randomPositionOnMicrotubule(cutLength, cutSeg);

    // sever the MT now !
    if (cutSeg != NULL)
        cutSeg->mt->sever(cutSeg,cutLength);
    else
        cerr << "ERROR (non-fatal): severing position outside total MT length. Ignoring event.\n";

    return;
}

void System::randomPositionOnMicrotubule(double& randomPos, Segment*& randomSeg)
{
    int ridx(0);
    // select position at which severing will take place
    double cutLength = randomGen.randDblExc(totalLength);
    // first, select relevant region. Go from both sides to reduce time.
    if (cutLength < 0.5*totalLength)
    {
        ridx = 0;
        while (true)
        {
            geometry->regions[ridx]->updateRegionLength();
            if (cutLength < geometry->regions[ridx]->totalLength)
                break;
            cutLength -= geometry->regions[ridx]->totalLength;
            ridx++;
            if (ridx == geometry->regions.size())
            {
                // apparently, cutLength was *just* too long...
                cutLength -= ZERO_CUTOFF;
                #ifdef DBG_ACID_TEST
                if (cutLength > geometry->regions[ridx-1]->totalLength)
                {
                    cerr << "ERROR: cannot locate intersection location\n";
                    exit(-1);
                }
                #endif
                break;
            }
        }
    }
    else
    {
        ridx = geometry->regions.size()-1;
        cutLength = totalLength - cutLength;
        while (true)
        {
            geometry->regions[ridx]->updateRegionLength();
            if (cutLength < geometry->regions[ridx]->totalLength)
                break;
            cutLength -= geometry->regions[ridx]->totalLength;
            ridx--;
            if (ridx == -1)
            {
                // apparently, cutLength was *just* too long...
                cutLength -= ZERO_CUTOFF;

                #ifdef DBG_ACID_TEST
                if (cutLength > geometry->regions[0]->totalLength)
                {
                    cerr << "ERROR: cannot locate intersection location\n";
                    exit(-1);
                }
                #endif
                break;
            }
        }
        //and wrap it back
        cutLength = geometry->regions[ridx]->totalLength - cutLength;
    }

    // Now that the area has been selected, cycle over trajectories and their segments
    Trajectory* trptr(NULL);
    list<Segment*>::iterator seg;
    double temp(0.0);
    bool stopFlag(false);
    // Again, let the search direction depend on the expected position
    if (cutLength < 0.5*(geometry->regions[ridx]->totalLength))
    {
        // search in the forward direction. Start with the first trajectory
        trptr = geometry->regions[ridx]->trajectories.first();
        while (true)
        {
            // get the first segment on this trajectory
            seg = trptr->segments.begin();
            // and loop over all segments
            while (seg != trptr->segments.end())
            {
                // make sure the lengths are accurate
                (*seg)->mt->updateLength();
                // break when the correct segment is found
                if (cutLength < (temp = (*seg)->length()))
                {
                    stopFlag = true;
                    break;
                }
                else
                {
                    // if not, continue
                    cutLength -= temp;
                    seg++;
                }
            }

            if (stopFlag)
                break;

            // if we ran out of trajectories, something is wrong. This could happen due to numerical inaccuracies.
            if (trptr->next() == NULL)
            {
                cerr << "Rare event: cutting on the edge\n";
                seg--;
                cutLength -= ZERO_CUTOFF;
                break;
            }
            trptr = trptr->next();
        }
    }
    else
    {
        // flip the length and trajectory search orders
        cutLength = geometry->regions[ridx]->totalLength - cutLength;
        trptr = geometry->regions[ridx]->trajectories.last();
        while (true)
        {
            // also, start with the last segment (one step beyond, in this case)
            seg = (trptr->segments.end());
            do
            {
                seg--;
                // make sure that the segment length is correct
                (*seg)->mt->updateLength();
                // break if the requested position lies on this segment
                if (cutLength < (temp = (**seg).length()))
                {
                    stopFlag = true;
                    break;
                }
                else
                    cutLength -= temp;
            }
            while (seg != trptr->segments.begin());

            if (stopFlag)
                break;

            if (trptr->previous() == NULL)
            {
                cerr << "Rare event: cutting on the edge\n";
                seg++;
                cutLength -= ZERO_CUTOFF;
                break;
            }
            trptr = trptr->previous();
        }
        cutLength = (*seg)->length() - cutLength;
    }

    randomSeg = *seg;
    randomPos = cutLength;
    return;
}

#ifdef CROSS_SEV
#ifndef NO_INLINE
inline
#endif

// TODO: verify this function
void System::handleSeveringAtCrossEvent()
{
    Microtubule* cutMT(NULL);
    IntersectionItr cutIS;
    OccupiedIntersection* cutOccIS(NULL);
    Segment* cutSeg(NULL);
    int cutSegNumber(0);
    list<Segment*>::iterator testSeg;
    list<Segment*>::iterator finalTestSeg;

    cutOccIS = OccupiedIntersectionList.ElementAddress(randomGen.randInt(OccupiedIntersectionList.size()-1));
    if (p.crossSeveringTop || (randomGen.randInt(1)==0))
        // if MT (bundle) on top should be cut, select that one, or (if random) pick it with 50% chance.
        // Note that no effort is made to ensure proportional selection for bundles of different thickness.
        cutIS = cutOccIS->intersectionToCut;
    else
        cutIS = cutOccIS->intersectionToCut->second.mirror;
    double diffAngle = cutIS->second.otherTrajectory->base.region->intersectionAngle(
                           cutIS->second.otherTrajectory,
                           cutIS->second.mirror->second.otherTrajectory);
    if (diffAngle > 0.5*PI)
        diffAngle = PI - diffAngle;
    // don't do anything if the angle between MTs is smaller than the cutoff angle
    if (diffAngle <= (p.crossSeveringStartAngle)*PI/180.)
        return;

    /*
    IMPORTANT NOTE: DO NOT REMOVE OCCUPIED INTERSECTION HERE! IT WILL BE REMOVED AUTOMATICALLY IN A LATER EVENT
    cutIS->second.occupancy--;		// do not decrease: will happen very soon in a new event
    {
    	removeOccupiedIntersection(cutIS->second);
    }
    */
    // find the segment to be cut
    cutSegNumber = randomGen.randInt(cutIS->second.occupancy - 1)  ; // nummers 0 - occupancy -1 .
    if (cutSegNumber < cutIS->second.occupancy/2)	// search from the beginning.
        // In principle, all searches could start from the beginning. Also starting from the end should save (on average) half of the time
    {
        testSeg = cutIS->second.mirror->second.otherTrajectory->segments.begin();
        finalTestSeg = cutIS->second.mirror->second.otherTrajectory->segments.end();
        cutSegNumber++;
        while (testSeg != finalTestSeg) // NOTE: test should not be necessary "while(true)"
        {
            if ((**testSeg).crossesIntersection(cutIS))
            {
                cutSegNumber--;
                if (cutSegNumber == 0)
                {
                    cutSeg = *testSeg;
                    break;
                }
            }
            testSeg++;
        }
    }
    else	// search from the end
    {
        testSeg = cutIS->second.mirror->second.otherTrajectory->segments.end();
        finalTestSeg = cutIS->second.mirror->second.otherTrajectory->segments.begin();
        while (testSeg != finalTestSeg)
        {
            testSeg--;
            if ((**testSeg).crossesIntersection(cutIS))
            {
                cutSegNumber++;
                if (cutSegNumber == cutIS->second.occupancy)
                {
                    cutSeg = *testSeg;
                    break;
                }
            }
        }
    }
    if(cutSeg != NULL)
    {
        cutMT = cutSeg->mt;
        cutMT->updateLength();
        if (cutMT != NULL)
            cutMT->severAtCross(cutIS, cutSeg);
        else
            cerr << "ERROR (non-fatal): severing position outside total MT length. Ignoring event.\n";
    }
    else
        cerr << "ERROR (non-fatal): Did not find the right segment. Ignoring event.\n";


    return;
}
#endif

#ifndef NO_INLINE
inline
#endif

void System::handleNucleationEvent()
{
    // if there is a finite tubulin pool and vPlus is essentially the same as vTM : leading to moving microtubules with zero length, causing ordering problems on collisions
    if ((p.restrictedPool != 0) && (p.vTM/(p.vPlus*(1. - totalLength/(p.poolDensity*geometry->area))) > 0.9999))
    return;

    // no nucleation at forbidden zone
    TrajectoryVector tv;

    // get random nucleation point (vector)
    SurfaceVector sv = geometry->randomSurfaceVector();

    // forbiden zone and parameter modified : no microtubule created
    if (sv.region->faceTag == -1)
        return;

    Segment* nucSeg(NULL);
    double pos(0.0),posOnSeg(0.0),overrideAngle(0.0),rand(0.0),rand2(0.0),preTheta(0.0),cosTheta(0.0),temp(0.0);

    switch (p.nucleationType)
    {
	// select isotropic nucleation
    	case nuc_isotropic:
        	overrideAngle = sv.angle;
        	break;
 	default:
  		break;
    }

    sv.angle = overrideAngle;
    tv = geometry->createTrajectory(sv);
    Microtubule* mt = growing_mts.create(this, tv);

    return;
}

void System::determineStochasticEvent()
{
   // for(regionIndex=0 ; oldtr->base.region != regions[regionIndex] ; regionIndex++);


    double nextEventInterval(0.0);
    double invTotalRate(0.0);
    double typeSelector(0.0);
    double logRandom(0.0);
    double temp(0.0);

    if ((growing_mts.size() == 0) || (p.restrictedPool == 0))
    {
        double constRate(0.0);
        double slopeRate(0.0);

        constRate = p.kNuc*geometry->area + p.kRes*shrinking_mts.size() +  p.kSev*totalLength;

        for(int dTag =0;dTag <=p.faceNumber; dTag++)
	constRate+=p.kCat*p.RegionKcatMultiplier[dTag]*growingTipsReg[dTag];

        #ifdef CROSS_SEV
        constRate += p.kCross*OccupiedIntersectionList.size();
        #endif

        slopeRate = p.kSev*((p.vPlus - p.vTM)*growing_mts.size() + (p.vMin - p.vTM)*shrinking_mts.size());

	// generate a random number larger than 0 with an exponential distribution
        logRandom = -log(randomGen.randDblExc());

        if (abs(slopeRate) < ZERO_CUTOFF)
        {
            nextEventInterval = logRandom / constRate;
        }
        else
        {
            if ((temp = pow(constRate,2) + 2*slopeRate*logRandom) < 0)

		// choose a time that is always larger than the disappearance time of all MTs
                nextEventInterval = VERY_LARGE;

            else
                nextEventInterval = (-constRate + sqrt(temp))/slopeRate;
        }

	// it will be negative for the case with no stochastic events
        invTotalRate = 1.0/(constRate + slopeRate*nextEventInterval);
    }

    else
    {
        double invLMax = 1.0/(p.poolDensity*geometry->area);
        double alpha = (p.vPlus - p.vTM)*growing_mts.size() + (p.vMin - p.vTM)*shrinking_mts.size();
        double beta = p.vPlus*growing_mts.size()*invLMax;
        double L_inf = alpha/beta;

	double baseRate = p.kNuc*geometry->area+p.kRes*shrinking_mts.size();
	for(int dTag =0;dTag <= p.faceNumber; dTag++)
	baseRate+=p.kCat*p.RegionKcatMultiplier[dTag]*growingTipsReg[dTag];

        #ifdef CROSS_SEV
        baseRate += p.kCross*OccupiedIntersectionList.size();
        #endif

        double x = p.kSev*(totalLength - L_inf)/(baseRate +  p.kSev*L_inf);
        double y = -beta*log(randomGen.randDblExc())/(baseRate +  p.kSev*L_inf);

        nextEventInterval = (-x+y+ LambertW(x*exp(x-y)))/beta;

        double totalRate = baseRate + p.kSev*(L_inf + (totalLength - L_inf)*exp(-beta*nextEventInterval));
        invTotalRate = 1.0/totalRate;

    }

    #ifdef DBG_ACID_TEST
    if ((nextEventInterval < 0) || (nextEventInterval != nextEventInterval ))
    {
        cerr << "Invalid deterministic event time calculated. Exiting\n";
        exit(-1);
    }
    #endif

    nextStochasticEventTime = systemTime + nextEventInterval;

    // even if  invTotalRate is negative! due to the timing, the event will never be executed.
    typeSelector = randomGen();

    double typeSelectRegSource(0.0);
    for(int dTag =0;dTag <= p.faceNumber; dTag++)
    typeSelectRegSource+=p.kCat*p.RegionKcatMultiplier[dTag]*growingTipsReg[dTag];

    if(typeSelector < p.kNuc*geometry->area*invTotalRate)
        nextStochasticEventType = nucleation;

    // select catastrophe event
    else if (typeSelector < (p.kNuc*geometry->area + typeSelectRegSource)*invTotalRate)
        nextStochasticEventType = catastrophe;

    // select rescue event
    else if (typeSelector < (p.kNuc*geometry->area + typeSelectRegSource + p.kRes*shrinking_mts.size())*invTotalRate)
        nextStochasticEventType = rescue;

    // select severing event
    #ifdef CROSS_SEV
      else if (typeSelector < (p.kNuc*geometry->area + typeSelectRegSource + p.kRes*shrinking_mts.size() + p.kCross*OccupiedIntersectionList.size())*invTotalRate)
        nextStochasticEventType = severingAtCross;
    #endif

    // select katanin induced event
    else
        nextStochasticEventType = katanin;

    #ifdef DBG_EVENT
    cerr << "determined next stochastic event. type=" << nextStochasticEventType << ", time=" << nextStochasticEventTime << "\n";
    #endif

    return;
}

void System::handleGlobalEvent(DeterministicEvent& event)
{
    // better to be updated
    updateAll();

    switch (event.global_type)
    {
    	case status:
        	if (!integrityCheck())
        	{
            		cerr << "Integrity check failed. Exiting.\n";
            		exit(-1);
        	}

        	performMeasurement();
        	nextStatusEventTime = systemTime + p.measurementInterval;
        	timeQueue.pushGlobal(nextStatusEventTime, status);
        	break;

    	case snapshot:
		// declaring variable within case: so the local blocking
        	{
			OrderParameters order;

        		geometry->getOrderParameters(order);

			// write a movie frame
        		movieFile << order.Rdirector[0] << "," << order.Rdirector[1] << "," << order.Rdirector[2] << "#";
        		movieFile << geometry->objectCM(0,0) << "," << geometry->objectCM(1,0) << "," << geometry->objectCM(2,0) << "#";
        		movieFile << geometry->objectPA(0,0) << "," << geometry->objectPA(1,0) << "," << geometry->objectPA(2,0) << "*";

			geometry->outputSnapshot(movieFile);

			// write a local order parameter heatmap frame
			geometry->outputOrderHeatMap(heatMapFile,order.localOrder,order.Sv);

        		nextSnapshotEventTime = systemTime + p.movieFrameInterval;
        		timeQueue.pushGlobal(nextSnapshotEventTime, snapshot);
		}
        	break;

    	case stop:
                {
			OrderParameters order;
        		geometry->getOrderParameters(order);

                        // division plane: polygon edge-list & area
                        Vector3d normalPPB(order.Rdirector[0],order.Rdirector[1],order.Rdirector[2]);
                        geometry->areaPPB = intersectingPolygon(geometry->ppb,geometry->globalVertex,geometry->regions,geometry->ppbEdgeList,normalPPB,geometry->nucleousPosition,p.geometry);

			if(p.showOutput ==1)
			cout << "Division plane area: "<< geometry->areaPPB <<endl;

        		stopSignal = true;
		}
        	break;

	case parameter_change:
		{cout << "Changing parameters. Reading from file [" << p.newParameterFile << "]\n";
		flushAndReload(false);
		string oldDir = p.outputDir;
		if (p.reinitialize(p.newParameterFile.c_str()))
		{
			if (p.outputDir != oldDir)
			{
				cout << "Switching directories.\n";
				writeMeasurementsToFile(0);
				closeFiles();
				initializeOutput();
			}
			p.writeToFile();

			if (p.movieEnabled)
			{
				if (nextSnapshotEventTime < ZERO_CUTOFF)
				{
					timeQueue.pushGlobal(p.movieFrameInterval, snapshot);
					nextSnapshotEventTime = p.movieFrameInterval;
				}
			}
			else
				nextSnapshotEventTime = 0;

			if (p.newParameterReadInterval > ZERO_CUTOFF)
			{
				timeQueue.pushGlobal(p.newParameterReadInterval, parameter_change);
				nextParameterEventTime = p.newParameterReadInterval;
			}

		}
		else
		{
			cerr << "Failed to read parameters from file. Exiting.\n";
			emergencyBreak();
			exit(-1);
		}

		timeQueue.pushGlobal(p.newParameterReadInterval, parameter_change);
		nextParameterEventTime = p.newParameterReadInterval;
		break;}

	case signalPPB:
        	cout << " Starting formation of PPB ...\n";
    		p.writeToFile();
		p.PPBkNucFraction = 0.8;
		cout << p.PPBkNucFraction <<" : "<< p.PPB << endl;
		OrderParameters order;
		geometry->getOrderParameters(order);
        	Vector3d normalPPB(order.Rdirector[0],order.Rdirector[1],order.Rdirector[2]);
        	geometry->areaPPB = intersectingPolygon(geometry->ppb,geometry->globalVertex,geometry->regions,geometry->ppbEdgeList,normalPPB,geometry->nucleousPosition,p.geometry);

		for(int dTag = 0; dTag <= p.faceNumber; dTag++)
		growingTipsReg[dTag] = 0;

		for(int eno = 0; eno < geometry->regions.size(); eno++)
		{
			if(geometry->regions[eno]->polyIntersectMark==1)
			{
				geometry->regions[eno]->faceTag = 1;

				for(int tip =0; tip< geometry->regions[eno]->growingPlusTipList.size();tip++)
				growingTipsReg[1]+=1;

				geometry->patchArea[1]+=geometry->regions[eno]->area;
                        	geometry->RegionsIndex[1].push_back(eno);
			}

			else
			{
				geometry->regions[eno]->faceTag = 2;

				for(int tip =0; tip< geometry->regions[eno]->growingPlusTipList.size();tip++)
				growingTipsReg[2]+=1;

				geometry->patchArea[2]+=geometry->regions[eno]->area;
                        	geometry->RegionsIndex[2].push_back(eno);
			}
		}
		timeQueue.pushGlobal(VERY_LARGE,signalPPB);
        	nextPPBEventTime = VERY_LARGE;

        	break;
    }

    return;
}

EventDescriptorIndex System::registerEventDescriptor(EventDescriptor* ei)
{
    // maybe use a pool of numbers instead of an increasing number....this could wrap around
    EventDescriptorIndex index = EventDescriptorID++;
    EventDescriptorMap[index] = ei;

    if (EventDescriptorID == numeric_limits<EventDescriptorIndex>::max() - 10)
        timeQueue.pushGlobal(systemTime, stop);

    return index;
}

void System::unregisterEventDescriptor(EventDescriptorIndex ei)
{
    EventDescriptorMap.erase(EventDescriptorMap.find(ei));
    return;
}

EventDescriptor* System::getEventDescriptor(EventDescriptorIndex idx)
{
    map<EventDescriptorIndex, EventDescriptor*>::iterator itr = EventDescriptorMap.find(idx);
    if (itr == EventDescriptorMap.end())
        return NULL;
    else
        return itr->second;
}

#ifdef CROSS_SEV
void System::removeOccupiedIntersection(Intersection& is)
{
    OccupiedIntersectionList.RemoveElement(is.occupiedListPtr->index);
    is.occupiedListPtr = NULL;
    is.mirror->second.occupiedListPtr = NULL;
    return;
}

void System::addOccupiedIntersection(IntersectionItr is)
{
    is->second.occupiedListPtr = OccupiedIntersectionList.create(is);
    is->second.mirror->second.occupiedListPtr = is->second.occupiedListPtr;
    return;
}
#endif

void System::makeBinomialTable(void)
{
    int i(0),j(0);
    double temp(0.0);
    double factorial[MAXBINOM];
    temp=1;
    factorial[0] = 1;
    for (i=1 ; i<MAXBINOM ; i++)
    {
        temp *= i;
        factorial[i] = temp;
    }
    for (i=0 ; i<MAXBINOM ; i++)
        for (j=0 ; j<=i ; j++)
            binomialTable[i][j] = factorial[i]/(factorial[j]*factorial[i-j]);
    return;
}


