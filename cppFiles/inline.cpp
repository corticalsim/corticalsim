#ifndef NO_INLINE
inline
#endif

void System::collisionProbabilities(double angle, double& Pcat, double& Pzip)
{
	switch (p.interactionType)
	{
		case int_zipFirst:
                        
			// angle < magic angle 
			if(angle < PI*p.magicAngle/180.0)
			{
				Pcat = 0.0;

				// zippering enabled 
				if(p.zipperingEnabled)
					Pzip = p.zipFraction;
				else
					Pzip=0.0;
				break;
			}


			Pzip=0.0;

                        // angle > p.catStartAngle and induced catastrophes enabled
			if (p.catastrophesEnabled && (angle > p.catStartAngle/180.0*PI))
			{
				// proportional catastrophe enabled 
				if(p.proportionalCatastrophes)
					Pcat = (angle - p.catStartAngle/180.0*PI)/(PI/2. - p.catStartAngle/180.0*PI)*p.inducedCatastropheFraction;

				// proportional catastrophe disabled
				else
					Pcat = p.inducedCatastropheFraction;
			}
			else
				Pcat = 0.0;

		break;

		case int_catFirst:

			// angle > p.catStartAngle and induced catastrophes enabled 
			if (p.catastrophesEnabled && (angle > p.catStartAngle/180.0*PI))
			{
				// proportional catastrophe enabled 
				if(p.proportionalCatastrophes)
					Pcat = (angle - p.catStartAngle/180.0*PI)/(PI/2. - p.catStartAngle/180.0*PI)*p.inducedCatastropheFraction;

				// proportional catastrophe disabled
				else
					Pcat = p.inducedCatastropheFraction;
			}
			else
				Pcat = 0.0;

			// angle < magic angle and zippering enabled 
			if(p.zipperingEnabled && (angle < p.magicAngle/180.0*PI))
				Pzip = 1.0 - Pcat;
			else
				Pzip = 0.0;
		break;
	}

	#ifdef DBG_INTERACTION_TEST
	cout << "InteractionTest\t" << angle << "\t" << Pzip << "\t" << Pcat << "\t" << 1-Pzip-Pcat << "\n";
	#endif

	return;
}

#ifndef NO_INLINE
inline
#endif

double System::multiPcross(int Npar, int Ncoll, double xSingle, double zSingle)
{
	double totalSum = 0;

	int i,j;

	for (i=0 ; i<= Npar ; i++)
		for (j=max(0,i-Ncoll); j<=max(0,i-1); j++)
			totalSum += (Npar + 1 - i)*binomialTable[Ncoll][i-j]*binomialTable[max(0,i-1)][j]*pow(xSingle,2*j+Ncoll-i)*pow(zSingle,2*(i-j));
	return totalSum/(Npar+1);
}

#ifndef NO_INLINE
inline
#endif

double System::multiPzip(int nPar, int nColl, double xSingle, double zSingle)
{
	double totalSum = 0;

	int i,j,k;

	for (i=0 ; i<= nPar ; i++)
		for (j=0 ; j <= nColl-1 ; j++)
			for (k=max(0,i-j); k<=i; k++)
				totalSum += binomialTable[j][i-k]*binomialTable[i][k]*pow(xSingle,2*k+j-i)*pow(zSingle,2*(i-k)+1);
	return totalSum/(nPar+1);
}

#ifndef NO_INLINE
inline
#endif

CollisionType System::collisionResult(double angle, int nParallel, int nCross)	// range: 0..PI/2
{
	#ifdef DBG_ASSERT
	if ((angle <0) || (angle > PI/2))
	{
		cout << "ERROR: angle out of range [0..PI/2] : " << angle << "\n";
		exit(-2);
	}
	#endif

	#ifdef DBG_ASSERT
	{
		if (nCross == 0)
		{
			cerr << "ERROR: crossover collision passed to collisionResult()\n";
			return ct_crossover;
		}
	}
	#endif

	double Pcat;
	double Pzip;
	double Psum;

        // get a random number 
	double rand = randomGen.randDblExc();

        // get the available Pcat and Pzip
	collisionProbabilities(angle, Pcat, Pzip);

	switch(p.bundleType)
	{
		case bdl_noZip:

			if (nParallel > 0)
				if (rand < Pcat)
					return ct_inducedCatastrophe;
		// break omitted on purpose!

		case bdl_sticky:

			if (nParallel > 0)
				return ct_crossover;
		// break omitted on purpose!

		case bdl_simple:

			// induced catastrophe selected 
			if (rand < Pcat)
				return ct_inducedCatastrophe;

			// zippering selected
			else if (rand < Pcat + Pzip)
				return ct_zipper;

			// cross over slected
			else
				return ct_crossover;

		// though never reached
		//break;	

		case bdl_Ncollision:

                        // add up total probabilities
			Pcat *= nCross;
			Pzip *= nCross;
			Psum = Pcat + Pzip;

                        // now normalize the probabilities 
			if (Psum > 1)
			{
				Pcat /= Psum;
				Pzip /= Pzip;
			}

			if (rand < Pcat)
				return ct_inducedCatastrophe;

			else if (rand < Pcat + Pzip)
				return ct_zipper;

			else
				return ct_crossover;

		// though never reached
		//break;	

		case bdl_multiCollision:

			double Pcross = 1.0 - Pcat - Pzip;
			double multiPx = multiPcross(min(nParallel,MAXBINOM-1), min(nCross,MAXBINOM-1), Pcross, Pzip);
			double multiPz = multiPzip(min(nParallel,MAXBINOM-1), min(nCross,MAXBINOM-1), Pcross, Pzip);

			if (rand < multiPx)
				return ct_crossover;

			else if (rand < multiPx + multiPz)
				return ct_zipper;

			else
				return ct_inducedCatastrophe;

		// though never reached
		// break;	

	}

	// should never occur
	return ct_crossover;	
}

#ifndef NO_INLINE
inline
#endif

double Microtubule::length()
{
	double temp = 0.0;
	Segment* seg = segments.first();

	// calculate total length of a MT 
	while (seg != NULL)
	{
		temp += seg->length();
		seg = seg->next();
	}
	return temp;
}

#ifndef NO_INLINE
inline
#endif

void Microtubule::updateLength(bool forceUpdate)
{
        // can go wrong on position cache refresh, so need a better solution
	if ((previousUpdateTag == system->currentTimeTag) && (!forceUpdate))
		return;

	#ifdef DBG_ACID_TEST
	if (segments.size() == 0)
		cerr << "ERROR: no segments left. address=" << this << "\n";
	if (plus.event.queue->progression(previousUpdateTag) < -ZERO_CUTOFF)
		cerr << "alarm\n";
	if (minus.event.queue->progression(previousUpdateTag) < -ZERO_CUTOFF)
		cerr << "alarm\n";
	#endif

        // move the position of the active (first and last) segment end
	segments.last()->end += plus.event.queue->progression(previousUpdateTag)*plus.dir*plus.velocity;
	segments.first()->start += minus.event.queue->progression(previousUpdateTag)*minus.dir*minus.velocity;

	previousUpdateTag = system->currentTimeTag;

	#ifdef DBG_ACID_TEST
	if (segments.size() == 1)
	{
		if ((segments.first()->end - segments.first()->start)*segments.first()->dir < -ZERO_CUTOFF)
		{
			cerr << "ERROR: negative segment length after updating. [mt address=" << this << ", length=" << (segments.first()->end - segments.first()->start)*segments.first()->dir << "]\n";
			cerr << "comparison value:" << -ZERO_CUTOFF << "\n";
			cerr << "previous update " << previousUpdateTag  << " current " << system->currentTimeTag << " progression=" << system->timeQueue.progression(previousUpdateTag) << "\n";
			exit(-1);
		}
	}
	#endif

	return;
}

#ifndef NO_INLINE
inline
#endif

void Trajectory::conditionalRemove()
{
	// if there is no segments and no tips, associted with a trajectory, then remove it rigth way
	if (segments.empty() && (notificationList.size()==0))
	base.region->removeTrajectory(this); 
	
	return;
}

#ifndef NO_INLINE
inline
#endif

void Region::updateRegionLength(bool forceUpdate)
{
	// can go wrong on position cache refresh, a cleaner solution required
	if ((previousUpdateTag == geometry->system->currentTimeTag) && (!forceUpdate))
		return;

	totalLength += geometry->system->timeQueue.progression(previousUpdateTag)
		*(geometry->system->p.vMin*(shrinkingPlusTipList.size()) - geometry->system->p.vTM*(minusTipList.size()) )
			+ geometry->system->vPlusQueue.progression(previousUpdateTag)*(geometry->system->p.vPlus*(growingPlusTipList.size()));

	previousUpdateTag = geometry->system->currentTimeTag;

	return;
}

