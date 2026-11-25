#include "microtubule.hpp"
#include "region.hpp"
#include "system.hpp"

double Microtubule::length()
{
    double temp = 0.0;
    Segment* seg = segments.first();

    // calculate total length of a MT
    while (seg != nullptr)
    {
        temp += seg->length();
        seg = seg->next();
    }
    return temp;
}

void Microtubule::update_length(bool force_update)
{
    // can go wrong on position cache refresh, so need a better solution
    if ((previousUpdateTag == system->currentTimeTag) && (!force_update))
    {
        return;
    }

#ifdef DBG_ACID_TEST
    if (segments.size() == 0)
    {
        cerr << "ERROR: no segments left. address=" << this << "\n";
    }
    if (plus.event.queue->progression(previousUpdateTag) < -ZERO_CUTOFF)
    {
        cerr << "alarm\n";
    }
    if (minus.event.queue->progression(previousUpdateTag) < -ZERO_CUTOFF)
    {
        cerr << "alarm\n";
    }
#endif

    // move the position of the active (first and last) segment end
    segments.last()->end
    += plus.event.queue->progression(previousUpdateTag) * static_cast<double>(plus.dir) * plus.velocity;
    segments.first()->start
    += minus.event.queue->progression(previousUpdateTag) * static_cast<double>(minus.dir) * minus.velocity;

    previousUpdateTag = system->currentTimeTag;

#ifdef DBG_ACID_TEST
    if (segments.size() == 1)
    {
        if ((segments.first()->end - segments.first()->start) * segments.first()->dir < -ZERO_CUTOFF)
        {
            cerr << "ERROR: negative segment length after updating. [mt address=" << this
                 << ", length=" << (segments.first()->end - segments.first()->start) * segments.first()->dir << "]\n";
            cerr << "comparison value:" << -ZERO_CUTOFF << "\n";
            cerr << "previous update " << previousUpdateTag << " current " << system->currentTimeTag
                 << " progression=" << system->timeQueue.progression(previousUpdateTag) << "\n";
            exit(-1);
        }
    }
#endif

    return;
}

Microtubule::Microtubule(System* s, TrajectoryVector tv, bool initialize):
    system(s),
    plus(this, tv, s->p.vPlus, &(s->vPlusQueue), s->p.vPlus),

    // reverse direction
    minus(this, tv.flipped(), -s->p.vTM, &(s->timeQueue), -s->p.vTM),
    disappearEvent(this, &(s->timeQueue), -s->p.vMin + s->p.vTM),
    type(MT_GROWING),
    nucleationTime(s->systemTime + s->systemTimeOffset),
    previousUpdateTag(s->currentTimeTag)
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule created.\n";
#endif

    if (initialize)
    {
        // create a single segment MT
        segments.create(this, tv);

        // intialize the plus tip of MT
        plus.initialize();

        // initialize the minus tip of MT
        minus.initialize();

        // microtubule is disappearing (segment length finishing)
        set_disappear_event();
    }
    return;
}

Microtubule::~Microtubule()
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule destroyed. address=" << this << "\n";
#endif

#ifdef DBG_ASSERT
    if ((segments.first()->start - segments.last()->end) > ZERO_CUTOFF)
    {
        cerr << "DBG/ASSERT: Microtubule::~Microtubule: endpoints do not coincide/\n";
        exit(-1);
    }
#endif

    // un-linking the MT from the trajectory(ies)
    plus.unlink_from_trajectory();
    minus.unlink_from_trajectory();

    return;
}

void Microtubule::set_disappear_event()
{
    // push to disappear event: (plus.velocity - minus.velocity < 0)
    if ((segments.size() == 1) && (type == MT_SHRINKING))
    {
        disappearEvent.push_on_queue(segments.first()->length(), EV_DISAPPEAR);
    }

    // flush disappear event queue
    else
    {
        disappearEvent.clear();
    }
    return;
}

void Microtubule::handle_event(const EventDescriptor* ei)
{
    // better update length before processing next event
    update_length();

#ifdef DBG_ACID_TEST
    if (!integrityCheck())
    {
        cerr << "failed before event\n";
        exit(-1);
    }
#endif

    switch (ei->type)
    {
        // wall collision event encountered
        case EV_WALL:
            wall();
            break;

            // MT-MT collision event encountered
        case EV_COLLISION:
            collision();
            break;

        // back track event encountered
        case EV_BACKTRACK:
            if (ei == &(plus.event))
            {
                backtrack(&plus);
            }
            else
            {
                backtrack(&minus);
            }
            break;

            // end of segment event encountered
        case EV_END_OF_SEGMENT:
            if (ei == &(plus.event))
            {
                end_of_segment(&plus);
            }
            else
            {
                end_of_segment(&minus);
            }
            break;

        // disappear event encountered
        case EV_DISAPPEAR:
            harakiri();
            break;

        case EV_NONE:
            break;
    };

#ifdef DBG_ACID_TEST
    if ((ei->type != ev_disappear) && (!integrityCheck()))
    {
        cerr << "failed after event\n";
        exit(-1);
    }
#endif

    return;
}

void Microtubule::catastrophe()
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule::catastrophe() called.\n";
#endif

#ifdef DBG_EXTRA_CHECK
    if (type != mt_growing)
    {
        cout << "ERROR: catastrophe called for non-growing MT\n";
        exit(-1);
    }
#endif

    // update MT length
    update_length();

    // growing MT will become shrinking MT, hence unregister the tip from growing-plus-tip list of this region
    plus.trajectory->base.region->unregister_from_region(plus.regionTag, T_PLUS, MT_GROWING);

    // now, this is shrinking MT
    type = MT_SHRINKING;

    // growing MT will become shrinking MT, hence register the tip on shrinking-plus-tip list of this region
    plus.regionTag = plus.trajectory->base.region->register_on_region(&plus, T_PLUS, MT_SHRINKING);

    // replace plus velocity by minus velocity
    plus.velocity = system->p.vMin;

    //
    plus.event.reinitialize(&(system->timeQueue), system->p.vMin);

    // import the growing MT as shrinking MT
    system->shrinking_mts.import(system->growing_mts, this);

    if (system->p.vMin < 0)
    {
        // velocity has changed sign, select a new event
        plus.advance_intersection();
    }

    // schedule next event, (if this is triggered by a collision, and vMin > 0, increase occupancy)
    plus.determine_event();

    // check possibility of a disappear event
    set_disappear_event();

    return;
}

void Microtubule::rescue()
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule::rescue() called.\n";
#endif

#ifdef DBG_EXTRA_CHECK
    if (type != mt_shrinking)
    {
        cout << "ERROR: rescue called for non-shrinking MT\n";
        exit(-1);
    }
#endif

    update_length();

    // shrinking MT will become growing  MT, hence unregister the tip from shrinking-plus-tip list of this region
    plus.trajectory->base.region->unregister_from_region(plus.regionTag, T_PLUS, MT_SHRINKING);

    // now, this is growing MT
    type = MT_GROWING;

    // shrinking MT will become growing MT, hence register the tip on growing-plus-tip list of this region
    plus.regionTag = plus.trajectory->base.region->register_on_region(&plus, T_PLUS, MT_GROWING);

    // replace minue velocity by plus velocity
    plus.velocity = system->p.vPlus;

    plus.event.reinitialize(&(system->vPlusQueue), system->p.vPlus);

    // import the shrinking MT as growing MT
    system->growing_mts.import(system->shrinking_mts, this);

    if (system->p.vMin < 0)
    {
        // velocity has changed sign, need to select a new event
        plus.advance_intersection();
    }

    // schedule next event
    plus.determine_event();

    // check possibility of a disappear event
    set_disappear_event();

    return;
}

void Microtubule::sever(Segment* cut_seg, double cut_pos)
{
    // sever MT at given position
    split_segment_at_traj_pos(cut_seg->start + static_cast<double>(cut_seg->dir) * cut_pos, cut_seg);
    system->totalLengthSeveringCount++;
    return;
}

void Microtubule::sever_at_cross(IntersectionItr is, Segment* cut_seg)
{
    // WARNING: take care to cut the segment right *after* the intersection. Not doing so can lead to problems
    // when treadmilling is disabled. The microtubule end is then located *exactly* on the intersection and can be cut
    // again and again, leading to problems.

    // this time, cutPos and cutSeg are already known.
    double split_pos(is->first);
    // shift the cut ever so slightly away from the intersection, in the direction of the plus end.
    if (cut_seg->dir == ::FORWARD)
    {
        split_pos += ZERO_CUTOFF;
    }
    else
    {
        split_pos -= ZERO_CUTOFF;
    }

    // only cut if the severing position is comfortably located within the segment interval
    if ((cut_seg->start - split_pos) * (cut_seg->end - split_pos) < ZERO_CUTOFF)
    {
        split_segment_at_traj_pos(split_pos, cut_seg);
        system->totalIntersectionSeveringCount++;
    }
    return;
}

void Microtubule::split_segment_at_traj_pos(double cut_pos, Segment* cut_seg)
{
    // used in sever and severAtCross
    Microtubule* new_mt(nullptr);
    TrajectoryVector tv;

    tv.trajectory = cut_seg->trajectory;
    tv.dir = cut_seg->dir;
    tv.pos = cut_pos;
    if (type == MT_GROWING)
    {
        new_mt = system->growing_mts.create(system, tv, false);
    }
    else
    {
        new_mt = system->shrinking_mts.create(system, tv, false);
    }

    new_mt->nucleationTime = nucleationTime;
    new_mt->type = type;

    // transfer segments
    new_mt->segments.create(new_mt, tv);
    new_mt->segments.first()->end = cut_seg->end;
    new_mt->segments.first()->endItr = cut_seg->endItr;
    new_mt->segments.first()->nucleationTime = cut_seg->nucleationTime;
    new_mt->segments.import_set(segments, cut_seg->next(), segments.last());
    Segment* temp_seg = new_mt->segments.first()->next();
    while (temp_seg != nullptr)
    {
        temp_seg->mt = new_mt;
        temp_seg = temp_seg->next();
    }

    // now update the plus tip trajectory/position
    new_mt->plus.trajectory = plus.trajectory;
    new_mt->plus.velocity = plus.velocity;
    new_mt->plus.event.reinitialize(plus.event.queue, plus.velocity);
    new_mt->plus.dir = plus.dir;
    new_mt->plus.nextCollision = plus.nextCollision;
    new_mt->plus.notificationTag = new_mt->plus.trajectory->register_for_notifications(&(new_mt->plus));
    new_mt->plus.regionTag = new_mt->plus.trajectory->base.region->register_on_region(&(new_mt->plus), T_PLUS, type);

    // next, update the minus end of the new MT
    new_mt->segments.first()->start = tv.pos;
    new_mt->minus.notificationTag = new_mt->minus.trajectory->register_for_notifications(&(new_mt->minus));
    new_mt->minus.regionTag = new_mt->minus.trajectory->base.region->register_on_region(&(new_mt->minus), T_MINUS, type);
    new_mt->minus.locate_intersection();

    new_mt->plus.determine_event();
    new_mt->minus.determine_event();
    new_mt->set_disappear_event();

    // now, finish up the base microtubule
    segments.last()->end = cut_pos;

    // always unregister from region (just to be sure)
    plus.trajectory->base.region->unregister_from_region(plus.regionTag, T_PLUS, type);
    if (type == MT_GROWING)
    {
        type = MT_SHRINKING;
        plus.velocity = system->p.vMin;
        plus.event.reinitialize(&(system->timeQueue), system->p.vMin);
        system->shrinking_mts.import(system->growing_mts, this);
    }

    if (plus.trajectory != tv.trajectory)
    {
        plus.trajectory->unregister_for_notifications(plus.notificationTag);
        plus.trajectory = tv.trajectory;
        plus.notificationTag = plus.trajectory->register_for_notifications(&(plus));
    }
    //...and register on region again
    plus.regionTag = plus.trajectory->base.region->register_on_region(&plus, T_PLUS, MT_SHRINKING);
    plus.dir = tv.dir;
    plus.nextCollision = new_mt->minus.nextCollision;
    if (system->p.vMin * new_mt->minus.velocity >= 0)
    {
        // both tips move in opposite directions [this should generally be the case]
        plus.advance_intersection();
    }

    plus.determine_event();
    if (segments.size() == 1)
    {
        minus.determine_event();
    }
    set_disappear_event();

    return;
}

void Microtubule::translate_position_m_t2_segment(double& cut_pos, Segment*& cut_seg)
{
    cut_seg = segments.first();
    while ((cut_seg != nullptr) && (cut_pos > cut_seg->length()))
    {
        cut_pos -= cut_seg->length();
        cut_seg = cut_seg->next();
    }
    if (cut_pos > cut_seg->length())
    {
        if (cut_pos - cut_seg->length() < ZERO_CUTOFF)
        {
            cut_pos = cut_seg->length() - ZERO_CUTOFF;
        }
        else
        {
            cerr << "ERROR: position beyond the end of the microtubule.\n";
            exit(-1);
        }
    }
    cut_pos = cut_seg->start + static_cast<double>(cut_seg->dir) * cut_pos;
    return;
}

void Microtubule::wall()
{
    // only for plus end

    TrajectoryVector tv;

    // set position of the plus end equal to the event position, to mitigate the effect of cumulating addition errors
    segments.last()->end = plus.nextEventPos;
    tv = plus.trajectory->next_trajectory(plus.dir);

    // next trajectory on 'forbidden zone', leads catastrophe
    if (tv.trajectory->base.region->faceTag == -1)
    {
        tv.trajectory->conditional_remove();
        catastrophe();

        // no wall crossing
        return;
    }

    // edge catastrophe enabled
    if (system->p.edgeCatastropheEnabled)
    {
        double cos_angle = (plus.dir == ::FORWARD) ? plus.trajectory->nextTrCosAngle : plus.trajectory->prevTrCosAngle;

        // MT is moving non-parallel to the edege
        if ((1.0 - cos_angle) > ZERO_CUTOFF)
        {
            double p_cat(0.0);

            // for sharp edged cube
            if (system->p.geometry == "cubeReal")
            {
                if (abs(plus.trajectory->base.region->faceTag - tv.trajectory->base.region->faceTag) > 0)
                {
                    if (abs(plus.trajectory->base.region->faceTag - tv.trajectory->base.region->faceTag) > 1)
                    {
                        p_cat = system->p.pCatSpecialEdgeMax;
                    }

                    else
                    {
                        p_cat = system->p.pCatRegularEdgeMax;
                    }
                }
            }

            // other
            else
            {
                p_cat = (plus.dir == ::FORWARD) ? plus.trajectory->nextTrpCat : plus.trajectory->prevTrpCat;
            }

            // when smooth catastrophe is enabled
            if (system->p.edgeCatastropheSmooth)
            {
                p_cat *= (1.0 - cos_angle);
            }

            // apply catastrophe
            if (system->randomGen.rand() <= p_cat)
            {
                tv.trajectory->conditional_remove();
                catastrophe();

                // no wall crossing
                return;
            }
        }
    }

    // update boundarty crossing count
    system->boundaryCrossingCount++;

#ifdef DBG_ASSERT
    if (segments.last()->endItr != plus.nextCollision)
    {
        cerr << "DBG/ASSERT: nextCollision does not match with endItr upon entering wall event.\n";
    }
#endif

    // wall crossing
    segments.last()->endItr = plus.nextCollision;

    // create a new segment
    segments.create(this, tv);
    segments.last()->startItr = tv.dir == ::FORWARD ? tv.trajectory->wall_begin() : tv.trajectory->wall_end();

    // switch trajectory
    plus.switch_trajectory(
    tv.trajectory, tv.dir, tv.dir == ::FORWARD ? tv.trajectory->wall_begin() : tv.trajectory->wall_end());

    if (segments.size() == 2)
    {
        // MT with two segment of zero (=0) length, consider possibilty of a disappear event
        minus.determine_event();
        set_disappear_event();
    }

    return;
}

void Microtubule::collision()
{
    // set position of the plus end equal to the event position, to mitigate the effect of cumulating addition errors
    segments.last()->end = plus.nextEventPos;
    CollisionType type;
    Direction dir;

    // no MT(s) passing throught the next collision point, resulting a simple crossover
    if (plus.nextCollision->second.mirror->second.occupancy == 0)
    {
        type = CT_CROSSOVER;
    }

    // MT(s) passing through the next collision point
    else
    {
        double diff_angle
        = plus.trajectory->base.region->intersection_angle(plus.trajectory, plus.nextCollision->second.otherTrajectory);

        if (diff_angle > 0.5 * PI)
        {
            dir = (plus.dir == ::FORWARD) ? BACKWARD : ::FORWARD;
            diff_angle = PI - diff_angle;
        }
        else
        {
            dir = plus.dir;
        }

        // find out type of collision event to be occured
        type = system->collision_result(
        diff_angle, plus.nextCollision->second.occupancy, plus.nextCollision->second.mirror->second.occupancy);
    }

    switch (type)
    {
            // zippering
        case CT_ZIPPER:
            system->totalZipperCount++;
            zipper(dir);
            break;

            // crossover
        case CT_CROSSOVER:
            system->totalCrossoverCount++;
            crossover();
            break;

            // induced catastrophe
        case CT_INDUCED_CATASTROPHE:
            system->totalInducedCatastropheCount++;
            catastrophe();
            break;

#ifdef DBG_EXTRA_CHECK
        default:
            cerr << "Unknown collision result type (" << type << ") used.\n";
            exit(-1);
#endif
    };

    return;
}

void Microtubule::zipper(Direction dir)
{
#ifdef DBG_ASSERT
    if (abs(plus.nextCollision->first - segments.last()->end) > ZERO_CUTOFF)
    {
        cout << "DBG/ASSERT: ERROR: end of segment does not match zippering location.[discrepancy: "
             << abs(plus.nextCollision->first - segments.last()->end) << "]\n";
        exit(-1);
    }
#endif

    // get a trajectory vetor on the new trajectory to be associted after zippering
    TrajectoryVector tv(plus.nextCollision->second.mirror->first, dir, plus.nextCollision->second.otherTrajectory);

    // link the end of the last segment (containing the plus tip) to the next plus end collision
    segments.last()->endItr = plus.nextCollision;

    // create a new zero (=0) length segment on new position and direction (of course on new trajectory)
    segments.create(this, tv);

    // link the bottom of the last segment (containing the new minus tip) to the next plus end collision mirror
    segments.last()->startItr = plus.nextCollision->second.mirror;

    // now do the switching of trajectory from incoming to outgoing
    plus.switch_trajectory(plus.nextCollision->second.otherTrajectory, dir, plus.nextCollision->second.mirror);

    if (segments.size() == 2)
    {
        // number of segment is two (=2), HERE it means two zero length segment, prepare to schedule a disappear event
        minus.determine_event();
        set_disappear_event();
    }

    return;
}

void Microtubule::end_of_segment(MTTip* tip)
{
#ifdef DBG_ASSERT
    if (segments.size() == 1)
    {
        cout << "DBG/ASSERT: ERROR: end-of-segment event called for single-segmented microtubule.\n";
        // exit(-1);
    }
#endif

    Direction new_dir;
    Segment* new_seg(nullptr);
    Segment* kill_seg(nullptr);

    // select plus tip associated segment
    if (tip == &plus)
    {
        kill_seg = segments.last();
        new_seg = segments.last()->previous();
        new_dir = new_seg->dir;
    }

    // select minus tip associated segment
    else
    {
        kill_seg = segments.first();
        new_seg = segments.first()->next();
        if (new_seg->dir == ::FORWARD)
        {
            new_dir = BACKWARD;
        }
        else
        {
            new_dir = ::FORWARD;
        }
    }

#ifdef DBG_ASSERT
    if (abs(killSeg->end - killSeg->start) > ZERO_CUTOFF)
    {
        cout << "DEBUG/ASSERT: ERROR: non-zero segment length on end of segment event\n";
        exit(-1);
    }
#endif

    // kill the selected segment
    segments.remove(kill_seg);

    // tip switching to a new trajectory (in a different region) after encountering wall-begin/wall-end
    if ((tip->nextCollision == tip->trajectory->wall_begin()) || (tip->nextCollision == tip->trajectory->wall_end()))
    {
        tip->switch_trajectory(new_seg->trajectory,
                              new_dir,
                              new_dir == ::FORWARD ? new_seg->trajectory->wall_end() : new_seg->trajectory->wall_begin());
        system->boundaryCrossingCount--;
    }

    // tip switching to a new trajectory (in the same region)
    else
    {
        tip->switch_trajectory(tip->nextCollision->second.otherTrajectory, new_dir, tip->nextCollision->second.mirror);
    }

    // MT length zero (=0) encountered
    if (segments.size() == 1)
    {
        // current tip is a plus tip, so determine the minus tip associated event
        if (tip == &(plus))
        {
            minus.determine_event();
        }

        // current tip is a minus tip, so determine the plus tip associated event
        else
        {
            plus.determine_event();
        }

        // schedule a disappear event
        set_disappear_event();
    }

    return;
}

void Microtubule::backtrack(MTTip* tip)
{
    // asociated with depolymerizing plus/minus end

    if (tip == &(plus))
    {
        segments.last()->end = plus.nextEventPos;
    }
    else
    {
        segments.first()->start = minus.nextEventPos;
    }

    // MT is back tracked from this site, resulting a decrease in occupancy
    tip->nextCollision->second.occupancy--;

#ifdef CROSS_SEV
    if ((tip->nextCollision->second.occupancy == 0) && (tip->nextCollision->second.mirror->second.occupancy > 0))
    {
        system->remove_occupied_intersection(tip->nextCollision->second);
    }
#endif

#ifdef DBG_ASSERT
    if (tip->nextCollision->second.occupancy < 0)
    {
        cerr << "DBG/ASSERT: ERROR: negative intersection occupancy after ";
        if (tip == &(plus))
        {
            cerr << "plus";
        }
        else
        {
            cerr << "minus";
        }
        cerr << " end backtracking.\n";

        cout << "number of equal elements: " << tip->trajectory->intersections.count(tip->nextCollision->first) << "\n";
        ;

        double temp = tip->nextCollision->first;

        do
        {
            cout << "pos=" << tip->mt->segments.last()->end << ", " << tip->mt->segments.first()->start
                 << ", occupancy=" << tip->nextCollision->second.occupancy << "  ";
            tip->advanceIntersection();
            tip->determineEvent();
            cout << "step (size=" << tip->nextCollision->first - temp << ")\n";
        } while (tip->event.type != ev_end_of_segment);
        cout << "pos=" << tip->mt->segments.last()->end << ", " << tip->mt->segments.first()->start << "  \n";

        exit(-1);
    }
#endif

    // schedule next event
    tip->advance_intersection();
    tip->determine_event();

    return;
}

void Microtubule::crossover()
{
#ifdef DBG_ASSERT
    if (plus.nextCollision->second.occupancy < 0)
    {
        cerr << "DBG/ASSERT: ERROR: negative intersection occupancy before crossover.\n";
        exit(-1);
    }
#endif

#ifdef CROSS_SEV
    if ((plus.nextCollision->second.occupancy == 0) && plus.nextCollision->second.mirror->second.occupancy > 0)
    {
        system->add_occupied_intersection(plus.nextCollision);
    }
#endif

    // skip any (active) event in the current intersection point, instead move the tip to the next intersection point
    plus.nextCollision->second.occupancy++;
    plus.advance_intersection();
    plus.determine_event();

    return;
}

void Microtubule::harakiri()
{
    // MT disappear event encountered, remove the MT
    if (type == MT_GROWING)
    {
        system->growing_mts.remove(this);
    }
    else
    {
        system->shrinking_mts.remove(this);
    }

    return;
}

bool Microtubule::integrity_check()
{
    bool valid(true);
    double max_deviation(0.0);

    Segment* seg(nullptr);
    IntersectionItr is;
    size_t count(0);

    // check connections between end-of segment iterators and match with positions (should be exact)
    if (segments.size() == 0)
    {
        // MT should contain at least one segment
        valid = false;
        cerr << "Microtubule contains no segments.\n";
    }

    else if (segments.size() != 1)
    {
        if (segments.size() == 0)
        {
            // MT should contain at least one segment
            cerr << "microtubule without segments!\n";
            exit(-1);
        }

        seg = segments.first();
        if (seg == nullptr)
        {
            // pointer to segment list can not be null
            cerr << "null pointer received\n";
            exit(-1);
        }

        // segment (some where) at the middle of a trajectory
        if ((seg->endItr != seg->trajectory->wall_end()) && (seg->endItr != seg->trajectory->wall_begin()))
        {
            if (seg->endItr->second.mirror != seg->next()->startItr)
            {
                // next-start-iterator should match with the end-second-mirror-iterator
                valid = false;
                cerr << "reciprocity violated\n";
            }

            if (abs(seg->end - seg->endItr->first) > ZERO_CUTOFF)
            {
                // segment end location does not coincide with the first intersection site
                valid = false;
                cerr << "Location of segment end does not coincide with intersection position\n";
            }

            max_deviation = max(max_deviation, abs(seg->end - seg->endItr->first));

            if (abs(seg->next()->start - seg->next()->startItr->first) > ZERO_CUTOFF)
            {
                // (next) segment end location does not coincide with the first intersection site
                valid = false;
                cerr << "Location of segment end does not coincide with intersection position\n";
            }
            max_deviation = max(max_deviation, abs(seg->next()->start - seg->next()->startItr->first));
        }

        // segment at the end of a trajectory
        else
        {
            // positions should match with wall end location
            if (seg->endItr == seg->trajectory->wall_end())
            {
                if (abs(seg->end - seg->trajectory->length) > ZERO_CUTOFF)
                {
                    // segment end location does not coincide with the first intersection site
                    valid = false;
                    cerr << "Segment end location does not coincide with end of trajectory\n";
                }
                max_deviation = max(max_deviation, abs(seg->end - seg->trajectory->length));
            }

            // positions should match with wall begin location
            else
            {
                if (abs(seg->end - 0) > ZERO_CUTOFF)
                {
                    // segment end location does not coincide with the first intersection site
                    valid = false;
                    cerr << "Segment end location does not coincide with end of trajectory\n";
                }
                max_deviation = max(max_deviation, abs(seg->end - 0));
            }

            if (seg->next()->startItr == seg->next()->trajectory->wall_end())
            {
                if (abs(seg->next()->start - seg->next()->trajectory->length) > ZERO_CUTOFF)
                {
                    valid = false;
                    cerr << "Segment start location does not coincide with end of trajectory\n";
                    cerr << "start: " << seg->next()->start << ", length: " << seg->next()->trajectory->length
                         << ", diff=" << abs(seg->next()->start - seg->next()->trajectory->length) << "\n";
                }
                max_deviation = max(max_deviation, abs(seg->next()->start - seg->next()->trajectory->length));
            }
            else
            {
                if (abs(seg->next()->start - 0) > ZERO_CUTOFF)
                {
                    valid = false;
                    cerr << "Segment start location does not coincide with start of trajectory\n";
                    cerr << "start: " << seg->next()->start << ", length: " << seg->next()->trajectory->length
                         << ", diff=" << abs(seg->next()->start) << "\n";
                }
                max_deviation = max(max_deviation, abs(seg->next()->start - 0));
            }
        }
        seg = seg->next();
    }

    if (!valid)
    {
        cout << "Broken segment connections\n";
    }

    // check begin and end points
    if (static_cast<double>(plus.dir) * plus.velocity > 0)
    {
        if (plus.nextCollision == plus.trajectory->wall_end())
        {
            // segment end location does not coincide with the first intersection site
            if (segments.last()->end - plus.trajectory->length > ZERO_CUTOFF)
            {
                valid = false;
            }
        }

        else
        {
            // segment end location does not coincide with the first intersection site
            if (segments.last()->end - plus.nextCollision->first > ZERO_CUTOFF)
            {
                valid = false;
            }
        }
    }

    else
    {
        if (plus.nextCollision == plus.trajectory->wall_begin())
        {
            // segment end location does not coincide with the first intersection site
            if (segments.last()->end < -ZERO_CUTOFF)
            {
                valid = false;
            }
        }
        else
        {
            // segment end location does not coincide with the first intersection site
            if (segments.last()->end - plus.nextCollision->first < -ZERO_CUTOFF)
            {
                valid = false;
            }
        }
    }

    if (!valid)
    {
        cerr << "ERROR: Problem with ordering of plus ends and next collisions.\n";
    }

    // minus end velocity can be zero in the absence of treadmilling, that case skip this check to avoid numerical
    // directionality problems
    if (static_cast<double>(minus.dir) * minus.velocity > ZERO_CUTOFF)
    {
        if (minus.nextCollision == minus.trajectory->wall_end())
        {
            if (segments.first()->start - minus.trajectory->length > ZERO_CUTOFF)
            {
                valid = false;
            }
        }
        else
        {
            if (segments.first()->start - minus.nextCollision->first > ZERO_CUTOFF)
            {
                valid = false;
            }
        }
    }
    else if (static_cast<double>(minus.dir) * minus.velocity < -ZERO_CUTOFF)
    {
        if (minus.nextCollision == minus.trajectory->wall_begin())
        {
            if (segments.first()->start < -ZERO_CUTOFF)
            {
                valid = false;
            }
        }
        else
        {
            if (segments.first()->start - minus.nextCollision->first < -ZERO_CUTOFF)
            {
                valid = false;
            }
        }
    }

    if (!valid)
    {
        cerr << "ERROR: Problem with ordering of minus ends and next collisions.\n";
    }

    // check validity of intersections by running them to the end of the trajectory
    seg = segments.first();
    while (seg != nullptr)
    {
        is = seg->startItr;
        count = 0;
        while ((is != seg->trajectory->wall_end()) && (count != seg->trajectory->intersections.size()))
        {
            is++;
            count++;
        }
        if (is != seg->trajectory->wall_end())
        {
            valid = false;
            cerr << "Error: segment start iterator does not belong to correct trajectory\n";
            cerr << "Number of segments in MT: " << segments.size() << "\n";
            ;
        }
        is = seg->endItr;
        count = 0;
        while ((is != seg->trajectory->wall_end()) && (count != seg->trajectory->intersections.size()))
        {
            is++;
            count++;
        }
        if (is != seg->trajectory->wall_end())
        {
            valid = false;
            cerr << "Error: segment end iterator does not belong to correct trajectory\n";
            cerr << "Number of segments in MT: " << segments.size() << "\n";
            ;
        }
        seg = seg->next();
    }

    // check the tip's tipmap pointer (round-trip)
    if (system->get_event_descriptor(plus.event.index) != &plus.event)
    {
        valid = false;
    }
    if (system->get_event_descriptor(minus.event.index) != &minus.event)
    {
        valid = false;
    }
    if (system->get_event_descriptor(disappearEvent.index) != &disappearEvent)
    {
        valid = false;
    }

    // round-trip check of the segment-trajectory pairs
    seg = segments.first();
    while (seg != nullptr)
    {
        if (*(seg->trajectoryTag) != seg)
        {
            valid = false;
        }

        seg = seg->next();
    }

    // round-trip check of tip-trajectory pairs
    if (*(plus.notificationTag) != &plus)
    {
        valid = false;
    }
    if (*(minus.notificationTag) != &minus)
    {
        valid = false;
    }

    // check orientation of segments / start-end, and tips
    if (segments.first()->dir == minus.dir)
    {
        valid = false;
    }
    if (segments.last()->dir != plus.dir)
    {
        valid = false;
    }

    seg = segments.first();
    while (seg != nullptr)
    {
        if (seg->dir == ::FORWARD)
        {
            if (seg->end - seg->start < -ZERO_CUTOFF)
            {
                cerr << "false directionality in positions. diff=" << seg->end - seg->start << "\n";
                valid = false;
            }
            is = seg->startItr;
            count = 0;
            while ((is != seg->endItr) && (count != seg->trajectory->intersections.size()))
            {
                is++;
                count++;
            }
            if (is != seg->endItr)
            {
                cerr << "false directionality in intersection iterators\n";
                valid = false;
            }
        }
        else
        {
            if (seg->end - seg->start > ZERO_CUTOFF)
            {
                cerr << "false directionality in positions. diff=" << seg->end - seg->start << "\n";
                valid = false;
            }
            is = seg->startItr;
            count = 0;
            while ((is != seg->endItr) && (count != seg->trajectory->intersections.size()))
            {
                is--;
                count++;
            }
            if (is != seg->endItr)
            {
                cerr << "false directionality in intersection iterators\n";
                valid = false;
            }
        }

        seg = seg->next();
    }

    // check whether segment end iterators are close to tip->nextCollision

    return valid;
}
