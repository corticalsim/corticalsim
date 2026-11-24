#include "microtubule.hpp"
#include "region.hpp"
#include "system.hpp"

double Microtubule::length()
{
    double temp = 0.0;
    Segment* seg = m_segments.first();

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
    if ((m_previous_update_tag == m_system->m_current_time_tag) && (!force_update))
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
    m_segments.last()->m_end
    += m_plus.m_event.m_queue->progression(m_previous_update_tag) * static_cast<double>(m_plus.m_dir) * m_plus.m_velocity;
    m_segments.first()->m_start
    += m_minus.m_event.m_queue->progression(m_previous_update_tag) * static_cast<double>(m_minus.m_dir) * m_minus.m_velocity;

    m_previous_update_tag = m_system->m_current_time_tag;

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
    m_system(s),
    m_plus(this, tv, s->m_p.m_v_plus, &(s->m_v_plus_queue), s->m_p.m_v_plus),

    // reverse direction
    m_minus(this, tv.flipped(), -s->m_p.m_v_tm, &(s->m_time_queue), -s->m_p.m_v_tm),
    m_disappear_event(this, &(s->m_time_queue), -s->m_p.m_v_min + s->m_p.m_v_tm),
    m_type(MT_GROWING),
    m_nucleation_time(s->m_system_time + s->m_system_time_offset),
    m_previous_update_tag(s->m_current_time_tag)
{
#ifdef DBG_MTS
    cout << "DBG/MTS: Microtubule created.\n";
#endif

    if (initialize)
    {
        // create a single segment MT
        m_segments.create(this, tv);

        // intialize the plus tip of MT
        m_plus.initialize();

        // initialize the minus tip of MT
        m_minus.initialize();

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
    m_plus.unlink_from_trajectory();
    m_minus.unlink_from_trajectory();

    return;
}

void Microtubule::set_disappear_event()
{
    // push to disappear event: (plus.velocity - minus.velocity < 0)
    if ((m_segments.size() == 1) && (m_type == MT_SHRINKING))
    {
        m_disappear_event.push_on_queue(m_segments.first()->length(), EV_DISAPPEAR);
    }

    // flush disappear event queue
    else
    {
        m_disappear_event.clear();
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

    switch (ei->m_type)
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
            if (ei == &(m_plus.m_event))
            {
                backtrack(&m_plus);
            }
            else
            {
                backtrack(&m_minus);
            }
            break;

            // end of segment event encountered
        case EV_END_OF_SEGMENT:
            if (ei == &(m_plus.m_event))
            {
                end_of_segment(&m_plus);
            }
            else
            {
                end_of_segment(&m_minus);
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
    m_plus.m_trajectory->m_base.m_region->unregister_from_region(m_plus.m_region_tag, T_PLUS, MT_GROWING);

    // now, this is shrinking MT
    m_type = MT_SHRINKING;

    // growing MT will become shrinking MT, hence register the tip on shrinking-plus-tip list of this region
    m_plus.m_region_tag = m_plus.m_trajectory->m_base.m_region->register_on_region(&m_plus, T_PLUS, MT_SHRINKING);

    // replace plus velocity by minus velocity
    m_plus.m_velocity = m_system->m_p.m_v_min;

    //
    m_plus.m_event.reinitialize(&(m_system->m_time_queue), m_system->m_p.m_v_min);

    // import the growing MT as shrinking MT
    m_system->m_shrinking_mts.import(m_system->m_growing_mts, this);

    if (m_system->m_p.m_v_min < 0)
    {
        // velocity has changed sign, select a new event
        m_plus.advance_intersection();
    }

    // schedule next event, (if this is triggered by a collision, and vMin > 0, increase occupancy)
    m_plus.determine_event();

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
    m_plus.m_trajectory->m_base.m_region->unregister_from_region(m_plus.m_region_tag, T_PLUS, MT_SHRINKING);

    // now, this is growing MT
    m_type = MT_GROWING;

    // shrinking MT will become growing MT, hence register the tip on growing-plus-tip list of this region
    m_plus.m_region_tag = m_plus.m_trajectory->m_base.m_region->register_on_region(&m_plus, T_PLUS, MT_GROWING);

    // replace minue velocity by plus velocity
    m_plus.m_velocity = m_system->m_p.m_v_plus;

    m_plus.m_event.reinitialize(&(m_system->m_v_plus_queue), m_system->m_p.m_v_plus);

    // import the shrinking MT as growing MT
    m_system->m_growing_mts.import(m_system->m_shrinking_mts, this);

    if (m_system->m_p.m_v_min < 0)
    {
        // velocity has changed sign, need to select a new event
        m_plus.advance_intersection();
    }

    // schedule next event
    m_plus.determine_event();

    // check possibility of a disappear event
    set_disappear_event();

    return;
}

void Microtubule::sever(Segment* cut_seg, double cut_pos)
{
    // sever MT at given position
    split_segment_at_traj_pos(cut_seg->m_start + static_cast<double>(cut_seg->m_dir) * cut_pos, cut_seg);
    m_system->m_total_length_severing_count++;
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
    if (cut_seg->m_dir == ::FORWARD)
    {
        split_pos += ZERO_CUTOFF;
    }
    else
    {
        split_pos -= ZERO_CUTOFF;
    }

    // only cut if the severing position is comfortably located within the segment interval
    if ((cut_seg->m_start - split_pos) * (cut_seg->m_end - split_pos) < ZERO_CUTOFF)
    {
        split_segment_at_traj_pos(split_pos, cut_seg);
        m_system->m_total_intersection_severing_count++;
    }
    return;
}

void Microtubule::split_segment_at_traj_pos(double cut_pos, Segment* cut_seg)
{
    // used in sever and severAtCross
    Microtubule* new_mt(nullptr);
    TrajectoryVector tv;

    tv.m_trajectory = cut_seg->m_trajectory;
    tv.m_dir = cut_seg->m_dir;
    tv.m_pos = cut_pos;
    if (m_type == MT_GROWING)
    {
        new_mt = m_system->m_growing_mts.create(m_system, tv, false);
    }
    else
    {
        new_mt = m_system->m_shrinking_mts.create(m_system, tv, false);
    }

    new_mt->m_nucleation_time = m_nucleation_time;
    new_mt->m_type = m_type;

    // transfer segments
    new_mt->m_segments.create(new_mt, tv);
    new_mt->m_segments.first()->m_end = cut_seg->m_end;
    new_mt->m_segments.first()->m_end_itr = cut_seg->m_end_itr;
    new_mt->m_segments.first()->m_nucleation_time = cut_seg->m_nucleation_time;
    new_mt->m_segments.import_set(m_segments, cut_seg->next(), m_segments.last());
    Segment* temp_seg = new_mt->m_segments.first()->next();
    while (temp_seg != nullptr)
    {
        temp_seg->m_mt = new_mt;
        temp_seg = temp_seg->next();
    }

    // now update the plus tip trajectory/position
    new_mt->m_plus.m_trajectory = m_plus.m_trajectory;
    new_mt->m_plus.m_velocity = m_plus.m_velocity;
    new_mt->m_plus.m_event.reinitialize(m_plus.m_event.m_queue, m_plus.m_velocity);
    new_mt->m_plus.m_dir = m_plus.m_dir;
    new_mt->m_plus.m_next_collision = m_plus.m_next_collision;
    new_mt->m_plus.m_notification_tag = new_mt->m_plus.m_trajectory->register_for_notifications(&(new_mt->m_plus));
    new_mt->m_plus.m_region_tag = new_mt->m_plus.m_trajectory->m_base.m_region->register_on_region(&(new_mt->m_plus), T_PLUS, m_type);

    // next, update the minus end of the new MT
    new_mt->m_segments.first()->m_start = tv.m_pos;
    new_mt->m_minus.m_notification_tag = new_mt->m_minus.m_trajectory->register_for_notifications(&(new_mt->m_minus));
    new_mt->m_minus.m_region_tag = new_mt->m_minus.m_trajectory->m_base.m_region->register_on_region(&(new_mt->m_minus), T_MINUS, m_type);
    new_mt->m_minus.locate_intersection();

    new_mt->m_plus.determine_event();
    new_mt->m_minus.determine_event();
    new_mt->set_disappear_event();

    // now, finish up the base microtubule
    m_segments.last()->m_end = cut_pos;

    // always unregister from region (just to be sure)
    m_plus.m_trajectory->m_base.m_region->unregister_from_region(m_plus.m_region_tag, T_PLUS, m_type);
    if (m_type == MT_GROWING)
    {
        m_type = MT_SHRINKING;
        m_plus.m_velocity = m_system->m_p.m_v_min;
        m_plus.m_event.reinitialize(&(m_system->m_time_queue), m_system->m_p.m_v_min);
        m_system->m_shrinking_mts.import(m_system->m_growing_mts, this);
    }

    if (m_plus.m_trajectory != tv.m_trajectory)
    {
        m_plus.m_trajectory->unregister_for_notifications(m_plus.m_notification_tag);
        m_plus.m_trajectory = tv.m_trajectory;
        m_plus.m_notification_tag = m_plus.m_trajectory->register_for_notifications(&(m_plus));
    }
    //...and register on region again
    m_plus.m_region_tag = m_plus.m_trajectory->m_base.m_region->register_on_region(&m_plus, T_PLUS, MT_SHRINKING);
    m_plus.m_dir = tv.m_dir;
    m_plus.m_next_collision = new_mt->m_minus.m_next_collision;
    if (m_system->m_p.m_v_min * new_mt->m_minus.m_velocity >= 0)
    {
        // both tips move in opposite directions [this should generally be the case]
        m_plus.advance_intersection();
    }

    m_plus.determine_event();
    if (m_segments.size() == 1)
    {
        m_minus.determine_event();
    }
    set_disappear_event();

    return;
}

void Microtubule::translate_position_m_t2_segment(double& cut_pos, Segment*& cut_seg)
{
    cut_seg = m_segments.first();
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
    cut_pos = cut_seg->m_start + static_cast<double>(cut_seg->m_dir) * cut_pos;
    return;
}

void Microtubule::wall()
{
    // only for plus end

    TrajectoryVector tv;

    // set position of the plus end equal to the event position, to mitigate the effect of cumulating addition errors
    m_segments.last()->m_end = m_plus.m_next_event_pos;
    tv = m_plus.m_trajectory->next_trajectory(m_plus.m_dir);

    // next trajectory on 'forbidden zone', leads catastrophe
    if (tv.m_trajectory->m_base.m_region->m_face_tag == -1)
    {
        tv.m_trajectory->conditional_remove();
        catastrophe();

        // no wall crossing
        return;
    }

    // edge catastrophe enabled
    if (m_system->m_p.m_edge_catastrophe_enabled)
    {
        double cos_angle = (m_plus.m_dir == ::FORWARD) ? m_plus.m_trajectory->m_next_tr_cos_angle : m_plus.m_trajectory->m_prev_tr_cos_angle;

        // MT is moving non-parallel to the edege
        if ((1.0 - cos_angle) > ZERO_CUTOFF)
        {
            double p_cat(0.0);

            // for sharp edged cube
            if (m_system->m_p.m_geometry == "cubeReal")
            {
                if (abs(m_plus.m_trajectory->m_base.m_region->m_face_tag - tv.m_trajectory->m_base.m_region->m_face_tag) > 0)
                {
                    if (abs(m_plus.m_trajectory->m_base.m_region->m_face_tag - tv.m_trajectory->m_base.m_region->m_face_tag) > 1)
                    {
                        p_cat = m_system->m_p.m_p_cat_special_edge_max;
                    }

                    else
                    {
                        p_cat = m_system->m_p.m_p_cat_regular_edge_max;
                    }
                }
            }

            // other
            else
            {
                p_cat = (m_plus.m_dir == ::FORWARD) ? m_plus.m_trajectory->m_next_trp_cat : m_plus.m_trajectory->m_prev_trp_cat;
            }

            // when smooth catastrophe is enabled
            if (m_system->m_p.m_edge_catastrophe_smooth)
            {
                p_cat *= (1.0 - cos_angle);
            }

            // apply catastrophe
            if (m_system->m_random_gen.rand() <= p_cat)
            {
                tv.m_trajectory->conditional_remove();
                catastrophe();

                // no wall crossing
                return;
            }
        }
    }

    // update boundarty crossing count
    m_system->m_boundary_crossing_count++;

#ifdef DBG_ASSERT
    if (segments.last()->endItr != plus.nextCollision)
    {
        cerr << "DBG/ASSERT: nextCollision does not match with endItr upon entering wall event.\n";
    }
#endif

    // wall crossing
    m_segments.last()->m_end_itr = m_plus.m_next_collision;

    // create a new segment
    m_segments.create(this, tv);
    m_segments.last()->m_start_itr = tv.m_dir == ::FORWARD ? tv.m_trajectory->wall_begin() : tv.m_trajectory->wall_end();

    // switch trajectory
    m_plus.switch_trajectory(
    tv.m_trajectory, tv.m_dir, tv.m_dir == ::FORWARD ? tv.m_trajectory->wall_begin() : tv.m_trajectory->wall_end());

    if (m_segments.size() == 2)
    {
        // MT with two segment of zero (=0) length, consider possibilty of a disappear event
        m_minus.determine_event();
        set_disappear_event();
    }

    return;
}

void Microtubule::collision()
{
    // set position of the plus end equal to the event position, to mitigate the effect of cumulating addition errors
    m_segments.last()->m_end = m_plus.m_next_event_pos;
    CollisionType type;
    Direction dir;

    // no MT(s) passing throught the next collision point, resulting a simple crossover
    if (m_plus.m_next_collision->second.m_mirror->second.m_occupancy == 0)
    {
        type = CT_CROSSOVER;
    }

    // MT(s) passing through the next collision point
    else
    {
        double diff_angle
        = m_plus.m_trajectory->m_base.m_region->intersection_angle(m_plus.m_trajectory, m_plus.m_next_collision->second.m_other_trajectory);

        if (diff_angle > 0.5 * PI)
        {
            dir = (m_plus.m_dir == ::FORWARD) ? BACKWARD : ::FORWARD;
            diff_angle = PI - diff_angle;
        }
        else
        {
            dir = m_plus.m_dir;
        }

        // find out type of collision event to be occured
        type = m_system->collision_result(
        diff_angle, m_plus.m_next_collision->second.m_occupancy, m_plus.m_next_collision->second.m_mirror->second.m_occupancy);
    }

    switch (type)
    {
            // zippering
        case CT_ZIPPER:
            m_system->m_total_zipper_count++;
            zipper(dir);
            break;

            // crossover
        case CT_CROSSOVER:
            m_system->m_total_crossover_count++;
            crossover();
            break;

            // induced catastrophe
        case CT_INDUCED_CATASTROPHE:
            m_system->m_total_induced_catastrophe_count++;
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
    TrajectoryVector tv(m_plus.m_next_collision->second.m_mirror->first, dir, m_plus.m_next_collision->second.m_other_trajectory);

    // link the end of the last segment (containing the plus tip) to the next plus end collision
    m_segments.last()->m_end_itr = m_plus.m_next_collision;

    // create a new zero (=0) length segment on new position and direction (of course on new trajectory)
    m_segments.create(this, tv);

    // link the bottom of the last segment (containing the new minus tip) to the next plus end collision mirror
    m_segments.last()->m_start_itr = m_plus.m_next_collision->second.m_mirror;

    // now do the switching of trajectory from incoming to outgoing
    m_plus.switch_trajectory(m_plus.m_next_collision->second.m_other_trajectory, dir, m_plus.m_next_collision->second.m_mirror);

    if (m_segments.size() == 2)
    {
        // number of segment is two (=2), HERE it means two zero length segment, prepare to schedule a disappear event
        m_minus.determine_event();
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
    if (tip == &m_plus)
    {
        kill_seg = m_segments.last();
        new_seg = m_segments.last()->previous();
        new_dir = new_seg->m_dir;
    }

    // select minus tip associated segment
    else
    {
        kill_seg = m_segments.first();
        new_seg = m_segments.first()->next();
        if (new_seg->m_dir == ::FORWARD)
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
    m_segments.remove(kill_seg);

    // tip switching to a new trajectory (in a different region) after encountering wall-begin/wall-end
    if ((tip->m_next_collision == tip->m_trajectory->wall_begin()) || (tip->m_next_collision == tip->m_trajectory->wall_end()))
    {
        tip->switch_trajectory(new_seg->m_trajectory,
                              new_dir,
                              new_dir == ::FORWARD ? new_seg->m_trajectory->wall_end() : new_seg->m_trajectory->wall_begin());
        m_system->m_boundary_crossing_count--;
    }

    // tip switching to a new trajectory (in the same region)
    else
    {
        tip->switch_trajectory(tip->m_next_collision->second.m_other_trajectory, new_dir, tip->m_next_collision->second.m_mirror);
    }

    // MT length zero (=0) encountered
    if (m_segments.size() == 1)
    {
        // current tip is a plus tip, so determine the minus tip associated event
        if (tip == &(m_plus))
        {
            m_minus.determine_event();
        }

        // current tip is a minus tip, so determine the plus tip associated event
        else
        {
            m_plus.determine_event();
        }

        // schedule a disappear event
        set_disappear_event();
    }

    return;
}

void Microtubule::backtrack(MTTip* tip)
{
    // asociated with depolymerizing plus/minus end

    if (tip == &(m_plus))
    {
        m_segments.last()->m_end = m_plus.m_next_event_pos;
    }
    else
    {
        m_segments.first()->m_start = m_minus.m_next_event_pos;
    }

    // MT is back tracked from this site, resulting a decrease in occupancy
    tip->m_next_collision->second.m_occupancy--;

#ifdef CROSS_SEV
    if ((tip->m_next_collision->second.m_occupancy == 0) && (tip->m_next_collision->second.m_mirror->second.m_occupancy > 0))
    {
        m_system->remove_occupied_intersection(tip->m_next_collision->second);
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
    if ((m_plus.m_next_collision->second.m_occupancy == 0) && m_plus.m_next_collision->second.m_mirror->second.m_occupancy > 0)
    {
        m_system->add_occupied_intersection(m_plus.m_next_collision);
    }
#endif

    // skip any (active) event in the current intersection point, instead move the tip to the next intersection point
    m_plus.m_next_collision->second.m_occupancy++;
    m_plus.advance_intersection();
    m_plus.determine_event();

    return;
}

void Microtubule::harakiri()
{
    // MT disappear event encountered, remove the MT
    if (m_type == MT_GROWING)
    {
        m_system->m_growing_mts.remove(this);
    }
    else
    {
        m_system->m_shrinking_mts.remove(this);
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
    if (m_segments.size() == 0)
    {
        // MT should contain at least one segment
        valid = false;
        cerr << "Microtubule contains no segments.\n";
    }

    else if (m_segments.size() != 1)
    {
        if (m_segments.size() == 0)
        {
            // MT should contain at least one segment
            cerr << "microtubule without segments!\n";
            exit(-1);
        }

        seg = m_segments.first();
        if (seg == nullptr)
        {
            // pointer to segment list can not be null
            cerr << "null pointer received\n";
            exit(-1);
        }

        // segment (some where) at the middle of a trajectory
        if ((seg->m_end_itr != seg->m_trajectory->wall_end()) && (seg->m_end_itr != seg->m_trajectory->wall_begin()))
        {
            if (seg->m_end_itr->second.m_mirror != seg->next()->m_start_itr)
            {
                // next-start-iterator should match with the end-second-mirror-iterator
                valid = false;
                cerr << "reciprocity violated\n";
            }

            if (abs(seg->m_end - seg->m_end_itr->first) > ZERO_CUTOFF)
            {
                // segment end location does not coincide with the first intersection site
                valid = false;
                cerr << "Location of segment end does not coincide with intersection position\n";
            }

            max_deviation = max(max_deviation, abs(seg->m_end - seg->m_end_itr->first));

            if (abs(seg->next()->m_start - seg->next()->m_start_itr->first) > ZERO_CUTOFF)
            {
                // (next) segment end location does not coincide with the first intersection site
                valid = false;
                cerr << "Location of segment end does not coincide with intersection position\n";
            }
            max_deviation = max(max_deviation, abs(seg->next()->m_start - seg->next()->m_start_itr->first));
        }

        // segment at the end of a trajectory
        else
        {
            // positions should match with wall end location
            if (seg->m_end_itr == seg->m_trajectory->wall_end())
            {
                if (abs(seg->m_end - seg->m_trajectory->m_length) > ZERO_CUTOFF)
                {
                    // segment end location does not coincide with the first intersection site
                    valid = false;
                    cerr << "Segment end location does not coincide with end of trajectory\n";
                }
                max_deviation = max(max_deviation, abs(seg->m_end - seg->m_trajectory->m_length));
            }

            // positions should match with wall begin location
            else
            {
                if (abs(seg->m_end - 0) > ZERO_CUTOFF)
                {
                    // segment end location does not coincide with the first intersection site
                    valid = false;
                    cerr << "Segment end location does not coincide with end of trajectory\n";
                }
                max_deviation = max(max_deviation, abs(seg->m_end - 0));
            }

            if (seg->next()->m_start_itr == seg->next()->m_trajectory->wall_end())
            {
                if (abs(seg->next()->m_start - seg->next()->m_trajectory->m_length) > ZERO_CUTOFF)
                {
                    valid = false;
                    cerr << "Segment start location does not coincide with end of trajectory\n";
                    cerr << "start: " << seg->next()->m_start << ", length: " << seg->next()->m_trajectory->m_length
                         << ", diff=" << abs(seg->next()->m_start - seg->next()->m_trajectory->m_length) << "\n";
                }
                max_deviation = max(max_deviation, abs(seg->next()->m_start - seg->next()->m_trajectory->m_length));
            }
            else
            {
                if (abs(seg->next()->m_start - 0) > ZERO_CUTOFF)
                {
                    valid = false;
                    cerr << "Segment start location does not coincide with start of trajectory\n";
                    cerr << "start: " << seg->next()->m_start << ", length: " << seg->next()->m_trajectory->m_length
                         << ", diff=" << abs(seg->next()->m_start) << "\n";
                }
                max_deviation = max(max_deviation, abs(seg->next()->m_start - 0));
            }
        }
        seg = seg->next();
    }

    if (!valid)
    {
        cout << "Broken segment connections\n";
    }

    // check begin and end points
    if (static_cast<double>(m_plus.m_dir) * m_plus.m_velocity > 0)
    {
        if (m_plus.m_next_collision == m_plus.m_trajectory->wall_end())
        {
            // segment end location does not coincide with the first intersection site
            if (m_segments.last()->m_end - m_plus.m_trajectory->m_length > ZERO_CUTOFF)
            {
                valid = false;
            }
        }

        else
        {
            // segment end location does not coincide with the first intersection site
            if (m_segments.last()->m_end - m_plus.m_next_collision->first > ZERO_CUTOFF)
            {
                valid = false;
            }
        }
    }

    else
    {
        if (m_plus.m_next_collision == m_plus.m_trajectory->wall_begin())
        {
            // segment end location does not coincide with the first intersection site
            if (m_segments.last()->m_end < -ZERO_CUTOFF)
            {
                valid = false;
            }
        }
        else
        {
            // segment end location does not coincide with the first intersection site
            if (m_segments.last()->m_end - m_plus.m_next_collision->first < -ZERO_CUTOFF)
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
    if (static_cast<double>(m_minus.m_dir) * m_minus.m_velocity > ZERO_CUTOFF)
    {
        if (m_minus.m_next_collision == m_minus.m_trajectory->wall_end())
        {
            if (m_segments.first()->m_start - m_minus.m_trajectory->m_length > ZERO_CUTOFF)
            {
                valid = false;
            }
        }
        else
        {
            if (m_segments.first()->m_start - m_minus.m_next_collision->first > ZERO_CUTOFF)
            {
                valid = false;
            }
        }
    }
    else if (static_cast<double>(m_minus.m_dir) * m_minus.m_velocity < -ZERO_CUTOFF)
    {
        if (m_minus.m_next_collision == m_minus.m_trajectory->wall_begin())
        {
            if (m_segments.first()->m_start < -ZERO_CUTOFF)
            {
                valid = false;
            }
        }
        else
        {
            if (m_segments.first()->m_start - m_minus.m_next_collision->first < -ZERO_CUTOFF)
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
    seg = m_segments.first();
    while (seg != nullptr)
    {
        is = seg->m_start_itr;
        count = 0;
        while ((is != seg->m_trajectory->wall_end()) && (count != seg->m_trajectory->m_intersections.size()))
        {
            is++;
            count++;
        }
        if (is != seg->m_trajectory->wall_end())
        {
            valid = false;
            cerr << "Error: segment start iterator does not belong to correct trajectory\n";
            cerr << "Number of segments in MT: " << m_segments.size() << "\n";
            ;
        }
        is = seg->m_end_itr;
        count = 0;
        while ((is != seg->m_trajectory->wall_end()) && (count != seg->m_trajectory->m_intersections.size()))
        {
            is++;
            count++;
        }
        if (is != seg->m_trajectory->wall_end())
        {
            valid = false;
            cerr << "Error: segment end iterator does not belong to correct trajectory\n";
            cerr << "Number of segments in MT: " << m_segments.size() << "\n";
            ;
        }
        seg = seg->next();
    }

    // check the tip's tipmap pointer (round-trip)
    if (m_system->get_event_descriptor(m_plus.m_event.m_index) != &m_plus.m_event)
    {
        valid = false;
    }
    if (m_system->get_event_descriptor(m_minus.m_event.m_index) != &m_minus.m_event)
    {
        valid = false;
    }
    if (m_system->get_event_descriptor(m_disappear_event.m_index) != &m_disappear_event)
    {
        valid = false;
    }

    // round-trip check of the segment-trajectory pairs
    seg = m_segments.first();
    while (seg != nullptr)
    {
        if (*(seg->m_trajectory_tag) != seg)
        {
            valid = false;
        }

        seg = seg->next();
    }

    // round-trip check of tip-trajectory pairs
    if (*(m_plus.m_notification_tag) != &m_plus)
    {
        valid = false;
    }
    if (*(m_minus.m_notification_tag) != &m_minus)
    {
        valid = false;
    }

    // check orientation of segments / start-end, and tips
    if (m_segments.first()->m_dir == m_minus.m_dir)
    {
        valid = false;
    }
    if (m_segments.last()->m_dir != m_plus.m_dir)
    {
        valid = false;
    }

    seg = m_segments.first();
    while (seg != nullptr)
    {
        if (seg->m_dir == ::FORWARD)
        {
            if (seg->m_end - seg->m_start < -ZERO_CUTOFF)
            {
                cerr << "false directionality in positions. diff=" << seg->m_end - seg->m_start << "\n";
                valid = false;
            }
            is = seg->m_start_itr;
            count = 0;
            while ((is != seg->m_end_itr) && (count != seg->m_trajectory->m_intersections.size()))
            {
                is++;
                count++;
            }
            if (is != seg->m_end_itr)
            {
                cerr << "false directionality in intersection iterators\n";
                valid = false;
            }
        }
        else
        {
            if (seg->m_end - seg->m_start > ZERO_CUTOFF)
            {
                cerr << "false directionality in positions. diff=" << seg->m_end - seg->m_start << "\n";
                valid = false;
            }
            is = seg->m_start_itr;
            count = 0;
            while ((is != seg->m_end_itr) && (count != seg->m_trajectory->m_intersections.size()))
            {
                is--;
                count++;
            }
            if (is != seg->m_end_itr)
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
