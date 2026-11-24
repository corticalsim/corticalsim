// #include <csignal>
#include "system.hpp"

using namespace std;

System* s1;

int main(int argc, char* argv[])
{
    cout << fixed;
    cout << setprecision(3);
    time_t start_time = time(nullptr);

    s1 = new System(argv[1]);

    // run the simulation (else show the mesh statistics)
    s1->run(s1->m_p.m_stop_time, s1->m_p.m_output_dir, s1->m_p.m_movie_dir);

    if (s1->m_p.m_show_output == 1)
    {
        cout << "corticalSim v" << PROGRAM_VERSION << endl;
        cout << "Total running time: " << static_cast<int>(time(nullptr)) - start_time << " seconds.\n";
        cout << "stochastic events: " << s1->m_total_s_event_count << ", deterministic events: " << s1->m_total_valid_d_event_count
             << ".\n";
#ifndef CROSS_SEV
        cout << "No severing at corners possible. If desired, #define CROSS_SEV.\n";
#endif
    }

    delete s1;
    return 0;
}
