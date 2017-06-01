/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma warning(disable:981)
//#pragma warning(disable:177)
#pragma warning(disable:522)
#pragma warning(disable:383)

#include <csignal>
#include "MersenneTwister.h"
#include "corticalSim.h"

using namespace std;

System* s1;
System* s2;

//Simon: disabled for use in VC++ 2010
//#pragma warning(disable:1418)
//#pragma warning(disable:1419) // let us define external functions in the primary source file

void exitNow(int sig)
{
	if (s1)
		s1->emergencyBreak();
	if (s2)
		s2->emergencyBreak();
	exit(sig);
}

void leave2(int sig) {
	cerr << "SIGINT received. Exiting.\n";
	exitNow(sig);
}
void leave3(int sig) {
	cerr << "SIGQUIT received. Exiting.\n";
	exitNow(sig);
}
void leave11(int sig) {
	cerr << "SIGSEGV received. Exiting.\n";
	exitNow(sig);
}
void leave15(int sig) {
	cerr << "SIGTERM received. Exiting.\n";
	exitNow(sig);
}



int main(int argc, char* argv[])
{

#ifdef _WIN32
#define SIGQUIT 3
#endif
	signal(SIGINT, leave2);
	signal(SIGQUIT, leave3);
	signal(SIGSEGV, leave11);
	signal(SIGTERM, leave15);

	cout << "corticalSim v" << PROGRAM_VERSION << "\n\n";
	
    time_t startTime = time(0);

	s1 = new System(argv[1]);

	s1->run(s1->p.stopTime);

	cout << "\nTotal running time: " << static_cast<int>(time(0)) - startTime << " seconds.\n";
	cout << "stochastic events: " << s1->totalSEventCount << ", deterministic events: " << s1->totalValidDEventCount << ".\n";


	delete s1;
	return 0;
}



