#include <csignal>
#include "corticalSimReal.h"

using namespace std;

System* s1;
System* s2;


int main(int argc, char* argv[])
{
	cout << fixed;
	cout << setprecision(3);
   	time_t startTime = time(0);
	
	s1 = new System(argv[1]);

	// run the simulation (else show the mesh statistics)
	s1->run(s1->p.stopTime,s1->p.outputDir,s1->p.movieDir);

	if(s1->p.showOutput == 1)
	{
	    	cout << "corticalSim v" << PROGRAM_VERSION << endl;
	    	cout << "Total running time: " << static_cast<int>(time(0)) - startTime << " seconds.\n";
	    	cout << "stochastic events: " << s1->totalSEventCount << ", deterministic events: " << s1->totalValidDEventCount << ".\n";
	    	#ifndef CROSS_SEV
	    	cout << "No severing at corners possible. If desired, #define CROSS_SEV.\n";
	    	#endif
	}

	delete s1;
	return 0;
}



