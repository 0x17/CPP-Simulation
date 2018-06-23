#include "Helpers.h"
#include "Runner.h"
#include <thread>

using namespace std;

int main(int argc, const char **argv) {
	Runner::commandLine(Helpers::extractArguments(argc, argv));
	//Runner::benchmark("Instances");
	//Runner::runOptimizers();
	//Runner::alphaPsiVariations();
	return 0;
}


