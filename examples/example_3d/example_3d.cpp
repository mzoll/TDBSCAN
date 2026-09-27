#include "tdbscan_algo/tdbscan_algo.h"

#include <iostream>
#include <cstdlib>

#include "tdbscan_algo/common_defs.h"

#include "example_3d/blib.h"
#include "example_3d/connectors.h"




// define time as a pure positive unit

using namespace std;
using namespace tdbscan;
using namespace ex3d;

TDBScan_Algo<Blib3d> construct_algo(const double distance_lim, const double time_lim) {
  auto limcon = new LimitingConnector(distance_lim, time_lim);

	TDBScan_Algo<Blib3d>::TDBScan_ParameterSet params;

	params.multiplicity=4;
	params.multiplicityTimeWindow=20;

	return {params, limcon};
}


double rand_ord() {return rand() * 100-50.;};
Position3d rand_pos() {return Position3d(rand_ord(), rand_ord(), rand_ord());};
ScalarTime_t rand_time() {return ScalarTime_t(rand() % 10000);};

// std::set<Blib4d, Blib4d::TimeOrder> construct_blibs() {
std::set<Blib3d> construct_blibs() {
	int many_blibs = 1000;

	std::set<Blib3d> blibs;
	for (int i = 0; i < many_blibs; i++) {
		blibs.insert(Blib3d(rand_pos(), rand_time()));
	}

	return blibs;
}

using namespace std;
using namespace tdbscan;

int main(int argc, char **argv) {
	auto my_algo = construct_algo({4.},{4});

	auto blibs = construct_blibs();

	auto result = my_algo.Process(blibs);

	return 0;
}
