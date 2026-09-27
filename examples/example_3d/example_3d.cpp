//
// Created by mzoll on 01/08/2026.
//

#include <cstdlib>

#include "tdbscan_algo/common_defs.h"

#include "blib.h"
#include "connectors.h"

#include "tdbscan_algo/tdbscan_algo.h"

using namespace std;
using namespace tdbscan;
using namespace ex3d;

TDBScan_Algo<Blib3dWithTrace> construct_algo(const double distance_lim, const double time_lim) {
  auto limcon = new LimitingConnector(distance_lim, time_lim);

	TDBScan_Algo<Blib3dWithTrace>::TDBScan_ParameterSet params;

	params.multiplicity=4;
	params.multiplicityTimeWindow=2.;
  params.emergenceTimeWindow=2.;//internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.earlyMergeOverlapRatio= 1.; //ZeroOne: its a ratio
  params.lateMergeOverlapRatio= 1.; //ZeroOne: its a ratio

	return TDBScan_Algo<Blib3dWithTrace>(params, limcon);
}


double rand_ord() {return rand() * 100-50.;};
Position3d rand_pos() {return Position3d(rand_ord(), rand_ord(), rand_ord());};
ScalarTime_t rand_time() {return ScalarTime_t(rand() % 10000);};

// std::set<Blib4d, Blib4d::TimeOrder> construct_blibs() {
std::set<Blib3dWithTrace> construct_blibs() {
	int many_blibs = 1000;

	std::set<Blib3dWithTrace> blibs;
	for (int i = 0; i < many_blibs; i++) {
		blibs.insert(Blib3dWithTrace(rand_pos(), rand_time(), Blib3dWithTrace::NOISE));
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
