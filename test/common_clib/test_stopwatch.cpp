//
// Created by marcel on 31.3.21.
//

#include <gtest/gtest.h>

#include "../../lib/tdbscan_algo/include/external/common_clib/stopwatch.h"
//#include "external/common_clib/stopwatch.h"

using namespace std;
using namespace common_clib;


TEST(stopwatch, functionality) {
  Stopwatch swatch("MyStopwatch", Stopwatch<>::policy::defer);
  EXPECT_NO_THROW(swatch.start());
  EXPECT_NO_THROW(swatch.lap());
  EXPECT_NO_THROW(swatch.lap());
  EXPECT_NO_THROW(auto _ = swatch.time());
  EXPECT_NO_THROW(swatch.stop());

  EXPECT_NO_THROW(swatch.reset());

  EXPECT_NO_THROW(swatch.start());
  EXPECT_NO_THROW(auto _ = swatch.time());
  EXPECT_NO_THROW(swatch.pause());
  EXPECT_NO_THROW(auto _ = swatch.time());
  EXPECT_NO_THROW(swatch.restart());
  EXPECT_NO_THROW(auto _ = swatch.time());
  EXPECT_NO_THROW(swatch.stop());

  EXPECT_ANY_THROW(swatch.stop());  // twice stop
  EXPECT_ANY_THROW(swatch.start());  // start after stop without calling reset()
}

TEST(stopwatch, reports) {
  Stopwatch swatch("MyReportWatch", Stopwatch<>::policy::start);
  EXPECT_NO_THROW(swatch.lap_report_elapsed("First Lap"));
  EXPECT_NO_THROW(swatch.lap_report_elapsed("Second Lap"));
  EXPECT_NO_THROW(swatch.report_elapsed("Time Report"));
  EXPECT_NO_THROW(swatch.stop_report_elapsed("Stoppoint"));
  EXPECT_NO_THROW(swatch.report_elapsed("Stoppoint"));
}
