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
  EXPECT_NO_THROW(swatch.time());
  EXPECT_NO_THROW(swatch.stop());

  EXPECT_NO_THROW(swatch.start());
  EXPECT_NO_THROW(swatch.restart());
  EXPECT_NO_THROW(swatch.time());
  EXPECT_NO_THROW(swatch.stop());
}

TEST(stopwatch, reports) {
  Stopwatch swatch("MyReportWatch", Stopwatch<>::policy::start);
  EXPECT_NO_THROW(swatch.lap_report("First Lap"));
  EXPECT_NO_THROW(swatch.lap_report("Second Lap"));
  EXPECT_NO_THROW(swatch.time_report("Time Report"));
  EXPECT_NO_THROW(swatch.stop_report("Stoppoint"));
}
