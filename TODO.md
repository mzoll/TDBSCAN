Features:
-----
* Discrete Variable Ordinate `IntPos1d` / `DiscreteBlib1d`
* Copy Stopwatch class from common_clib:
  * Add stopwatch in code and log_INFO/DEBUG(?!) some statistics at the end of each algo-run.
  * Mention performance in README.
* write an algo class that could function as a continuously running service:
  * (need to learn about that!) what is a meaningful input and output interface.
* write and algo class that especially deals with discrete Ordinate Blib classes, such as segmented detectors.
* review if some of the legacy code-classes need to be ported or reused: for example
  * Geometry Hash maps
  * configuration / serialization
  * Punsh card maps. 
* fix logging: from `std:cout` to `boost::logging` or similar. Maybe conditionally build against boost.
* pipelining:
  * setup pre-commit hooks for git.
  * setup github actions.
    * clang tidy
    * clang format
    * ...

Fixes:
---
* Rethink how the earl merge process are done, and optimize:
  * Withhold adding current blib to active and newly established clusters before merge process.
* Reevaluate the late merge process, if it fulfills its purpose.
* Reevaluate the AlgoParameter class.
* Reevaluate how the ``sync_time`` is evaluated and can be meaningfully used: See also Feature 'continuously running service'. 
* see if unit-test in core-library are complete.
* Create a meaningful scenario for `example3d`, mostly focusing on the signal source.


README / Tutorial /Docs
------

* Finish README:
  * Describe ``Connectors``:
    * Guiding principles on high level.
    * Construction in code
  * Describe Parameters to the algorithm: THATS GONNA BE HUGE TASK
* Break up README into parts of a structured documentation.
* Demonstrate a Blib that is using boost/datetime or a unix/iso-datetime for the time ordinate and how seemslessly it integrates
* More documentation to classes and functions, especially in code.


Examples
--------

* Locate / ask for a raw IceCube 86 detector sample for demonstration purposes,
* Possibly other detectors of similar build.
* Locate other useful data samples for demonstration purposes.
* 