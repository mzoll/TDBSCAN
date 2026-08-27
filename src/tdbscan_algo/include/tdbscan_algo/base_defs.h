//
// Created by netsu on 07/08/2026.
//

#ifndef TDBSCAN_BASE_DEFS_H
#define TDBSCAN_BASE_DEFS_H

namespace tdbscan {
    /** make a typedef for what is the notion of Time;
    * On first principles time is a continuous monotonic increasing variable
    */
    typedef double Time_t;

    /// forward declare of Positional value Prototype
    class Position_t {};

    typedef double Distance_t;

    typedef double Timediff_t;


}


#endif //TDBSCAN_BASE_DEFS_H
