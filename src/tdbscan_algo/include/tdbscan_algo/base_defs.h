//
// Created by netsu on 07/08/2026.
//

#ifndef TDBSCAN_BASE_DEFS_H
#define TDBSCAN_BASE_DEFS_H

namespace tdbscan {
    typedef double Timediff_t;
    typedef double Time_t;

    // /** make a typedef for what is the notion of Time;
    // * On first principles time is a continuous monotonic increasing variable
    // */
    // class Time_t {
    // public:
    //     virtual
    //     bool operator<(const Time_t &rhs) const = 0;
    //     virtual
    //     bool operator==(const Time_t &rhs) const = 0;
    //     virtual
    //     Timediff_t operator-(const Time_t &rhs) const = 0;
    //
    //     static double min();  //needs to be implemented by subclass
    //
    //     static double max();  //needs to be implemented by subclass
    // };


    /// forward declare of Positional value Prototype

    typedef double Distance_t;

    class Position_t {
    public:
        virtual ~Position_t() = default;

        [[nodiscard]]
        bool operator<(const Position_t& other) const;

        [[nodiscard]]
        bool operator==(const Position_t& other) const;
    };

}


#endif //TDBSCAN_BASE_DEFS_H
