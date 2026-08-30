//
// Created by netsu on 07/08/2026.
//

#ifndef TDBSCAN_BASE_DEFS_H
#define TDBSCAN_BASE_DEFS_H

namespace tdbscan {
    // typedef double Timediff_t;
    // typedef double Time_t;

    /** Prototype for a time-like object
    * On first principles time is a continuous monotonic increasing variable,
    * that adheres a strict order principle
    */
    class Time_t {
    public:
        /// this needs to be defined by the Final Class!
        typedef void* Timediff_t;
    public:
        bool operator<(const Time_t &rhs) const;
        bool operator==(const Time_t &rhs) const;
        Timediff_t operator-(const Time_t &rhs) const;

        static void* min();  //needs to be implemented by subclass

        static void* max();  //needs to be implemented by subclass
    };


    /// forward declare of Positional value Prototype
    class Ordinate_t {
    public:
        /// this needs to be defined by the Final Class!
        typedef void* Distance_t;
    public:
        [[nodiscard]]
        bool operator<(const Ordinate_t& rhs) const;

        [[nodiscard]]
        bool operator==(const Ordinate_t& rhs) const;

        [[nodiscard]]
        Distance_t magnitude() const;
        [[nodiscard]]
        Distance_t distance( const Ordinate_t& other) const;
    };
}


#endif //TDBSCAN_BASE_DEFS_H
