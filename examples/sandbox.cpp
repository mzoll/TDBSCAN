//
// Created by mzoll on 31/08/2026.
//

#include <set>
#include <iostream>
#include <list>

#include "tdbscan_algo/common_defs.h"

// using namespace tdbscan_algo;
using namespace std;

class AbsBase_t {
  virtual int getData() const = 0;
};

class Base_t {
public:
  const int data_{42};
  int getData() const {return data_;};
};

class Derived_t : public Base_t {
public:
  const int more_data{0};
};



template <class t_type>
class Functional_Base {
  const int prop_;
public:
  Functional_Base(int some_int) : prop_(some_int) {};
  virtual ~Functional_Base() = default;
  virtual bool eval(t_type a) const = 0;
};

//-----------------------


class Functional_A__ : public Functional_Base<Base_t> {
public:
  Functional_A__(int some_int) : Functional_Base(some_int) {};
  [[nodiscard]] bool eval(const Base_t a) const override { return a.data_ == 42;};
};

class Functional_A : public Functional_A__, public Functional_Base<Derived_t> {
public:
  Functional_A(int some_int) : Functional_A__(some_int), Functional_Base<Derived_t>(some_int) {};
  bool eval(const Derived_t a) const {return Functional_A__::eval(a);};
};

class Functional_B__ : public Functional_Base<Base_t> {
public:
  Functional_B__(int some_int) : Functional_Base(some_int) {};
  [[nodiscard]] bool eval(const Base_t a) const override { return a.data_ % 2 == 0;};
};

class Functional_B : public Functional_B__, public Functional_Base<Derived_t> {
public:
  Functional_B(int some_int) : Functional_B__(some_int), Functional_Base<Derived_t>(some_int) {};
  bool eval(const Derived_t a) const {return Functional_B__::eval(a);};
};




//---------------------------------------------

// template <>
// class Functional_Base<Base_t> {
// public:
//   bool eval(const Base_t a) const { return a.data_ == 42;};
// };
//
// template <>
// class Functional_Base<Derived_t> : public Functional_Base<Base_t> {};
//
// using Functional = Functional_Base<Derived_t>;
//
//
//



template <class t_type>
class Host {
public:
  const Functional_Base<t_type>* fkt_ptr_;
  explicit Host(const Functional_Base<t_type>* fkt_ptr) : fkt_ptr_(fkt_ptr) {};

  t_type eval(const t_type a) const {
    if (! fkt_ptr_->eval(a))
      throw 1;
    return a;
  };
};


int main() {
  auto fkt = new Functional_A(1);
  Derived_t d;

  auto h = Host<Derived_t>(fkt);

  return h.eval(d).more_data;
}
