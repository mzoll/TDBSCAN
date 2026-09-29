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
  int getData() const { return data_; };
};

class Derived_t : public Base_t {
public:
  const int more_data{0};
};


template<class t_type>
class Functional_Base {
protected:
  const bool hidden_prop_;

public:
  Functional_Base(bool invert = false) : hidden_prop_(invert) {};

  virtual ~Functional_Base() = default;

  [[nodiscard]] virtual bool eval(t_type a) const = 0;
};


class Functional_A_ : public Functional_Base<Base_t> {
public:
  Functional_A_() : Functional_Base(false) {};
  [[nodiscard]] bool eval(const Base_t a) const override { return a.data_ == 42; };
};

class Functional_A : public Functional_A_, public Functional_Base<Derived_t> {
public:
  Functional_A() : Functional_A_(), Functional_Base<Derived_t>(Functional_A_::hidden_prop_) {};
  bool eval(const Derived_t a) const { return Functional_A_::eval(static_cast<Base_t>(a)); };
};

class Functional_B_ : public Functional_Base<Base_t> {
public:
  Functional_B_() : Functional_Base(false) {};
  [[nodiscard]] bool eval(const Base_t a) const override { return hidden_prop_ ^ a.data_ % 2 == 0; };
};

class Functional_B : public Functional_B_, public Functional_Base<Derived_t> {
public:
  Functional_B() : Functional_B_(), Functional_Base<Derived_t>(Functional_B_::hidden_prop_) {};
  bool eval(const Derived_t a) const { return Functional_B_::eval(static_cast<Base_t>(a)); };
};

//-----------------------------------------------------------------

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


template<class t_type>
class Host {
public:
  const Functional_Base<t_type> *fkt_ptr_;
  explicit Host(const Functional_Base<t_type> *fkt_ptr) : fkt_ptr_(fkt_ptr) {};

  t_type eval(const t_type a) const {
    if (!fkt_ptr_->eval(a))
      throw 1;
    return a;
  };
};


int main() {
  auto fkt = new Functional_A();
  Derived_t d;

  auto h = Host<Derived_t>(fkt);

  //return h.eval(d).more_data;

  return true ^ (1 == 1);
}
