//
// Created by marcel on 6.7.21.
//

#ifndef COMMONCLIB__THREADSAFE__INTERRUPT_H
#define COMMONCLIB__THREADSAFE__INTERRUPT_H

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <exception>

namespace common_clib::threadsafe {

/// we need something that behaves like an exception
struct interrupt_exception : std::exception {};

/**
 * A materialized handle that can be copied and everything
 * @tparam Tclass this block works for that class
 *
 * @example
 *
 */
template <class Tclass>
class Interrupt {
 public:  // ctor
  /// constructor
  explicit Interrupt(Tclass& reference);

 public:  // methods
  /// block the waiting on the callers to the semaphore
  void trigger() noexcept;

  /// check the toggle, is the block has already been triggered
  [[nodiscard]] bool is_triggered() const noexcept {
    return triggered_.load();
  }

  /// reset the toggle, so that the Interrupt can be reused
  void reset() noexcept;

 private:  // state
  /// the target of this blocking Interrupt
  Tclass* referencing_;
  /// a toggle
  std::atomic<bool> triggered_;
};

/// a protoclass for something that is interruptable
class Interruptable {
  mutable std::mutex mtx_;
  std::atomic<long> external_interrupt_set_counter_{0L};
  std::atomic<bool> internal_interrupt_set_{false};
public:
  virtual ~Interruptable() noexcept;
  std::condition_variable cond;  //make stuff wait on this cond-var
protected:
  virtual void when_setting_interrupt() noexcept {};
  virtual void when_resetting_interrupt() noexcept {};
protected:
  void set_internal_interrupt() noexcept;
  void reset_internal_interrupt() noexcept;
private:
  void set_one_external_interrupt() noexcept;
  void reset_one_external_interrupt() noexcept;

  friend Interrupt<Interruptable>;
public:
  bool is_interrupted() const noexcept;

  /// obtain the Interrupt handle, that when triggered throws the `interrupt_exception` resurfacing all threads from
  /// any waiting.
  Interrupt<Interruptable> get_interrupt() noexcept {
    return Interrupt(*this);
  };
};

}  // namespace common_clib::threadsafe

// ========================
// ===== DEFINITIONS =====
// ========================

namespace common_clib::threadsafe {

template <class Tclass>
Interrupt<Tclass>::Interrupt(Tclass& reference)
  : referencing_(&reference), triggered_(false){};

template <class Tclass>
void Interrupt<Tclass>::reset() noexcept {
  if (!triggered_.load())  // this interrupt was never triggered to begin with
    return;
  referencing_->reset_one_external_interrupt();
  triggered_.store(false);
}

template <class Tclass>
void Interrupt<Tclass>::trigger() noexcept {
  assert(referencing_);
  if (triggered_.load()) return;
  referencing_->set_one_external_interrupt();
  triggered_.store(true);
}

inline Interruptable::~Interruptable() noexcept {
  std::lock_guard lock(mtx_);
  cond.notify_all();
};

inline
void Interruptable::set_internal_interrupt() noexcept {
  std::lock_guard lock(mtx_);
  internal_interrupt_set_.store(true);
  when_setting_interrupt();
  cond.notify_all();
};

inline
void Interruptable::reset_internal_interrupt() noexcept {
  std::lock_guard lock(mtx_);
  internal_interrupt_set_.store(false);
  when_resetting_interrupt();
};

inline
void Interruptable::set_one_external_interrupt() noexcept {
  std::lock_guard lock(mtx_);
  ++external_interrupt_set_counter_;
  when_setting_interrupt();
  cond.notify_all();
};

inline
void Interruptable::reset_one_external_interrupt() noexcept {
  std::lock_guard lock(mtx_);
  when_resetting_interrupt();
  --external_interrupt_set_counter_;
};

inline
bool Interruptable::is_interrupted() const noexcept {
  return (external_interrupt_set_counter_.load() + static_cast<long>(internal_interrupt_set_.load()) != 0);
}

}  // namespace common_clib::threadsafe

#endif  // COMMONCLIB__THREADSAFE__INTERRUPT_H
