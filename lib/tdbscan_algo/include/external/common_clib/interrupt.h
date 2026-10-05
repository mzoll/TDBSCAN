//
// Created by marcel on 6.7.21.
//

#ifndef COMMONCLIB__THREADSAFE__INTERRUPT_H
#define COMMONCLIB__THREADSAFE__INTERRUPT_H

#include <atomic>
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

}  // namespace common_clib::threadsafe

#endif  // COMMONCLIB__THREADSAFE__INTERRUPT_H
