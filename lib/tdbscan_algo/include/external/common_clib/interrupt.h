//
// Created by marcel on 6.7.21.
//

#ifndef COMMONCLIB__THREADSAFE__INTERRUPT_H
#define COMMONCLIB__THREADSAFE__INTERRUPT_H

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <string>
#include <exception>

namespace common_clib::threadsafe {

/// we need something that behaves like an exception
struct interrupt_exception : std::exception {
private:
  const std::string message_;
public:
  // Constructor accepting const char*
  interrupt_exception(const char* msg) :
  message_(msg) {}

  const char* what() const noexcept {
    return message_.c_str();
  }
};


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
  /// trigger an interrupt, which is a wake on the conditional-variable
  void trigger() noexcept;


  /// block the waiting on the callers to the semaphore
  void block() noexcept;


  /// check the toggle, is the block has already been triggered
  [[nodiscard]] bool is_triggered() const noexcept {
    return blocking_.load();
  }

  /// reset the toggle, so that the Interrupt can be reused
  void unblock() noexcept;

 private:  // state
  /// the target of this blocking Interrupt
  Tclass* referencing_;
  /// a toggle
  std::atomic<bool> blocking_;
};

/// a protoclass for something that is interruptable
class Interruptable {
public:
  /// we need something that behaves like an exception
  struct interrupt_set : interrupt_exception {
    interrupt_set(const char* msg) : interrupt_exception(msg) {};
  };

  /// we need something that behaves like an exception
  struct block_set : interrupt_exception {
    block_set(const char* msg) : interrupt_exception(msg) {};
  };
private:
  mutable std::mutex mtx_;
  std::atomic<unsigned long> external_block_set_counter_{0L};
  std::atomic<bool> internal_block_set_{false};
public:
  virtual ~Interruptable() noexcept;
protected:
  std::condition_variable cond;  //make stuff wait on this cond-var
protected:
  virtual void when_setting_block() noexcept {};
  virtual void when_resetting_block() noexcept {};
protected:
  virtual void interrupt() noexcept;
protected:
  void set_internal_block() noexcept;
  void reset_internal_block() noexcept;
private:
  void set_one_external_block() noexcept;
  void reset_one_external_block() noexcept;

  friend Interrupt<Interruptable>;
public:
  bool is_blocked() const noexcept;

  /// obtain the Interrupt handle, that when triggered throws the `interrupt` resurfacing all threads from
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
  : referencing_(&reference), blocking_(false){};

template <class Tclass>
void Interrupt<Tclass>::trigger() noexcept {
  assert(referencing_);
  if (blocking_.load())
    return;
  referencing_->interrupt();
  blocking_.store(true);
}

template <class Tclass>
void Interrupt<Tclass>::block() noexcept {
  assert(referencing_);
  if (blocking_.load())
    return;
  referencing_->set_one_external_block();
  blocking_.store(true);
}

template <class Tclass>
void Interrupt<Tclass>::unblock() noexcept {
  if (!blocking_.load())  // this interrupt was never triggered to begin with
    return;
  referencing_->reset_one_external_block();
  blocking_.store(false);
}

inline
Interruptable::~Interruptable() noexcept {
  std::lock_guard lock(mtx_);
  cond.notify_all();
};

inline
void
Interruptable::interrupt() noexcept {
  cond.notify_all();
};

inline
void Interruptable::set_internal_block() noexcept {
  std::lock_guard lock(mtx_);
  internal_block_set_.store(true);
  when_setting_block();
};

inline
void Interruptable::reset_internal_block() noexcept {
  std::lock_guard lock(mtx_);
  internal_block_set_.store(false);
  when_resetting_block();
};

inline
void Interruptable::set_one_external_block() noexcept {
  std::lock_guard lock(mtx_);
  ++external_block_set_counter_;
  when_setting_block();
};

inline
void Interruptable::reset_one_external_block() noexcept {
  std::lock_guard lock(mtx_);
  when_resetting_block();
  --external_block_set_counter_;
};

inline
bool Interruptable::is_blocked() const noexcept {
  return (external_block_set_counter_.load() + static_cast<long>(internal_block_set_.load()) != 0);
}

}  // namespace common_clib::threadsafe

#endif  // COMMONCLIB__THREADSAFE__INTERRUPT_H
