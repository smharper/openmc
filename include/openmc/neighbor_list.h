#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <cstdint>
#include <forward_list>
#include <mutex>

#include "openmc/openmp_interface.h"

namespace openmc {

// Forward declare the neighbor list iterator type.
template <typename T_value, typename T_prefix_iter, typename T_suffix_iter>
class NeighborListIter;

//==============================================================================
//! A threadsafe, dynamic container for listing neighboring cells.
//
//! This container is a reduced interface for a linked list with an added OpenMP
//! lock for write operations.  It allows for threadsafe dynamic growth; any
//! number of threads can safely read data without locks or reference counting.
//==============================================================================

class NeighborList
{
public:
  using value_type = int32_t;
  using prefix_iter = std::vector<value_type>::const_iterator;
  using suffix_iter = std::forward_list<value_type>::const_iterator;
  using const_iterator = NeighborListIter<value_type, prefix_iter, suffix_iter>;

  // Attempt to add an element.
  //
  // If the relevant OpenMP lock is currently owned by another thread, this
  // function will return without actually modifying the data.  It has been
  // found that returning the transport calculation and possibly re-adding the
  // element later is slightly faster than waiting on the lock to be released.
  void push_back(int new_elem)
  {
    // Try to acquire the lock.
    std::unique_lock<OpenMPMutex> lock(mutex_, std::try_to_lock);
    if (lock) {
      // It is possible another thread already added this element to the list
      // while this thread was searching for a cell so make sure the given
      // element isn't a duplicate before adding it.
      if (std::find(suffix_.cbegin(), suffix_.cend(), new_elem)
          == suffix_.cend()) {
        // Find the end of the list and add the the new element there.
        if (!suffix_.empty()) {
          auto it1 = suffix_.cbegin();
          auto it2 = ++suffix_.cbegin();
          while (it2 != suffix_.cend()) it1 = it2++;
          suffix_.insert_after(it1, new_elem);
        } else {
          suffix_.push_front(new_elem);
        }
      }
    }
  }

  // Move data from the linked-list suffix to the consecutive vector prefix.
  //
  // The consecutive data slightly improves runtime (likely due to cache
  // locality).  Note that this function is not guaranteed threadsafe---the
  // caller is responsible for making sure only one thread at a time calls this
  // function.
  void make_consecutive()
  {
    while (!suffix_.empty()) {
      prefix_.push_back(suffix_.front());
      suffix_.pop_front();
    }
  }

  const_iterator cbegin() const;
  const_iterator cend() const;


private:
  std::vector<value_type> prefix_;
  std::forward_list<value_type> suffix_;
  OpenMPMutex mutex_;

  friend class NeighborListIter<value_type, prefix_iter, suffix_iter>;
};

//==============================================================================

template <typename T_value, typename T_prefix_iter, typename T_suffix_iter>
class NeighborListIter
{
public:
  // Construct from a prefix iterator.
  NeighborListIter(const NeighborList* nl, T_prefix_iter it)
  {
    // If we were given an iterator to the end of the prefix, immediately switch
    // over to suffix mode.
    base_ = nl;
    if (it != base_->prefix_.cend()) {
      in_prefix_ = true;
      prefix_iter_ = it;
    } else {
      in_prefix_ = false;
      suffix_iter_ = base_->suffix_.cbegin();
    }
  }

  // Construct from a suffix iterator.
  NeighborListIter(const NeighborList* nl, T_suffix_iter it)
  {
    in_prefix_ = false;
    suffix_iter_ = it;
  }

  T_value operator*()
  {
    if (in_prefix_) {
      return *prefix_iter_;
    } else {
      return *suffix_iter_;
    }
  }

  bool operator==(const NeighborListIter& other)
  {
    if (in_prefix_ != other.in_prefix_) return false;

    if (in_prefix_) {
      return prefix_iter_ == other.prefix_iter_;
    } else {
      return suffix_iter_ == other.suffix_iter_;
    }
  }

  bool operator!=(const NeighborListIter& other)
  {return !(*this == other);}

  NeighborListIter& operator++()
  {
    if (in_prefix_) {
      // We are in the prefix so increment the prefix iterator.
      ++prefix_iter_;

      // If we've reached the end of the prefix, switch to the suffix iterator.
      if (prefix_iter_ == base_->prefix_.cend()) {
        in_prefix_ = false;
        suffix_iter_ = base_->suffix_.cbegin();
      }

    } else {
      // We are in the suffix so increment the suffix iterator.
      ++suffix_iter_;
    }
    return *this;
  }

private:
  const NeighborList* base_;

  // This type essentially wraps the implementation for two different external
  // iterators.  A union is used to contain one iterator or the other, and the
  // in_prefix_ flag indicates which version of that union is valid.
  union {
    T_prefix_iter prefix_iter_;
    T_suffix_iter suffix_iter_;
  };
  bool in_prefix_;
};

//==============================================================================

inline NeighborList::const_iterator
NeighborList::cbegin() const
{
  return const_iterator(this, prefix_.cbegin());
}

inline NeighborList::const_iterator
NeighborList::cend() const
{
  return const_iterator(this, suffix_.cend());
} 

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
