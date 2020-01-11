#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <cstdint>
#include <forward_list>
#include <iterator> // for forward_iterator_tag
#include <mutex>

#include "openmc/openmp_interface.h"

namespace openmc {

//==============================================================================
//! A list node type that will be chained to implement `NeighborList`.
//==============================================================================

struct NeighborListNode
{
  explicit NeighborListNode(int32_t v)
    : value(v)
  {}

  NeighborListNode* next {nullptr};
  int32_t value;
};

//==============================================================================
//! An iterator that traverses elements of a `NeighborList`.
//==============================================================================

class NeighborListIter
{
public:
  using value_type = int32_t;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = int32_t;
  using iterator_category = std::forward_iterator_tag;

  NeighborListIter(const NeighborListNode* n, difference_type position,
    difference_type length)
    : node_(n),
    position_(position),
    length_(length)
  {}

  reference operator*()
  {return node_->value;}

  bool operator==(const NeighborListIter& other)
  {return node_ == other.node_;}

  bool operator!=(const NeighborListIter& other)
  {return node_ != other.node_;}

  NeighborListIter& operator++()
  {
    // Only read node_->next values for the first `length_-1` nodes. Subsequent
    // nodes may have invalid next pointers due to concurrent read/writes.
    if (position_ < length_ - 1) {
      node_ = node_->next;
      ++position_;
    } else {
      node_ = nullptr;
    }
    return *this;
  }

private:
  const NeighborListNode* node_;
  difference_type position_;
  difference_type length_;
};

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
  using const_iterator = NeighborListIter;

  NeighborList() = default;

  ~NeighborList()
  {
    // It is assumed this function is called in a serial section of the code.
    NeighborListNode* this_node = head_;
    while (this_node) {
      NeighborListNode* next_node = this_node->next;
      delete this_node;
      this_node = next_node;
    }
  }

  const_iterator cbegin() const
  {
    // Atomically read the list length to prevent split reads.  If the length
    // is zero, explicitly pass nullptr to the iterator.  (head_ should be a
    // nullptr, but it may be an invalid value due to a concurrent write.)
    // head_ will be valid for all lengths > 0.
    const_iterator::difference_type length;
    #pragma omp atomic read
    length = length_;
    if (length == 0) {
      return const_iterator(nullptr, 0, 0);
    } else {
      return const_iterator(head_, 0, length_);
    }
  }

  const_iterator cend() const
  {return const_iterator(nullptr, 0, 0);}

  //! Attempt to add an element.
  //
  //! If the relevant OpenMP lock is currently owned by another thread, this
  //! function will return without actually modifying the data.  It has been
  //! found that returning the transport calculation and possibly re-adding the
  //! element later is slightly faster than waiting on the lock to be released.
  void push_back(value_type new_elem)
  {
    // Update the list if a lock is available.
    std::unique_lock<OpenMPMutex> lock(write_mutex_, std::try_to_lock);
    if (lock) {
      // It is possible another thread already added this element to the list
      // while this thread was searching for a cell so make sure the given
      // element isn't a duplicate before adding it.
      if (std::find(cbegin(), cend(), new_elem) == cend()) {
        // Find the end of the list and add the the new element there.
        if (auto* tail = head_) {
          while (tail->next) {
            tail = tail->next;
          }
          tail->next = new NeighborListNode(new_elem);
        } else {
          head_ = new NeighborListNode(new_elem);
        }
        // Use a flush to ensure the appropriate node pointer is updated and
        // visible to all threads before adjusting the length.  This will
        // prevent iterators from trying to read an invalid pointer.  Update
        // the length atomically.
        #pragma omp flush
        #pragma omp atomic update
        ++length_;
      }
    }
  }

private:
  NeighborListNode* head_ {nullptr};
  OpenMPMutex write_mutex_;
  const_iterator::difference_type length_ {0};
};

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
