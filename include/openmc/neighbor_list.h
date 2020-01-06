#ifndef OPENMC_NEIGHBOR_LIST_H
#define OPENMC_NEIGHBOR_LIST_H

#include <algorithm>
#include <atomic>
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

  std::atomic<NeighborListNode*> next {nullptr};
  int32_t value;
};

//==============================================================================
//! An iterator that traverses elements of a `NeighborList`.
//==============================================================================

class NeighborListIter
{
public:
  using value_type = int32_t;
  using reference = const int32_t&;
  using pointer = const int32_t*;
  using difference_type = ptrdiff_t;
  using iterator_category = std::forward_iterator_tag;

  NeighborListIter(const NeighborListNode* n)
    : node_(n)
  {}

  reference operator*()
  {return node_->value;}

  bool operator==(const NeighborListIter& other)
  {return node_ == other.node_;}

  bool operator!=(const NeighborListIter& other)
  {return node_ != other.node_;}

  NeighborListIter& operator++()
  {
    node_ = node_->next.load(std::memory_order_relaxed);
    return *this;
  }

private:
  const NeighborListNode* node_;
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
    NeighborListNode* this_node = head_;
    while (this_node) {
      NeighborListNode* next_node = this_node->next;
      delete this_node;
      this_node = next_node;
    }
  }

  const_iterator cbegin() const
  {return const_iterator(head_.load(std::memory_order_relaxed));}

  const_iterator cend() const
  {return const_iterator(nullptr);}

  //! Attempt to add an element.
  //
  //! If the relevant OpenMP lock is currently owned by another thread, this
  //! function will return without actually modifying the data.  It has been
  //! found that returning the transport calculation and possibly re-adding the
  //! element later is slightly faster than waiting on the lock to be released.
  void push_back(value_type new_elem)
  {
    // Try to acquire the lock.
    std::unique_lock<OpenMPMutex> lock(write_mutex_, std::try_to_lock);
    if (lock) {
      // It is possible another thread already added this element to the list
      // while this thread was searching for a cell so make sure the given
      // element isn't a duplicate before adding it.
      if (std::find(cbegin(), cend(), new_elem) == cend()) {
        // Find the end of the list and add the the new element there.
        auto* head = head_.load(std::memory_order_relaxed);
        if (head) {
          auto* tail = head;
          auto* next = tail->next.load(std::memory_order_relaxed);
          while (next) {
            tail = next;
            next = tail->next.load(std::memory_order_relaxed);
          }
          auto* new_tail = new NeighborListNode(new_elem);
          tail->next.store(new_tail, std::memory_order_relaxed);
        } else {
          head_.store(new NeighborListNode(new_elem),
                      std::memory_order_relaxed);
        }
      }
    }
  }

private:
  std::atomic<NeighborListNode*> head_ {nullptr};
  OpenMPMutex write_mutex_;
};

} // namespace openmc
#endif // OPENMC_NEIGHBOR_LIST_H
