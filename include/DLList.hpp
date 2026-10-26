#ifndef DLLIST_HPP
#define DLLIST_HPP
#include <iostream>
#include <cstddef>
#include <iterator>
#include <list>
#include <array>
#include <memory>
#include <memory_resource>
#include <type_traits>
#include <vector>

template <class T>
using DLIter = typename std::pmr::list<T>::iterator;

template <class T>
using DLConstIter = typename std::pmr::list<T>::const_iterator;

// An auxiliary class that makes any derived class
// become self-aware of its location in a list.
template <class T>
class DLItem
{
    // This is a thin shim over an iterator that makes any instance
    // of a class deriving from this one *self-aware*, meaning that
    // it can erase itself from the list in O(1) time without having
    // to iterate over the list.
    // It also provides some convenience methods, such as prev() and next()
    // to extract the preceding and following elements from the list.

    DLIter<T> self;

  public:

    DLIter<T> prev()
    {
        // Return the previous element in the list.
        return std::prev(self);
    }

    DLIter<T> next()
    {
        // Return the next element in the list.
        return std::next(self);
    }

    DLIter<T> get_self()
    {
        // Return an iterator.
        return self;
    }

    void set_self(const DLIter<T>& s)
    {
        // Store a new iterator.
        self = s;
    }
};

template <class T, size_t granularity = 256>
class DLList
{
    // Restrict this class to classes derived from DLItem<T>.
    static_assert(std::is_base_of_v<DLItem<T>, T>, "This class can be used only with classes that derive from DLItem.");

    std::array<T, granularity> buffer;
    std::pmr::monotonic_buffer_resource mbr{ buffer.data(), buffer.size() };
    std::pmr::synchronized_pool_resource pool{ &mbr };
    std::pmr::list<T> list_{ &pool };

  public:

    ~DLList() { removeAll(); }

    inline size_t size() { return list_.size(); }

    inline DLIter<T> begin() { return list_.begin(); }

    inline DLIter<T> cbegin() { return list_.cbegin(); }

    inline DLIter<T> end() { return list_.end(); }

    inline DLIter<T> cend() { return list_.cend(); }

    inline DLIter<T> next(DLIter<T> elem) { return std::next(elem); }

    inline DLIter<T> next(T* elem) { return std::next(elem->get_self()); }

    inline DLIter<T> prev(DLIter<T> elem) { return std::prev(elem); }

    inline DLIter<T> prev(T* elem) { return std::prev(elem->get_self()); }

    inline T* first() { return list_.empty() ? nullptr : std::addressof(*begin()); }

    inline T* last() { return list_.empty() ? nullptr : std::addressof(*prev(end())); }

    template <class... Params>
    T* create(Params&&... params)
    {
        // Pass any parameters that T's constructors accept.
        const auto it{ list_.emplace(list_.end(), std::forward<Params>(params)...) };
        it->set_self(it);
        return std::addressof(*it);
    }

    void remove(T* element)
    {
        // Removes a single element.
        // A check is performed to see if there are any elements in the list.
        // Returns the end iterator if the removed element was the last one,
        // otherwise it returns the element *after* the one that was removed.
        if (element == nullptr)
        {
            return;
        }

        DLIter<T> it{ element->get_self() };
        if (it != end())
        {
            list_.erase(it);
        }
    }

    void removeAll()
    {
        // TODO: Work out if it's necessary to clean up the `self` attributes.
        list_.clear();
        pool.release();
    }

    T* import(DLList<T>& other, T* element)
    {
        // Import an element from another list at the end of this list.
        if (element == nullptr)
        {
            return element;
        }
        return import(other, element, list_.end());
    }

    T* import(DLList<T>& other, T* element, DLIter<T> pos)
    {
        if (element == nullptr)
        {
            return element;
        }

        // Import an element from another list at the specified position.
        // It must be copied because we are using custom memory pools,
        // so the allocators don't compare as equal.
        // Returns a pointer to the inserted element.
        const auto it{ list_.insert(pos, *element) };
        other.list_.erase(element->get_self());
        it->set_self(it);

        return std::addressof(*it);
    }

    void importSet(DLList<T>& other, T* first, T* last)
    {
        // Import a sequence of elements from another list at the end of this list.
        return importSet(other, first, last, list_.end());
    }

    void importSet(DLList<T>& other, T* first, T* last, DLIter<T> pos)
    {
        // Import a sequence of elements from another list at the specified position.
        const auto _first{ first->get_self() };
        const auto _last{ last->next() };
        size_t dist{ std::distance(_first, _last) };
        
        auto it{ list_.insert(pos, _first, _last) };
        other.list_.erase(_first, _last);

        // Update the selfies :).
        for (size_t i{ 0 }; i < dist; ++i)
        {
            it->set_self(it);
            it = it->next();
        }
    }
};

#endif // DLLIST_HPP
