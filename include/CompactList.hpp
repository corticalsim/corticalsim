#ifndef COMPACTLIST_H
#define COMPACTLIST_H
#include <boost/pool/singleton_pool.hpp>
#include <new>
#include <iostream>

using namespace std;

// Properties:
//- data remains at the same memory address
//- overhead vector assures quick access by index
//- indices in vector change to keep a closed row
// CompactListItems contain their index number; this is updated when changed.
template <class T, int AllocatorGranularity = 32>
class CompactList
{
  private:

    vector<T*> m_elements;

    // Memory management: choose to place items (T) together in blocks of
    // AllocatorGranularity. IMPORTANT NOTE: Not tread-safe!
    struct PoolTag
    {
    };

    typedef boost::singleton_pool<PoolTag,
                                  sizeof(T),
                                  boost::default_user_allocator_new_delete,
                                  boost::details::pool::null_mutex,
                                  AllocatorGranularity>
    MemPool;
    T* insert_element(T* new_element);
    void remove(T*);
    void remove_all();

  public:

    CompactList():
        m_elements()
    {
        m_elements.reserve(AllocatorGranularity);
    };

    ~CompactList();

    T* element_address(int idx) { return m_elements[idx]; };

    void remove_element(int idx);

    int size() { return m_elements.size(); };

    T& operator[](int idx) { return *(m_elements[idx]); };

    T* create();
    template <class Q>
    T* create(Q);
    template <class Q, class R>
    T* create(Q, R);
    template <class Q, class R, class S>
    T* create(Q, R, S);
};

template <class T, int I>
CompactList<T, I>::~CompactList()
{
    remove_all();
    return;
}

template <class T, int I>
void CompactList<T, I>::remove_all()
// private function (called by ~)
{
    for (int j = m_elements.size() - 1; j >= 0; j--)
    {
        remove(m_elements[j]);
    }
    m_elements.pop_back();
}

template <class T, int I>
inline void CompactList<T, I>::remove(T* element)
// private function (called by removeAll (~) and RemoveElement)
// Removes a single element. No checking is performed on the pointer.
{
    element->~T();
    MemPool::free(element);
    return;
}

template <class T, int I>
inline T* CompactList<T, I>::create()
{
    T* new_element = new (MemPool::malloc()) T();
    return insert_element(new_element);
}

template <class T, int I>
template <class Q>
inline T* CompactList<T, I>::create(Q par1)
{
    T* new_element = new (MemPool::malloc()) T(par1);
    return insert_element(new_element);
}

template <class T, int I>
template <class Q, class R>
inline T* CompactList<T, I>::create(Q par1, R par2)
{
    T* new_element = new (MemPool::malloc()) T(par1, par2);
    return insert_element(new_element);
}

template <class T, int I>
template <class Q, class R, class S>
inline T* CompactList<T, I>::create(Q par1, R par2, S par3)
{
    T* new_element = new (MemPool::malloc()) T(par1, par2, par3);
    return insert_element(new_element);
}

template <class T, int I>
inline T* CompactList<T, I>::insert_element(T* new_element)
// private function!
{
    new_element->m_index = m_elements.size();
    m_elements.push_back(new_element);
    return new_element;
}

template <class T, int I>
inline void CompactList<T, I>::remove_element(int idx)
// check (only the outer!!) can eliminate?
{
    int size = m_elements.size();
    if (idx < size)
    {
        remove(m_elements[idx]);
        // if idx was not the last element, move last element to the gap
        if (idx != size - 1)
        {
            m_elements[idx] = m_elements.back();
            m_elements[idx]->m_index = idx;
        }
        m_elements.pop_back();
    }
    else
    {
        cerr << "ERROR CompactList: trying to remove non-existing element. "
                "(index = "
             << idx << ", size = " << size << ")\n";
    }
    return;
}

template <class T>
class CompactListItem
{
  public:

    CompactListItem():
        m_index(0) {};
    int m_index;

  protected:

    virtual ~CompactListItem() {};
    template <class Q, int I>
    friend class CompactList;
};

#endif
// COMPACTLIST_H
