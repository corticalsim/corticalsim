#ifndef COMPACTLIST_H_
#define COMPACTLIST_H_
#include <boost/pool/singleton_pool.hpp>
#include <new>

using namespace std;

//Properties:
//- data remains at the same memory address
//- overhead vector assures quick access by index
//- indices in vector change to keep a closed row
//CompactListItems contain their index number; this is updated when changed.
template<class T, int AllocatorGranularity = 32>
class CompactList
{
private:
    vector < T* > elements;

    //Memory management: choose to place items (T) together in blocks of AllocatorGranularity.
    //IMPORTANT NOTE: Not tread-safe!
    struct PoolTag {};
    typedef boost::singleton_pool<PoolTag, sizeof(T), boost::default_user_allocator_new_delete, boost::details::pool::null_mutex, AllocatorGranularity> MemPool;
    T* InsertElement(T* newElement);
    void remove(T*);
    void removeAll();

	public:
    CompactList() : elements()
    {
        elements.reserve(AllocatorGranularity);
    } ;
    ~CompactList() ;
    T* ElementAddress(int idx)
    {
        return elements[idx];
    };
    void RemoveElement(int idx);
    int size()
    {
        return elements.size();
    } ;
    T& operator[] (int idx)
    {
        return *(elements[idx]);
    };
    T* create();
    template<class Q> T* create(Q);
    template<class Q, class R> T* create(Q,R);
    template<class Q, class R, class S> T* create(Q,R,S);
};

template<class T, int i>
CompactList<T, i>::~CompactList()
{
    removeAll();
    return;
}

template<class T, int i>
void CompactList<T, i>::removeAll()
//private function (called by ~)
{
    for (int j=elements.size() - 1; j >= 0 ; j--)
    {
        remove(elements[j]);
    }
    elements.pop_back(); 
}

template<class T, int i>
inline void CompactList<T, i>::remove(T* element)
//private function (called by removeAll (~) and RemoveElement)
//Removes a single element. No checking is performed on the pointer.
{
    element->~T();
    MemPool::free(element);
    return;
}

template<class T,int i>
inline T* CompactList<T,i>::create()
{
    T* newElement = new (MemPool::malloc()) T();
    return InsertElement(newElement);
}

template<class T, int i> template<class Q>
inline T* CompactList<T,i>::create(Q par1)
{
    T* newElement = new (MemPool::malloc()) T(par1);
    return InsertElement(newElement);
}

template<class T, int i> template<class Q, class R>
inline T* CompactList<T, i>::create(Q par1, R par2)
{
    T* newElement = new (MemPool::malloc()) T(par1,par2);
    return InsertElement(newElement);
}

template<class T, int i> template<class Q, class R, class S>
inline T* CompactList<T, i>::create(Q par1, R par2, S par3)
{
    T* newElement = new (MemPool::malloc()) T(par1,par2,par3);
    return InsertElement(newElement);
}

template<class T, int i>
inline T* CompactList<T,i>::InsertElement( T* newElement)
// private function!
{
    newElement->index = elements.size();
    elements.push_back(newElement);
    return newElement;
}

template<class T, int i>
inline void CompactList<T,i>::RemoveElement(int idx)
// check (only the outer!!) can eliminate?
{
    int size = elements.size();
    if (idx < size)
    {
        remove(elements[idx]);
        // if idx was not the last element, move last element to the gap
        if (idx != size -1)	
        {
            elements[idx] = elements.back();
            elements[idx]->index = idx;
        }
        elements.pop_back();
    }
    else
        cerr << "ERROR CompactList: trying to remove non-existing element. (index = " << idx << ", size = " << size << ")\n";
    return;
}

template<class T>
class CompactListItem
{
public:
    CompactListItem() : index(0) {} ;
    int index;
protected:
    virtual ~CompactListItem() {};
    template<class Q, int i> friend class CompactList;
};

#endif 
//COMPACTLIST_H_
