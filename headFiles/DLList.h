#ifndef DLLIST_H_
#define DLLIST_H_
#include <boost/pool/singleton_pool.hpp>
#include <new>

using namespace std;

//DLList Class v1.0
//Simon Tindemans, 1/4/2008
//The DLList class provides a doubly-linked list that gets around a number of problems
//with the STL list class:
//1) The STL class always copies the elements on insertion. DLList initializes them 'on the spot',
//making it suitable for non-copyable classes.
//2) The STL class does not allow for conversion from data pointers to list element iterators. The
//DLList class makes this possible, leading to easier deletion (no need to store/retrieve iterators) and
//list traversal from within objects.
//3) When transferring elements (single-element splice), the STL list invalidates iterators. DLList
//allows one to transfer elements without invalidating pointers.

//The DLList list can be constructed using two types of list elements.
//First, the DLBaseItem class should be used as a base class for the individual list
//items. This has the advantage that the pointers to the list container elements and
//the individual objects are the same.

//The second method is using DLContainerItem objects. These objects contain an object
//of the type that is to be stored. This has the advantage that the objects don't need
//to be compiled for list use.

//Lists of objects of type T are declared as follows:
//DLList<T> name;  					[for DLBaseItem-derived classes]
//DLList<DLContainerItem<T> > name;	[for other classes]

template<class T, int AllocatorGranularity = 32>
class DLList
{
private:
    int numElements;
    //This class keeps a pointer to both the first and the last element.
    //This is not strictly necessary (first is enough), but this has the advantage
    //that elements can be inserted in the back and iterated over from the front of
    //the set. In the absence of deletions, this leads to a traversal of the elements
    //in the order of creation - and hopefully in the order of memory allocation.
    T* firstElement;
    T* lastElement;
    //definitions for the Boost memory pool
    //Currently, it uses a *non*-thread safe allocation mechanism (for speed) and a user-defined granularity
    struct PoolTag {};
    typedef boost::singleton_pool<PoolTag, sizeof(T), boost::default_user_allocator_new_delete, boost::details::pool::null_mutex, AllocatorGranularity> MemPool;
    void unlinkElement(T*);
    T* linkElement(T*);
public:
    DLList() : numElements(0), firstElement(NULL), lastElement(NULL) {};
    ~DLList();
    inline int size()
    {
        return numElements;
    }
    inline T* first()
    {
        return firstElement;
    }
    inline T* last()
    {
        return lastElement;
    }
    inline T* next(T* el)
    {
        return el->nextElement;
    }
    inline T* previous(T* el)
    {
        return el->previousElement;
    }
    T* create();
    template<class Q> T* create(Q);
    template<class Q, class R> T* create(Q,R);
    template<class Q, class R, class S> T* create(Q,R,S);
    void remove(T*);
    void removeAll();
    void import(DLList<T, AllocatorGranularity>&, T*);
    void importSet(DLList<T, AllocatorGranularity>&, T*, T*);
};

template<class T, int i>
DLList<T, i>::~DLList()
{
    removeAll();
    return;
}

template<class T, int i>
void DLList<T, i>::removeAll()
{
    while (firstElement != NULL)
        remove(firstElement);
}

template<class T, int i>
inline T* DLList<T, i>::linkElement(T* newElement)
//Inserts a new element into the linked list.
//Insertion takes place at the back of the list for memory-traversal speed considerations.
{
    numElements++;
    newElement->previousElement = lastElement;
    newElement->nextElement = NULL;
    if (lastElement != NULL)
        lastElement->nextElement = newElement;
    else // if list was empty
        firstElement = newElement;
    lastElement = newElement;
    return newElement;
}

template<class T, int i>
inline void DLList<T, i>::unlinkElement(T* element)
//Removes an element from the linked list.
{
    numElements--;
    if (element == firstElement)
        firstElement = element->nextElement;
    if (element == lastElement)
        lastElement = element->previousElement;
    if (element->nextElement != NULL)
        element->nextElement->previousElement = element->previousElement;
    if (element->previousElement != NULL)
        element->previousElement->nextElement = element->nextElement;
    return;
}

template<class T, int i>
inline void DLList<T, i>::remove(T* element)
//Removes a single element. No checking is performed on the pointer.
{
    unlinkElement(element);
    element->~T();
    MemPool::free(element);
    return;
}

/// The following functions create a new element from the heap and insert it into the linked list.
//These differ only in the number of arguments that is passed to the constructor.
template<class T,int i>
inline T* DLList<T,i>::create()
{
    T* newElement = new (MemPool::malloc()) T();
    return linkElement(newElement);
}

template<class T, int i> template<class Q>
inline T* DLList<T,i>::create(Q par1)
{
    T* newElement = new (MemPool::malloc()) T(par1);
    return linkElement(newElement);
}

template<class T, int i> template<class Q, class R>
inline T* DLList<T, i>::create(Q par1, R par2)
{
    T* newElement = new (MemPool::malloc()) T(par1,par2);
    return linkElement(newElement);
}

template<class T, int i> template<class Q, class R, class S>
inline T* DLList<T, i>::create(Q par1, R par2, S par3)
{
    T* newElement = new (MemPool::malloc()) T(par1,par2,par3);
    return linkElement(newElement);
}

template<class T, int i>
//This functions transfers 'element' from 'oldList' to this list.
inline void DLList<T, i>::import(DLList<T, i>& oldList, T* element)
{
    oldList.unlinkElement(element);
    linkElement(element);
    return;
}

template<class T, int i>
inline void DLList<T, i>::importSet(DLList<T, i>& oldList, T* firstElement, T* lastElement)
//This functions transfers a number of elements from 'oldList' to this list.
//WARNING: not efficient at all!
{
    T* current = firstElement;
    T* temp;
    while (current != NULL)
    {
        temp = current;
        current = current->next();
        import(oldList, temp);
        if (temp == lastElement)
            break;
    }
    return;
}

template<class T>
//this list item should be used as a base class [type 1 - see above]
class DLBaseItem
{
public:
    template<class X, int i> friend class DLList;
private:
    T* nextElement;
    T* previousElement;
protected:
    // make sure that this cannot be deleted...
    virtual ~DLBaseItem() {};
public:
    inline T* next()
    {
        return nextElement;
    }
    inline T* previous()
    {
        return previousElement;
    }
};

template<class T>
//This is a stand-alone version of the list item. It contains a data object of type T
//and functions to convert between the pointers to the data and pointers to the element
//itself
class DLContainerItem
{
    friend class DLList<DLContainerItem<T> >;
private:
    DLContainerItem<T>* nextElement;
    DLContainerItem<T>* previousElement;
protected:
    // make sure that items can only be deleted by DLList<DLContainerItem<T> > objects
    virtual ~DLContainerItem() {};
public:
    T data;
    inline DLContainerItem<T>* next()
    {
        return nextElement;
    }
    inline DLContainerItem<T>* previous()
    {
        return previousElement;
    }
    inline T* dataPtr()
    {
        return &data;
    }
    // very ugly hack, but it works [tested with G++ & icc]. For a data object with an unknown location in a known list, call
    // list->first()->containerPtr(data*);
    // for some reason, I need to convert to void* and char*... This should have been done in a better way...
    inline DLContainerItem<T>* containerPtr(T* dPtr)
    {
        return static_cast<DLContainerItem<T>*>(
                   static_cast<void*>(
                       static_cast<short*>(static_cast<void*>(dPtr)) - (static_cast<short*>(static_cast<void*>(&data)) - static_cast<short*>(static_cast<void*>(this)))
                   )
               );
    }
    DLContainerItem() {}
    template<class Q> inline DLContainerItem(Q par1) : data(par1) {}
    template<class Q, class R> inline DLContainerItem(Q par1,R par2) : data(par1,par2) {}
    template<class Q, class R, class S> inline DLContainerItem(Q par1,R par2,S par3) : data(par1,par2,par3) {}
};

#endif //DLLIST_H_
