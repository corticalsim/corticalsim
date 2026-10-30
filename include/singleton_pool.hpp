#ifndef INCLUDE_SINGLETON_POOL_HPP
#define INCLUDE_SINGLETON_POOL_HPP

#include <memory_resource>

namespace crts
{

    template <typename ElementType>
    struct SingletonPool
    {
      private:

        std::pmr::monotonic_buffer_resource mbr;
        std::pmr::polymorphic_allocator<ElementType> pa;
        SingletonPool():
            pa{ &mbr } {};

      public:

        static SingletonPool<ElementType>& instance()
        {
            static SingletonPool<ElementType> instance{};
            return instance;
        }

        static ElementType* malloc() { return SingletonPool<ElementType>::instance().pa.allocate(sizeof(ElementType)); }

        static void free(ElementType* /*element*/) {};
    };

} // namespace crts
#endif // INCLUDE_SINGLETON_POOL_HPP
