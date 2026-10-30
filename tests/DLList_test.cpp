#include <cstdint>
#include <gtest/gtest.h>
#include <scoped_allocator>
#include "DLList.hpp"
#include <tuple>

namespace
{

    struct DummyTestType: DLBaseItem<DummyTestType>
    {
        using PlayloadType = std::uint64_t;
        PlayloadType payload1;
        PlayloadType payload2;
        PlayloadType payload3;

        DummyTestType():
            payload1{ 71UL },
            payload2{ 72UL },
            payload3{ 73UL }
        {
        }

        explicit DummyTestType(const PlayloadType& p):
            payload1{ p },
            payload2{ 72UL },
            payload3{ 73UL }
        {
        }

        DummyTestType(const PlayloadType& p1, const PlayloadType& p2):
            payload1{ p1 },
            payload2{ p2 },
            payload3{ 73UL }
        {
        }

        DummyTestType(const PlayloadType& p1, const PlayloadType& p2, const PlayloadType& p3):
            payload1{ p1 },
            payload2{ p2 },
            payload3{ p3 }
        {
        }

        bool operator==(const DummyTestType& other) const
        {
            return this->payload1 == other.payload1 && this->payload2 == other.payload2
                   && this->payload3 == other.payload3;
        }

        DummyTestType(const DummyTestType&) = default;
        DummyTestType(DummyTestType&&) = delete;
        DummyTestType& operator=(const DummyTestType&) = default;
        DummyTestType& operator=(DummyTestType&&) = delete;
        ~DummyTestType() override = default;
    };

} // namespace

TEST(DLBaseItemTest, base_item_can_not_be_created_outside_of_the_DLList)
{
    static_assert(!std::is_destructible_v<DLBaseItem<DummyTestType>>,
                  "DLBaseItem can not be constructed outside of DLList, since it's destructor is protected");
}

TEST(DLListTest, empty_DLList_can_be_created)
{
    // When creating a DLList with the default constructor
    DLList<DummyTestType> dl_list{};
    // Then it's initial size is 0
    EXPECT_EQ(dl_list.size(), 0);
}

TEST(DLListTest, size_reports_correct_size)
{
    // Given a newly created DLList
    DLList<DummyTestType> dl_list{};

    // When we create several new elements
    constexpr std::size_t num_new_elements{ 17 };

    for (std::size_t i = 0; i < num_new_elements; ++i)
    {
        std::ignore = dl_list.create();
    }

    // Then the list size corresponds to the number of new elements
    EXPECT_EQ(dl_list.size(), num_new_elements);
}

TEST(DLListTest, first_returns_the_pointer_to_the_first_crteated_element)
{
    // Given a newly created DLList, to which several items were added
    DLList<DummyTestType> dl_list{};
    auto* first_el = dl_list.create();
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();

    // When we call first on the list
    auto* retreived_first_element = dl_list.first();

    // Then we retreave the same first element that was created in the beginning
    EXPECT_EQ(first_el, retreived_first_element);
}

TEST(DLListTest, last_returns_the_pointer_to_the_last_crteated_element)
{
    // Given a newly created DLList, to which several items were added
    DLList<DummyTestType> dl_list{};
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();
    auto* last_el = dl_list.create();

    // When we call last on the list
    auto* retreived_last_element = dl_list.last();

    // Then we retreave the same last element that was created in the beginning
    EXPECT_EQ(last_el, retreived_last_element);
}

TEST(DLListTest, next_returns_the_pointer_to_the_next_crteated_element)
{
    // Given a newly created DLList, to which several items were added
    DLList<DummyTestType> dl_list{};
    auto* first_el = dl_list.create();
    auto* next_el = dl_list.create();
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();

    // When we call next on the first element
    auto* retreived_next_element = first_el->next();

    // Then we retreave the same next element that was created in the beginning
    EXPECT_EQ(next_el, retreived_next_element);
}

TEST(DLListTest, calling_last_on_the_last_element_returns_a_null_pointer)
{
    // Given a newly created DLList, to which several items were added
    DLList<DummyTestType> dl_list{};
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();
    auto* last_el = dl_list.create();

    // When we call next on the last element
    auto* one_past_last = last_el->next();

    // Then we get a nullptr back
    EXPECT_EQ(one_past_last, nullptr);
}

TEST(DLListTest, calling_first_on_the_first_element_returns_a_null_pointer)
{
    // Given a newly created DLList, to which several items were added
    DLList<DummyTestType> dl_list{};
    auto* first_el = dl_list.create();
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();
    std::ignore = dl_list.create();

    // When we call previous on the first elemenk
    auto* one_before_first = first_el->previous();

    // Then we get a nullptr back
    EXPECT_EQ(one_before_first, nullptr);
}

TEST(DLListTest, element_created_onto_the_list_with_argumentless_create_is_default_constructed)
{
    // Given a  DLList
    DLList<DummyTestType> dl_list{};

    // When we add an element with create() without giving create arguments
    auto* new_el = dl_list.create();

    // Then the created element will be default constructed from the underlying data type
    auto expected = DummyTestType{};
    EXPECT_EQ(*new_el, expected);
}

TEST(DLListTest, element_created_onto_the_list_with1_argument_create_uses_one_argument_constructor)
{
    // Given a  DLList
    DLList<DummyTestType> dl_list{};

    // When we add an element with create(payload1)
    auto payload1 = 99UL;
    auto expected = DummyTestType{ payload1 };
    auto* new_el = dl_list.create(payload1);

    // Then the created element will be constructed using the one argument constructor of the onderlying data type
    EXPECT_EQ(*new_el, expected);
}

TEST(DLListTest, element_created_onto_the_list_with2_argument_create_uses_one_argument_constructor)
{
    // Given a  DLList
    DLList<DummyTestType> dl_list{};

    // When we add an element with create(payload1, payload2)
    auto payload1 = 99UL;
    auto payload2 = 77UL;
    auto expected = DummyTestType{ payload1, payload2 };
    auto* new_el = dl_list.create(payload1, payload2);

    // Then the created element will be constructed using the two argument constructor of the onderlying data type
    EXPECT_EQ(*new_el, expected);
}

TEST(DLListTest, element_created_onto_the_list_with3_argument_create_uses_one_argument_constructor)
{
    // Given a  DLList
    DLList<DummyTestType> dl_list{};

    // When we add an element with create(payload1, payload2, payload3)
    auto payload1 = 99UL;
    auto payload2 = 93UL;
    auto payload3 = 11UL;
    auto expected = DummyTestType{ payload1, payload2, payload3 };
    auto* new_el = dl_list.create(payload1, payload2, payload3);

    // Then the created element will be constructed using the two argument constructor of the onderlying data type
    EXPECT_EQ(*new_el, expected);
}

TEST(DLListTest, element_created_onto_the_list_with_explicit_create_has_the_provided_value)
{
    // Given a DLList and data compatible with the list element type
    DLList<DummyTestType> dl_list{};
    DummyTestType element_data{ 11 };

    // When we add an element with create(element)
    auto* new_el = dl_list.create(element_data);

    // Then the created element contain the same data as the original
    EXPECT_EQ(new_el->payload1, element_data.payload1);
}
