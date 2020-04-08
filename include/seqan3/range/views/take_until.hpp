// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::take_until and seqan3::views::take_until_or_throw.
 */

#pragma once

#include <seqan3/core/type_traits/iterator.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/type_traits/transformation_trait_or.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/detail.hpp>
#include <seqan3/range/detail/inherited_iterator_base.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/type_traits>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

// ============================================================================
//  view_take_until
// ============================================================================

/*!\brief The type returned by seqan3::views::take_until and seqan3::views::take_until_or_throw.
 * \tparam urng_t    The type of the underlying range, must model std::ranges::view.
 * \tparam fun_t     Type of the callable that will be evaluated on every member; must model
 *                   std::invocable with seqan3::reference_t<urng_t> as argument and return `bool`.
 * \tparam or_throw  Whether to throw an exception when the input is exhausted before the end of line is reached.
 * \implements std::ranges::view
 * \implements std::ranges::random_access_range
 * \ingroup views
 *
 * \details
 *
 * Note that most members of this class are generated by ranges::view_interface which is not yet documented here.
 */
template <std::ranges::view urng_t, typename fun_t, bool or_throw, bool and_consume>
class view_take_until : public std::ranges::view_interface<view_take_until<urng_t, fun_t, or_throw, and_consume>>
{
private:

    static_assert(std::invocable<fun_t, reference_t<urng_t>>,
                  "The functor type for views::take_until must model std::invocable<fun_t, reference_t<urng_t>>.");
    static_assert(std::boolean<std::result_of_t<fun_t&&(reference_t<urng_t>)>>,
                  "The functor type for views::take_until must return std::boolean.");

    //!\brief The underlying range.
    urng_t urange;

    //!\brief The functor.
    ranges::semiregular_t<fun_t> fun;

    //!\brief Whether this view is const_iterable or not.
    static constexpr bool const_iterable = const_iterable_range<urng_t> &&
                                           std::regular_invocable<fun_t, reference_t<urng_t>>;

    //!\brief The iterator type inherits from the underlying type, but overwrites several operators.
    //!\tparam rng_t Should be `urng_t` for defining #iterator and `urng_t const` for defining #const_iterator.
    template <typename rng_t>
    class iterator_type : public inherited_iterator_base<iterator_type<rng_t>, std::ranges::iterator_t<rng_t>>
    {
    private:
        //!\brief The iterator type of the underlying range.
        using base_base_t = std::ranges::iterator_t<rng_t>;
        //!\brief The CRTP wrapper type.
        using base_t      = inherited_iterator_base<iterator_type, std::ranges::iterator_t<rng_t>>;

        //!\brief The sentinel type is identical to that of the underlying range.
        using sentinel_type = std::ranges::sentinel_t<rng_t>;

        //!\brief Auxiliary type.
        using fun_ref_t = std::conditional_t<std::is_const_v<rng_t>,
                                             std::remove_reference_t<fun_t> const &,
                                             std::remove_reference_t<fun_t> &>;
        //!\brief Reference to the functor stored in the view.
        ranges::semiregular_t<fun_ref_t> fun;

    public:
        /*!\name Constructors, destructor and assignment
         * \brief Exceptions specification is implicitly inherited.
         * \{
         */
        constexpr iterator_type()                                      = default; //!< Defaulted.
        constexpr iterator_type(iterator_type const & rhs)             = default; //!< Defaulted.
        constexpr iterator_type(iterator_type && rhs)                  = default; //!< Defaulted.
        constexpr iterator_type & operator=(iterator_type const & rhs) = default; //!< Defaulted.
        constexpr iterator_type & operator=(iterator_type && rhs)      = default; //!< Defaulted.
        ~iterator_type()                                               = default; //!< Defaulted.

        //!\brief Constructor that delegates to the CRTP layer.
        iterator_type(base_base_t it) noexcept(noexcept(base_t{it})) :
            base_t{std::move(it)}
        {}

        //!\brief Constructor that delegates to the CRTP layer and initialises the callable.
        iterator_type(base_base_t it,
                      fun_ref_t _fun,
                      sentinel_type /*only used by the consuming iterator*/) noexcept(noexcept(base_t{it})) :
            base_t{std::move(it)}, fun{_fun}
        {}
        //!\}

        /*!\name Associated types
         * \brief All are derived from the base_base_t.
         * \{
         */

        //!\brief The difference type.
        using difference_type       = typename std::iterator_traits<base_base_t>::difference_type;
        //!\brief The value type.
        using value_type            = typename std::iterator_traits<base_base_t>::value_type;
        //!\brief The reference type.
        using reference             = typename std::iterator_traits<base_base_t>::reference;
        //!\brief The pointer type.
        using pointer               = typename std::iterator_traits<base_base_t>::pointer;
        //!\brief The iterator category tag.
        using iterator_category     = iterator_tag_t<base_base_t>;
        //!\}

        /*!\name Comparison operators
         * \brief We define comparison against self and against the sentinel.
         * \{
         */
        //!\brief Delegate comparison to base_base_t.
        bool operator==(iterator_type const & rhs) const
            noexcept(noexcept(std::declval<base_base_t &>() == std::declval<base_base_t &>()))
        //!\cond
            requires std::forward_iterator<base_base_t>
        //!\endcond
        {
            return *this->this_to_base() == *rhs.this_to_base();
        }

        //!\brief Evaluate functor, possibly throw.
        bool operator==(sentinel_type const & rhs) const
            noexcept(!or_throw &&
                     noexcept(std::declval<base_base_t &>() == std::declval<sentinel_type &>()) &&
                     noexcept(fun(std::declval<reference>())))
        {
            if (*this->this_to_base() == rhs) // [[unlikely]]
            {
                if constexpr (or_throw)
                    throw unexpected_end_of_input{"Reached end of input before functor evaluated to true."};
                else
                    return true;
            }

            return fun(**this);
        }

        //!\brief Switch lhs and rhs for comparison.
        friend bool operator==(sentinel_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(rhs == lhs))
        {
            return rhs == lhs;
        }

        //!\brief Switch lhs and rhs for comparison.
        bool operator!=(sentinel_type const & rhs) const
            noexcept(noexcept(std::declval<iterator_type &>() == rhs))
        {
            return !(*this == rhs);
        }

        //!\brief Delegate comparison to base_base_t.
        bool operator!=(iterator_type const & rhs) const
            noexcept(noexcept(std::declval<iterator_type &>() == rhs))
        //!\cond
            requires std::forward_iterator<base_base_t>
        //!\endcond
        {
            return !(*this == rhs);
        }

        //!\brief Switch lhs and rhs for comparison.
        friend bool operator!=(sentinel_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(rhs != lhs))
        {
            return rhs != lhs;
        }
        //!\}
    }; // class iterator_type

    //!\brief Special iterator type used when consuming behaviour is selected.
    //!\tparam rng_t Should be `urng_t` for defining #iterator and `urng_t const` for defining #const_iterator.
    template <typename rng_t>
    class iterator_type_consume_input : public inherited_iterator_base<iterator_type_consume_input<rng_t>, std::ranges::iterator_t<rng_t>>
    {
    private:
        //!\brief The iterator type of the underlying range.
        using base_base_t = std::ranges::iterator_t<rng_t>;
        //!\brief The CRTP wrapper type.
        using base_t      = inherited_iterator_base<iterator_type_consume_input, std::ranges::iterator_t<rng_t>>;

        //!\brief Auxiliary type.
        using fun_ref_t = std::conditional_t<std::is_const_v<rng_t>,
                                             std::remove_reference_t<fun_t> const &,
                                             std::remove_reference_t<fun_t> &>;
        //!\brief Reference to the functor stored in the view.
        ranges::semiregular_t<fun_ref_t> fun;

        //!\brief The sentinel type is identical to that of the underlying range.
        using sentinel_type = std::ranges::sentinel_t<rng_t>;

        //!\brief Whether this iterator has reached the end (cache is only used on pure input ranges).
        sentinel_type stored_end;

        //!\brief Whether the end was reached by evaluating the functor.
        bool at_end_gracefully = false;

    public:
        /*!\name Constructors, destructor and assignment
         * \brief Exceptions specification is implicitly inherited.
         * \{
         */
        //!\brief Defaulted.
        constexpr iterator_type_consume_input()                                                    = default;
        //!\brief Defaulted.
        constexpr iterator_type_consume_input(iterator_type_consume_input const & rhs)             = default;
        //!\brief Defaulted.
        constexpr iterator_type_consume_input(iterator_type_consume_input && rhs)                  = default;
        //!\brief Defaulted.
        constexpr iterator_type_consume_input & operator=(iterator_type_consume_input const & rhs) = default;
        //!\brief Defaulted.
        constexpr iterator_type_consume_input & operator=(iterator_type_consume_input && rhs)      = default;
        //!\brief Defaulted.
        ~iterator_type_consume_input()                                                             = default;

        //!\brief Constructor that delegates to the CRTP layer and initialises the callable.
        iterator_type_consume_input(base_base_t it,
                                    fun_ref_t _fun,
                                    sentinel_type sen) noexcept(noexcept(base_t{it})) :
            base_t{std::move(it)}, fun{_fun}, stored_end{std::move(sen)}
        {
            if ((*this->this_to_base() != stored_end) && fun(**this))
            {
                at_end_gracefully = true;
                ++(*this);
            }
        }
        //!\}

        /*!\name Associated types
         * \brief All are derived from the base_base_t.
         * \{
         */
        using difference_type       = typename std::iterator_traits<base_base_t>::difference_type;  //!< From base.
        using value_type            = typename std::iterator_traits<base_base_t>::value_type;       //!< From base.
        using reference             = typename std::iterator_traits<base_base_t>::reference;        //!< From base.
        using pointer               = typename std::iterator_traits<base_base_t>::pointer;          //!< From base.
        using iterator_category     = std::input_iterator_tag;                                      //!< Always input.
        //!\}

        /*!\name Arithmetic operators
         * \brief seqan3::detail::inherited_iterator_base operators are used unless specialised here.
         * \{
         */
        //!\brief Override pre-increment to implement consuming behaviour.
        iterator_type_consume_input & operator++()
            noexcept(noexcept(++std::declval<base_t &>()) &&
                     noexcept(std::declval<base_base_t &>() != std::declval<sentinel_type &>()) &&
                     noexcept(fun(std::declval<reference>())))
        {
            base_t::operator++();

            while ((*this->this_to_base() != stored_end) && fun(**this))
            {
                at_end_gracefully = true;
                base_t::operator++();
            }

            return *this;
        }

        //!\brief Post-increment implemented via pre-increment.
        iterator_type_consume_input operator++(int)
            noexcept(noexcept(++std::declval<iterator_type_consume_input &>()) &&
                     std::is_nothrow_copy_constructible_v<iterator_type_consume_input>)
        {
            iterator_type_consume_input cpy{*this};
            ++(*this);
            return cpy;
        }
        //!\}

        /*!\name Comparison operators
         * \brief We define comparison against self and against the sentinel.
         * \{
         */
        //!\brief Return the saved at_end state.
        bool operator==(sentinel_type const & rhs) const
            noexcept(!or_throw &&
                     noexcept(std::declval<base_base_t &>() != std::declval<sentinel_type &>()) &&
                     noexcept(fun(std::declval<reference>())))
        {
            if (at_end_gracefully)
                return true;

            if (*this->this_to_base() == rhs)
            {
                if constexpr (or_throw)
                    throw unexpected_end_of_input{"Reached end of input before functor evaluated to true."};
                else
                    return true;
            }

            return fun(**this);
        }

        //!\brief Return the saved at_end state.
        friend bool operator==(sentinel_type const & lhs, iterator_type_consume_input const & rhs)
            noexcept(noexcept(rhs == lhs))
        {
            return rhs == lhs;
        }

        //!\brief Return the saved at_end state.
        bool operator!=(sentinel_type const & rhs) const
            noexcept(noexcept(std::declval<iterator_type_consume_input &>() == rhs))
        {
            return !(*this == rhs);
        }

        //!\brief Return the saved at_end state.
        friend bool operator!=(sentinel_type const & lhs, iterator_type_consume_input const & rhs)
            noexcept(noexcept(rhs != lhs))
        {
            return rhs != lhs;
        }
        //!\}
    }; // class iterator_type_consume_input

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The reference_type.
    using reference         = reference_t<urng_t>;
    //!\brief The const_reference type is equal to the reference type if the underlying range is const-iterable.
    using const_reference   = detail::transformation_trait_or_t<seqan3::reference<urng_t const>, void>;
    //!\brief The value_type (which equals the reference_type with any references removed).
    using value_type        = value_type_t<urng_t>;
    //!\brief The size_type is void, because this range is never sized.
    using size_type         = void;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::ranges::range_difference_t<urng_t>;
    //!\brief The iterator type of this view (a random access iterator).
    using iterator          = std::conditional_t<and_consume && !std::ranges::forward_range<urng_t>,
                                                iterator_type_consume_input<urng_t>,
                                                iterator_type<urng_t>>;

    //!\brief The const_iterator type is equal to the iterator type if the underlying range is const-iterable.
    using const_iterator    = std::conditional_t<and_consume && !std::ranges::forward_range<urng_t>,
                                                 void,
                                                 detail::transformation_trait_or_t<
                                                    std::type_identity<iterator_type<urng_t const>>, void>>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    view_take_until()                                                  = default; //!< Defaulted.
    constexpr view_take_until(view_take_until const & rhs)             = default; //!< Defaulted.
    constexpr view_take_until(view_take_until && rhs)                  = default; //!< Defaulted.
    constexpr view_take_until & operator=(view_take_until const & rhs) = default; //!< Defaulted.
    constexpr view_take_until & operator=(view_take_until && rhs)      = default; //!< Defaulted.
    ~view_take_until()                                                 = default; //!< Defaulted.

    /*!\brief Construct from another range.
     * \param[in] _urange The underlying range.
     * \param[in] _fun    The functor that acts as termination criterium.
     */
    view_take_until(urng_t _urange, fun_t _fun)
        : urange{std::move(_urange)}, fun{std::forward<fun_t>(_fun)}
    {}

    /*!\brief Construct from another viewable_range.
     * \tparam rng_t      Type of the passed range; `urng_t` must be constructible from this.
     * \param[in] _urange The underlying range.
     * \param[in] _fun    The functor that acts as termination criterium.
     */
    template <std::ranges::viewable_range rng_t>
    //!\cond
        requires std::constructible_from<rng_t, std::ranges::all_view<rng_t>>
    //!\endcond
    view_take_until(rng_t && _urange, fun_t _fun)
        : view_take_until{std::views::all(std::forward<rng_t>(_urange)), std::move(_fun)}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to seqan3::views::take_until::end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return {std::ranges::begin(urange), static_cast<fun_t &>(fun), std::ranges::end(urange)};
    }

    //!\copydoc begin()
    const_iterator begin() const noexcept
        requires const_iterable
    {
        return {std::ranges::begin(urange), static_cast<fun_t const &>(fun), std::ranges::end(urange)};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
        requires const_iterable
    {
        return {std::ranges::begin(urange), static_cast<fun_t const &>(fun), std::ranges::end(urange)};
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto end() noexcept
    {
        return std::ranges::end(urange);
    }

    //!\copydoc end()
    auto end() const noexcept
        requires const_iterable
    {
        return std::ranges::cend(urange);
    }

    //!\copydoc end()
    auto cend() const noexcept
        requires const_iterable
    {
        return std::ranges::cend(urange);
    }
    //!\}
};

//!\brief Type deduction guide that strips references.
//!\relates seqan3::detail::view_take_until
template <typename urng_t, typename fun_t, bool or_throw = false, bool and_consume = false>
view_take_until(urng_t &&, fun_t) -> view_take_until<std::ranges::all_view<urng_t>, fun_t, or_throw, and_consume>;

// ============================================================================
//  take_until_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for views::take_until and views::take_until_or_throw.
 * \tparam or_throw Whether to throw an exception when the input is exhausted before the end of line is reached.
 */
template <bool or_throw, bool and_consume>
struct take_until_fn
{
    //!\brief Store the arguments and return a range adaptor closure object.
    template <typename fun_t>
    constexpr auto operator()(fun_t && fun) const
    {
        return adaptor_from_functor{*this, std::forward<fun_t>(fun)};
    }

    /*!\brief Call the view's constructor with the given parameters.
     * \tparam    urng_t The underlying range type; must model std::ranges::view.
     * \tparam    fun_t  The type of the callable; concept checks done in class.
     * \param[in] urange The underlying range.
     * \param[in] fun    The callable that will be evaluated on every element.
     * \returns An instance of seqan3::detail::view_take_until.
     */
    template <std::ranges::viewable_range urng_t, typename fun_t>
    constexpr auto operator()(urng_t && urange, fun_t && fun) const
    {
        return view_take_until<std::ranges::all_view<urng_t>, fun_t, or_throw, and_consume>
        {
            std::views::all(std::forward<urng_t>(urange)),
            std::forward<fun_t>(fun)
        };
    }
};

} // namespace seqan3::detail

// ============================================================================
//  views::take_until (adaptor instance definition)
// ============================================================================

namespace seqan3::views
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns elements from the underlying range until the functor evaluates to
 *                      true (or the end of the underlying range is reached).
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam fun_t        The type of the functor; must model std::invocable with seqan3::reference_t<urng_t>
 *                      and return a type convertible to `bool`.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] fun       The functor.
 * \returns             All elements of the underlying range up until (but excluding) the element that evaluates the
 *                      functor to true.
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/take_until.hpp}
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       |                                       | *preserved*¹                                       |
 * | std::ranges::bidirectional_range |                                       | *preserved*¹                                       |
 * | std::ranges::random_access_range |                                       | *preserved*¹                                       |
 * | std::ranges::contiguous_range    |                                       | *preserved*¹                                       |
 * |                                  |                                       |                                                    |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | *lost*                                             |
 * | std::ranges::common_range        |                                       | *lost*                                             |
 * | std::ranges::output_range        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_range     |                                       | *preserved*¹                                       |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   |                                       | seqan3::reference_t<urng_t>                        |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ¹ The marked properties are only *preserved* if the specified functor models
 * `std::regular_invocable<fun_t, reference_t<urng_t>`, i.e. applying the functor doesn't change the functor.
 * If the functor only models `std::invocable` and not `std::regular_invocable` these concepts are *lost*.
 *
 * Throwing: `seqan3::views::take_until_or_throw` and `seqan3::views::take_until_or_throw_and_consume` throw an exception
 * if the end of the underlying range is reached before their own termination criterium is met. This is useful
 * if you want a "strict" evaluation of the functor.
 *
 * Consuming: `seqan3::views::take_until_and_consume` and `seqan3::views::take_until_or_throw_and_consume` behave the same
 * as their non-consuming counter-parts if the underlying range models at least `std::forward_range`.
 * If, however, the underlying range is a pure `std::input_range`, the view will keep moving the underlying
 * iterator forward as long as the termination criterium holds and the underlying range is not at end.
 * This is useful for string tokenisation among other things.
 *
 * ### Example
 *
 * \include test/snippet/range/views/take_until.cpp
 *
 * \hideinitializer
 */
inline auto constexpr take_until = detail::take_until_fn<false, false>{};

// ============================================================================
//  views::take_until_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns elements from the underlying range until the functor evaluates to true
 *        (**throws** if the end of the underlying range is reached).
 * \throws seqan3::unexpected_end_of_input If the underlying range contains no element that satisfies the functor.
 * \ingroup views
 *
 * \copydetails seqan3::views::take_until
 * \hideinitializer
 */
inline auto constexpr take_until_or_throw = detail::take_until_fn<true, false>{};

// ============================================================================
//  views::take_until_and_consume (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns elements from the underlying range until the functor evaluates to true
 *        (or the end of the underlying range is reached; consumes end in single-pass ranges).
 * \throws seqan3::unexpected_end_of_input If the underlying range contains no element that satisfies the functor.
 * \ingroup views
 *
 * \copydetails seqan3::views::take_until
 * \hideinitializer
 */
inline auto constexpr take_until_and_consume = detail::take_until_fn<false, true>{};

// ============================================================================
//  views::take_until_or_throw_and_consume (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns elements from the underlying range until the functor evaluates to true
 *        (**throws** if the end of the underlying range is reached; consumes end in single-pass ranges).
 * \throws seqan3::unexpected_end_of_input If the underlying range contains no element that satisfies the functor.
 * \ingroup views
 *
 * \copydetails seqan3::views::take_until
 * \hideinitializer
 */
inline auto constexpr take_until_or_throw_and_consume = detail::take_until_fn<true, true>{};

//!\}

} // namespace seqan3::views
