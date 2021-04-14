// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include <seqan3/std/iterator>

#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>

//NOTE(h-2): This class is extensively tested via *_file_input. This is just a minimal test.

struct fake_file_t : std::basic_istream<char>
{
    using base = std::basic_istream<char>;

    using base::base;

    using iterator = seqan3::detail::in_file_iterator<fake_file_t>;
    using value_type = char;
    using reference = const char&;
    using size_type = size_t;
    using difference_type = std::ptrdiff_t;
    using stream_ptr_t = std::unique_ptr<std::basic_istream<char>,
                                         std::function<void(std::basic_istream<char>*)>>;

    bool at_end{false};
    char record_buffer{};
    stream_ptr_t secondary_stream{};
    stream_ptr_t primary_stream{};
    std::streampos position_buffer{};

    void read_next_record()
    {
        position_buffer = secondary_stream->tellg();

        // at end if we could not read further
        if (secondary_stream->eof())
        {
            at_end = true;
            return;
        }

        secondary_stream->get(record_buffer);
    }

    iterator begin()
    {
        secondary_stream->seekg(0);
        this->read_next_record();
        return {*this};
    }

    static void stream_deleter_default(std::basic_istream<char> * ptr) { delete ptr; }

    fake_file_t(std::istringstream in) :
        primary_stream{new std::istringstream{std::move(in)}, stream_deleter_default}
    {
        secondary_stream = seqan3::detail::make_secondary_istream(*primary_stream);
    }
};

TEST(in_file_iterator, concepts)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    EXPECT_TRUE((std::input_iterator<it_t>));
}

TEST(in_file_iterator, member_types)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;
    EXPECT_TRUE((std::is_same_v<typename it_t::value_type,
                                char>));
    EXPECT_TRUE((std::is_same_v<typename it_t::reference,
                                const char &>));
    EXPECT_TRUE((std::is_same_v<typename it_t::const_reference,
                                const char &>));
    EXPECT_TRUE((std::is_same_v<typename it_t::difference_type,
                                std::ptrdiff_t>));
    EXPECT_TRUE((std::is_same_v<typename it_t::size_type,
                                size_t>));
    EXPECT_TRUE((std::is_same_v<typename it_t::iterator_category,
                                std::input_iterator_tag>));
}

TEST(in_file_iterator, operations)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    fake_file_t f{std::istringstream{"helloworld"}};

    // construct
    it_t it = f.begin();

    // deref
    EXPECT_EQ(*it, 'h');

    // pre-inc
    EXPECT_EQ(*(++it), 'e');

    // post-inc
    it++;

    // deref
    EXPECT_EQ(*it, 'l');
}

TEST(in_file_iterator, comparison)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    fake_file_t f{std::istringstream{"helloworld"}};
    it_t it = f.begin();

    // not at end
    EXPECT_FALSE(it == std::default_sentinel);

    // consume the entire range
    ++it; ++it; ++it; ++it; ++it; ++it; ++it; ++it; ++it; ++it; ++it;

    // at end
    EXPECT_TRUE(it == std::default_sentinel);
}

TEST(in_file_iterator, file_position)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    fake_file_t f{std::istringstream{"helloworld"}};
    it_t it = f.begin();
    auto beginning = it.file_position();

    // Go to the 6th character (w) and store it.
    ++it; ++it; ++it; ++it; ++it;
    EXPECT_EQ(*it, 'w');
    auto w_position = it.file_position();

    // Go back to the beginning.
    it.seek_to(beginning);
    EXPECT_EQ(*it, 'h');

    // Go directly to the w.
    it.seek_to(w_position);
    EXPECT_EQ(*it, 'w');
}
