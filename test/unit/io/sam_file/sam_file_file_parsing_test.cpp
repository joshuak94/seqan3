// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream/byte.hpp>
#include <seqan3/core/debug_stream/optional.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/debug_stream/variant.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/input_format_concept.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/output_format_concept.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;
using seqan3::operator""_phred42;
using seqan3::operator""_tag;

struct verbose_file_data : public ::testing::Test
{
    verbose_file_data()
    {
        ref_sequences = std::vector<seqan3::dna5_vector>{ref_seq};
        ref_ids = std::vector<std::string>{ref_id};
        header = seqan3::sam_file_header{ref_ids};
        header.ref_id_info.emplace_back(ref_seq.size(), "");
        header.ref_dict[header.ref_ids()[0]] = 0; // set up header which is otherwise done on file level

        tag_dicts[0]["NM"_tag] = -7;
        tag_dicts[0]["AS"_tag] = 2;
        tag_dicts[0]["CC"_tag] = 300;
        tag_dicts[0]["cc"_tag] = -300;
        tag_dicts[0]["aa"_tag] = 'c';
        tag_dicts[0]["ff"_tag] = 3.1f;
        tag_dicts[0]["zz"_tag] = "str";
        tag_dicts[1]["bc"_tag] = std::vector<int8_t>{-3};
        tag_dicts[1]["bC"_tag] = std::vector<uint8_t>{3u, 200u};
        tag_dicts[1]["bs"_tag] = std::vector<int16_t>{-3, 200, -300};
        tag_dicts[1]["bS"_tag] = std::vector<uint16_t>{300u, 40u, 500u};
        tag_dicts[1]["bi"_tag] = std::vector<int32_t>{-3, 200, -66000};
        tag_dicts[1]["bI"_tag] = std::vector<uint32_t>{294967296u};
        tag_dicts[1]["bf"_tag] = std::vector<float>{3.5f, 0.1f, 43.8f};
        tag_dicts[1]["bH"_tag] = std::vector<std::byte>{std::byte{0x1A}, std::byte{0xE3}, std::byte{0x01}};
    }

    std::vector<seqan3::dna5_vector> seqs
    {
        "ACGT"_dna5,
        "AGGCTGNAG"_dna5,
        "GGAGTATA"_dna5
    };

    std::vector<std::string> ids
    {
        "read1",
        "read2",
        "read3"
    };

    std::vector<std::vector<seqan3::phred42>> quals
    {
        { "!##$"_phred42 },
        { "!##$&'()*"_phred42 },
        { "!!*+,-./"_phred42 },
    };

    std::vector<int32_t> offsets
    {
        1,
        0,
        1
    };

    seqan3::dna5_vector ref_seq = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;

    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped1 = {'A'_dna5, 'C'_dna5, 'T'_dna5, seqan3::gap{}};
    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped2 = {'C'_dna5, 'T'_dna5, 'G'_dna5, 'A'_dna5,
                                                                 'T'_dna5, 'C'_dna5, 'G'_dna5, 'A'_dna5, 'G'_dna5};
    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped3 = {'T'_dna5, seqan3::gap{}, 'G'_dna5, seqan3::gap{},
                                                                 'A'_dna5, seqan3::gap{}, 'T'_dna5, 'C'_dna5};

    std::string ref_id = "ref";

    std::vector<int32_t> ref_offsets
    {
        0,
        1,
        2
    };

    std::vector<std::pair<std::vector<seqan3::gapped<seqan3::dna5>>,
                          std::vector<seqan3::gapped<seqan3::dna5>>>> alignments
    {
        {ref_seq_gapped1, std::vector<seqan3::gapped<seqan3::dna5>>{'C'_dna5, seqan3::gap{}, 'G'_dna5, 'T'_dna5}},
        {ref_seq_gapped2, std::vector<seqan3::gapped<seqan3::dna5>>{'A'_dna5, 'G'_dna5, 'G'_dna5, 'C'_dna5, 'T'_dna5,
                                                                    'G'_dna5, 'N'_dna5, seqan3::gap{}, 'A'_dna5}},
        {ref_seq_gapped3, std::vector<seqan3::gapped<seqan3::dna5>>{'G'_dna5, seqan3::gap{}, 'A'_dna5, 'G'_dna5,
                                                                    'T'_dna5, 'A'_dna5, seqan3::gap{}, 'T'_dna5}}
    };

    std::vector<seqan3::sam_flag> flags
    {
        seqan3::sam_flag{41u},
        seqan3::sam_flag{42u},
        seqan3::sam_flag{43u}
    };

    std::vector<uint8_t> mapqs
    {
        61u,
        62u,
        63u
    };

    std::vector<std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>> mates
    {
        {0, 9, 300},
        {0, 9, 300},
        {0, 9, 300}
    };

    std::vector<seqan3::sam_tag_dictionary> tag_dicts
    {
        seqan3::sam_tag_dictionary{},
        seqan3::sam_tag_dictionary{},
        seqan3::sam_tag_dictionary{}
    };

    std::vector<seqan3::dna5_vector> ref_sequences{};
    std::vector<std::string> ref_ids{};
    seqan3::sam_file_header<std::vector<std::string>> header{};

    std::vector<std::streampos> file_positions{};
};

TEST_F(verbose_file_data, read_in_bam)
{
    std::filesystem::path in_path{CURDIR"simple_three_verbose_reads.bam"};
    seqan3::sam_file_input fin{in_path};

    size_t i{0};
    auto rec = fin.begin();
    for (; rec != fin.end(); ++rec)
    {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        EXPECT_EQ(seqan3::get<seqan3::field::seq>(*rec), this->seqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::id>(*rec), this->ids[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::qual>(*rec), this->quals[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::offset>(*rec), this->offsets[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::ref_id>(*rec), 0);
        EXPECT_EQ(*seqan3::get<seqan3::field::ref_offset>(*rec), this->ref_offsets[i]);
        // EXPECT_RANGE_EQ(std::get<0>(seqan3::get<seqan3::field::alignment>(*rec)), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>(seqan3::get<seqan3::field::alignment>(*rec)), std::get<1>(this->alignments[i]));
        EXPECT_EQ(seqan3::get<seqan3::field::flag>(*rec), this->flags[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mapq>(*rec), this->mapqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mate>(*rec), this->mates[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::tags>(*rec), this->tag_dicts[i]);
#pragma GCC diagnostic pop

        EXPECT_EQ((*rec).sequence(), this->seqs[i]);
        EXPECT_EQ((*rec).id(), this->ids[i]);
        EXPECT_EQ((*rec).base_qualities(), this->quals[i]);
        EXPECT_EQ((*rec).sequence_position(), this->offsets[i]);
        EXPECT_EQ((*rec).reference_id(), 0);
        EXPECT_EQ(*(*rec).reference_position(), this->ref_offsets[i]);
        // EXPECT_RANGE_EQ(std::get<0>((*rec).alignment()), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>((*rec).alignment()), std::get<1>(this->alignments[i]));
        EXPECT_EQ((*rec).flag(), this->flags[i]);
        EXPECT_EQ((*rec).mapping_quality(), this->mapqs[i]);
        EXPECT_EQ((*rec).mate_reference_id(), std::get<0>(this->mates[i]));
        EXPECT_EQ((*rec).mate_position(), std::get<1>(this->mates[i]));
        EXPECT_EQ((*rec).template_length(), std::get<2>(this->mates[i]));
        EXPECT_EQ((*rec).tags(), this->tag_dicts[i]);

        this->file_positions.push_back(rec.file_position());
        ++i;
    }
    --i;

    for (; i > 0; --i)
    {
        rec.seek_to(this->file_positions[i]);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        EXPECT_EQ(seqan3::get<seqan3::field::seq>(*rec), this->seqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::id>(*rec), this->ids[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::qual>(*rec), this->quals[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::offset>(*rec), this->offsets[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::ref_id>(*rec), 0);
        EXPECT_EQ(*seqan3::get<seqan3::field::ref_offset>(*rec), this->ref_offsets[i]);
        // EXPECT_RANGE_EQ(std::get<0>(seqan3::get<seqan3::field::alignment>(*rec)), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>(seqan3::get<seqan3::field::alignment>(*rec)), std::get<1>(this->alignments[i]));
        EXPECT_EQ(seqan3::get<seqan3::field::flag>(*rec), this->flags[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mapq>(*rec), this->mapqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mate>(*rec), this->mates[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::tags>(*rec), this->tag_dicts[i]);
#pragma GCC diagnostic pop

        EXPECT_EQ((*rec).sequence(), this->seqs[i]);
        EXPECT_EQ((*rec).id(), this->ids[i]);
        EXPECT_EQ((*rec).base_qualities(), this->quals[i]);
        EXPECT_EQ((*rec).sequence_position(), this->offsets[i]);
        EXPECT_EQ((*rec).reference_id(), 0);
        EXPECT_EQ(*(*rec).reference_position(), this->ref_offsets[i]);
        // EXPECT_RANGE_EQ(std::get<0>((*rec).alignment()), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>((*rec).alignment()), std::get<1>(this->alignments[i]));
        EXPECT_EQ((*rec).flag(), this->flags[i]);
        EXPECT_EQ((*rec).mapping_quality(), this->mapqs[i]);
        EXPECT_EQ((*rec).mate_reference_id(), std::get<0>(this->mates[i]));
        EXPECT_EQ((*rec).mate_position(), std::get<1>(this->mates[i]));
        EXPECT_EQ((*rec).template_length(), std::get<2>(this->mates[i]));
        EXPECT_EQ((*rec).tags(), this->tag_dicts[i]);

        EXPECT_EQ(rec.file_position(), this->file_positions[i]);
    }
}

TEST_F(verbose_file_data, read_in_sam)
{

    std::filesystem::path in_path{CURDIR"simple_three_verbose_reads.sam"};
    seqan3::sam_file_input fin{in_path};

    size_t i{0};
    auto rec = fin.begin();
    for (; rec != fin.end(); ++rec)
    {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        EXPECT_EQ(seqan3::get<seqan3::field::seq>(*rec), this->seqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::id>(*rec), this->ids[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::qual>(*rec), this->quals[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::offset>(*rec), this->offsets[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::ref_id>(*rec), 0);
        EXPECT_EQ(*seqan3::get<seqan3::field::ref_offset>(*rec), this->ref_offsets[i]);
        // EXPECT_RANGE_EQ(std::get<0>(seqan3::get<seqan3::field::alignment>(*rec)), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>(seqan3::get<seqan3::field::alignment>(*rec)), std::get<1>(this->alignments[i]));
        EXPECT_EQ(seqan3::get<seqan3::field::flag>(*rec), this->flags[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mapq>(*rec), this->mapqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mate>(*rec), this->mates[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::tags>(*rec), this->tag_dicts[i]);
#pragma GCC diagnostic pop

        EXPECT_EQ((*rec).sequence(), this->seqs[i]);
        EXPECT_EQ((*rec).id(), this->ids[i]);
        EXPECT_EQ((*rec).base_qualities(), this->quals[i]);
        EXPECT_EQ((*rec).sequence_position(), this->offsets[i]);
        EXPECT_EQ((*rec).reference_id(), 0);
        EXPECT_EQ(*(*rec).reference_position(), this->ref_offsets[i]);
        // EXPECT_RANGE_EQ(std::get<0>((*rec).alignment()), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>((*rec).alignment()), std::get<1>(this->alignments[i]));
        EXPECT_EQ((*rec).flag(), this->flags[i]);
        EXPECT_EQ((*rec).mapping_quality(), this->mapqs[i]);
        EXPECT_EQ((*rec).mate_reference_id(), std::get<0>(this->mates[i]));
        EXPECT_EQ((*rec).mate_position(), std::get<1>(this->mates[i]));
        EXPECT_EQ((*rec).template_length(), std::get<2>(this->mates[i]));
        EXPECT_EQ((*rec).tags(), this->tag_dicts[i]);

        this->file_positions.push_back(rec.file_position());
        ++i;
    }
    --i;

    for (; i > 0; --i)
    {
        rec.seek_to(this->file_positions[i]);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        EXPECT_EQ(seqan3::get<seqan3::field::seq>(*rec), this->seqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::id>(*rec), this->ids[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::qual>(*rec), this->quals[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::offset>(*rec), this->offsets[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::ref_id>(*rec), 0);
        EXPECT_EQ(*seqan3::get<seqan3::field::ref_offset>(*rec), this->ref_offsets[i]);
        // EXPECT_RANGE_EQ(std::get<0>(seqan3::get<seqan3::field::alignment>(*rec)), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>(seqan3::get<seqan3::field::alignment>(*rec)), std::get<1>(this->alignments[i]));
        EXPECT_EQ(seqan3::get<seqan3::field::flag>(*rec), this->flags[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mapq>(*rec), this->mapqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mate>(*rec), this->mates[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::tags>(*rec), this->tag_dicts[i]);
#pragma GCC diagnostic pop

        EXPECT_EQ((*rec).sequence(), this->seqs[i]);
        EXPECT_EQ((*rec).id(), this->ids[i]);
        EXPECT_EQ((*rec).base_qualities(), this->quals[i]);
        EXPECT_EQ((*rec).sequence_position(), this->offsets[i]);
        EXPECT_EQ((*rec).reference_id(), 0);
        EXPECT_EQ(*(*rec).reference_position(), this->ref_offsets[i]);
        // EXPECT_RANGE_EQ(std::get<0>((*rec).alignment()), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>((*rec).alignment()), std::get<1>(this->alignments[i]));
        EXPECT_EQ((*rec).flag(), this->flags[i]);
        EXPECT_EQ((*rec).mapping_quality(), this->mapqs[i]);
        EXPECT_EQ((*rec).mate_reference_id(), std::get<0>(this->mates[i]));
        EXPECT_EQ((*rec).mate_position(), std::get<1>(this->mates[i]));
        EXPECT_EQ((*rec).template_length(), std::get<2>(this->mates[i]));
        EXPECT_EQ((*rec).tags(), this->tag_dicts[i]);

        EXPECT_EQ(rec.file_position(), this->file_positions[i]);
    }
}
