// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Meta-Header for the argument parser module.
 */

#pragma once

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/auxiliary.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>


/*!\defgroup argument_parser Argument Parser
 * \brief The Argument Parser Module
 *
 * ## The Argument Parser Class
 *
 * \copydetails seqan3::argument_parser
 *
 * ## Parsing Command Line Arguments
 *
 * \copydetails seqan3::argument_parser::parse
 *
 * ## Argument Validation
 *
 * The SeqAn 3 Argument Parser offers a validation mechanism for (positional_)options
 * via callables. You can pass any functor that fulfills the seqan3::validator_concept
 * and takes the value passed to the add_(positional_)option function call as
 * a parameter. We provide some commonly used functor that might come in handy:
 *
 * - seqan3::regex_validator
 * - seqan3::value_list_validator
 * - seqan3::arithmetic_range_validator
 * - seqan3::file_ext_validator
 */
