// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_char_literal.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5 letter1{'A'_dna5};
    auto letter2 = 'A'_dna5;
}
//![main]
