#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main

#include "src/fasta.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

using namespace fasta;

BOOST_AUTO_TEST_SUITE(test_fasta)


BOOST_AUTO_TEST_CASE(test_parse_fasta)
{
    FastaData fd = parse_fasta("tests/TAU_SPECI.fasta.7c", 7, false, 0);
    BOOST_CHECK_EQUAL(fd.sequences.size(), 19);
    BOOST_CHECK_EQUAL(fd.probabilities.size(), 97);
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_nonexistent)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/nonexistent.fasta.7c", 7, false, 0),
        std::invalid_argument
    );
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_invalid_probabilities_format)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/wrong_probs_format.fasta.7c", 7, false, 0),
        std::runtime_error
    );
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_invalid_codons)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/invalid_codon.fasta.7c", 7, false,  0),
        std::runtime_error
    );
}

BOOST_AUTO_TEST_CASE(test_make_sequence)
{
  Sequence s = make_sequence("desc", "AaaaaaBbbbbbCccccc", 6);
  BOOST_CHECK_EQUAL(s.residues.size(), 3);
}

BOOST_AUTO_TEST_CASE(test_sequence_to_string)
{
  Sequence s = make_sequence("desc", "LSKAL", 1);
  std::string string_seq = sequence_to_string(s);
  std::string expected = "LSKAL";
  BOOST_CHECK_EQUAL(string_seq, expected);
}


BOOST_AUTO_TEST_SUITE_END()
