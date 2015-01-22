#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "outfile.h"

#include <boost/test/auto_unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>


BOOST_AUTO_TEST_SUITE(test_outfile)


namespace fs = boost::filesystem;
BOOST_AUTO_TEST_CASE(test_write_encoded_alignment)
{
  int codon_length = 6;
  fasta::Sequence s = fasta::make_sequence(">desc", "AAAAAALAAAAACAAAAA",
                                           codon_length);
  fasta::SequenceList sequences = {s};
  seq_data::SequenceData sequence_data;
  sequence_data.sequences = sequences;
  std::string filename_prefix = "tests/test_out_encoded.fasta";
  outfile::write_encoded_alignment(sequences, sequence_data, filename_prefix);
  fasta::FastaData fd = fasta::parse_fasta(filename_prefix + "_al", codon_length);
  fs::path p(filename_prefix + "_al");
  BOOST_CHECK(fs::exists(p));
  BOOST_CHECK_EQUAL(fd.sequences.size(), 1);
  for (size_t i = 0; i < s.residues.size(); ++i) {
    BOOST_CHECK_EQUAL(s.residues[i].codon, fd.sequences[0].residues[i].codon);
  }
  BOOST_CHECK_EQUAL(s.description, fd.sequences[0].description);
}


BOOST_AUTO_TEST_CASE(test_write_decoded_alignment)
{
  int codon_length = 6;
  fasta::Sequence s = fasta::make_sequence(">desc", "AAAAAALAAAAACAAAAA",
                                           codon_length);
  fasta::Sequence s_out = fasta::make_sequence(">desc", "ALC", 1);
  fasta::SequenceList sequences = {s};
  seq_data::SequenceData sequence_data;
  sequence_data.sequences = sequences;
  std::string filename_prefix = "tests/test_out_decoded.fasta";
  outfile::write_decoded_alignment(sequences, sequence_data, filename_prefix);
  fasta::FastaData fd = fasta::parse_fasta(filename_prefix + "_al", 1);
  fs::path p(filename_prefix + "_al");
  BOOST_CHECK(fs::exists(p));
  BOOST_CHECK_EQUAL(fd.sequences.size(), 1);
  for (size_t i = 0; i < s_out.residues.size(); ++i) {
    BOOST_CHECK_EQUAL(s_out.residues[i].codon,
                      fd.sequences[0].residues[i].codon);
  }
  BOOST_CHECK_EQUAL(s_out.description, fd.sequences[0].description);
}



BOOST_AUTO_TEST_SUITE_END()
