#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "f_config.h"
#include "msa.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE(test_msa)

BOOST_AUTO_TEST_CASE(test_run_msa)
{
  int codon_length = 7;
  fasta::Sequence s1;
  // AKLCAKL
  s1 = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
                                 "SAAAZAA"
                                 "YAAAAAAAAAAAAAKAAAAAA"
                                 "LAAAAAA", codon_length);
  fasta::Sequence s2;
  // AKLAKL
  s2 = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
                                 "YAAAZAA"
                                 "AAAAAAAKAAAAAALAAAAAA"
                                 "RAAAAAA", codon_length);
  FeatureNamesList feature_list = {"ptm_phosph0", "ptm_phosph1",
                                   "ptm_phosph2", "ptm_phosph3",
                                   "ptm_phosphP", "ptm_acet0",
                                   "ptm_acet1", "ptm_acet2",
                                   "ptm_acet3", "ptm_Nglyc0",
                                   "ptm_Nglyc1", "ptm_Nglyc2",
                                   "ptm_Nglyc3", "ptm_amid0",
                                   "ptm_amid1", "ptm_amid2",
                                   "ptm_amid3", "ptm_hydroxy0",
                                   "ptm_hydroxy1", "ptm_hydroxy2",
                                   "ptm_hydroxy3", "ptm_methyl0",
                                   "ptm_methyl1", "ptm_methyl2",
                                   "ptm_methyl3", "ptm_Oglyc0",
                                   "ptm_Oglyc1", "ptm_Oglyc2",
                                   "ptm_Oglyc3"}; 
  std::map<std::string, double> probabilities;
  f_config::FeatureSettingsMap f_set;
  fasta::SequenceList query_seq_list = {s1};
  fasta::SequenceList sequences = {s1, s2};
  int domain_modifier = 4;
  int motif_modifier = 3;
  int ptm_modifier = 10;
  double gap_open_pen = 5;
  double gap_ext_pen = 1;
  double end_pen = 1;
  bool one_round = false;
  seq_data::SequenceData sequence_data;
  sequence_data.sequences = sequences;
  sequence_data.feature_list = feature_list;
  std::vector<fasta::SequenceList> alignment;
  alignment = msa::run_msa(sequence_data, f_set, gap_open_pen, gap_ext_pen,
                           end_pen, domain_modifier, motif_modifier, 
                           ptm_modifier, codon_length, one_round);
  fasta::Sequence e_s1;
  // AKLCAKL
  e_s1 = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
                                   "SAAAZAA"
                                   "YAAAAAAAAAAAAAKAAAAAA"
                                   "LAAAAAA", codon_length);
  fasta::Sequence e_s2;
  // AKLAKL
  e_s2 = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
                                   "YAAAZAA"
                                   "-AAAAAAAAAAAAAKAAAAAA"
                                   "LAAAAAA", codon_length);
  fasta::Sequence e_s1_lower;
  // AKLCAKL
  e_s1_lower = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
                                         "SAAAZAA"
                                         "YAAAAAAAAAAAAAKAAAAAA"
                                         "LAAAAAA", codon_length);
  fasta::Sequence e_s2_lower;
  // AKLAKL
  e_s2_lower = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
                                         "YAAAZAA"
                                         "-AAAAAAAAAAAAAKAAAAAA"
                                         "lAAAAAA", codon_length);
  std::vector<fasta::SequenceList> expected_alignment = {{e_s1, e_s2},
                                                         {e_s1_lower,
                                                          e_s2_lower}};

  BOOST_CHECK_EQUAL(alignment[0][0].residues.size(),
                    expected_alignment[0][0].residues.size());
  BOOST_CHECK_EQUAL(alignment[0][1].residues.size(),
                    expected_alignment[0][1].residues.size());
  BOOST_CHECK_EQUAL(alignment[1][0].residues.size(),
                    expected_alignment[1][0].residues.size());
  BOOST_CHECK_EQUAL(alignment[1][1].residues.size(),
                    expected_alignment[1][1].residues.size());

  for (size_t i = 0; i < alignment[0][0].residues.size(); ++i) {
    BOOST_CHECK_EQUAL(alignment[0][0].residues[i].codon,
                      expected_alignment[0][0].residues[i].codon);
  }
  for (size_t i = 0; i < alignment[0][1].residues.size(); ++i) {
    BOOST_CHECK_EQUAL(alignment[0][1].residues[i].codon,
                      expected_alignment[0][1].residues[i].codon);
  }
  for (size_t i = 0; i < alignment[1][0].residues.size(); ++i) {
    BOOST_CHECK_EQUAL(alignment[1][0].residues[i].codon,
                      expected_alignment[1][0].residues[i].codon);
  }
  for (size_t i = 0; i < alignment[1][1].residues.size(); ++i) {
    BOOST_CHECK_EQUAL(alignment[1][1].residues[i].codon,
                      expected_alignment[1][1].residues[i].codon);
  }
}

BOOST_AUTO_TEST_SUITE_END()
