#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "features_profile.h"
#include "f_config.h"
#include "profile.h"
#include "scoring_matrix.h"
#include "types.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE(test_scoring_matrix)

BOOST_AUTO_TEST_CASE(test_backtrace_algorithm_path)
{
  int codon_length = 7;
  fasta::Sequence s1;
  // AKLCAKL
  s1 = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
                                 "CAAAAAAAAAAAAAKAAAAAA"
                                 "LAAAAAA", codon_length);
  fasta::Sequence s2;
  // AKLAKL
  s2 = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
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
  ProfileMap profile = create_score_profile(query_seq_list);
  int domain_modifier = 4;
  int motif_modifier = 3;
  int ptm_modifier = 10;
  FeaturesProfile f_profile(feature_list, domain_modifier,
                            ptm_modifier, motif_modifier,
                            probabilities);
  f_profile.update_scores(query_seq_list, f_set);
  double gap_open_pen = 5;
  double gap_ext_pen = 1;
  double end_pen = 1;
  ScoringMatrix scores(s1.residues.size(), s2.residues.size(), gap_open_pen,
                       end_pen, gap_ext_pen);
  scores.calculate_scores(s2, profile, f_profile, codon_length);
  fasta::SequenceList alignment;
  alignment = scores.backtrace_alignment_path(s2, 
                                              profile, f_profile,
                                              codon_length);
  fasta::Sequence e_s1;
  // AKLCAKL
  e_s1 = fasta::make_sequence("d", "AAAAAAAAAAAAAAAAAAAAA"
                                   "AAAAAAAAAAAAAAAAAAAAA"
                                   "AAAAAAA-AAAAAA", codon_length);
  fasta::Sequence e_s2;
  // AKLAKL
  e_s2 = fasta::make_sequence("d", "AAAAAAAKAAAAAALAAAAAA"
                                   "-AAAAAAAAAAAAAKAAAAAA"
                                   "LAAAAAARAAAAAA", codon_length);
  fasta::SequenceList expected_alignment = {e_s1, e_s2};
  BOOST_CHECK_EQUAL(alignment[0].residues.size(),
                    expected_alignment[0].residues.size());
  BOOST_CHECK_EQUAL(alignment[1].residues.size(),
                    expected_alignment[1].residues.size());
  for (size_t i = 0; i < alignment[0].residues.size(); ++i) {
    BOOST_CHECK_EQUAL(alignment[0].residues[i].codon,
                      expected_alignment[0].residues[i].codon);
  }
  for (size_t i = 0; i < alignment[1].residues.size(); ++i) {
    BOOST_CHECK_EQUAL(alignment[1].residues[i].codon,
                      expected_alignment[1].residues[i].codon);
  }


}

BOOST_AUTO_TEST_SUITE_END()
