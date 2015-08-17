#define BOOST_TEST_DYN_LINK

#include "src/fasta.h"
#include "src/feature_scores.h"
#include "src/f_config.h"
#include "src/profile.h"
#include "src/scoring_matrix.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE(test_scoring_matrix)

BOOST_AUTO_TEST_CASE(test_calculate_scores) {
  fasta::Sequence s1;
  fasta::Sequence s2;
  int codon_length = 1;
  s1 = fasta::make_sequence("d", "ASLKSLKPT", codon_length);
  s2 = fasta::make_sequence("d", "ASLRP", codon_length);
  fasta::SequenceList query_seq_list = {s1};
  fasta::SequenceList sequences = {s1, s2};
  std::string sbst_mat = "BLOSUM";
  profile::ProfileMap profile = profile::create_score_profile(query_seq_list,
                                                              sbst_mat);
  int domain_modifier = 4;
  int motif_modifier = 3;
  int ptm_modifier = 10;
  std::map<std::string, double> probabilities;
  f_config::FeatureSettingsMap f_set;
  std::vector<std::string> feature_list;
  FeatureScores f_profile(feature_list, domain_modifier,
                          ptm_modifier, motif_modifier,
                          probabilities);
  std::vector<double> identities(query_seq_list.size(), 1.0);
  bool fade_out = false;
  f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
  double gap_open_pen = -5;
  double gap_ext_pen = -1;
  double end_pen = -1;
  bool no_feat = false;
  ScoringMatrix scores(s1.residues.size(), s2.residues.size(), gap_open_pen,
                       end_pen, gap_ext_pen, no_feat);
  scores.calculate_scores(s2, profile, f_profile, codon_length);
  SingleScoringMatrix matrix_V = scores.get_V_matrix();
  SingleScoringMatrix expected_V = {{0, -10000000, -10000000, -10000000,
                                     -10000000, -10000000},
                                    {-10000000,  4,  0, -3, -4, -5},
                                    {-10000000,  0,  8, -2, -3, -4},
                                    {-10000000, -3, -2, 12,  1, -1},
                                    {-10000000, -4, -2,  1, 14,  6},
                                    {-10000000, -3,  1,  0,  6, 13},
                                    {-10000000, -6, -5,  5,  4,  6},
                                    {-10000000, -7, -5, -2,  7,  7},
                                    {-10000000, -8, -7, -4,  2, 14},
                                    {-10000000, -8, -6, -3,  2,  5}};

  for (size_t i = 0; i < matrix_V.size(); ++i) {
    BOOST_CHECK_EQUAL_COLLECTIONS(matrix_V[i].begin(), matrix_V[i].end(),
                                  expected_V[i].begin(), expected_V[i].end());
  }

  // AND THE OTHER WAY ROUND (s2 becomes the profile)
  query_seq_list = {s2};
  profile = profile::create_score_profile(query_seq_list, sbst_mat);
  f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
  ScoringMatrix scores2(s2.residues.size(), s1.residues.size(), gap_open_pen,
                        end_pen, gap_ext_pen, no_feat);
  scores2.calculate_scores(s1, profile, f_profile, codon_length);
  matrix_V = scores2.get_V_matrix();

  expected_V = {{0, -10000000, -10000000, -10000000, 
                 -10000000, -10000000, -10000000,
                 -10000000, -10000000, -10000000},
                {-10000000,  4,  0, -3, -4, -3, -6, -7, -8, -8},
                {-10000000,  0,  8, -2, -2,  1, -5, -5, -7, -6}, 
                {-10000000, -3, -2, 12,  1,  0,  5, -2, -4, -3},
                {-10000000, -4, -3,  1, 14,  6,  4,  7,  2,  2},
                {-10000000, -5, -4, -1,  6, 13,  6,  7, 14,  5}};
  for (size_t i = 0; i < matrix_V.size(); ++i) {
    BOOST_CHECK_EQUAL_COLLECTIONS(matrix_V[i].begin(), matrix_V[i].end(),
                                  expected_V[i].begin(), expected_V[i].end());
  }
}
BOOST_AUTO_TEST_CASE(test_backtrace_alignment_path) {
  fasta::Sequence s1;
  fasta::Sequence s2;
  int codon_length = 1;
  s1 = fasta::make_sequence("d", "ASLKSLKPT", codon_length);
  s2 = fasta::make_sequence("d", "ASLRP", codon_length);
  fasta::SequenceList query_seq_list = {s1};
  fasta::SequenceList sequences = {s1, s2};
  int domain_modifier = 4;
  int motif_modifier = 3;
  int ptm_modifier = 10;
  std::map<std::string, double> probabilities;
  f_config::FeatureSettingsMap f_set;
  std::vector<std::string> feature_list;
  FeatureScores f_profile(feature_list, domain_modifier,
                          ptm_modifier, motif_modifier,
                          probabilities);
  std::vector<double> identities(query_seq_list.size(), 1.0);
  bool fade_out = false;
  f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
  double gap_open_pen = -5;
  double gap_ext_pen = -1;
  double end_pen = -1;
  const bool no_feat = false;
  fasta::Sequence e_s1;
  fasta::Sequence e_s2;
  e_s1 = fasta::make_sequence("d", "AAAAAAAAA", codon_length);
  e_s2 = fasta::make_sequence("d", "A---SLRP-", codon_length);
  std::string sbst_mat = "BLOSUM";
  profile::ProfileMap profile = profile::create_score_profile(query_seq_list,
                                                              sbst_mat);
  ScoringMatrix scores(s1.residues.size(), s2.residues.size(), gap_open_pen,
                       end_pen, gap_ext_pen, no_feat);
  scores.calculate_scores(s2, profile, f_profile, codon_length);
  fasta::SequenceList result = scores.backtrace_alignment_path(s2, profile,
                                                               f_profile,
                                                               codon_length);
  BOOST_CHECK_EQUAL(e_s1.residues.size(), result[0].residues.size());
  BOOST_CHECK_EQUAL(e_s2.residues.size(), result[1].residues.size());
  for (size_t i = 0; i < e_s1.residues.size(); ++i) {
    BOOST_CHECK_EQUAL(e_s1.residues[i].codon, result[0].residues[i].codon);
    BOOST_CHECK_EQUAL(e_s2.residues[i].codon, result[1].residues[i].codon);
  }
  query_seq_list = {s2};
  f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
  e_s1 = fasta::make_sequence("d", "A---AAAA-", codon_length);
  e_s2 = fasta::make_sequence("d", "ASLKSLKPT", codon_length);
  profile = profile::create_score_profile(query_seq_list, sbst_mat);
  ScoringMatrix scores2(s2.residues.size(), s1.residues.size(), gap_open_pen,
                        end_pen, gap_ext_pen, no_feat);
  scores2.calculate_scores(s1, profile, f_profile, codon_length);
  result = scores2.backtrace_alignment_path(s1, profile,
                                            f_profile,
                                            codon_length);
  BOOST_CHECK_EQUAL(e_s1.residues.size(), result[0].residues.size());
  BOOST_CHECK_EQUAL(e_s2.residues.size(), result[1].residues.size());
  for (size_t i = 0; i < e_s1.residues.size(); ++i) {
    BOOST_CHECK_EQUAL(e_s1.residues[i].codon, result[0].residues[i].codon);
    BOOST_CHECK_EQUAL(e_s2.residues[i].codon, result[1].residues[i].codon);
  }
  gap_ext_pen = -0.1;
  end_pen = -0.01;
  gap_open_pen = -1;
  s1 = fasta::make_sequence("d", "RRRDDRR", codon_length);
  s2 = fasta::make_sequence("d", "RRRWWRR", codon_length);
  query_seq_list = {s1};
  f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
  profile = profile::create_score_profile(query_seq_list, sbst_mat);
  ScoringMatrix scores3(s1.residues.size(), s2.residues.size(), gap_open_pen,
                        end_pen, gap_ext_pen, no_feat);
  scores3.calculate_scores(s2, profile, f_profile, codon_length);
  result = scores3.backtrace_alignment_path(s2, profile,
                                            f_profile,
                                            codon_length);
  e_s1 = fasta::make_sequence("d", "AAAAA--AA", codon_length);
  e_s2 = fasta::make_sequence("d", "RRR--WWRR", codon_length);
  for (size_t i = 0; i < e_s1.residues.size(); ++i) {
    BOOST_CHECK_EQUAL(e_s1.residues[i].codon, result[0].residues[i].codon);
    BOOST_CHECK_EQUAL(e_s2.residues[i].codon, result[1].residues[i].codon);
  }
  query_seq_list = {s2};
  f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
  profile = profile::create_score_profile(query_seq_list, sbst_mat);
  scores3 = ScoringMatrix(s2.residues.size(), s1.residues.size(), gap_open_pen,
                         end_pen, gap_ext_pen, no_feat);
  scores3.calculate_scores(s1, profile, f_profile, codon_length);
  result = scores3.backtrace_alignment_path(s1, profile,
                                            f_profile,
                                            codon_length);
  e_s1 = fasta::make_sequence("d", "AAAAA--AA", codon_length);
  e_s2 = fasta::make_sequence("d", "RRR--DDRR", codon_length);
  for (size_t i = 0; i < e_s1.residues.size(); ++i) {
    BOOST_CHECK_EQUAL(e_s1.residues[i].codon, result[0].residues[i].codon);
    BOOST_CHECK_EQUAL(e_s2.residues[i].codon, result[1].residues[i].codon);
  }
  codon_length = 7;
  feature_list = {"ptm_phosph0", "ptm_phosph1",
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
                  "ptm_Oglyc3", "motif_aa"};
  s1 = fasta::make_sequence("d", "TAAAZAATAAAZAARAAAAAARAAAAAARAAAAAADAAAAAA"
                                 "DAAAAAARAAAAAARAAAAAA", codon_length);
  s2 = fasta::make_sequence("d", "RAAAZAARAAAAaaRAAAAAAWAAAAAAWAAAAAARAAAAAA"
                                 "RAAAAAA", codon_length);
  query_seq_list = {s1};
  f_profile = FeatureScores(feature_list, domain_modifier,
                            ptm_modifier, motif_modifier,
                            probabilities);
  f_profile.update_scores(query_seq_list, f_set, identities, fade_out);
  profile = profile::create_score_profile(query_seq_list, sbst_mat);
  scores3 = ScoringMatrix(s2.residues.size(), s1.residues.size(), gap_open_pen,
                         end_pen, gap_ext_pen, no_feat);
  scores3.calculate_scores(s2, profile, f_profile, codon_length);
  result = scores3.backtrace_alignment_path(s2, profile,
                                            f_profile,
                                            codon_length);
  e_s1 = fasta::make_sequence("d", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA-AAAAAA"
                                   "-AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                                   "AAAAAAA", codon_length);
  e_s2 = fasta::make_sequence("d", "-AAAAAARAAAZAARAAAAaaRAAAAAAWAAAAAA"
                                   "WAAAAAARAAAAAA-AAAAAA-AAAAAARAAAAAA"
                                   "-AAAAAA", codon_length);
  for (size_t i = 0; i < e_s1.residues.size(); ++i) {
    BOOST_CHECK_EQUAL(e_s1.residues[i].codon, result[0].residues[i].codon);
    BOOST_CHECK_EQUAL(e_s2.residues[i].codon, result[1].residues[i].codon);
  }

}

BOOST_AUTO_TEST_SUITE_END()
