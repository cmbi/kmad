#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "feature_scores.h"
#include "f_config.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE(test_feature_scores)

BOOST_AUTO_TEST_CASE(test_update_scores)
{
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
                                   "ptm_Oglyc3", "motif_aa", 
                                   "motif_ab", "motif_ac",
                                   "domain_aa", "domain_ac",
                                   "USR_feature1", "USR_feature2"};
  
  // SEQUENCE S1
  fasta::Residue r1_1("AAAAdaa", std::vector<std::string>({"ptm_phosphP",
                                                          "motif_aa"}));
  fasta::Residue r1_2("MAAAAAA", std::vector<std::string>({}));
  fasta::Residue r1_3("EAaadaa", std::vector<std::string>({"domain_aa"}));
  fasta::Residue r1_4("LAAAAAA", std::vector<std::string>({}));
  // SEQUENCE S2
  fasta::Residue r2_1("AAAAdaa", std::vector<std::string>({"ptm_phosphP",
                                                          "motif_aa"}));
  fasta::Residue r2_2("EAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Residue r2_3("EAaaZAA", std::vector<std::string>({"ptm_phosph0",
                                                           "domain_aa"}));
  fasta::Residue r2_4("KAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  // SEQUENCE S3
  fasta::Residue r3_1("AAAAZaa", std::vector<std::string>({"ptm_phosph0",
                                                           "motif_aa"}));
  fasta::Residue r3_2("MAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Residue r3_3("EAacZAA", std::vector<std::string>({"ptm_phosph0",
                                                           "domain_ac"}));
  fasta::Residue r3_4("LAaaZaa", std::vector<std::string>({"ptm_phosph0",
                                                           "motif_aa",
                                                           "domain_aa"}));
  // SEQUENCE S4
  fasta::Residue r4_1("MAAAZab", std::vector<std::string>({"ptm_phosph0",
                                                           "motif_ab"}));
  fasta::Residue r4_2("MAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Residue r4_3("KAaaZAA", std::vector<std::string>({"ptm_phosph0",
                                                           "domain_aa"}));
  fasta::Residue r4_4("LAAAAaa", std::vector<std::string>({"motif_aa"}));
  
  // SEQUENCE S5
  fasta::Residue r5_1("AAAAAac", std::vector<std::string>({"motif_ac"}));
  fasta::Residue r5_2("MAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Residue r5_3("AAAAAAA", std::vector<std::string>({"USR_feature2"}));
  fasta::Residue r5_4("LAAAAAA", std::vector<std::string>({"USR_feature1"}));

  fasta::Sequence s1 = fasta::make_sequence({r1_1, r1_2, r1_3, r1_4});
  fasta::Sequence s2 = fasta::make_sequence({r2_1, r2_2, r2_3, r2_4});
  fasta::Sequence s3 = fasta::make_sequence({r3_1, r3_2, r3_3, r3_4});
  fasta::Sequence s4 = fasta::make_sequence({r4_1, r4_2, r4_3, r4_4});
  fasta::Sequence s5 = fasta::make_sequence({r5_1, r5_2, r5_3, r5_4});

  fasta::SequenceList sequences = {s1, s2, s3, s4, s5};

  std::map<std::string, double> probs = {{"motif_aa", 1.0},
                                         {"motif_ab", 0.5},
                                         {"motif_ac", 0.8}};

  FeatureScores f_profile(feature_list, 4, 10, 3, probs);
  f_config::FeatureSettings settings1;
  settings1.add_score = 1;
  settings1.subtract_score = 0;
  settings1.add_features = {"USR_feature1"};
  f_config::FeatureSettings settings2;
  settings2.add_score = 1;
  settings2.subtract_score = 1;
  settings2.add_features = {"USR_feature2"};
  settings2.subtract_features = {"USR_feature1"};
  f_config::FeatureSettingsMap s_map; 
  s_map = {{"USR_feature1", settings1},
           {"USR_feature2", settings2}};
  f_profile.update_scores(sequences, s_map);
  std::map<std::string, Scores> p = f_profile.get_scores();
  std::map<std::string, Scores> expected_scores;
  expected_scores = {{"ptm_phosph0", {5.2, 8, 6, 4}},
                     {"ptm_phosph1", {4.68, 7.2, 5.4, 3.6}},
                     {"ptm_phosph2", {4.16, 6.4, 4.8, 3.2}},
                     {"ptm_phosph3", {3.64, 5.6, 4.2, 2.8}},
                     {"ptm_phosphP", {1.56, 2.4, 1.8, 1.2}},
                     {"ptm_acet0", {0, 0, 0, 0}},
                     {"ptm_acet1", {0, 0, 0, 0}},
                     {"ptm_acet2", {0, 0, 0, 0}},
                     {"ptm_acet3", {0, 0, 0, 0}},
                     {"ptm_Nglyc0", {0, 0, 0, 0}},
                     {"ptm_Nglyc1", {0, 0, 0, 0}},
                     {"ptm_Nglyc2", {0, 0, 0, 0}},
                     {"ptm_Nglyc3", {0, 0, 0, 0}},
                     {"ptm_amid0", {0, 0, 0, 0}},
                     {"ptm_amid1", {0, 0, 0, 0}},
                     {"ptm_amid2", {0, 0, 0, 0}},
                     {"ptm_amid3", {0, 0, 0, 0}},
                     {"ptm_hydroxy0", {0, 0, 0, 0}},
                     {"ptm_hydroxy1", {0, 0, 0, 0}},
                     {"ptm_hydroxy2", {0, 0, 0, 0}},
                     {"ptm_hydroxy3", {0, 0, 0, 0}},
                     {"ptm_methyl0", {0, 0, 0, 0}},
                     {"ptm_methyl1", {0, 0, 0, 0}},
                     {"ptm_methyl2", {0, 0, 0, 0}},
                     {"ptm_methyl3", {0, 0, 0, 0}},
                     {"ptm_Oglyc0", {0, 0, 0, 0}},
                     {"ptm_Oglyc1", {0, 0, 0, 0}},
                     {"ptm_Oglyc2", {0, 0, 0, 0}},
                     {"ptm_Oglyc3", {0, 0, 0, 0}},
                     {"motif_aa", {1.8, 0, 0, 1.2}},
                     {"motif_ab", {0.3, 0, 0, 0}},
                     {"motif_ac", {0.48, 0, 0, 0}},
                     {"domain_aa", {0, 0, 1.6, 0.8}},
                     {"domain_ac", {0, 0, -1.6, -0.8}},
                     {"USR_feature2", {0, 0,  0.2, -0.2}},
                     {"USR_feature1", {0, 0, 0, 0.2}}};
  for (auto& feat: feature_list){
    for (size_t i = 0; i < p[feat].size(); ++i) {
      BOOST_CHECK(std::abs(p[feat][i] - expected_scores[feat][i]) < 0.001);
    }
  }
  feature_list = {"ptm_phosph9", "motif_aa", "domain_aa"};
  f_profile = FeatureScores(feature_list, 4, 10, 3, probs);
  // SEQUENCE S1
  r1_1 = fasta::Residue ("AAAAdaa", std::vector<std::string>({"ptm_phosph9",
                                                              "motif_aa"}));
  r1_2 = fasta::Residue ("MAAAAAA", std::vector<std::string>({}));
  r1_3 = fasta::Residue ("EAaadaa", std::vector<std::string>({"domain_aa"}));
  r1_4 = fasta::Residue ("LAAAAAA", std::vector<std::string>({}));
  s1 = fasta::make_sequence({r1_1, r1_2, r1_3, r1_4});
  sequences = {s1};
  BOOST_CHECK_THROW(f_profile.update_scores(sequences, s_map),
                    std::invalid_argument);

}

BOOST_AUTO_TEST_SUITE_END()