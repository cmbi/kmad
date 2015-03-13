#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "f_config.h"
#include "seq_data.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <string>


BOOST_AUTO_TEST_SUITE(test_kman_suite)

BOOST_AUTO_TEST_CASE(test_seq_data)
{

  f_config::FeatureSettingsMap test_map;
  f_config::FeatureSettings test_settings;
  f_config::FeaturePositions positions1;
  positions1.seq_no = 0;
  positions1.positions = {1, 2};
  f_config::FeaturePositions positions2;
  positions2.seq_no = 1;
  positions2.positions = {0};
  test_settings.positions = {positions1, positions2};
  test_map["USR_feature1"] = test_settings;


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

  fasta::Sequence s1 = fasta::make_sequence({r1_1, r1_2, r1_3, r1_4});
  fasta::Sequence s2 = fasta::make_sequence({r2_1, r2_2, r2_3, r2_4});

  fasta::SequenceList sequences = {s1, s2};
  fasta::FastaData test_data;
  test_data.sequences = sequences;

  // SEQUENCE S1
  fasta::Residue e1_1("AAAAdaa", std::vector<std::string>({"ptm_phosphP",
                                                          "motif_aa"}));
  fasta::Residue e1_2("MAAAAAA", std::vector<std::string>({"USR_feature1"}));
  fasta::Residue e1_3("EAaadaa", std::vector<std::string>({"domain_aa",
                                                           "USR_feature1"}));
  fasta::Residue e1_4("LAAAAAA", std::vector<std::string>({}));
  // SEQUENCE S2
  fasta::Residue e2_1("AAAAdaa", std::vector<std::string>({"ptm_phosphP",
                                                          "motif_aa",
                                                          "USR_feature1"}));
  fasta::Residue e2_2("EAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Residue e2_3("EAaaZAA", std::vector<std::string>({"ptm_phosph0",
                                                           "domain_aa"}));
  fasta::Residue e2_4("KAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Sequence exp_s1 = fasta::make_sequence({e1_1, e1_2, e1_3, e1_4});
  fasta::Sequence exp_s2 = fasta::make_sequence({e2_1, e2_2, e2_3, e2_4});
  fasta::SequenceList expected_sequences = {exp_s1, exp_s2};

  FeatureNamesList expected_feature_list = {"ptm_phosph0", "ptm_phosph1",
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
                                            "USR_feature1", "domain_aa"}; 

  bool gapped = true;
  seq_data::SequenceData test_result = seq_data::process_fasta_data(test_data,
                                                                    test_map,
                                                                    gapped);

  BOOST_CHECK_EQUAL_COLLECTIONS(expected_feature_list.begin(), 
                                expected_feature_list.end(),
                                test_result.feature_list.begin(),
                                test_result.feature_list.end());
  
  for (size_t i = 0; i < test_result.sequences.size(); ++i) {
    for (size_t j = 0; j < test_result.sequences[i].residues.size(); ++j) {
      BOOST_CHECK_EQUAL_COLLECTIONS(
          expected_sequences[i].residues[j].features.begin(),
          expected_sequences[i].residues[j].features.end(),
          test_result.sequences[i].residues[j].features.begin(),
          test_result.sequences[i].residues[j].features.end());
    }
  }
}


BOOST_AUTO_TEST_CASE(test_remove_gaps) {
  fasta::Sequence s1 = fasta::make_sequence("ABAD-C", 1);
  fasta::Sequence s2 = fasta::make_sequence("ATSSAC", 1);
  fasta::Sequence s3 = fasta::make_sequence("--BAA-", 1);
  fasta::SequenceList s = {s1, s2 , s3};
  s = seq_data::remove_gaps(s);
  std::vector<std::string> result;
  for (auto& seq : s) {
    result.push_back(fasta::make_string(seq));
  }
  expected = {"ABADC", "ATSSAC", "BAA"};
  BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
                                expected.begin(), expected.end());
  s1 = fasta::make_sequence("ABADC", 1);
  s2 = fasta::make_sequence("ATSSAC", 1);
  s3 = fasta::make_sequence("BAA", 1);
  s = {s1, s2, s3};
  s = seq_data::remove_gaps(s);
  result.clear()
  for (auto& seq : s) {
    result.push_back(fasta::make_string(seq));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
                                expected.begin(), expected.end());
}


BOOST_AUTO_TEST_SUITE_END()
