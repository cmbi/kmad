#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main

#include "src/fasta.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

using namespace fasta;

BOOST_AUTO_TEST_SUITE(test_fasta)


BOOST_AUTO_TEST_CASE(test_parse_fasta)
{
    FastaData fd = parse_fasta("tests/TAU_SPECI.fasta.7c", 7);
    BOOST_CHECK_EQUAL(fd.sequences.size(), 19);
    BOOST_CHECK_EQUAL(
        sequence_to_string(fd.sequences[0]),
        "MAEPRQEFDTAEDHAEGYALLQDQEGEHGLKASPLQTPADDGPEEPVSETSDAKSTPTAEDVT"
        "APLVDERTPGEQAATQPPTDIPEGTTAEEAGIGDTPNMEDQAAGHVTQARMVSKGKEGTGSED"
        "RKAKGADSKTGTKIATPRGTAPPGQKGTANATRIPAKTTPSPKTPPGTGEPAKSGDRSGYSSP"
        "GSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSTKSRLQTAPVPMPDLKNVRSKIGST"
        "ENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLG"
        "NIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHG"
        "AEIVYKSPVVSGDTSPRHLSNVSSTGSINMVDSPQLATLADEVSASLAKQGL"
    );
    BOOST_CHECK_EQUAL(
        sequence_to_string(fd.sequences[1]),
        "MAEPRQEFDVMEDHAQGDYTLQDHEGDMEPGLKESPLQTPADDGSEEPGSETSDAKSTPTAEA"
        "EEAGIGDTSNLEDQAAGHVTQARMVSKGKDGTGPDDKKAKGADGKPGTKIATPRGAAPPGQKG"
        "QANATRIPAKTTPTPKTSPGTGESGKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVA"
        "VVRTPPKSPSAAKSRLQAAPGPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSK"
        "CGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSK"
        "IGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGS"
        "IDMVDSPQLATLADEVSASLAKQGL"
    );
    BOOST_CHECK_EQUAL(fd.probabilities.size(), 97);
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_nonexistent)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/nonexistent.fasta.7c", 7), std::invalid_argument
    );
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_invalid_probabilities_format)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/wrong_probs_format.fasta.7c", 7), std::runtime_error
    );
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_invalid_codons)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/invalid_codon.fasta.7c", 7), std::runtime_error
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

BOOST_AUTO_TEST_CASE(test_get_conf_data)
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


  //// SEQUENCE S1
  std::string sequence1_str = "AAAAdaaMAAAAAAEAaadaaLAAAAAA";
  std::string sequence2_str = "AAAAdaaEAAAZAAEAaaZAAKAAAZAA";
  fasta::Sequence s1 = fasta::make_sequence("d", sequence1_str, 7);
  fasta::Sequence s2 = fasta::make_sequence("d", sequence2_str, 7);
  fasta::SequenceList sequences = {s1, s2};
  fasta::FastaData test_data;
  test_data.sequences = sequences;

  std::string exp_sequence1_str = "AAAAdaaMAAAAAAEAaadaaLAAAAAA";
  std::string exp_sequence2_str = "AAAAdaaEAAAZAAEAaaZAAKAAAZAA";
  fasta::Sequence exp_s1 = fasta::make_sequence("d", exp_sequence1_str, 7);
  fasta::Sequence exp_s2 = fasta::make_sequence("d", exp_sequence2_str, 7);
  exp_s1.residues[1].features.push_back("USR_feature1");
  exp_s1.residues[2].features.push_back("USR_feature1");
  exp_s2.residues[0].features.push_back("USR_feature1");
  fasta::SequenceList expected_sequences = {exp_s1, exp_s2};

  FeatureNamesList expected_feature_list = {
      "p_phosph0", "p_phosph1", "p_phosph2", "p_phosph3", "p_phosphP",
      "p_acet0", "p_acet1", "p_acet2", "p_acet3", "p_Nglyc0", "p_Nglyc1",
      "p_Nglyc2", "p_Nglyc3", "p_amid0", "p_amid1", "p_amid2", "p_amid3",
      "p_hydroxy0", "p_hydroxy1", "p_hydroxy2", "p_hydroxy3", "p_methyl0",
      "p_methyl1", "p_methyl2", "p_methyl3", "p_Oglyc0", "p_Oglyc1",
      "p_Oglyc2", "p_Oglyc3", "p_cys_bridge0", "s_a_helix", "s_turn",
      "s_b_ladder", "s_b_bridge", "s_310_helix", "s_pi_helix", "s_b_ladder",
      "m_aa", "USR_feature1", "d_aa"};

  bool gapped = true;
  auto test_result = fasta::get_conf_data(test_data, test_map, gapped);

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
  fasta::Sequence s1 = fasta::make_sequence("", "ABAD-C", 1);
  fasta::Sequence s2 = fasta::make_sequence("", "ATSSAC", 1);
  fasta::Sequence s3 = fasta::make_sequence("", "--BAA-", 1);
  fasta::SequenceList s = {s1, s2 , s3};
  s = fasta::remove_gaps(s);
  std::vector<std::string> result;
  for (auto& seq : s) {
    result.push_back(fasta::sequence_to_string(seq));
  }
  std::vector<std::string> expected = {"ABADC", "ATSSAC", "BAA"};
  BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
                                expected.begin(), expected.end());
  s1 = fasta::make_sequence("", "ABADC", 1);
  s2 = fasta::make_sequence("", "ATSSAC", 1);
  s3 = fasta::make_sequence("", "BAA", 1);
  s = {s1, s2, s3};
  s = fasta::remove_gaps(s);
  result.clear();
  for (auto& seq : s) {
    result.push_back(fasta::sequence_to_string(seq));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
                                expected.begin(), expected.end());
}

BOOST_AUTO_TEST_CASE(test_assign_feature_by_pattern) {
  fasta::Sequence s1 = fasta::make_sequence("", "ABADB-C", 1);
  fasta::Sequence s2 = fasta::make_sequence("", "ATSS-AC", 1);
  fasta::Sequence s3 = fasta::make_sequence("", "--B-AA-", 1);
  fasta::SequenceList s = {s1, s2 , s3};
  fasta::SequenceList e = {s1, s2 , s3};
  std::string pattern = "BA";
  std::string feat_name = "testfeat";
  e[0].residues[1].features.push_back(feat_name);
  e[0].residues[2].features.push_back(feat_name);
  e[2].residues[2].features.push_back(feat_name);
  e[2].residues[4].features.push_back(feat_name);
  fasta::assign_feature_by_pattern(s, pattern, feat_name);
  for (size_t i = 0; i < s.size(); ++i) {
    for (size_t j = 0; j < s[i].residues.size(); ++j) {
      BOOST_CHECK_EQUAL_COLLECTIONS(
            s[i].residues[j].features.begin(), s[i].residues[j].features.end(),
            e[i].residues[j].features.begin(), e[i].residues[j].features.end());
    }
  }
}



BOOST_AUTO_TEST_SUITE_END()
