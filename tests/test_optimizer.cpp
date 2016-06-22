#define BOOST_TEST_DYN_LINK

#include "src/fasta.h"
#include "src/optimizer.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

namespace sbst = substitution_matrix;
BOOST_AUTO_TEST_SUITE(test_optimizer)

BOOST_AUTO_TEST_CASE(test_find_gap_end)
{
  fasta::Sequence s = fasta::make_sequence("", "WF---F", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_end(s, 2), 4);

  s = fasta::make_sequence("", "WF---", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_end(s, 2), 5);

  s = fasta::make_sequence("", "---WF", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_end(s, 0), 2);

}


BOOST_AUTO_TEST_CASE(test_find_gap_start)
{
  fasta::Sequence s = fasta::make_sequence("", "WF---WF", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_start(s, 4), 2);

  s = fasta::make_sequence("", "---WF", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_start(s, 2), -1);

  s = fasta::make_sequence("", "WF---", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_start(s, 4), 2);
}

BOOST_AUTO_TEST_CASE(test_single_move_score)
{
  std::vector<fasta::SequenceList> alignment = {{
      fasta::make_sequence("", "AAE--AA", 1),
      fasta::make_sequence("", "AAE-EAA", 1),
      fasta::make_sequence("", "AA--EAA", 1)
  }};

  const sbst::SimilarityScoresMap* sim_scores = &sbst::BLOSUM;
  size_t seq_no = 0;
  int position = 2;
  std::string side = "left";
  double domain = 0;
  double ptm = 0;
  double motif = 0;

  optimizer::MoveData m = optimizer::single_move_score(alignment, seq_no,
                                                       position, side,
                                                       sim_scores, domain,
                                                       motif, ptm);
  BOOST_CHECK_EQUAL(m.score_gain, 5);
  seq_no = 2;
  position = 4;
  side = "right";
  m = optimizer::single_move_score(alignment, seq_no,
                                   position, side,
                                   sim_scores, domain, motif, ptm);

  BOOST_CHECK_EQUAL(m.score_gain, 5);
  alignment = {{fasta::make_sequence("", "AAE--", 1),
                fasta::make_sequence("", "AAE-E", 1),
                fasta::make_sequence("", "AA--E", 1)}};
  seq_no = 0;
  position = 2;
  side = "left";
  m = optimizer::single_move_score(alignment, seq_no,
                                   position, side,
                                   sim_scores, domain, motif, ptm);
  BOOST_CHECK_EQUAL(m.score_gain, -5);

  alignment = {{fasta::make_sequence("", "E--AA", 1),
                fasta::make_sequence("", "E-EAA", 1),
                fasta::make_sequence("", "--EAA", 1)}};
  seq_no = 2;
  position = 2;
  side = "right";
  m = optimizer::single_move_score(alignment, seq_no,
                                   position, side,
                                   sim_scores, domain, motif, ptm);
  BOOST_CHECK_EQUAL(m.score_gain, -105);
}


BOOST_AUTO_TEST_CASE(test_calculate_move_scores) 
{
  std::vector<fasta::SequenceList> alignment = {{
    fasta::make_sequence("", "AAE--AA", 1),
    fasta::make_sequence("", "AAE-EAA", 1),
    fasta::make_sequence("", "AA--EAA", 1)}};

  std::string sbst_mat = "BLOSUM";
  double domain = 0;
  double ptm = 0;
  double motif = 0;
  std::vector<optimizer::MoveData> m = optimizer::calculate_move_scores(
      alignment, domain, motif, ptm, sbst_mat);
  std::vector<optimizer::MoveData> expected = {
    optimizer::MoveData(0, 2, 4, 5),
    optimizer::MoveData(2, 4, 2, 5)};

  for (size_t i = 0; i < expected.size(); ++i) {
    BOOST_CHECK_EQUAL(expected[i].seq_number, m[i].seq_number);
    BOOST_CHECK_EQUAL(expected[i].old_position, m[i].old_position);
    BOOST_CHECK_EQUAL(expected[i].new_position, m[i].new_position);
    BOOST_CHECK_EQUAL(expected[i].score_gain, m[i].score_gain);
  }
}


BOOST_AUTO_TEST_CASE(test_filter_move_data)
{
  std::vector<optimizer::MoveData> move_data = {
    optimizer::MoveData(0, 2, 4, 5),
    optimizer::MoveData(2, 4, 2, 5)};

  optimizer::filter_move_data(move_data);
  std::vector<optimizer::MoveData> expected = {
    optimizer::MoveData(0, 2, 4, 5)};
  BOOST_CHECK_EQUAL(expected.size(), move_data.size());
  for (size_t i = 0; i < expected.size(); ++i) {
    BOOST_CHECK_EQUAL(expected[i].seq_number, move_data[i].seq_number);
    BOOST_CHECK_EQUAL(expected[i].old_position, move_data[i].old_position);
    BOOST_CHECK_EQUAL(expected[i].new_position, move_data[i].new_position);
    BOOST_CHECK_EQUAL(expected[i].score_gain, move_data[i].score_gain);
  }
}


BOOST_AUTO_TEST_CASE(test_remove_residues)
{
  std::vector<fasta::SequenceList> alignment = {
    {fasta::make_sequence("", "AAE--AA", 1),
     fasta::make_sequence("", "AAE-EAA", 1),
     fasta::make_sequence("", "AA--EAA", 1)},
    {fasta::make_sequence("", "AAE--AA", 1),
     fasta::make_sequence("", "AAE-EAA", 1),
     fasta::make_sequence("", "AA--EAA", 1)}};
  std::vector<optimizer::MoveData> move_data = {
    optimizer::MoveData(0, 2, 4, 5)};

  alignment = optimizer::remove_residues(alignment, move_data);

  std::vector<std::string> result;
  for (auto& item : alignment) {
    for (auto& seq : item) {
      result.push_back(fasta::sequence_to_string(seq));
    }
  }

  std::vector<std::string> expected = {"AA--EAA",
                                       "AAE-EAA",
                                       "AA--EAA",
                                       "AA--EAA",
                                       "AAE-EAA",
                                       "AA--EAA"};
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                result.begin(), result.end());
}


BOOST_AUTO_TEST_CASE(test_optimize_alignment)
{
  std::vector<fasta::SequenceList> alignment = {
    {fasta::make_sequence("", "AAE--AA", 1),
     fasta::make_sequence("", "AAE-EAA", 1),
     fasta::make_sequence("", "AA--EAA", 1)},
    {fasta::make_sequence("", "AAE--AA", 1),
     fasta::make_sequence("", "AAE-EAA", 1),
     fasta::make_sequence("", "AA--EAA", 1)}};


  std::string sbst_mat = "BLOSUM";
  double domain = 0;
  double motif = 0;
  double ptm = 0;
  alignment = optimizer::optimize_alignment(alignment, domain, motif, ptm,
                                            sbst_mat);

  std::vector<std::string> result;
  for (auto& item : alignment) {
    for (auto& seq : item) {
      result.push_back(fasta::sequence_to_string(seq));
    }
  }

  std::vector<std::string> expected = {"AA-EAA",
                                       "AAEEAA",
                                       "AA-EAA",
                                       "AA-EAA",
                                       "AAEEAA",
                                       "AA-EAA"};
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                result.begin(), result.end());
}

BOOST_AUTO_TEST_CASE(test_score_ptm) {
  fasta::Residue res1("A", {"p_phosph0"});
  fasta::Residue res2("A", {"p_acet0"});
  double p_modifier = 10;
  BOOST_CHECK_EQUAL(optimizer::score_ptm(res1, res2, p_modifier), 0);
  res1 = fasta::Residue("A", {"p_phosph0"});
  res2 = fasta::Residue("A", {"p_phosph2"});
  BOOST_CHECK_EQUAL(optimizer::score_ptm(res1, res2, p_modifier), 8);
  res1 = fasta::Residue("A", {"p_phosph0"});
  res2 = fasta::Residue("A", {});
  BOOST_CHECK_EQUAL(optimizer::score_ptm(res1, res2, p_modifier), 0);
}

BOOST_AUTO_TEST_CASE(test_score_domain) {
  fasta::Residue res1("A", {"d_aa"});
  fasta::Residue res2("A", {"d_aa"});
  double d_modifier = 10;
  BOOST_CHECK_EQUAL(optimizer::score_domain(res1, res2, d_modifier),
                    d_modifier);
  res1 = fasta::Residue("A", {"d_aa"});
  res2 = fasta::Residue("A", {"d_ab"});
  BOOST_CHECK_EQUAL(optimizer::score_domain(res1, res2, d_modifier),
                    - d_modifier);
  res1 = fasta::Residue("A", {"d_aa"});
  res2 = fasta::Residue("A", {});
  BOOST_CHECK_EQUAL(optimizer::score_domain(res1, res2, d_modifier), 0);
}

BOOST_AUTO_TEST_CASE(test_score_motif) {
  fasta::Residue res1("A", {"m_aa"});
  fasta::Residue res2("A", {"m_aa"});
  double m_modifier = 10;
  BOOST_CHECK_EQUAL(optimizer::score_motif(res1, res2, m_modifier),
                    m_modifier);
  res1 = fasta::Residue("A", {"m_aa"});
  res2 = fasta::Residue("A", {"m_ab"});
  BOOST_CHECK_EQUAL(optimizer::score_motif(res1, res2, m_modifier), 0);
  res1 = fasta::Residue("A", {"m_aa"});
  res2 = fasta::Residue("A", {});
  BOOST_CHECK_EQUAL(optimizer::score_motif(res1, res2, m_modifier), 0);
}


BOOST_AUTO_TEST_SUITE_END()
