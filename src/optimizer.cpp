#include "msa.h"
#include "optimizer.h"
#include "profile.h"

#include <boost/filesystem.hpp>
#include <algorithm>
#include <set>

namespace sbst = substitution_matrix;


std::vector<fasta::SequenceList> optimizer::optimize_alignment(
    const std::vector<fasta::SequenceList>& alignment,
    double domain_modifier, double motif_modifier, double ptm_modifier,
    const std::string& sbst_mat) {
  std::vector<optimizer::MoveData> m = optimizer::calculate_move_scores(
      alignment, domain_modifier, motif_modifier, ptm_modifier, sbst_mat);
  optimizer::filter_move_data(m);
  std::vector<fasta::SequenceList> new_alignment;
  new_alignment = optimizer::move_residues(alignment, m);
  new_alignment = msa::remove_gapcolumns(new_alignment);
  return new_alignment;
}


std::vector<optimizer::MoveData> optimizer::calculate_move_scores(
    const std::vector<fasta::SequenceList>& alignment,
    double domain_modifier, double motif_modifier, double ptm_modifier,
    const std::string& sbst_mat)
  {
  size_t alignment_length = alignment[0][0].residues.size();
  std::vector<optimizer::MoveData> move_data;
  std::string side;
  const sbst::SimilarityScoresMap* sim_scores;
  if (sbst_mat == "BLOSUM") {
    sim_scores = &sbst::BLOSUM;
  } else {
    sim_scores = &sbst::DISORDER;
  }
  optimizer::MoveData left;
  optimizer::MoveData right;
  for (size_t i  = 0; i < alignment[0].size(); ++i) {
    bool in_gap = false;
    for (size_t j = 0; j < alignment_length; ++j) {
      if (alignment[0][i].residues[j].codon[0] == '-' && !in_gap && j > 0
            && alignment[0][i].residues[j - 1].codon[0] != '-') {
          in_gap = true;
          side = "left";
          left = optimizer::single_move_score(alignment, i, j - 1, side,
                                              sim_scores);
      } else if (in_gap && alignment[0][i].residues[j].codon[0] != '-') {
          side = "right";
          right = optimizer::single_move_score(alignment, i, j, side,
                                               sim_scores);
          in_gap = false;
          optimizer::MoveData winner = (left.score_gain < right.score_gain) ? right : left;
          if (winner.score_gain > 0 ) {
            move_data.push_back(winner);
          }
      }
    }
  }
  return move_data;
}


bool optimizer::reverse_sort(int i, int j) {
  return j > i;
}


void optimizer::filter_move_data(
    std::vector<optimizer::MoveData>& move_data) {
  std::vector<int> del_indexes;
  for (size_t i = 0; i < move_data.size(); ++i) {
    for (size_t j = i; j < move_data.size(); ++j) {
      if (i != j && move_data[i].old_position == move_data[j].new_position) {
         if (move_data[i].score_gain > move_data[j].score_gain
             || (i > j
                 && move_data[i].score_gain == move_data[j].score_gain)) {
          del_indexes.push_back(j);
         } else {
          del_indexes.push_back(j);
         }
      }
    }
  }
  // remove duplicates
  std::set<int> s(del_indexes.begin(), del_indexes.end());
  del_indexes.assign(s.begin(), s.end());
  // sort in reverse order - for the sake of easy deletion
  std::sort(del_indexes.begin(), del_indexes.end(), optimizer::reverse_sort);

  for (auto& i : del_indexes) {
    move_data.erase(move_data.begin() + i);
  }
}


std::vector<fasta::SequenceList> optimizer::move_residues(
    const std::vector<fasta::SequenceList>& alignment,
    const std::vector<optimizer::MoveData>& move_data) {

  assert(alignment.size() == 2);
  assert(alignment[0].size() == alignment[1].size());

  std::vector<fasta::SequenceList> new_alignment = alignment;
  for (auto& i : move_data) {
    fasta::Residue tmp = new_alignment[0][i.seq_number].residues[i.new_position];
    new_alignment[0][i.seq_number].residues[i.new_position] = 
      new_alignment[0][i.seq_number].residues[i.old_position];
    new_alignment[0][i.seq_number].residues[i.old_position] = tmp;

    tmp = new_alignment[1][i.seq_number].residues[i.new_position];
    new_alignment[1][i.seq_number].residues[i.new_position] = 
      new_alignment[1][i.seq_number].residues[i.old_position];
    new_alignment[1][i.seq_number].residues[i.old_position] = tmp;
  }
  return new_alignment;
}


optimizer::MoveData optimizer::single_move_score(
    const std::vector<fasta::SequenceList>& alignment,
    size_t seq_number, int position, const std::string& side,
    const sbst::SimilarityScoresMap* sim_scores) {
  double pre_score = 0;
  char res1 = alignment[0][seq_number].residues[position].codon[0];
  int index = (std::find(sbst::ALPHABET.begin(), sbst::ALPHABET.end(), res1)
               - sbst::ALPHABET.begin());

  for (size_t i = 0; i < alignment[0].size(); ++i) {
    char res2 = alignment[0][i].residues[position].codon[0];
    if (i != seq_number && res2 != '-') {
      pre_score += sim_scores->at(res2)[index];
    }
  }

  double post_score = 0;
  int position2 = 0;
  if (side == "left") {
    // position2 - position of the last character of the gap
    position2 = optimizer::find_gap_end(alignment[0][seq_number],
                                        position + 1);
  } else {
    // position2 - position of the first character of the gap
    position2 = optimizer::find_gap_start(alignment[0][seq_number],
                                          position - 1);
  }
  if (position2 != -1
      && (unsigned)position2 != alignment[0][seq_number].residues.size()) {
    for (size_t i = 0; i < alignment[0].size(); ++i) {
      char res2 = alignment[0][i].residues[position2].codon[0];
      if (i != seq_number && res2 != '-') {
        post_score += sim_scores->at(res2)[index];
      }
    }
  }
  optimizer::MoveData m(seq_number, position, position2,
                        post_score - pre_score);
  return m;
}


int optimizer::find_gap_end(const fasta::Sequence& seq, int start) {
  int gap_end = seq.residues.size();
  bool not_found = true;
  for(size_t i = start; not_found && i < seq.residues.size(); i++) {
    if (seq.residues[i].codon[0] != '-') {
      not_found = false;  
      gap_end = i - 1;
    }
  }
  return gap_end;
}


int optimizer::find_gap_start(const fasta::Sequence& seq, int gap_end) {
  int gap_start = -1;
  bool not_found = true;
  for(size_t i = gap_end; not_found && i > 0; --i) {
    if (seq.residues[i].codon[0] != '-') {
      not_found = false;  
      gap_start = i + 1;
    }
  }
  return gap_start;
}
