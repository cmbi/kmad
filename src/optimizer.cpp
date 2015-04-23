#include "optimizer.h"
#include "profile.h"

#include <algorithm>

namespace sbst = substitution_matrix;

// TODO: implement
std::vector<fasta::SequenceList> optimizer::optimize_alignment(
    const std::vector<fasta::SequenceList>& alignment,
    double domain_modifier, double motif_modifier, double ptm_modifier,
    const std::string& sbst_mat) {
  std::vector<fasta::SequenceList> new_alignment;
  return new_alignment;
}

// TODO: implement
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
      if (alignment[0][i].residues[j].codon[0] == '-' && !in_gap) {
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

// TODO: implement
std::vector<fasta::SequenceList> optimizer::remove_residues(
    const std::vector<fasta::SequenceList>& alignment,
    const std::vector<optimizer::MoveData>& move_data) {
  std::vector<fasta::SequenceList> new_alignment;
  return new_alignment;
}

// TODO: test it
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
                                            position);
  } else {
    // position2 - position of the first character of the gap
    position2 = optimizer::find_gap_start(alignment[0][seq_number],
                                              position);
  }
  for (size_t i = 0; i < alignment[0].size(); ++i) {
    char res2 = alignment[0][i].residues[position2].codon[0];
    if (i != seq_number && res2 != '-') {
      post_score += sim_scores->at(res2)[index];
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
