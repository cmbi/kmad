#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "fasta.h"
#include "substitution_matrix.h"
#include <vector>

namespace sbst = substitution_matrix;

namespace optimizer {

  struct MoveData {
   MoveData(int seq_number, int old_position, int new_position,
       double score_gain) 
      : seq_number(seq_number),
        old_position(old_position),
        new_position(new_position),
        score_gain(score_gain) {}
   MoveData() {}
   int seq_number; 
   int old_position;
   int new_position;
   double score_gain;
  };

  ///
  /// takes an alignment and parameters; returns an optimized alignment
  ///
  std::vector<fasta::SequenceList> optimize_alignment(
      const std::vector<fasta::SequenceList>& alignment,
      double domain_modifier, double motif_modifier, double ptm_modifier,
      const std::string& sbst_mat);

  ///
  /// finds residues to move; returns their positions and scores for removing
  /// {{ sequence_number, old_position, new_position, score_gain }}
  ///
  std::vector<MoveData> calculate_move_scores(
      const std::vector<fasta::SequenceList>& alignment,
      double domain_modifier, double motif_modifier, double ptm_modifier,
      const std::string& sbst_mat);

  ///
  /// takes an alignment and residues to remove; returns an alignment with
  /// removed residues
  ///
  std::vector<fasta::SequenceList> remove_residues(
      const std::vector<fasta::SequenceList>& alignment,
      const std::vector<MoveData>& move_data);

  ///
  /// 
  ///
  void filter_move_data(std::vector<MoveData>& move_data);
  bool reverse_sort(int i, int j);


  MoveData single_move_score(
      const std::vector<fasta::SequenceList>& alignment,
      size_t seq_number, int position, const std::string& side,
      const sbst::SimilarityScoresMap* sim_scores, double domain_modifier,
      double motif_modifier, double ptm_modifier);

  int find_gap_end(const fasta::Sequence& seq, int start);
  int find_gap_start(const fasta::Sequence& seq, int start);


  double get_two_res_score(fasta::Residue res1, fasta::Residue res2,
                           int res1_index,
                           const sbst::SimilarityScoresMap* sim_scores,
                           double domain_modifier,
                           double motif_modifier, 
                           double ptm_modifier);
  
  double score_ptm(fasta::Residue res1, fasta::Residue res2,
                   double ptm_modifier);

  double score_motif(fasta::Residue res1, fasta::Residue res2,
                   double motif_modifier);

  double score_domain(fasta::Residue res1, fasta::Residue res2,
                      double domain_modifier);
}

#endif /* OPTIMIZER_H */
