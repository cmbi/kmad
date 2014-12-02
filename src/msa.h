#ifndef MSA_H
#define MSA_H

#include "types.h"

#include <iostream>
#include <vector>


class Sequences;
namespace msa{
  ///
  /// performs the full multiple sequence alignment, returns aligned sequences
  ///
  StringSequences run_msa(Sequences sequences, std::string conf_filename, 
                           double gap_open_pen, double gap_ext_pen, 
                           double end_pen, int domain_score, 
                           int motif_score, int phosph_score, int codon_length,
                           IDsList motifs_ids, ProbsList motifs_probs);
}

#endif /* MSA_H */
