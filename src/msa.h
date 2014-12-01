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
  string_sequences run_msa(Sequences sequences, std::string conf_filename, 
                           double gapPen, double gapExt, 
                           double endPenalty, double lcr_mod, int domainScore, 
                           int motifScore, int phosphScore, int codonLength,
                           ids_list motifs_ids, probs_list motifs_probs);
}

#endif /* MSA_H */
