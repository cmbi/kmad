#include "seq_data.h"

#include <iostream>


seq_data::SequenceData seq_data::process_fasta_data(
    const fasta::FastaData& fasta_data, 
    const f_config::FeatureSettingsMap& f_set) {
  seq_data::SequenceData s;
  s.probabilities = fasta_data.probabilities;
  s.sequences = fasta_data.sequences;

  for (auto feat_it = f_set.begin(); feat_it != f_set.end(); feat_it++) {
    for (auto& seq : feat_it->second.positions) {
      for (auto& pos : seq.positions) {
        s.sequences[seq.seq_no].residues[pos].features.push_back(
            feat_it->first);
      }
    }
  }
  return s;
}
