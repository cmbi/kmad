#include "seq_data.h"

#include <algorithm>
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

  s.feature_list = make_feature_list(fasta_data.sequences);
  return s;
}


std::vector<std::string> seq_data::make_feature_list(
    const fasta::SequenceList& sequences) {
  std::vector<std::string> feature_list;
  for (auto& seq : sequences) {
    for (auto& res : seq.residues) {
      for (auto& feat_name : res.features) {
        if (std::find(feature_list.begin(), feature_list.end(), feat_name) 
            == feature_list.end()) {
          feature_list.push_back(feat_name);
        }
      }
    }
  }
  return feature_list;
}
