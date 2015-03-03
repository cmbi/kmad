#include "seq_data.h"

#include <algorithm>


seq_data::SequenceData seq_data::process_fasta_data(
    const fasta::FastaData& fasta_data, 
    const f_config::FeatureSettingsMap& f_set) {
  seq_data::SequenceData s;
  s.probabilities = fasta_data.probabilities;
  s.sequences = fasta_data.sequences;

  for (auto feat_it = f_set.begin(); feat_it != f_set.end(); ++feat_it) {
    for (auto& seq : feat_it->second.positions) {
      for (auto& pos : seq.positions) {
        s.sequences[seq.seq_no].residues[pos].features.push_back(
            feat_it->first);
      }
    }
  }
  s.feature_list = make_feature_list(s.sequences);
  return s;
}


FeatureNamesList seq_data::make_feature_list(
    const fasta::SequenceList& sequences) {
  FeatureNamesList feature_list = {"ptm_phosph0", "ptm_phosph1",
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
                                   "ptm_Oglyc3"}; 
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
