#ifndef SEQDATA_H
#define SEQDATA_H

#include "fasta.h"
#include "f_config.h"


namespace seq_data {
  struct SequenceData {
    fasta::SequenceList sequences;
    std::map<std::string, double> probabilities;  
    std::vector<std::string> feature_list;
  };

  SequenceData process_fasta_data(const fasta::FastaData& fasta_data,
                                  const f_config::FeatureSettingsMap& f_set,
                                  bool gapped);
  fasta::SequenceList remove_gaps(const fasta::SequenceList& sequences);
  FeatureNamesList make_feature_list(const fasta::SequenceList& sequences);
  void assign_feature_by_pattern(fasta::SequenceList& sequences,
                                 const std::string& pattern,
                                 const std::string& feat_name);
  int find_real_pos(const std::string& sequence, int position);

  bool compare_alignments(const std::vector<fasta::SequenceList>& al1,
      const std::vector<fasta::SequenceList>& al2);
}

#endif /* SEQDATA_H */
