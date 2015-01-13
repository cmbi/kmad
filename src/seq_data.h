#ifndef SEQDATA_H
#define SEQDATA_H

#include "fasta.h"
#include "f_config.h"

#include <iostream>


namespace seq_data {
  struct SequenceData {
    fasta::SequenceList sequences;
    std::map<std::string, double> probabilities;  
  };

  SequenceData process_fasta_data(const fasta::FastaData& fasta_data,
                                  const f_config::FeatureSettingsMap& f_set);
}

#endif /* SEQDATA_H */
