#include "fasta.h"

#include <iostream>
#include <unordered_map>
#include <vector>

namespace feature_analysis {
  typedef std::unordered_map<std::string, std::vector<std::string>> CodesMap;
  
  struct ConsensusResidue {
    std::string amino_acid;
    std::string domain;
    std::string ptm;
    std::string motif;
    std::vector<std::string> user_defined_features;
  };
  typedef std::vector<ConsensusResidue> ConsensusSequence;

  CodesMap parse_mapfile(std::string filename);
  ConsensusSequence analyze_alignment(
      CodesMap codes_map,
      std::vector<fasta::SequenceList> alignment,
      double conservation_cutoff);
  void write_consensus_to_file(ConsensusSequence cons_seq,
                               std::string out_cons_filename);


}
