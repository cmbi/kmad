#include "outfile.h"

#include <iostream>
#include <fstream>
#include <sstream>


void outfile::write_encoded_alignment(
    const fasta::SequenceList& sequences,
    const seq_data::SequenceData& sequence_data,
    const std::string& filename_prefix) {
  std::stringstream sstr;
  sstr << filename_prefix << "_al";
  std::ofstream output_file(sstr.str().c_str(), std::ios::out);
  for (size_t i = 0; i < sequences.size(); ++i) {
    output_file << sequence_data.sequences[i].description << "\n";
    for (auto& res : sequences[i].residues) {
      output_file << res.codon;
    }
    output_file << "\n";
  }
}


void outfile::write_decoded_alignment(
    const fasta::SequenceList& sequences,
    const seq_data::SequenceData& sequence_data,
    const std::string& filename_prefix) {
  std::stringstream sstr;
  sstr << filename_prefix << "_al";
  std::ofstream output_file(sstr.str().c_str(), std::ios::out);
  for (size_t i = 0; i < sequences.size(); ++i) {
    output_file << sequence_data.sequences[i].description << "\n";
    for (auto& res : sequences[i].residues) {
      output_file << res.codon[0];
    }
    output_file << "\n";
  }
}
