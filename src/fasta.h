#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <unordered_map>
#include <vector>

/// 
///  Parses a fasta or semi-fasta(encoded) file to a vec of Sequence structs, 
/// 
namespace fasta {

  ///
  /// Holds a single residue codon, and a list of its features
  ///
  struct Residue {
    Residue(std::string codon, std::vector<std::string> features) 
      : codon(codon),
        features(features) {}
    Residue(std::string codon) : codon(codon) {}
    Residue(){};

    std::string codon;
    std::vector<std::string> features;
  };
  ///
  /// hodls description which is later used as a fasta header and a list of
  /// residues
  ///
  struct Sequence {
    std::string description;
    std::vector<Residue> residues;
  };

  typedef std::vector<fasta::Sequence> SequenceList;

  ///
  /// holds all the data from the (semi)fasta file - sequences with their headers
  /// and a list of motif probabilities
  ///
  struct FastaData {
    SequenceList sequences;
    std::unordered_map<std::string, double> probabilities;
  };

  ///
  /// parses a 'fasta' file, returns a Sequences object,
  /// writes motif ids and probabilities to ids and probs
  ///
  FastaData parse_fasta(std::string filename, int codonLength,
                        bool refine, int refine_seq);

  Sequence make_sequence(const std::string& description,
                         const std::string& codons, int codon_length);
  Sequence make_sequence(const std::vector<Residue>&);
  Sequence make_sequence(unsigned long sequence_length,
                         const fasta::Residue residue);
  Sequence make_sequence(const std::string sequence_string, int codon_length);
  void extend_sequence(Sequence& seq, const std::string sequence_string,
                           int codon_length);
  ///
  /// converts a Sequence struct to a string
  ///
  std::string make_string(const Sequence seq);
  Residue make_residue(const std::string& codon);
  ///
  /// check if all sequence lengths are equal (for the refinement mode)
  ///
  bool check_length(SequenceList sequences, int refine_seq);
}


#endif /* FASTA_H */
