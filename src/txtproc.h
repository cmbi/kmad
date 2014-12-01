#ifndef TXTPROC_H
#define TXTPROC_H

#include "sequences.h"

#include <iostream>
#include <string>
#include <vector>


class FeaturesProfile;
class Sequences;
namespace txtProc{
  ///
  /// converts a string to double
  ///
	double convertStringToDouble(std::string& s);
  ///
  /// converts a char to a string
  ///
	std::string charToString(char mychar);
  ///
  /// conactenates two chars to a string
  ///
	std::string charToString(char mychar1, char mychar2);
  ///
  /// writes a line of text from the istream to string t 
  /// (removes newline characters)
  ///
	std::istream& safeGetline(std::istream& is, std::string& t);
  ///
  /// splits a string by the delimiter, returns a vector of strings
  ///
	std::vector<std::string> split(const std::string& s, char delim);
  ///
  /// checks if character from the fasta file is on the list of accepetd
  /// characters
  ///
	bool acceptedChar(char my_char);
  ///
  /// reads a 'fasta' file, returns a Sequences object, 
  /// writes motif ids and probabilities to ids and probs
  ///
	Sequences read_fasta(std::string filename, int codonLength, 
                       IDsList* ids, ProbsList* probs);
  ///
  /// writes alignment to file (encoded)
  ///
	void writeAlignmentToFile(StringSequences& sequences, 
                            SeqNames& sequence_names,
                            std::string filename);
  ///
  /// writes alignment to a file in the regular fasta format (decoded)
  ///
	void writeAlignmentWithoutCodeToFile(StringSequences& sequences, 
                                       SeqNames& sequence_names, 
                                       std::string filename, int codon_length);
  ///
  /// reades configuration file, adds user defined features and their alignment
  /// rules
  ///
	void process_conf_file(std::string filename, FeaturesProfile& feat_profile, 
                         Sequences& sequences_aa);
  ///
  /// from a conf_file string creates a list of indexes of features to be
  /// scored
  ///
	FeaturesList unfold(std::string conf_string, 
                      FeatureNamesList& listOfFeatures);
}

#endif /* TXTPROC_H */
