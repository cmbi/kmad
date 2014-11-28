#ifndef TXTPROC_H
#define TXTPROC_H

#include "Sequences.h"

#include <iostream>
#include <string>
#include <vector>
class FeaturesProfile;
class Sequences;
namespace txtProc{
	double convertStringToDouble(std::string& s);
	std::string charToString(char mychar);
	std::string charToString(char mychar1, char mychar2);
	std::istream& safeGetline(std::istream& is, std::string& t);
	std::vector<std::string> split(const std::string& s, char delim);
	bool acceptedChar(char my_char);
	Sequences read_fasta(std::string filename, int codonLength, 
                       ids_list* ids, probs_list* probs);
	void writeAlignmentToFile(string_sequences& sequences, 
                            seqNames& sequence_names,
                            std::string filename);
	void writeAlignmentWithoutCodeToFile(string_sequences& sequences, 
                                       seqNames& sequence_names, 
                                       std::string filename, int codon_length);
	void process_conf_file(std::string filename, FeaturesProfile& feat_profile, 
                         Sequences& sequences_aa);
	featuresList unfold(std::string conf_string, 
                      featureNamesList& listOfFeatures);
}

#endif /* TXTPROC_H */
