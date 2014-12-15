#include "fasta.h"
#include "txtproc.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <iostream>
#include <iterator>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>


namespace fs = boost::filesystem;


typedef std::vector<std::string> InputLine;

Sequences fasta::parse_fasta(std::string filename,
                             int codonLength,
                             IDsList* ids,
                             ProbsList* probs) {
  CodonSeqWithNamesList resultSequences;
  fs::path p(filename);
  if (!fs::exists(p)) {
    throw std::invalid_argument("File not found: " + filename);
  } else {
    std::string fastaSymbol = ">";
    std::ifstream fastafile (filename.c_str());
    SeqNames newName;
    CodonSeq newSequence;
    CodonSeqWithName newEntry;
    bool sequences = true;
    int seqNo = -1;

    std::string line;
    while (std::getline(fastafile, line)) {
        std::string firstChar = line.substr(0, 1);
        if (line != std::string("## PROBABILITIES")) {
          if (sequences && firstChar == fastaSymbol) {
            seqNo++;
            resultSequences.push_back(newEntry);
            newName.push_back(line);
            resultSequences[seqNo].push_back(newName);
            newName.clear();
            resultSequences[seqNo].push_back(newSequence);
          } else if (sequences) {
            for (unsigned int i = 0; i < line.size();i++) {
              if (i % codonLength == 0) {
                std::string newResidue = "";
                //j for goes through all codon postions of this residue
                for (unsigned int j = i;j < i + codonLength; j++) {
                  // Use boost regular expression because compiler support for
                  // c++11 regular expressions is incomplete.
                  boost::regex re("\\w");

                  if (boost::regex_match(std::string(1, line[j]), re)) {
                    newResidue += line[j];
                  } else {
                    std::string msg = "Unexpected character: " + line[j];
                    throw std::runtime_error(msg);
                  }
                }
                resultSequences[seqNo][1].push_back(newResidue);
              }
            }
          } else{
            //else means we're already in the motifs probs section
            std::istringstream iss(line);
            InputLine motif ((std::istream_iterator<std::string>(iss)),
                              std::istream_iterator<std::string>());
            if (motif.size() == 2 && motif[0] != "motif") {
              ids->push_back(motif[0]);
              probs->push_back(std::stod(motif[1]));
            }
          }
        } else {
          sequences = false;
        }
      }
    fastafile.close();
  }
  Sequences sequences_out(resultSequences);
  return sequences_out;
}
