#include "fasta.h"
#include "misc.h"
#include "txtproc.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <iterator>
#include <stdexcept>


typedef std::vector<std::string> InputLine;

Sequences fasta::parse_fasta(std::string filename,
                             int codonLength, 
                             IDsList* ids, 
                             ProbsList* probs) {
  CodonSeqWithNamesList resultSequences;
  if (!misc::CheckIfFileExists(&filename)) {
    throw std::runtime_error("Input file doesn't exist");
  } else {
    std::string fastaSymbol = ">";
    std::ifstream fastafile (filename.c_str());
    SeqNames newName;
    CodonSeq newSequence;
    CodonSeqWithName newEntry;
    bool sequences = true;
    int seqNo = -1;
    std::string line;
    while(!txtproc::SafeGetline(fastafile, line).eof()) {
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
                    if (txtproc::AcceptedChar(line[j])) {
                      newResidue += line[j];
                    } else {
                      std::cout << "I found a weird character (" << line[j] 
                                << std::endl;
                      std::exit(EXIT_FAILURE);
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
            if (motif.size() == 2) {
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
  Sequences sequences(resultSequences);
  return sequences;
}
