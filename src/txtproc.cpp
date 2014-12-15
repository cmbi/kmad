#include "txtproc.h"
#include "vec_util.h"
#include "residue.h"
#include "sequences.h"
#include "features_profile.h"

#include <boost/algorithm/string.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <tuple>
#include <algorithm>
#include <iterator>

typedef std::vector<std::string> FeatDescriptor;
typedef std::vector<std::string> SplitLine;
namespace {
  static const AlphabetVec AcceptedCharacters = { 'a', 'b', 'c', 'd', 'e', 'f',
                                                  'g', 'h', 'i', 'j', 'k', 'l',
                                                  'm', 'n', 'o', 'p', 'q', 'r',
                                                  's', 't', 'u', 'v', 'w', 'x',
                                                  'y', 'z', 'A', 'B', 'C', 'D',
                                                  'E', 'F', 'G', 'H', 'I', 'J',
                                                  'K', 'L', 'M', 'N', 'O', 'P',
                                                  'Q', 'R', 'S', 'T', 'U', 'V',
                                                  'W', 'X', 'Y', 'Z', '0', '1',
                                                  '2', '3', '4', '5', '6', '7',
                                                  '8', '9'};
}


void txtproc::WriteAlignmentToFile(StringSequences& sequences,
                                   SeqNames& sequence_names,
                                   std::string filename) {
  std::stringstream sstr;
  sstr << filename << "_al";
  std::ofstream outputFile(sstr.str().c_str(), std::ios::out);
  for (unsigned int i = 0; i < sequences.size(); i++) {
    outputFile << sequence_names[i] << "\n" << sequences[i] << "\n";
  }
}


void txtproc::WriteAlignmentWithoutCodeToFile(StringSequences& sequences,
                                              SeqNames& sequence_names,
                                              std::string filename,
                                              int codon_length) {
  std::stringstream sstr;
  sstr << filename << "_al";
  std::ofstream outputFile(sstr.str().c_str(), std::ios::out);
  for (unsigned int i = 0; i < sequences.size(); i++) {
    outputFile << sequence_names[i] << "\n";
    std::string seq = "";
    for (unsigned int j = 0; j < sequences[i].size(); j+=codon_length) {
      seq += sequences[i][j];
    }
    outputFile << seq << std::endl;
  }
}


std::istream& txtproc::SafeGetline(std::istream& is, std::string& t) {
  t.clear();
  std::streambuf* sb = is.rdbuf();
    for (;;) {
          int c = sb->sbumpc();
    switch (c) {
      case '\n':
        return is;
      case '\r':
        if (sb->sgetc() == '\n')
        sb->sbumpc();
        return is;
      case EOF:
        // Also handle the case when the last line has no line ending
        if (t.empty())
          is.setstate(std::ios::eofbit);
        return is;
      default:
        t += (char)c;
    }
  }
}


bool txtproc::AcceptedChar(char my_char) {
  bool result = false;
  for (auto &acc_char : AcceptedCharacters) {
    if (acc_char == my_char) {
      result = true;
      break;
    }
  }
  return result;
}
