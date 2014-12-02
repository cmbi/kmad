#include "misc.h"
#include "txtproc.h"
#include "vec_util.h"
#include "residue.h"
#include "sequences.h"
#include "features_profile.h"

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
typedef std::vector<std::string> SplitFeatName;
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


std::vector<std::string> txtproc::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


void txtproc::writeAlignmentToFile(StringSequences& sequences,
                                   SeqNames& sequence_names, 
                                   std::string filename) {
  std::stringstream sstr;
  sstr << filename << "_al";
  std::ofstream outputFile(sstr.str().c_str(), std::ios::out);
  for (unsigned int i = 0; i < sequences.size(); i++) {
    outputFile << sequence_names[i] << "\n" << sequences[i] << "\n";
  }
}


void txtproc::writeAlignmentWithoutCodeToFile(StringSequences& sequences,
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


std::string txtproc::charToString(char mychar) {
  return std::string(1, mychar);
}


std::string txtproc::charToString(char mychar1, char mychar2) {
  std::string newstring = std::string(1, mychar1);
  newstring.push_back(mychar2);
  return newstring;
}


std::istream& txtproc::safeGetline(std::istream& is, std::string& t) {
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


bool txtproc::acceptedChar(char my_char) {
  bool result = false;
  for (auto &acc_char : AcceptedCharacters) {
    if (acc_char == my_char) {
      result = true;
      break;
    }
  }
  return result;
}


void txtproc::process_conf_file(std::string filename, 
                                FeaturesProfile& feat_profile, 
                                Sequences& sequences_aa) {
  std::ifstream conf_file(filename.c_str());
  std::string line;
  RuleTuplesList usr_feature_rules;
  DefaultRulesList feature_rules;
  bool features = true;
  std::string tag_usr = "## USER DEFINED";
  while(!safeGetline(conf_file, line).eof()) {
       if (features) {
          std::size_t found = line.find(tag_usr);
          if (found != std::string::npos) {
            features = false;
          }
       } else if (line[0] != '#' && features) {
         line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
         line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
         FeatDescriptor tmp_vector = split(line, ';');
         feature_rules.push_back(std::make_tuple(tmp_vector[0], tmp_vector[1], 
                                                 std::stoi(tmp_vector[2]), 
                                                 std::stoi(tmp_vector[3]), 
                                                 std::stoi(tmp_vector[4])));
       } else if (line[0] != '#') {
         line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
         line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
         FeatDescriptor tmp_vector = split(line, ';');
         usr_feature_rules.push_back(std::make_tuple(tmp_vector[0], 
                                                     tmp_vector[1], 
                                                     std::stoi(tmp_vector[2]), 
                                                     std::stoi(tmp_vector[3]), 
                                                     std::stoi(tmp_vector[4]),
                                                     std::stod(tmp_vector[5]),
                                                     std::stod(tmp_vector[6]),
                                                     std::stod(tmp_vector[7]),
                                                     std::stod(tmp_vector[8]),
                                                     tmp_vector[9],
                                                     tmp_vector[10]));
    }
  }
  feat_profile.add_usr_features(usr_feature_rules);
  sequences_aa.add_usr_features(usr_feature_rules);
  feat_profile.set_rules(usr_feature_rules);
}


FeaturesList txtproc::unfold(std::string conf_string, 
                             FeatureNamesList& list_of_features) {
  FeatureNamesList tmp_vector = split(conf_string, ',');
  FeaturesList out_vector;
  for (auto &item : tmp_vector) {
    if (split(item, '_').size() > 1) {            
      // this is a single feature entry, e.g. 'PF_A'
      std::string feat_name = std::string("USR_") + item;
      out_vector.push_back(vec_util::FindIndex(feat_name, list_of_features));
    } else if (split(item, '[').size() == 1) {            
      // this is an entry with only the tag specified (without any exceptions)
      for (unsigned int j = 0; j < list_of_features.size(); j++) {
        SplitFeatName singlefeat = split(list_of_features[j], '_');
        if (singlefeat.size() > 1 && singlefeat[1] == item) {
          out_vector.push_back(j);
        }
      }
    } else {
      //TAG with exceptions
      SplitFeatName tagfeat = split(item, '[');
      std::string tag = tagfeat[0];
      FeatureNamesList exceptions = split(split(tagfeat[1], ']')[0], '.');
      for (unsigned int j = 0; j < list_of_features.size(); j++) {
        SplitFeatName singlefeat = split(list_of_features[j], '_');
        if (singlefeat.size() > 1 && singlefeat[1] == tag && 
            !vec_util::CheckIfContains(exceptions, singlefeat[2])) {
          out_vector.push_back(j);
        }
      }
    }
  }
  return out_vector;
}
