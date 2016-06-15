#include "fasta.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <fstream>
#include <sstream>


namespace fs = boost::filesystem;
namespace {
  std::unordered_map<char, std::string> strct_code_map = {
          {'H', "s_a_helix"}, {'T', "s_turn"}, {'S', "s_b_ladder"},
          {'B', "s_b_bridge"}, {'G', "s_310_helix"}, {'I', "s_pi_helix"},
          {'E', "s_b_ladder"},
  };
  std::unordered_map<char, std::string> ptm_code_map = {
          {'N', "p_phosph0"}, {'O', "p_phosph1"}, {'P', "p_phosph2"},
          {'Q', "p_phosph3"}, {'B', "p_acet0"}, {'C', "p_acet1"},
          {'D', "p_acet2"}, {'E', "p_acet3"}, {'F', "p_Nglyc0"},
          {'G', "p_Nglyc1"}, {'H', "p_Nglyc2"}, {'I', "p_Nglyc3"},
          {'J', "p_amid0"}, {'K', "p_amid1"}, {'L', "p_amid2"},
          {'M', "p_amid3"}, {'R', "p_hydroxy0"}, {'S', "p_hydroxy1"},
          {'T', "p_hydroxy2"}, {'U', "p_hydroxy3"}, {'V', "p_methyl0"},
          {'W', "p_methyl1"}, {'X', "p_methyl2"}, {'Y', "p_methyl3"},
          {'Z', "p_Oglyc0"}, {'a', "p_Oglyc1"}, {'b', "p_Oglyc2"},
          {'c', "p_Oglyc3"}, {'d', "p_phosphP"}, {'s', "p_cys_bridge0"}};
}


fasta::FastaData fasta::parse_fasta(
        std::string filename, int codon_length, bool refine, int refine_seq) {

        fs::path p(filename);
        if (!fs::exists(p)) {
                throw std::invalid_argument("File not found: " + filename);
        }

  std::ifstream fastafile(filename.c_str());
  std::string line;
  bool in_sequence_section = true;
  fasta::FastaData fd;
  while (std::getline(fastafile, line)) {
    if (line == std::string("## PROBABILITIES")) {
      in_sequence_section = false;
      continue;
    }

    if (in_sequence_section && line.substr(0, 1) == ">") {
      std::string description = line;
      std::getline(fastafile, line);
      line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
      fd.sequences.push_back(fasta::make_sequence(description, line,
                                                  codon_length));
    } else if (in_sequence_section) {
      int last_index = fd.sequences.size() - 1;
      fasta::extend_sequence(fd.sequences[last_index], line, codon_length);

    } else {
      std::vector<std::string> result;
      boost::split(result, line, boost::is_any_of("\t "));

      if (result.size() != 2) {
        throw std::runtime_error("Invalid probability format: " + line);
      }
      fd.probabilities["m_" + result[0]] = std::stod(result[1]);
    }
  }
  fastafile.close();
  if (refine && !check_length(fd.sequences, refine_seq)) {
        throw std::runtime_error("In the 'refine' mode all sequences should "
                                 "have the same length");
  }

  return fd;
}


fasta::Sequence fasta::make_sequence(const std::string& description,
                                     const std::string& codons,
                                     int codon_length)
{
  fasta::Sequence s;
  s.description = description;

  for (unsigned i = 0; i < codons.size(); i += codon_length) {
    boost::regex re("[a-zA-Z0-9-]{" + std::to_string(codon_length) + "}");
    std::string codon = codons.substr(i, codon_length);
    if (!boost::regex_match(codon, re)) {
      throw std::runtime_error("Invalid codon: " + codon);
    }
    fasta::Residue r = fasta::make_residue(codon);
    s.residues.push_back(r);
  }
  return s;
}

fasta::Sequence fasta::make_sequence(const std::vector<Residue>& residue_list)
{
  fasta::Sequence s;
  s.residues = residue_list;
  return s;
}

fasta::Sequence fasta::make_sequence(unsigned long sequence_length,
                                     const fasta::Residue& residue) {
  fasta::Sequence s;
  s.residues = std::vector<fasta::Residue>(sequence_length, residue);
  return s;
}

fasta::Sequence fasta::make_sequence(const std::string& sequence_string,
                int codon_length) {
  fasta::Sequence s;
  for (size_t i = 0; i < sequence_string.size(); ++i) {
    if (i % codon_length == 0) {
      std::string codon = sequence_string.substr(i, codon_length);
      fasta::Residue res = fasta::make_residue(codon);
      s.residues.push_back(res);
    }
  }
  return s;
}
void fasta::extend_sequence(fasta::Sequence& seq,
                     const std::string& sequence_string, int codon_length) {
  fasta::Sequence tmp_seq = fasta::make_sequence(sequence_string,
      codon_length);
  for (auto& res : tmp_seq.residues) {
    seq.residues.push_back(res);
  }
}

std::string fasta::sequence_to_string(const fasta::Sequence& seq) {
  std::string result;
  for (auto& residue : seq.residues) {
    result.push_back(residue.codon[0]);
  }
  return result;
}

fasta::Residue fasta::make_residue(const std::string& codon) {
  std::vector<std::string> features;
  if (codon.size() >= 5 && codon[4] != 'A') {
    features.push_back(ptm_code_map.at(codon[4]));
  }
  if (codon.size() > 1
      && strct_code_map.find(codon[1]) != strct_code_map.end()) {
    features.push_back(strct_code_map.at(codon[1]));
  }
  if (codon.size() >= 7 && codon.substr(5, 2) != "AA") {
    features.push_back("m_" + codon.substr(5, 2));
  }
  if (codon.size() >= 4 && codon.substr(2, 2) != "AA") {
    features.push_back("d_" + codon.substr(2, 2));
  }
  return Residue(codon, features);
}

bool fasta::check_length(fasta::SequenceList sequences, int refine_seq) {
  if (refine_seq == 0) {
    refine_seq = sequences.size();
  }
  bool result = true;
  size_t prev_length = sequences[0].residues.size();
  int i = 1;
  while (result && i < refine_seq){
    size_t length = sequences[i].residues.size();
    if (length != prev_length) {
      result = false;
    }
    ++i;
  }
  return result;
}
