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


fasta::FastaData fasta::parse_fasta(std::string const& filename,
                                    int codon_length)
{
  fs::path p(filename);
  if (!fs::exists(p)) {
    throw std::invalid_argument("File not found: " + filename);
  }

  std::ifstream fastafile(filename.c_str());
  std::string line;
  std::string header;
  std::string seq_line;
  bool in_sequence_section = true;
  fasta::FastaData fd;
  while (std::getline(fastafile, line)) {
    if (line.substr(0, 1) == ">") {
        if (!seq_line.empty()) {
            auto sequence = fasta::make_sequence(
                header, seq_line, codon_length
            );
            fd.sequences.push_back(sequence);
            seq_line = "";
        }
        in_sequence_section = true;

        // Parse header
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        header = line;
        continue;
    }

    if (line.substr(0, 1) == "#") {
        assert(!seq_line.empty());

        auto sequence = fasta::make_sequence(header, seq_line, codon_length);
        fd.sequences.push_back(sequence);
        in_sequence_section = false;
        continue;
    }

    if (in_sequence_section) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        seq_line += line;
    } else {
        std::vector<std::string> result;
        boost::split(result, line, boost::is_any_of("\t "));

        if (result.size() != 2) {
            throw std::runtime_error("Invalid probability format: " + line);
        }
        fd.probabilities["m_" + result[0]] = std::stod(result[1]);
    }
  }

  // If the file doesn't have a probability section, add the last sequence that
  // was being processed.
  if (in_sequence_section) {
      auto sequence = fasta::make_sequence(header, seq_line, codon_length);
      fd.sequences.push_back(sequence);
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

bool fasta::check_length(fasta::SequenceList const& sequences, int limit) {
  if (limit == 0) {
    limit = sequences.size();
  }
  bool result = true;
  size_t prev_length = sequences[0].residues.size();
  int i = 1;
  while (result && i < limit) {
    size_t length = sequences[i].residues.size();
    if (length != prev_length) {
      result = false;
    }
    ++i;
  }
  return result;
}

fasta::FastaData fasta::get_conf_data(
    const fasta::FastaData& fasta_data,
    const f_config::FeatureSettingsMap& f_set, bool gapped) {
  fasta::FastaData f;
  f.probabilities = fasta_data.probabilities;
  f.sequences = fasta_data.sequences;
  if (!gapped) {
    f.sequences = remove_gaps(f.sequences);
  }
  for (auto feat_it = f_set.begin(); feat_it != f_set.end(); ++feat_it) {
    assign_feature_by_pattern(f.sequences, feat_it->second.pattern,
        feat_it->first);
    for (auto& seq : feat_it->second.positions) {
      if ((signed)f.sequences.size() > seq.seq_no && seq.seq_no >= 0) {
        for (auto& pos : seq.positions) {
          if ((signed)f.sequences[seq.seq_no].residues.size() > pos && pos >= 0) {
            f.sequences[seq.seq_no].residues[pos].features.push_back(
                feat_it->first);
          }
          else {
            std::cout << "Warning: feature positions should be in range: 1 - "
                          << "sequence length, feature " << feat_it->first
                          << " cannot be annotated at position " << pos
                          << " in sequence " << seq.seq_no << std::endl;
          }
        }
      }
      else {
        std::cout << "Warning: sequence numbers should be in range: 1 - "
                      << "number of sequences (" << f.sequences.size()
                      << "), feature " << feat_it->first
                      << " cannot be annotated in sequence "
                      << seq.seq_no << std::endl;
      }
    }
  }
  f.feature_list = make_feature_list(f.sequences);
  return f;
}

fasta::SequenceList fasta::remove_gaps(
    const fasta::SequenceList& sequences) {
  fasta::SequenceList s = sequences;
  for (auto& seq : s) {
    seq.residues.clear();
  }
  for (size_t i = 0; i < sequences.size(); ++i) {
    for (size_t j = 0; j < sequences[i].residues.size(); ++j) {
      if (sequences[i].residues[j].codon[0] != '-') {
        s[i].residues.push_back(sequences[i].residues[j]);
      }
    }
  }
  return s;
}


void fasta::assign_feature_by_pattern(fasta::SequenceList& sequences,
                                      const std::string& pattern,
                                      const std::string& feat_name)
{
  if (pattern.size() > 0) {
    boost::regex re(pattern);
    for (size_t i = 0; i < sequences.size(); ++i) {
      std::string seq = fasta::sequence_to_string(sequences[i]);
      std::string seq_nogaps = seq;
      seq_nogaps.erase(std::remove(seq_nogaps.begin(), seq_nogaps.end(), '-'),
          seq_nogaps.end());
      for(auto it = boost::sregex_iterator(seq_nogaps.begin(), seq_nogaps.end(),
            re);
              it != boost::sregex_iterator();
                   ++it)
      {
        int match_start = find_real_pos(seq, it->position());
        int match_end = find_real_pos(seq, match_start + it->str().size());
        for (int j = match_start; j < match_end; ++j) {
          if (sequences[i].residues[j].codon[0] != '-') {
            sequences[i].residues[j].features.push_back(feat_name);
          }
        }
      }
    }
 }
}

int fasta::find_real_pos(const std::string& sequence, int position) {
  int pos = 0;
  size_t i = 0;
  while (i < sequence.size() && pos < position) {
    if (sequence[i] != '-') {
      ++pos;
    }
    ++i;
  }
  return i;
}

FeatureNamesList fasta::make_feature_list(
    const fasta::SequenceList& sequences) {
  FeatureNamesList feature_list = {
          "p_phosph0", "p_phosph1", "p_phosph2", "p_phosph3", "p_phosphP",
          "p_acet0", "p_acet1", "p_acet2", "p_acet3", "p_Nglyc0", "p_Nglyc1",
          "p_Nglyc2", "p_Nglyc3", "p_amid0", "p_amid1", "p_amid2", "p_amid3",
          "p_hydroxy0", "p_hydroxy1", "p_hydroxy2", "p_hydroxy3", "p_methyl0",
          "p_methyl1", "p_methyl2", "p_methyl3", "p_Oglyc0", "p_Oglyc1",
          "p_Oglyc2", "p_Oglyc3", "p_cys_bridge0", "s_a_helix", "s_turn",
          "s_b_ladder", "s_b_bridge", "s_310_helix", "s_pi_helix",
          "s_b_ladder"};
  for (auto& seq : sequences) {
    for (auto& res : seq.residues) {
      for (auto& feat_name : res.features) {
        if (std::find(feature_list.begin(), feature_list.end(), feat_name)
            == feature_list.end()) {
          feature_list.push_back(feat_name);
        }
      }
    }
  }
  return feature_list;
}
