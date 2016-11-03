#include "f_config.h"

#include <algorithm>
#include <boost/regex.hpp>
#include <iostream>

namespace lcg = libconfig;

f_config::FeatureSettingsMap f_config::ConfParser::parse_conf_file(
    const std::string& filename) {
  lcg::Config cnfg;
  try
  {
    cnfg.readFile(filename.c_str());
  }
  catch (const lcg::FileIOException &fioex) {
    std::cerr << "I/O error while reading file." << std::endl;
    throw;
  }
  catch (const lcg::ParseException &pex) {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    throw;
  }
  f_config::RawFeatureSettingsMap raw_map = process_config(cnfg);
  return process_settings(raw_map);
}


f_config::RawFeatureSettingsMap f_config::ConfParser::process_config(
    const lcg::Config& cnfg) {
  f_config::RawFeatureSettingsMap feat_config;
  try
  {
    const lcg::Setting& root  = cnfg.getRoot();
    const lcg::Setting& features = root["feature_settings"]["usr_features"];
    int count = features.getLength();
    for (int i = 0; i < count; ++i) {
      RawFeatureSettings feat_set;
      const lcg::Setting& feature = features[i];
      std::string name;
      if (!feature.lookupValue("name", name))
        continue;
      feature.lookupValue("tag", feat_set.tag);
      bool found_add_score = feature.lookupValue("add_score",
                                                 feat_set.add_score);
      bool found_sbtrct_score = feature.lookupValue("subtract_score",
                                                    feat_set.subtract_score);
      if (!(found_add_score || found_sbtrct_score))
        continue;
      feature.lookupValue("pattern", feat_set.pattern);

      lcg::Setting& add_features_set = feature["add_features"];
      for (int j = 0; j < add_features_set.getLength(); ++j) {
        feat_set.add_features.push_back(add_features_set[j]);
      }
      lcg::Setting& add_tags_set = feature["add_tags"];
      for (int j = 0; j < add_tags_set.getLength(); ++j) {
        feat_set.add_tags.push_back(add_tags_set[j]);
      }
      lcg::Setting& add_exceptions_set = feature["add_exceptions"];
      for (int j = 0; j < add_exceptions_set.getLength(); ++j) {
        feat_set.add_exceptions.push_back(add_exceptions_set[j]);
      }
      lcg::Setting& subtract_features_set = feature["subtract_features"];
      for (int j = 0; j < subtract_features_set.getLength(); ++j) {
        feat_set.subtract_features.push_back(subtract_features_set[j]);
      }
      lcg::Setting& subtract_tags_set = feature["subtract_tags"];
      for (int j = 0; j < subtract_tags_set.getLength(); ++j) {
        feat_set.subtract_tags.push_back(subtract_tags_set[j]);
      }
      lcg::Setting& subtract_exceptions_set = feature["subtract_exceptions"];
      for (int j = 0; j < subtract_exceptions_set.getLength(); ++j) {
        feat_set.subtract_exceptions.push_back(subtract_exceptions_set[j]);
      }
      lcg::Setting& positions_set = feature["positions"];
      for (int j = 0; j <  positions_set.getLength(); ++j) {
        FeaturePositions feat_pos;
        positions_set[j].lookupValue("seq", feat_pos.seq_no);
          --feat_pos.seq_no;
        lcg::Setting& single_pos_set = positions_set[j]["pos"];
        for (int k = 0; k < single_pos_set.getLength(); ++k) {
          feat_pos.positions.push_back(single_pos_set[k]);
        }
        for (auto& pos : feat_pos.positions) {
          --pos;
        }
        feat_set.positions.push_back(feat_pos);
      }
      feat_config[name] = feat_set;
    }
  }
  catch(const lcg::SettingNotFoundException &nfex)
  {
    std::cerr << "Setting not found" << std::endl;
    throw;
  }
  return feat_config;
}


f_config::FeatureSettingsMap f_config::ConfParser::process_settings(
    const f_config::RawFeatureSettingsMap& raw_map) {
  f_config::FeatureSettingsMap processed_map;
  for (auto feat_it = raw_map.begin(); feat_it != raw_map.end(); ++feat_it) {
     FeatureSettings processed_settings;
     processed_settings.add_score = feat_it->second.add_score;
     processed_settings.subtract_score = feat_it->second.subtract_score;
     processed_settings.positions = feat_it->second.positions;
     processed_settings.pattern = feat_it->second.pattern;
     ///
     /// filter out the features from 'add_features' and 'subtract_features'
     /// that have no settings (are not in keys in the RawFetaureSettings map)
     ///
     for (auto& feat_name : feat_it->second.add_features) {
       if (raw_map.find(feat_name) != raw_map.end()) {
         processed_settings.add_features.push_back("USR_"+ feat_name);
       }
     }

     for (auto& feat_name : feat_it->second.subtract_features) {
       if (raw_map.find(feat_name) != raw_map.end()) {
         processed_settings.subtract_features.push_back("USR_"+ feat_name);
       }
     }
     ///
     /// add all features that contain the given tag, unless they are mentioned
     /// in the exceptions
     ///
     for (auto& feat_tag : feat_it->second.add_tags) {
       for (auto feat_it_j = raw_map.begin(); feat_it_j != raw_map.end();
            ++feat_it_j) {
         if (feat_it_j->second.tag == feat_tag
             && std::find(feat_it->second.add_exceptions.begin(),
                          feat_it->second.add_exceptions.end(),
                          feat_it_j->first)
                == feat_it->second.add_exceptions.end()
             && std::find(processed_settings.add_features.begin(),
                          processed_settings.add_features.end(),
                          "USR_" + feat_it_j->first)
                == processed_settings.add_features.end()) {
           processed_settings.add_features.push_back(
               "USR_" + feat_it_j->first);
         }
       }
     }

     for (auto& feat_tag : feat_it->second.subtract_tags) {
       for (auto feat_it_j = raw_map.begin(); feat_it_j != raw_map.end();
            ++feat_it_j) {
         if (feat_it_j->second.tag == feat_tag
             && std::find(feat_it->second.subtract_exceptions.begin(),
                          feat_it->second.subtract_exceptions.end(),
                          feat_it_j->first)
                == feat_it->second.subtract_exceptions.end()
             && std::find(processed_settings.subtract_features.begin(),
                          processed_settings.subtract_features.end(),
                          "USR_" + feat_it_j->first)
                == processed_settings.subtract_features.end()) {
           processed_settings.subtract_features.push_back(
               "USR_" + feat_it_j->first);
         }
       }
     }
     ///
     /// add the feature to the new map
     ///
     processed_map["USR_" + feat_it->first] = processed_settings;
  }
  return processed_map;
}

fasta::FastaData f_config::get_conf_data(
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

void f_config::assign_feature_by_pattern(fasta::SequenceList& sequences,
                                      const std::string& pattern,
                                      const std::string& feat_name)
{
  if (pattern.size() > 0) {
    boost::regex re(pattern);
    for (size_t i = 0; i < sequences.size(); ++i) {
      //std::string seq = fasta::sequence_to_string(sequences[i]);
      auto seq = fasta::sequence_to_string(sequences[i]);
      auto seq_nogaps = seq;
      seq_nogaps.erase(std::remove(seq_nogaps.begin(), seq_nogaps.end(), '-'),
          seq_nogaps.end());
      for(auto it = boost::sregex_iterator(seq_nogaps.begin(), seq_nogaps.end(),
            re);
              it != boost::sregex_iterator();
                   ++it)
      {
        auto match_start = find_real_pos(seq, it->position());
        auto match_end = find_real_pos(seq, match_start + it->str().size());
        for (int j = match_start; j < match_end; ++j) {
          if (sequences[i].residues[j].codon[0] != '-') {
            sequences[i].residues[j].features.push_back(feat_name);
          }
        }
      }
    }
 }
}

int f_config::find_real_pos(const std::string& sequence, int position) {
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

FeatureNamesList f_config::make_feature_list(
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
