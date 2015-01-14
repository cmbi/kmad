#include "features_profile.h"

#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>


FeaturesProfile::FeaturesProfile(std::vector<std::string> features,
    int domain_modifier, int ptm_modifier, int motif_modifier,
    std::map<std::string, double> motif_probabilities)
: m_domain_modifier(domain_modifier),
  m_ptm_modifier(ptm_modifier),
  m_motif_modifier(motif_modifier),
  m_features(features),
  m_motif_probabilities(motif_probabilities)
{}


void FeaturesProfile::update_scores(const fasta::SequenceList& sequences, 
    const f_config::FeatureSettingsMap& f_set) {
  m_occurences = update_occurences(sequences);
  // convert occurences to probabilities
  for (auto& occ: m_occurences) {
    size_t i = 0;
    for (auto& v: occ.second) {
      occ.second[i] = v / sequences.size();
      i++;
    }
  }
  // create an empty score map 
  std::map<std::string, Scores> scores;
  for (auto& f: m_features) {
    scores[f] = std::vector<double>(sequences[0].residues.size(), 0);
  }
  for (size_t i = 0; i < sequences[0].residues.size(); i++) {
    for (auto &feat : m_features){
      if (feat.substr(0,3) == "ptm") {
        scores[feat][i] = score_ptm(i, feat);
      } else if (feat.substr(0,6) == "domain") {
        scores[feat][i] = score_domain(i, feat);
      } else if (feat.substr(0,5) == "motif") {
        scores[feat][i] = score_motif(i, feat);
      } 
      else if (feat.substr(0,3) == "USR") {
        scores[feat][i] = score_usr_feature(i, feat, f_set.at(feat));
      }
    }
  }
  m_scores = scores;
}


std::map<std::string, Occurences> FeaturesProfile::update_occurences(
    const fasta::SequenceList& sequences) {
  std::map<std::string, Occurences> p;
  for (auto& f: m_features) {
    p[f] = std::vector<double>(sequences[0].residues.size(), 0);
  }
  for (size_t i = 0; i < sequences[0].residues.size(); i++) {
    for (size_t j = 0; j < sequences.size(); j++) {
      for (auto& f : sequences[j].residues[i].features) {
        assert(p.find(f) != p.end());
        p[f][i] += 1.0;
      }
    }
  }
  return p;
}


double FeaturesProfile::score_ptm(unsigned long position, 
                                  std::string ptm_name) {
  double result  = 0;
  std::string ptm_type = ptm_name;
  //pop back last character to get just the ptm type
  ptm_type.pop_back();
  // level of annotation - last character of feature's name
  char ptm_level = ptm_name.back();
  double ptm_score = 0.;
  // first set ptm_score based on annotation level of the query ptm
  if (ptm_level == '0') {
    ptm_score = 1.0;
  } else if (ptm_level == '1') {
    ptm_score = 0.9;
  } else if (ptm_level == '2') {
    ptm_score = 0.8;
  } else if (ptm_level == '3') {
    ptm_score = 0.7;
  } else if (ptm_level == 'P') {
    ptm_score = 0.3;
  } else {
    std::string msg = "wrong annotation level on position "
                    + std::to_string(position)
                    + " ptmname: " + ptm_name;
    throw std::invalid_argument(msg);
  }
  for (auto feat_it = m_occurences.begin();
       feat_it != m_occurences.end(); feat_it++) {
    std::string i_name = feat_it->first;
    std::string i_type = i_name;
    // popping back last character, to get just the ptm type
    // (without its level of annotation)
    i_type.pop_back();
    if (i_type == ptm_type) {
      char i_level = i_name.back();
      if (i_level == '0') {
        result += feat_it->second[position];
      } else if (i_level == '1') {
        result += feat_it->second[position] * 0.9;
      } else if (i_level == '2') {
        result += feat_it->second[position] * 0.8;
      } else if (i_level == '3') {
        result += feat_it->second[position] * 0.7;
      } else if (i_level == 'P') {
        result += feat_it->second[position] * 0.3;
      }
    }
  }
  result = result * ptm_score * m_ptm_modifier;
  return result;
}


double FeaturesProfile::score_domain(unsigned long position, 
                                     std::string dom_name) {
  double result = 0;
  for (auto feat_it = m_occurences.begin(); 
       feat_it != m_occurences.end(); feat_it++) {
    if (feat_it->first == dom_name) {
      result += feat_it->second[position];
    } else if (feat_it->first.substr(0,6) == "domain") {
      result -= feat_it->second[position];
    }
  }
  return result * m_domain_modifier;
}


double FeaturesProfile::score_motif(unsigned long position,
                                    std::string feat_name) {
  double result = 0;
  bool not_found = true;
  for (auto feat_it = m_occurences.begin();
       feat_it != m_occurences.end() && not_found; feat_it++) {
    if (feat_it->first == feat_name) {
      result = feat_it->second[position] * m_motif_modifier \
               * m_motif_probabilities[feat_name];
      not_found = false;     
    }
  }
  return result;
}


double FeaturesProfile::score_usr_feature(unsigned long position,
                                          std::string feat_name,
                                          f_config::FeatureSettings settings) {
  double result = 0;
  // features that add scores
  for (auto& feature : settings.add_features) {
    result += m_occurences[feat_name][position] \
              * m_occurences[feature][position] \
              * settings.add_score;
  }
  // features that subtract scores
  for (auto& feature : settings.subtract_features) {
    result -= m_occurences[feat_name][position] \
              * m_occurences[feature][position] \
              * settings.subtract_score;
  }
  return result;
}


std::map<std::string, Scores> FeaturesProfile::get_scores() {
  return m_scores;
}

double FeaturesProfile::get_score(const std::string& feat_name,
                                        unsigned long position) const {
  return m_scores.at(feat_name)[position];
}
