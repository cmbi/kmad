#include "feature_scores.h"

#include <boost/algorithm/string.hpp>

#include <algorithm>


FeatureScores::FeatureScores(const std::vector<std::string>& features,
    double domain_modifier, double ptm_modifier, double motif_modifier,
    double strct_modifier, std::unordered_map<std::string, double> motif_probabilities)
: m_domain_modifier(domain_modifier),
  m_ptm_modifier(ptm_modifier),
  m_strct_modifier(strct_modifier),
  m_motif_modifier(motif_modifier),
  m_features(features),
  m_motif_probabilities(motif_probabilities)
{}


void FeatureScores::update_scores(const fasta::SequenceList& sequences,
    const f_config::FeatureSettingsMap& f_set,
    const std::vector<double>& identities,
    const bool fade_out) {
  m_occurences = update_occurences(sequences, identities, fade_out);
  // convert occurences to probabilities
  for (auto& occ: m_occurences) {
    size_t i = 0;
    for (auto& v: occ.second) {
      occ.second[i] = v / sequences.size();
      ++i;
    }
  }
  // create an empty score map
  std::unordered_map<std::string, Scores> scores;
  for (auto& f: m_features) {
    scores[f] = std::vector<double>(sequences[0].residues.size(), 0);
  }
  for (size_t i = 0; i < sequences[0].residues.size(); ++i) {
    for (auto &feat : m_features){
      if (feat.substr(0, 2) == "p_") {
        scores[feat][i] = score_ptm(i, feat);
      } else if (feat.substr(0, 2) == "d_") {
        scores[feat][i] = score_domain(i, feat);
      } else if (feat.substr(0, 2) == "m_") {
        scores[feat][i] = score_motif(i, feat);
      } else if (feat.substr(0, 2) == "s_") {
        scores[feat][i] = score_strct(i, feat);
      }
      else if (feat.substr(0, 3) == "USR") {
        scores[feat][i] = score_usr_feature(i, feat, f_set.at(feat));
      }
    }
  }
  m_scores = scores;
}


std::unordered_map<std::string, Occurences> FeatureScores::update_occurences(
    const fasta::SequenceList& sequences,
    const std::vector<double>& identities, const bool fade_out) {

  std::unordered_map<std::string, Occurences> p;
  for (auto& f: m_features) {
    p[f] = std::vector<double>(sequences[0].residues.size(), 0);
  }

  for (size_t i = 0; i < sequences[0].residues.size(); ++i) {
    for (size_t j = 0; j < sequences.size(); ++j) {
      for (auto& f : sequences[j].residues[i].features) {
        assert(p.find(f) != p.end());
        if (fade_out) {
          p[f][i] += 0.5 * (1 + identities[j]);
        } else {
          p[f][i] += 1.0;
        }
      }
    }
  }
  return p;
}


double FeatureScores::score_ptm(
                unsigned long position, const std::string& ptm_name) {
  double result  = 0;
  std::string ptm_type = ptm_name;
  //pop back last character to get just the ptm type
  ptm_type.pop_back();
  // level of annotation - last character of feature's name
  auto ptm_level = ptm_name.back();
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
  // get occurrence based score
  for (auto feat_it = m_occurences.begin();
       feat_it != m_occurences.end(); ++feat_it) {
    auto i_name = feat_it->first;
    // get just the ptm type without its level of annotation (last char)
    std::string i_type = i_name.substr(0, i_name.size() - 1);
    if (i_type == ptm_type) {
      auto i_level = i_name.back();
      double score = feat_it->second[position];
      if (i_level == '0') {
        result += score;
      } else if (i_level == '1') {
        result += score * 0.9;
      } else if (i_level == '2') {
        result += score * 0.8;
      } else if (i_level == '3') {
        result += score * 0.7;
      } else if (i_level == 'P') {
        result += score * 0.3;
      }
    }
  }
  return result * ptm_score * m_ptm_modifier;
}


double FeatureScores::score_strct(unsigned long position,
                                  const std::string& strct_name) {
  double result = 0;
  for (auto feat_it = m_occurences.begin();
       feat_it != m_occurences.end(); ++feat_it) {
    if (feat_it->first == strct_name) {
      result += feat_it->second[position];
    } else if (feat_it->first.substr(0, 2) == "s_") {
      result -= feat_it->second[position];
    }
  }
  return result * m_strct_modifier;
}


double FeatureScores::score_domain(
                unsigned long position, const std::string& dom_name) {
  double result = 0;
  for (auto feat_it = m_occurences.begin();
       feat_it != m_occurences.end(); ++feat_it) {
    if (feat_it->first == dom_name) {
      result += feat_it->second[position];
    } else if (feat_it->first.substr(0, 2) == "d_") {
      result -= feat_it->second[position];
    }
  }
  return result * m_domain_modifier;
}


double FeatureScores::score_motif(
                unsigned long position, const std::string& feat_name) {
  return m_occurences[feat_name][position] * m_motif_modifier \
           * m_motif_probabilities[feat_name];
}


double FeatureScores::score_usr_feature(
                unsigned long position, const std::string& feat_name,
                const f_config::FeatureSettings& settings) {
  double result = 0;
  // features that add scores
  for (auto& feature : settings.add_features) {
    result += m_occurences[feature][position] * settings.add_score;
  }
  // features that subtract scores
  for (auto& feature : settings.subtract_features) {
    result -= m_occurences[feature][position] * settings.subtract_score;
  }
  return result;
}


double FeatureScores::get_score(
                const std::string& feat_name, unsigned long position) const {
  return m_scores.at(feat_name)[position];
}

std::unordered_map<std::string, std::vector<double>> FeatureScores::get_scores() {
        return m_scores;
}
