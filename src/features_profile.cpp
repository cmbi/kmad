#include "features_profile.h"
#include "txtproc.h"
#include "vec_util.h"

#include <boost/range/numeric.hpp>
#include <boost/algorithm/string.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <tuple>


FeaturesProfileMap create_score_features_profile(
    const fasta::SequenceList& sequences, 
    const std::vector<std::string>& features, int ptm_modifier, 
    int domain_modifier, int motif_modifier, 
    std::map<std::string, double> probabilities) {
  FeaturesProfileMap p = create_features_profile(sequences, features);
  // convert occurences to probabilities
  for (auto& occ: p) {
    size_t i = 0;
    for (auto& v: occ.second) {
      occ.second[i] = v / sequences.size();
      i++;
    }
  }
  // create an empty score map 
  FeaturesProfileMap score_p;
  for (auto& f: features) {
    score_p[f] = std::vector<double>(sequences[0].residues.size(), 0);
  }
  for (size_t i = 0; i < sequences[0].residues.size(); i++) {
    for (auto &feat : features){
      if (feat.substr(0,3) == "ptm") {
        score_p[feat][i] = score_ptm(p, i, feat, ptm_modifier);
      } else if (feat.substr(0,6) == "domain") {
        score_p[feat][i] = score_domain(p, i, feat, domain_modifier);
      } else if (feat.substr(0,5) == "motif") {
        score_p[feat][i] = score_motif(p, i, feat, motif_modifier, 
                                       probabilities);
      } 
      // else if (feat.substr(0,3) == "USR") {
      //   score_p[feat][i] = score_usr_feature(p, i, feat);
      // }
    }
  }
  return p;
}


FeaturesProfileMap create_features_profile(
    const fasta::SequenceList& sequences,
    const std::vector<std::string>& features) {
  FeaturesProfileMap p;
  for (auto& f: features) {
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


namespace {
  // 0 - highest level of annotation, 3 - lowest, P - predicted
  FeatureNamesList list_of_features = {"ptm_phosph0", "ptm_phosph1",
                                       "ptm_phosph2", "ptm_phosph3",
                                       "ptm_phosphP", "ptm_acet0",
                                       "ptm_acet1", "ptm_acet2",
                                       "ptm_acet3", "ptm_Nglyc0",
                                       "ptm_Nglyc1", "ptm_Nglyc2",
                                       "ptm_Nglyc3", "ptm_amid0",
                                       "ptm_amid1", "ptm_amid2",
                                       "ptm_amid3", "ptm_hydroxy0",
                                       "ptm_hydroxy1", "ptm_hydroxy2",
                                       "ptm_hydroxy3", "ptm_methyl0",
                                       "ptm_methyl1", "ptm_methyl2",
                                       "ptm_methyl3", "ptm_Oglyc0",
                                       "ptm_Oglyc1", "ptm_Oglyc2",
                                       "ptm_Oglyc3", "domain_0",
                                       "motif_0", "lcr"};
  IndexList domain_indexes = {29};
  IndexList motif_indexes = {30};
  std::string nothing = "AA";
  std::string domain = "domain";
}


FeaturesProfile::FeaturesProfile(int dom, int phosph, int motif,
                                 std::map<std::string, double> probs)
: m_domain_score(dom),
  m_phosph_score(phosph),
  m_motif_score(motif)
{
  // TODO: Use probabilities directly
  for (const auto& kv: probs) {
    m_motifs_ids.push_back(kv.first);
    m_motifs_probs.push_back(kv.second);
  }
}


void FeaturesProfile::CreateProfile(const std::vector<fasta::Sequence>& alignment,
                                    int codon_length) {
  CountOccurences(alignment, codon_length);
  ProcessProfile();
}


void FeaturesProfile::ProcessProfile() {
  m_prf_matrix.clear();
  for (unsigned int i = 0; i < m_occurences_matrix.size(); i++) {
    std::string feat_name = list_of_features[i];
    ProfileMatrixRow feature_row(m_occurences_matrix[i].size());
    for (unsigned int j = 0; j < m_occurences_matrix[i].size(); j++) {
      if (feat_name.substr(0,3) == "ptm") {
        feature_row[j] = ScorePTMs(j, feat_name);
      } else if (feat_name.substr(0,6) == "domain") {
        feature_row[j] = ScoreDomains(j, feat_name);
      } else if (feat_name.substr(0,5) == "motif") {
        feature_row[j] = ScoreMotifs(j, feat_name);
      } else if (feat_name.substr(0,3) == "USR") {
        feature_row[j] = ScoreUsrFeatures(j, feat_name);
      }
    }
    m_prf_matrix.push_back(feature_row);
  }
}


// TODO: Implement
void FeaturesProfile::CountOccurences(const std::vector<fasta::Sequence>& alignment,
                                      int codon_length) {
  //m_occurences_matrix.clear();
  //std::string nothing = "AA";
  //int no_of_sequences = alignment.size();

  //for (unsigned int i = 0; i < alignment[0].size(); i++) {
    //ProfileMatrixColumn profile_column(list_of_features.size(),0);
    //int non_gaps = 0;
    //for (unsigned int j = 0; j < alignment.size(); j++) {
      //if (alignment[j][i].get_aa() != '-') {
        //FeatureNamesList features = alignment[j][i].get_features();
        //for (unsigned int k = 0; k < features.size(); k++) {
          //std::string feat_name = features[k];
          //if (feat_name != nothing) {
            //int feat_index = FindFeaturesIndex(feat_name);
            //profile_column[feat_index] += GetModifier(feat_name);
          //}
        //}
      //non_gaps++;
      //}
    //}
    //vec_util::DivideVectorByAScalar(profile_column, no_of_sequences);
    ////vec_util::DivideVectorByAScalar(profileColumn,nonGaps);
    //m_occurences_matrix.push_back(profile_column);
  //}
  //vec_util::TransposeVec(m_occurences_matrix);
}


double FeaturesProfile::GetMotifsProb(std::string& m_id) {
  double prob = 0;

  FeatNameSplit motif_name;
  boost::split(motif_name, m_id, boost::is_any_of("_"));
  std::string id_code = motif_name[1];

  assert(m_motifs_ids.size() == m_motifs_probs.size());
  for (unsigned int i = 0; i < m_motifs_ids.size(); i++) {
    if (m_motifs_ids[i] == id_code) {
      prob = m_motifs_probs[i];
      break;
    }
  }
  return prob;
}


void FeaturesProfile::get_score(unsigned int position, FeaturesList& features,
                                double& add_score) {
  for (auto &feat : features) {
       add_score += m_prf_matrix[feat][position];
  }
}


double FeaturesProfile::ScoreMotifs(unsigned int& position,
                                    std::string& feat_name) {
  int features_index = FindFeaturesIndex(feat_name);
  double result = 0;
  if (features_index == -1) {
    result = 0;
  } else {
    result = m_occurences_matrix[features_index][position];
  }
  return result;
}


int FeaturesProfile::FindFeaturesIndex(std::string& feat_name) {
  int features_index = -1;
  for (unsigned int i = 0; i < list_of_features.size(); i++) {
    if (feat_name == list_of_features[i]) {
      features_index = i;
      break;
    }
  }
  return features_index;
}


void FeaturesProfile::ExpandListOfFeatures(const std::vector<fasta::Sequence> sequences)
{
  // TODO: Implement
  // TODO: Don't mess with stuff in the anonymous namespace.

  //for (unsigned int i = 0; i < sequences.size(); i++) {
    //for (unsigned int j = 0; j < sequences[i].size(); j++) {
      //FeatureNamesList features = sequences[i][j].get_features();
      //for (unsigned int k = 0; k < features.size(); k++) {
        //std::string feature_k = features[k];
        ////check whether this domain is already in the list of features
        //if ((std::find(list_of_features.begin(), list_of_features.end(),
                      //feature_k) == list_of_features.end())
            //&& feature_k != nothing) {
          //list_of_features.push_back(feature_k);
          //if (feature_k.substr(0, 6) == "domain") {
            ////to look later for domains only in these positions (saves time)
            //domain_indexes.push_back(list_of_features.size() - 1);
          //} else if (feature_k.substr(0, 5) == "motif") {
            //// same as with domain_indexes
            //motif_indexes.push_back(list_of_features.size()-1);
          //}
        //}
      //}
    //}
  //}
}


double FeaturesProfile::ScoreDomains(unsigned int& position,
                                     std::string& dom_name) {
  double result = 0;
  for (unsigned int i = 0; i < domain_indexes.size(); i++) {
    int dom_index = domain_indexes[i];
    if (list_of_features[dom_index] == dom_name) {
      result += m_occurences_matrix[dom_index][position];
    } else {
      // subtract score for every other domain
      result -= m_occurences_matrix[dom_index][position];
    }
  }
  return result;
}


double score_domain(FeaturesProfileMap& p, unsigned long position,
                    std::string dom_name, int domain_modifier) {
  double result = 0;
  for (auto feat_it = p.begin(); feat_it != p.end(); feat_it++) {
    if (feat_it->first == dom_name) {
      result += feat_it->second[position];
    } else {
      result -= feat_it->second[position];
    }
  }
  return result * domain_modifier;
}



double score_motif(FeaturesProfileMap& p, unsigned long position,
                   std::string feat_name, int motif_modifier,
                    std::map<std::string, double> probs) {
  double result = 0;
  bool not_found = true;
  for (auto feat_it = p.begin(); feat_it != p.end() && not_found; feat_it++) {
    if (feat_it->first == feat_name) {
      result = feat_it->second[position] * motif_modifier * probs[feat_name];
      not_found = false;     
    }
  }
  // int features_index = FindFeaturesIndex(feat_name);
  // if (features_index == -1) {
  //   result = 0;
  // } else {
  //   result = m_occurences_matrix[features_index][position];
  // }
  return result;
}


double score_ptm(FeaturesProfileMap& p, unsigned long position,
                 std::string ptm_name, int ptm_modifier) {
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
  for (auto feat_it = p.begin(); feat_it != p.end(); feat_it++) {
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
  result = result * ptm_score * ptm_modifier;
  return result;
}


double FeaturesProfile::ScorePTMs(unsigned int& position,
                                  std::string& ptm_name) {
  double result  = 0;
  std::string ptm_type = ptm_name;
  //pop back last character to get just the ptm type
  ptm_type.pop_back();
  // level of annotation - last character of feature's name
  char ptm_level = ptm_name.back();
  double ptm_score;
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
  // now go through list of features to find in which rows in profile
  // are the features that we're gonna score for
  assert(m_occurences_matrix.size() == list_of_features.size());
  for (unsigned int i = 0; i < list_of_features.size(); i++) {
    std::string i_name = list_of_features[i];
    std::string i_type = i_name;
    // popping back last character, to get just the ptm type
    // (without its level of annotation)
    i_type.pop_back();
    if (i_type == ptm_type) {
      char i_level = i_name.back();
      if (i_level == '0') {
        result += m_occurences_matrix[i][position];
      } else if (i_level == '1') {
        result += m_occurences_matrix[i][position] * 0.9;
      } else if (i_level == '2') {
        result += m_occurences_matrix[i][position] * 0.8;
      } else if (i_level == '3') {
        result += m_occurences_matrix[i][position] * 0.7;
      } else if (i_level == 'P') {
        result += m_occurences_matrix[i][position] * 0.3;
      }
    }
  }
  result = result * ptm_score;
  return result;
}


double FeaturesProfile::ScoreUsrFeatures(unsigned int& position,
                                         std::string& feat_name) {
  //first find tuple(s) with rules for this feature
  double result = 0;
  for (auto &rule : m_rules) {
    if (std::get<0>(rule) == feat_name) {
      double add_tmp = std::get<2>(rule);
      //positions of increasing features
      FeaturesList incr_features = std::get<6>(rule);
      //go through features that increase the score
      for (auto &feat : incr_features) {
        double prf_score = m_occurences_matrix[feat][position];
        if (prf_score != 0) {
           result += add_tmp*prf_score;
        }
      }
      //the same for decreasing features
      add_tmp = std::get<4>(rule);
      FeaturesList decr_features = std::get<7>(rule);
      //go through features that increase the score
      for (auto &feat : decr_features) {
        double prf_score = m_occurences_matrix[feat][position];
        if (prf_score != 0) {
          result -= add_tmp*prf_score;
        }
      }
    }
  }
  return result;
}


double FeaturesProfile::GetModifier(std::string& feat_name) {

  FeatNameSplit feat_name_vec;
  boost::split(feat_name_vec, feat_name, boost::is_any_of("_"));
  std::string feat_code  = feat_name_vec[0];

  double modifier = 1;
  if (feat_code == "ptm") {
    modifier = m_phosph_score;
  } else if (feat_code == "domain") {
    modifier = m_domain_score;
  } else if (feat_code == "motif") {
    modifier = m_motif_score*GetMotifsProb(feat_name);
  }
  return modifier;
}
