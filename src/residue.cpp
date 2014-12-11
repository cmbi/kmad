#include "residue.h"
#include "txtproc.h"
#include <iostream>
#include <string>
#include <vector>


Residue::Residue(std::string codon)
: m_codon(codon) {
  m_aa = codon[0];
  CodonToFeatures();
}


Residue::Residue(char aa, int codon_length){
  m_codon = std::string(codon_length, 'A');
  m_codon[0] = aa;
  m_aa = aa;
}


Residue::Residue() {}


// adds features (in this order: ptm, domains, motifs)
void Residue::CodonToFeatures() {
  std::string nothing = "AA";
  std::string feat = "";
  if (m_codon.size() >= 5) {
    switch (m_codon[4]) {
      case 'N': feat = "ptm_phosph0";
        break;
      case 'O': feat = "ptm_phosph1";
        break;
      case 'P': feat = "ptm_phosph2";
        break;
      case 'Q': feat = "ptm_phosph3";
        break;
      case 'B': feat = "ptm_acet0";
        break;
      case 'C': feat = "ptm_acet1";
        break;
      case 'D': feat = "ptm_acet2";
        break;
      case 'E': feat = "ptm_acet3";
        break;
      case 'F': feat = "ptm_Nglyc0";
        break;
      case 'G': feat = "ptm_Nglyc1";
        break;
      case 'H': feat = "ptm_Nglyc2";
        break;
      case 'I': feat = "ptm_Nglyc3";
        break;
      case 'J': feat = "ptm_amid0";
        break;
      case 'K': feat = "ptm_amid1";
        break;
      case 'L': feat = "ptm_amid2";
        break;
      case 'M': feat = "ptm_amid3";
        break;
      case 'R': feat = "ptm_hydroxy0";
        break;
      case 'S': feat = "ptm_hydroxy1";
        break;
      case 'T': feat = "ptm_hydroxy2";
        break;
      case 'U': feat = "ptm_hydroxy3";
        break;
      case 'V': feat = "ptm_methyl0";
        break;
      case 'W': feat = "ptm_methyl1";
        break;
      case 'X': feat = "ptm_methyl2";
        break;
      case 'Y': feat = "ptm_methyl3";
        break;
      case 'Z': feat = "ptm_Oglyc0";
        break;
      case 'a': feat = "ptm_Oglyc1";
        break;
      case 'b': feat = "ptm_Oglyc2";
        break;
      case 'c': feat = "ptm_Oglyc3";
        break;
      case 'd': feat = "ptm_phosphP";  //predicted phosphorylation
        break;
      default: feat = nothing;
        break;
    }
    m_features.push_back(feat);
  }
  //DOMAIN
  if (m_codon.size() >= 4) {
    feat = std::string(1, m_codon[2]) + std::string(1, m_codon[3]);
    if (feat != nothing) {
      feat = "domain_" + feat;
      m_features.push_back(feat);
    }
  }
  //MOTIF
  if (m_codon.size() >= 7) {
    feat = std::string(1, m_codon[5]) + std::string(1, m_codon[6]);
    if (feat != nothing) {
      feat = "motif_" + feat;
      m_features.push_back(feat);
    }
  }
}


char Residue::get_aa() const{
  return m_aa;
}


std::string Residue::get_codon() const{
  return m_codon;
}


void Residue::change_to_lowercase() {
  m_aa = tolower(m_aa);
  m_codon[0] = tolower(m_codon[0]);
}


FeatureNamesList Residue::get_features() const{
  return m_features;
}


void Residue::add_feature(std::string new_feat) {
  m_features.push_back(new_feat);
}


FeaturesList Residue::get_feat_indexes() {
  return m_feature_indexes;
}


void Residue::set_feat_indexes(FeaturesList new_features) {
  m_feature_indexes = new_features;
}
