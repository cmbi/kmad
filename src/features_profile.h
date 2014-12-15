#ifndef FEATURESPROFILE_H
#define FEATURESPROFILE_H
#include "types.h"

#include <iostream>
#include <string>
#include <vector>
class Residue;
class FeaturesProfile{
public:
  ///
  /// constructor, creates empty profile
  /// @param dom domain score
  /// @param phosph PTM score
  /// @param motif motif score
  /// @param lcr LCR gap modifier
  /// @param m_ids list of motif ids
  /// @param m_probs list of motif probabilities
  ///
  FeaturesProfile(int dom, int phosph, int motif,
                  IDsList m_ids, ProbsList m_probs);
  ///
  /// extracts the score for features on a particular position from profile
  /// @param position position in the profile
  /// @param features list of features assigned to the residue
  /// @param add_score output score
  ///
  void get_score(unsigned int position, FeaturesList& features,
                double& add_score);
  ///
  /// creates a profile matrix from the given alignment
  ///
  void CreateProfile(const SequenceList& alignment,
                     int codon_length);
  void ProcessProfile();
  ///
  /// processes the matrix of feature occurences creating a matrix of scores
  /// for aligning features at certain positions
  ///
  void CountOccurences(const SequenceList& alignment,
                       int codon_length);
  ///
  /// Takes a list of sequences, finds motifs and domains in it and adds them
  /// to the list of features
  ///
  void ExpandListOfFeatures(const SequenceList&);

  /// Finds index of a certain feature in the profile
  int FindFeaturesIndex(std::string& feat_name);
private:
  double GetMotifsProb(std::string& m_id);
  double ScoreMotifs(unsigned int& position, std::string& feat_name);
  double ScoreDomains(unsigned int& position, std::string& dom_name);
  double ScorePTMs(unsigned int& position, std::string& ptm_name);
  double ScoreUsrFeatures(unsigned int& position, std::string& feat_name);
  ///
  /// returns score modifier for a given feature
  ///
  double GetModifier(std::string& feat_name);
  int m_domain_score, m_phosph_score, m_motif_score;
  IDsList m_motifs_ids;
  ProbsList m_motifs_probs;
  PrcRulesList m_rules;
  ProfileMatrix m_prf_matrix;
  ProfileMatrix m_occurences_matrix;
};

#endif /* FEATURESPROFILE_H */
