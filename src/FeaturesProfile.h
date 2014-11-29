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
  FeaturesProfile(int dom, int phosph, int motif, int lcr, 
                  ids_list m_ids, probs_list m_probs);
  ///
  /// extracts the score for features on a particular position from profile
  /// @param position position in the profile
  /// @param features list of features assigned to the residue
  /// @param add_score output score
  ///
	void getScore(unsigned int position, featuresList& features, 
                double& add_score);
	profile_matrix getMatrix();
  ///
  /// creates a profile matrix from the given alignment
  ///
  void createProfile(const sequenceList& alignment, 
                     const identitiesList& sequenceIdentityValues, 
                     bool weightsModeOn, int codon_length);
  void processProfile();
  ///
  /// processes the matrix of feature occurences creating a matrix of scores
  /// for aligning features at certain positions
  ///
	void countOccurences(const sequenceList& alignment, 
                       const identitiesList& sequenceIdentityValues, 
                       bool weightsModeOn, int codon_length);
  ///
  /// Takes a list of sequences, finds motifs and domains in it and adds them 
  /// to the list of features
  ///
	void expandListOfFeatures(const sequenceList&);
  ///
  /// Sets rules for aligning user defined features
  ///
  void setRules(rulesTuplesList& new_rules);
  ///
  /// Adds user defined features to the list of features
  ///
  void add_USR_features(rulesTuplesList& new_rules);
  /// Finds index of a certain feature in the profile
	int findFeaturesIndex(std::string& feat_name);
private:
	double get_motifs_prob(std::string& m_id);
	double score_motifs(unsigned int& position, std::string& feat_name);
	double score_domains(unsigned int& position, std::string& dom_name);
	double score_PTMs(unsigned int& position, std::string& ptm_name);
	double score_USR_features(unsigned int& position, std::string& feat_name); 
  ///
  /// returns score modifier for a given feature
  ///
	double get_modifier(std::string& feat_name);
	int m_domainScore, m_phosphScore, m_motifScore;
	ids_list m_motifs_ids;
	probs_list m_motifs_probs;
  prcRulesList m_rules;
	profile_matrix m_prfMatrix;	
  profile_matrix m_occurences_matrix;
};

#endif /* FEATURESPROFILE_H */
