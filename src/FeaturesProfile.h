#ifndef FEATURESPROFILE_H
#define FEATURESPROFILE_H
#include "types.h"
#include <iostream>
#include <string>
#include <vector>
class Residue;
class FeaturesProfile{
public:
  FeaturesProfile(int dom, int phosph, int motif, int lcr, 
                  ids_list m_ids, probs_list m_probs);
	//getters
	void getScore(unsigned int position, featuresList& features, 
                double& add_score);
	profile_matrix getMatrix();
	//setters
  void createProfile(const sequenceList& alignment, 
                     const identitiesList& sequenceIdentityValues, 
                     bool weightsModeOn, int codon_length);
  void processProfile();
	void countOccurences(const sequenceList& alignment, 
                       const identitiesList& sequenceIdentityValues, 
                       bool weightsModeOn, int codon_length);
	void expandListOfFeatures(const sequenceList&);
  void setRules(rulesTuplesList& new_rules);
  void add_USR_features(rulesTuplesList& new_rules);
	int findFeaturesIndex(std::string& feat_name);
private:
	double get_motifs_prob(std::string& m_id);
	double score_motifs(unsigned int& position, std::string& feat_name);
	double score_domains(unsigned int& position, std::string& dom_name);
	double score_PTMs(unsigned int& position, std::string& ptm_name);
	double score_USR_features(unsigned int& position, std::string& feat_name); 
	double get_modifier(std::string& feat_name);
	int m_domainScore, m_phosphScore, m_motifScore;
	ids_list m_motifs_ids;
	probs_list m_motifs_probs;
  prcRulesList m_rules;
	profile_matrix m_prfMatrix;	
  profile_matrix m_occurences_matrix;
};

#endif /* FEATURESPROFILE_H */
