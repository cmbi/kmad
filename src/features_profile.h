#ifndef FEATURESPROFILE_H
#define FEATURESPROFILE_H


#include "fasta.h"
#include "types.h"

#include <iostream>
#include <string>
#include <vector>


typedef std::map<std::string, std::vector<double>> FeaturesProfileMap;
FeaturesProfileMap create_features_profile(
    const fasta::SequenceList& sequences, 
    const std::vector<std::string>& features);
FeaturesProfileMap create_score_features_profile(
    const fasta::SequenceList& sequences,
    const std::vector<std::string>& features);
double score_ptm(FeaturesProfileMap& p, unsigned int position,
                 std::string feat_name, int ptm_modifier);
double score_domain(FeaturesProfileMap& p, unsigned int position,
                    std::string feat_name);
double score_motif(FeaturesProfileMap& p, unsigned int position,
                   std::string feat_name);
double score_usr_feature(FeaturesProfileMap& p, unsigned int position,
                         std::string feat_name);

class FeaturesProfile {
public:
  ///
  /// constructor, creates empty profile
  /// @param dom domain score
  /// @param phosph PTM score
  /// @param motif motif score
  /// @param lcr LCR gap modifier
  ///
  FeaturesProfile(int dom, int phosph, int motif,
                  std::map<std::string, double> probs);
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
  // TODO: How does this CreateProfile differ to that in 'profile.h'?
  void CreateProfile(const std::vector<fasta::Sequence>& alignment,
                     int codon_length);
  void ProcessProfile();
  ///
  /// processes the matrix of feature occurences creating a matrix of scores
  /// for aligning features at certain positions
  ///
  void CountOccurences(const std::vector<fasta::Sequence>& alignment, int codon_length);
  ///
  /// Takes a list of sequences, finds motifs and domains in it and adds them
  /// to the list of features
  ///
  void ExpandListOfFeatures(const std::vector<fasta::Sequence> sequences);

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
  std::vector<std::string> m_motifs_ids;
  std::vector<double> m_motifs_probs;
  PrcRulesList m_rules;
  ProfileMatrix m_prf_matrix;
  ProfileMatrix m_occurences_matrix;
};

#endif /* FEATURESPROFILE_H */
