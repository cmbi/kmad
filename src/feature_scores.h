#ifndef FEATURESPROFILE_H
#define FEATURESPROFILE_H

#include "fasta.h"
#include "f_config.h"


typedef std::vector<double> Occurences; 
typedef std::vector<double> Scores;

class FeatureScores {
  public:
    ///
    /// costructor
    ///
    FeatureScores(const std::vector<std::string> features,
                  double domain_modifier, double ptm_modifier, double motif_modifier,
                    std::map<std::string, double> motif_probabilities);
    ///
    /// updates the 'profile-like' feature scores based on an alignment
    ///
    void update_scores(
        const fasta::SequenceList& sequences,
        const f_config::FeatureSettingsMap& usr_feature_settings, 
        const std::vector<double>& identities,
        const bool fade_out);

    ///
    /// get the complete 'profile-like' matrix of feature scores
    ///
    std::map<std::string, Scores> get_scores();
    ///
    /// get a single score (for a particular feature on a particular position)
    ///
    double get_score(const std::string& feat_name,
                     unsigned long position) const;
  private:
    double m_domain_modifier;
    double m_ptm_modifier;
    double m_motif_modifier;
    std::vector<std::string> m_features;
    std::map<std::string, double> m_motif_probabilities;
    std::map<std::string, Scores> m_scores;
    std::map<std::string, Occurences> m_occurences;

    std::map<std::string, Occurences> update_occurences(
        const fasta::SequenceList& sequences,
        const std::vector<double>& identities, const bool fade_out);
    double score_ptm(unsigned long position, std::string feat_name);
    double score_domain(unsigned long position, std::string feat_name);
    double score_motif(unsigned long position, std::string feat_name);
    double score_usr_feature(unsigned long position, std::string feat_name,
                             f_config::FeatureSettings settings);
};

#endif /* FEATURESPROFILE_H */
