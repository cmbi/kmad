#ifndef FEATURESPROFILE_H
#define FEATURESPROFILE_H

#include "fasta.h"
#include "f_config.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

typedef std::vector<double> Occurences; 
typedef std::vector<double> Scores;

class FeaturesProfile {
  public:
    FeaturesProfile(const std::vector<std::string> features,
                    int domain_modifier, int ptm_modifier, int motif_modifier,
                    std::map<std::string, double> motif_probabilities);
    void update_scores(const fasta::SequenceList& sequences,
                       const f_config::FeatureSettingsMap&
                       usr_feature_settings);
    std::map<std::string, Scores> get_scores();
    double get_score(const std::string& feat_name,
                     unsigned long position) const;
  private:
    int m_domain_modifier;
    int m_ptm_modifier;
    int m_motif_modifier;
    std::vector<std::string> m_features;
    std::map<std::string, double> m_motif_probabilities;
    std::map<std::string, Scores> m_scores;
    std::map<std::string, Occurences> m_occurences;

    std::map<std::string, Occurences> update_occurences(
        const fasta::SequenceList& sequences);
    double score_ptm(unsigned long position, std::string feat_name);
    double score_domain(unsigned long position, std::string feat_name);
    double score_motif(unsigned long position, std::string feat_name);
    double score_usr_feature(unsigned long position, std::string feat_name,
                             f_config::FeatureSettings settings);
};

#endif /* FEATURESPROFILE_H */
