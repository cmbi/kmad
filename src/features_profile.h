#ifndef FEATURESPROFILE_H
#define FEATURESPROFILE_H

#include "fasta.h"
#include "f_config.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

typedef std::map<std::string, std::vector<double>> FeaturesProfileMap;

struct FeaturesProfile {
    int m_domain_modifier;
    int m_ptm_modifier;
    int m_motif_modifier;
    std::vector<std::string> m_features;
    std::map<std::string, double> m_motif_probabilities;
    FeaturesProfileMap m_scores;
    FeaturesProfileMap m_occurences;

    FeaturesProfile(const std::vector<std::string> features,
                    int domain_modifier, int ptm_modifier, int motif_modifier,
                    std::map<std::string, double> motif_probabilities);
    FeaturesProfileMap create_features_profile(
        const fasta::SequenceList& sequences);
    void create_score_features_profile(const fasta::SequenceList& sequences,
                                       const f_config::FeatureSettingsMap&
                                        usr_feature_settings);
    double score_ptm(unsigned long position, std::string feat_name);
    double score_domain(unsigned long position, std::string feat_name);
    double score_motif(unsigned long position, std::string feat_name);
    double score_usr_feature(unsigned long position, std::string feat_name,
                             f_config::FeatureSettings settings);
};

#endif /* FEATURESPROFILE_H */
