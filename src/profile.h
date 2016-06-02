#ifndef PROFILE_H
#define PROFILE_H

#include "fasta.h"


namespace profile {
  typedef std::unordered_map<char, std::vector<double>> ProfileMap;

  ProfileMap create_profile(const fasta::SequenceList& sequences);
  ProfileMap create_score_profile(const fasta::SequenceList& sequences,
                                  const std::string& sbst_mat);
  double get_score(const ProfileMap& p, int position, char aa);
}

#endif /* PROFILE_H */
