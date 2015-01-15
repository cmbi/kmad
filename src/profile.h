#ifndef PROFILE_H
#define PROFILE_H

#include "fasta.h"


namespace profile {
  typedef std::map<char, std::vector<double>> ProfileMap;

  ProfileMap create_profile(const fasta::SequenceList& sequences);
  ProfileMap create_score_profile(const fasta::SequenceList& sequences);
}

#endif /* PROFILE_H */
