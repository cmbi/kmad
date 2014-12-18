#ifndef PROFILE_H
#define PROFILE_H

#include "fasta.h"
#include "types.h"

#include <iostream>
#include <string>
#include <vector>


// TODO: namespace

typedef std::map<char, std::vector<double>> ProfileMap;

//ProfileMap create_profile(const fasta::SequenceList& sequences);
ProfileMap create_profile(const fasta::SequenceList& sequences);
ProfileMap create_score_profile(const fasta::SequenceList& sequences);


class Profile{
public:
  ///
  /// constructor; creates a Profile object with profile matrix mat
  ///
  Profile(ProfileMatrix mat);
  ///
  /// constructor; creates an empty profile
  ///
  Profile();
  ///
  /// builds a profile from the matrix of occurences and the substitution
  /// matrix with appropriate weights
  ///
  void ProcessProfile(std::vector<fasta::Sequence>& alignment);
  ///
  /// returns the profile matrix
  ///
  ProfileMatrix get_matrix() const;
  ///
  /// returns a score for amino acid aAcid on a certain position
  ///
  double get_element(int position, char aacid);
  ///
  /// returns a score for amino acid with index aAcidInt on a certain position
  ///
  double get_element(int aacid_index, int position);
private:
  void CreateProfile(std::vector<fasta::Sequence>& alignment);
  ProfileMatrix m_prf_matrix;
};

#endif /* PROFILE_H */
