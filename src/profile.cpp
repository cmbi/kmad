#include "profile.h"

#include <boost/range/numeric.hpp>

#include <algorithm>
#include <iostream>

typedef std::map<char, std::vector<double>> SimilarityScoresMap;
namespace {
  static const std::vector<char> ALPHABET = {'A','R','N','D','C','Q','E','G',
                                             'H','I','L','K','M','F','P','S',
                                             'T','W','Y','V'};
  static const SimilarityScoresMap SIM_SCORES = {
    {'A', { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0}},
    {'R', {-1,  5,  0, -2, -3,  1,  0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3}},
    {'N', {-2,  0,  6,  1, -3,  0,  0,  0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3}},
    {'D', {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3}},
    {'C', { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1}},
    {'Q', {-1,  1,  0,  0, -3,  5,  2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2}},
    {'E', {-1,  0,  0,  2, -4,  2,  5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2}},
    {'G', { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3}},
    {'H', {-2,  0,  1, -1, -3,  0,  0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3}},
    {'I', {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3}},
    {'L', {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1}},
    {'K', {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2}},
    {'M', {-1, -1, -2, -3, -1,  0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1}},
    {'F', {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1}},
    {'P', {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2}},
    {'S', { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2}},
    {'T', { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0}},
    {'W', {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3}},
    {'Y', {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1}},
    {'V', { 0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4}}};
}

ProfileMap create_score_profile(const fasta::SequenceList& sequences) {
  ProfileMap p = create_profile(sequences);

  // Convert the profile occurrences to probabilities.
  for (auto& occ: p) {
    size_t i = 0;
    for (auto& v: occ.second) {
      occ.second[i] = v / sequences.size();
      i++;
    }
  }
  for (unsigned i = 0; i < p['A'].size(); i++) {
    std::vector<double> score_column(20, 0); 
    for (auto &prob: p){
      std::vector<double> sbst_column = SIM_SCORES.at(prob.first);
      for (size_t k = 0; k < sbst_column.size(); ++k) {
        score_column[k] += sbst_column[k] * prob.second[i];
      }
    }
    // replace the old column (with probs) with the new column (with scores)
    for (size_t j = 0; j < ALPHABET.size(); ++j) {
      p[ALPHABET[j]][i] = score_column[j];
    }
  }
  return p;
}


ProfileMap create_profile(const fasta::SequenceList& sequences) {
  ProfileMap p;

  // Initialise profile map with all valid letters to a vector of the correct
  // size, with all values set to 0.
  for (auto& c: ALPHABET) {
    p[c] = std::vector<double>(sequences[0].residues.size(), 0);
  }

  for (size_t i = 0; i < sequences[0].residues.size(); i++) {
    for (size_t j = 0; j < sequences.size(); j++) {
      char amino_acid = sequences[j].residues[i].codon[0];

      if (amino_acid == '-') {
        continue;
      }

      if (amino_acid == 'B') {
        assert(p.find('D') != p.end() && p.find('N') != p.end());
        p['D'][i] += 0.5;
        p['N'][i] += 0.5;
      } else if (amino_acid == 'Z') {
        assert(p.find('E') != p.end() && p.find('G') != p.end());
        p['E'][i] += 0.5;
        p['G'][i] += 0.5;
      } else if (amino_acid == 'X') {
        for (auto& kv: p) {
          kv.second[i] += 0.05;
        }
      } else {
        assert(p.find(amino_acid) != p.end());
        p[amino_acid][i] += 1.0;
      }
    }
  }
  return p;
}
