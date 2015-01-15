#include "profile.h"

#include "substitution_matrix.h"

#include <boost/range/numeric.hpp>

#include <algorithm>
#include <iostream>


namespace {
  static const std::vector<char> ALPHABET = {'A','R','N','D','C','Q','E','G',
                                             'H','I','L','K','M','F','P','S',
                                             'T','W','Y','V'};
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
      std::vector<double> sbst_column = substitution_matrix::get_column(
          prob.first);
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
