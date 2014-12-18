#include "substitution_matrix.h"

#include "vec_util.h"

#include<boost/range/numeric.hpp>

#include <iostream>
#include <vector>


typedef std::vector<char> AlphaList;
typedef std::vector<SbstMatColumn> SbstMatColumnsList;


namespace {
  static const AlphabetVec Alphabet = {'A','R','N','D','C','Q','E','G',
                                       'H','I','L','K','M','F','P','S',
                                       'T','W','Y','V'};
  //BLOSUM62
  static const SbstMatrix SimScores = {
    { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0},
    {-1,  5,  0, -2, -3,  1,  0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3},
    {-2,  0,  6,  1, -3,  0,  0,  0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3},
    {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3},
    { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
    {-1,  1,  0,  0, -3,  5,  2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2},
    {-1,  0,  0,  2, -4,  2,  5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2},
    { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3},
    {-2,  0,  1, -1, -3,  0,  0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3},
    {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3},
    {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1},
    {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2},
    {-1, -1, -2, -3, -1,  0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1},
    {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1},
    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2},
    { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2},
    { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0},
    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3},
    {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1},
    { 0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4}};
}


int substitution_matrix::get_element(char char1, char char2) {
  int index1 = 0;
  int index2 = 0;
  for (unsigned int i = 0; i < Alphabet.size(); i++) {
    if (Alphabet[i] == char1) {
      index1 = i;
    }
    if (Alphabet[i] == char2) {
      index2 = i;
    }
  }
  return SimScores[index1][index2];
}


ProfileMatrix substitution_matrix::ConvertToProfileFormat(
    std::vector<fasta::Residue>& seq)
{
  ProfileMatrix result(seq.size());
  SbstMatColumnsList new_sbst_row;
  for (unsigned int i = 0; i < result.size(); i++) {
    if (seq[i].codon[0] == 'B') {
      new_sbst_row.clear();
      new_sbst_row.push_back(SimScores[2]);
      new_sbst_row.push_back(SimScores[3]);
      result[i] = vec_util::Average(new_sbst_row);
    } else if (seq[i].codon[0] == 'Z') {
      new_sbst_row.clear();
      new_sbst_row.push_back(SimScores[6]);
      new_sbst_row.push_back(SimScores[7]);
      result[i] = vec_util::Average(new_sbst_row);
    } else if (seq[i].codon[0] == 'X') {
      result[i] = vec_util::Average(SimScores);
    } else {
      int aAcidInt = FindAminoAcidsIndex(seq[i].codon[0]);
      SbstMatColumn sbst_column_int = SimScores[aAcidInt];
      //adds a column to the result(converted from int to double)
      result[i] = vec_util::ConvertIntVecToDoubleVec(sbst_column_int);
    }
  }
  vec_util::TransposeVec(result);
  return result;
}


int substitution_matrix::get_element(int i, int j) {
  return SimScores[i][j];
}


void substitution_matrix::get_column(unsigned int& column_no,
                                     SbstMatColumn& column_int) {
  column_int = SimScores[column_no];
}


std::vector<double> substitution_matrix::get_column(const char& aa) {
  int column_no = FindAminoAcidsIndex(aa);
  std::vector<double> column;
  for (auto& item : SimScores[column_no]) {
    column.push_back(double(item));
  }
  return column;
}


int substitution_matrix::FindAminoAcidsIndex(char aa) {
  int aacid_index = -1;
  for (unsigned int i = 0; i < Alphabet.size(); i++) {
    if (aa == Alphabet[i]) {
      aacid_index = i;
      break;
    }
  }
  return aacid_index;
}
