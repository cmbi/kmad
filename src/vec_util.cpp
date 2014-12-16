#include "vec_util.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>


void vec_util::TransposeVec(ProfileMatrix& vec) {
  ProfileMatrix new_vec;
  ProfileMatrixRow new_row;
  for (unsigned int i = 0; i < vec[0].size(); i++) {
    for (unsigned int j = 0; j < vec.size(); j++) {
      new_row.push_back(vec[j][i]);
    }
    new_vec.push_back(new_row);
    new_row.clear();
  }
  vec = new_vec;
}


void vec_util::DivideVectorByAScalar(ProfileMatrixRow& vec, int scalar) {
  ProfileMatrixRow result;
  for (auto &item : vec) {
    result.push_back(item / scalar);
  }
  vec = result;
}


void vec_util::DivideVectorByAScalar(ProfileMatrixRow& vec, double& scalar) {
  ProfileMatrixRow result;
  for (auto &item : vec) {
    result.push_back(item / scalar);
  }
  vec = result;
}


void vec_util::MultiplyVectorByAScalar(ProfileMatrixRow& vec, double& scalar) {
  ProfileMatrixRow result;
  for (auto &item : vec) {
    result.push_back(item * scalar);
  }
  vec = result;
}


ProfileMatrixRow vec_util::AddUp(Matrix2D& vec) {
  ProfileMatrixRow new_vec;
  for (unsigned int i = 0; i < vec[0].size(); i++) {
    double sum = 0;
    for (unsigned int j = 0; j < vec.size(); j++) {
      sum += vec[j][i];
    }
    new_vec.push_back(sum);
  }
  return new_vec;
}


ProfileMatrixColumn vec_util::ConvertIntVecToDoubleVec(SbstMatColumn& vec) {
  ProfileMatrixColumn result;
  for (auto &item : vec) {
    result.push_back(double(item));
  }
  return result;
}


std::vector<std::string> vec_util::Flatten(
    const std::vector<fasta::Sequence>& sequences)
{
  // TODO: Use functional tools
  std::vector<std::string> result;
  for (auto& sequence: sequences) {
    std::string new_seq = "";
    for (auto& residue: sequence.residues) {
      new_seq += residue.codon;
    }
    result.push_back(new_seq);
  }
  return result;
}


ProfileMatrixRow vec_util::Average(const SbstMatrixColumns& vec) {
  ProfileMatrixRow result;
  for (unsigned int i = 0; i < vec[0].size(); i++) {
    double sum = 0;
    for (unsigned int j = 0; j < vec.size(); j++) {
      sum += vec[j][i];
    }
    result.push_back(sum / vec.size());
  }
  return result;
}
