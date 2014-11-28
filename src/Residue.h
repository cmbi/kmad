#ifndef RESIDUE_H
#define RESIDUE_H

#include "types.h"
#include <iostream>
#include <vector>
class Residue{
	private:
		std::string m_codon;
    featuresList m_feature_indexes;
		featureNamesList m_features;
		char m_aa;
		void codon_to_features();
	public:
		Residue(std::string codon, featureNamesList additional_features);
		Residue();
		char getAA() const;
		char getAA();
		std::string getCodon() const;
		std::string getCodon();
		void setAA(char);
		void lowercase();
		void add_feature(std::string new_feat);
		featureNamesList getFeatures() const;
		featureNamesList getFeatures();
    featuresList getFeatIndexes();
    void setFeatIndexes(featuresList new_features);
};

#endif /* RESIDUE_H */
