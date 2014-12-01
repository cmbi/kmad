#ifndef RESIDUE_H
#define RESIDUE_H

#include "types.h"
#include <iostream>
#include <vector>
class Residue{
	public:
    ///
    /// constructor; creates a Residue object based on the given codon
    ///
		Residue(std::string codon);
    ///
    /// constructor; creates an empty Residue object
    ///
		Residue();
    ///
    /// returns the residue's amino acid code
    ///
		char getAA() const;
    ///
    /// returns residue's codon
    ///
		std::string getCodon() const;
    ///
    /// changes amino acid charactre to lowercase
    ///
		void lowercase();
    ///
    /// adds a new feature to the residue's list of features
    ///
		void add_feature(std::string new_feat);
    ///
    /// returns list of features assigned to the residue
    ///
		FeatureNamesList getFeatures() const;
    ///
    /// return a list of feature indexes (from the profile) assigned to the
    /// residue
    ///
    FeaturesList getFeatIndexes();
    ///
    /// sets feature indexes
    ///
    void setFeatIndexes(FeaturesList new_features);
	private:
		std::string m_codon;
    FeaturesList m_feature_indexes;
		FeatureNamesList m_features;
		char m_aa;
		void codon_to_features();
};

#endif /* RESIDUE_H */
