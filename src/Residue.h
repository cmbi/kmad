#ifndef RESIDUE_H
#define RESIDUE_H

#include <iostream>
#include <vector>
class Residue{
	private:
		std::string codon;
    std::vector<int> feature_indexes;
		std::vector<std::string> features;
		char aa;
		std::vector<std::string> codon_to_features(std::string);
	public:
		Residue(std::string, std::vector<std::string>);
		Residue();
		char getAA() const;
		char getAA();
		std::string getCodon() const;
		std::string getCodon();
		void setAA(char);
		void lowercase();
		void add_feature(std::string);
		std::vector<std::string> getFeatures() const;
		std::vector<std::string> getFeatures();
    std::vector<int> getFeatIndexes();
    void setFeatIndexes(std::vector<int>);
};

#endif /* RESIDUE_H */
