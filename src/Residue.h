#ifndef RESIDUE_H
#define RESIDUE_H

#include <iostream>
#include <vector>
class Residue{
	private:
		std::string m_codon;
    std::vector<int> m_feature_indexes;
		std::vector<std::string> m_features;
		char m_aa;
		//std::vector<std::string> codon_to_features();
		void codon_to_features();
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

typedef std::vector<Residue> sequence;
typedef std::vector<sequence> sequenceList;
#endif /* RESIDUE_H */
