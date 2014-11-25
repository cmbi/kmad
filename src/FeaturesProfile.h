#ifndef FEATURESPROFILE_H
#define FEATURESPROFILE_H

#include <iostream>
#include <string>
#include <vector>
class Residue;
class FeaturesProfile{
public:
	FeaturesProfile(int, int, int, int, std::vector<std::string>, std::vector<double>);
	//getters
	void getScore(unsigned int, std::vector<std::string>&, double&);
	//double getScore(int,std::string);
	double getGapMod(int, std::vector<std::string>);
	std::vector<std::vector<double> > getMatrix();
	//setters
	void setMatrix(std::vector<std::vector<double> >);
	void printProfile();
	void printOcc();
	void createProfile(const std::vector< std::vector<Residue> >&, 
                     const std::vector<double>&, bool, int);
  void processProfile();
	void countOccurences(const std::vector< std::vector<Residue> >&, 
                     const std::vector<double>&, bool, int);

	void expandListOfFeatures(const std::vector< std::vector< Residue> > &);
	void setRules(std::vector< std::tuple<std::string, std::string, int, int, 
                                        int, double, double, double, double, 
                                        std::string, std::string> >&);

	void add_USR_features(std::vector< std::tuple<std::string, std::string, int,
                                                int, int, double, double, 
                                                double, double, std::string,
                                                std::string> >&);
	void printFeatures();
private:
	double motifs_prob(std::string&);
	double score_motifs(unsigned int&, std::string&);
	double score_domains(unsigned int&, std::string&);
	double score_PTMs(unsigned int&, std::string&);
	double score_USR_features(unsigned int&, std::string&); 
	double modifier(std::string&);
	int findFeaturesIndex(std::string&);
	std::string name(std::string&,int&);
	int domainScore, phosphScore, motifScore;
	std::vector<std::string> motifs_ids;
	std::vector<double> motifs_probs;
	std::vector<std::tuple<std::string,int,double,double,double,double,std::vector<int>,std::vector<int> > > rules;
	std::vector< std::vector<double> > prfMatrix;	
  std::vector< std::vector<double> > occurences_matrix;
};

#endif /* FEATURESPROFILE_H */
