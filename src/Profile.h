#ifndef PROFILE_H
#define PROFILE_H

#include <iostream>
#include <string>
#include <vector>
class Residue;
class Profile{
public:
	Profile(std::vector< std::vector<double> >); //constructor
	Profile();
	void buildPseudoProfile(std::vector< std::vector< Residue > >&, const std::vector<double>&, bool);
	//getters/setters
	std::vector< std::vector<double> > getMatrix() const;
	void setMatrix(std::vector<std::vector<double> >);
	double getElement(int, char);
	double getElement(int, int);
	void printProfile(int,int);
	void printProfile();
private:
	//functions
	void createProfile(std::vector< std::vector<Residue> >&,const std::vector<double>&,bool);
	double countNonGaps(int);
	int getMaxDoubleValue(std::vector<double>);
	//variables
	std::vector< std::vector<double> > prfMatrix;
};

#endif /* PROFILE_H */
