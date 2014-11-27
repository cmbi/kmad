#ifndef PROFILE_H
#define PROFILE_H

#include "types.h"
#include <iostream>
#include <string>
#include <vector>
class Residue;
class Profile{
public:
	Profile(std::vector< std::vector<double> >); //constructor
	Profile();
	void buildPseudoProfile(std::vector< std::vector< Residue > >&, 
                          const std::vector<double>&, bool);
	//getters/setters
	std::vector< std::vector<double> > getMatrix() const;
	double getElement(int, char);
	double getElement(int, int);
private:
	//functions
	void createProfile(std::vector<std::vector<Residue>>&,
                     const std::vector<double>&, bool);

	//variables
	std::vector< std::vector<double> > m_prfMatrix;
};

#endif /* PROFILE_H */
