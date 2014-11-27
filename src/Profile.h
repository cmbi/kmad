#ifndef PROFILE_H
#define PROFILE_H

#include "types.h"
#include <iostream>
#include <string>
#include <vector>
class Residue;
class Profile{
public:
	Profile(profile_matrix mat); //constructor
	Profile();
	void buildPseudoProfile(sequenceList& alignment, 
                          const identitiesList& sequenceIdentityValues, 
                          bool weightsModeOn);
	//getters/setters
	profile_matrix getMatrix() const;
	double getElement(int position, char aAcid);
	double getElement(int aAcidInt, int position);
private:
	//functions
	void createProfile(sequenceList&,
                     const identitiesList& sequenceIdentityValues, 
                     bool weightsModeOn);
	//variables
	profile_matrix m_prfMatrix;
};

#endif /* PROFILE_H */
