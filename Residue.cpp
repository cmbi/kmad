#include "Residue.h"
#include "txtProc.h"
#include <iostream>
#include <string>
#include <vector>
Residue::Residue(std::string codon, std::vector<std::string> additional_features){
	this->codon = codon;
	aa = codon[0];
	features = codon_to_features(codon);
	for (int i = 0; i < additional_features.size(); i++){
		features.push_back(additional_features[i]);
	}
}
Residue::Residue(){}
// adds features (in this order: ptm, domains, motifs)
std::vector<std::string> Residue::codon_to_features(std::string){
	std::string nothing = "AA";
	std::string feat;
	if (codon.size() >= 5){
		switch (codon[4]){	
			case 'N': feat = "ptm_phosph0";
				break;
			case 'O': feat = "ptm_phosph1"; 
				break;
			case 'P': feat = "ptm_phosph2"; 
				break;
			case 'Q': feat = "ptm_phosph3"; 
				break;
			case 'B': feat = "ptm_acet0"; 
				break;
			case 'C': feat = "ptm_acet1"; 
				break;
			case 'D': feat = "ptm_acet2"; 
				break;
			case 'E': feat = "ptm_acet3"; 
				break;
			case 'F': feat = "ptm_Nglyc0"; 
				break;
			case 'G': feat = "ptm_Nglyc1"; 
				break;
			case 'H': feat = "ptm_Nglyc2"; 
				break;
			case 'I': feat = "ptm_Nglyc3"; 
				break;
			case 'J': feat = "ptm_amid0";
				break;
			case 'K': feat = "ptm_amid1";
				break;
			case 'L': feat = "ptm_amid2";
				break;
			case 'M': feat = "ptm_amid3";
				break;
			case 'R': feat = "ptm_hydroxy0";
				break;
			case 'S': feat = "ptm_hydroxy1";
				break;
			case 'T': feat = "ptm_hydroxy2";
				break;
			case 'U': feat = "ptm_hydroxy3";
				break;
			case 'V': feat = "ptm_methyl0";
				break;
			case 'W': feat = "ptm_methyl1";
				break;
			case 'X': feat = "ptm_methyl2";
				break;
			case 'Y': feat = "ptm_methyl3";
				break;
			case 'Z': feat = "ptm_Oglyc0";
				break;
			case 'a': feat = "ptm_Oglyc1";
				break;
			case 'b': feat = "ptm_Oglyc2";
				break;
			case 'c': feat = "ptm_Oglyc3";
				break;
			case 'd': feat = "ptm_phosphP";	//predicted phosphorylation
				break;
			default: feat = nothing;
				break;
		}
		features.push_back(feat);
	}
	else features.push_back(nothing);
	//DOMAIN
	if (codon.size() >= 4){
		feat = txtProc::charToString(codon[2],codon[3]);
		if (feat != nothing){
			feat = "domain_"+feat;
			features.push_back(feat);
		}
		else{
			features.push_back(nothing);
		}
	}
	else features.push_back(nothing);
	//MOTIF
	if (codon.size() >= 7){
		feat = txtProc::charToString(codon[5],codon[6]);
		if (feat != nothing){
			feat = "motif_"+feat;
			features.push_back(feat);
		}
		else{
			features.push_back(nothing);
		}
	}
	else features.push_back(nothing);
	return features;
}
char Residue::getAA() const{
	return aa;
}
char Residue::getAA() {
	return aa;
}
std::string Residue::getCodon() const{
	return codon;
}
std::string Residue::getCodon() {
	return codon;
}
void Residue::setAA(char new_aa){
	aa = new_aa;
}
//change amino acid char to lowercase
void Residue::lowercase(){
	aa = tolower(aa);
	codon[0] = tolower(codon[0]);
}
std::vector<std::string> Residue::getFeatures() const{
	return features;
}
std::vector<std::string> Residue::getFeatures() {
	return features;
}
void Residue::add_feature(std::string new_feat){
	features.push_back(new_feat);
}
