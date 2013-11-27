#include "Sequences.h"
#include "ScoringMatrix.h"
#include "Profile.h"
#include "txtProc.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
int main(int argc, char *argv[]){
	if (argc == 4){		
		int gapPen = txtProc::convertStringToInt(argv[2]);						//assign arguments
		std::string arg3 = argv[3];
		bool verboseMode;
		if (arg3=="-v"){
			verboseMode = true;
		}
		else{
			verboseMode = false;
		}
		Sequences rawSequences(txtProc::processFASTA(argv[1]));			//read data from file
		Profile prf;										//this prf will be useful for next rounds of alignments
		std::vector<std::string> multipleAlignment(rawSequences.performMSA(&prf,gapPen,verboseMode));	//create multiple sequence alignment	
		txtProc::writeAlignmentToFile(multipleAlignment,rawSequences.getSequences(),argv[1]);	//write multiple alignment to a file
	}
	else {
		std::cout << "MSA filename.fasta gapPenalty verboseMode" << "\n";
	}
	return 0;
}
