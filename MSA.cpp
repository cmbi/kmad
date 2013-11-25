#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include "Sequences.h"
#include "ScoringMatrix.h"
#include "Profile.h"
#include "UsefulStuff.h"
using namespace std;
// MAIN function
UsefulStuff util;
int main(int argc, char *argv[]){
	if (argc == 5){		
		util = UsefulStuff();
		int gapPen = util.convertStringToInt(argv[2]);						//assign arguments
		string arg4 = argv[4];
		bool verboseMode = false;
		if (arg4=="-v"){
			verboseMode = true;
		}
		int codonLength = util.convertStringToInt(argv[3]);
		Sequences rawSequences;
		if (codonLength>1){
			rawSequences = Sequences(Sequences(util.processFASTA(argv[1],codonLength)));			//read data from file
		}
		else{
			rawSequences = Sequences(Sequences(util.processFASTA(argv[1])));			//read data from file
		}
		rawSequences.printEncodedSequence(0);
		/*Profile prf;										//this prf will be useful for next rounds of alignments
		vector<string> multipleAlignment(rawSequences.performMSA(&prf,gapPen,verboseMode));	//create multiple sequence alignment	
		util.writeAlignmentToFile(multipleAlignment,rawSequences.getSequences(),argv[1]);	//write multiple alignment to a file
		*/
	}
	else {
		cout << "MSA filename.fasta gapPenalty codonLength(>=1) verboseMode(-nv/-v)" << endl;
	}
	return 0;
}
