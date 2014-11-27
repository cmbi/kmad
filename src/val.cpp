#include "val.h"
#include <iostream>
#include <vector>
#include <string>


/* function maxValuesFromVector
finds n max values from vector<int>  vectorA */
std::vector< std::vector<int> > val::nMaxValues(std::vector<int>& vectorA, int n){
	int i = 0;
	std::vector< std::vector<int> > result;
	std::vector<int> newEntry;
	if ((signed)vectorA.size() < n){
		n = vectorA.size();
	}

	while (i < n){
		newEntry.push_back(vectorA[i]);
		newEntry.push_back(i);
		result.push_back(newEntry);
		newEntry.clear();
		i++;
	}

	for (unsigned int i = n; i < vectorA.size(); i++){
		for (int j = 0; j < n; j++){
			if (vectorA[i] > result[j][0] ){
				int min = 100000;
				int minIndex = -1;
				for (int k = 0; k < n; k++){
					if (min > result[k][0]){
						min = result[k][0];
						minIndex = k;	
					}	
				}
				result[minIndex][0] = vectorA[i];
				result[minIndex][1] = i;
				break;
			}
		}
	}
	return result;
}
