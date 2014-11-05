#ifndef FINDVAL_H
#define FINDVAL_H

#include <iostream>
#include <vector>
#include <string>
namespace findVal{
	double maxValueDoubles(double,double,double); 		
	std::vector< std::vector<int> > nMaxValues(std::vector<int>&, int);
	int getMaxDoubleValuesIndex(std::vector<double>&);
}

#endif /* FINDVAL_H */
