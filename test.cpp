#include <iostream>
#include <vector>
using namespace std;

int main(){
	const vector<int> a;
	a = {1,2};
	const vector< vector<int> > hello={{1,2},{2,3}};
	cout << hello.at(1).at(1);	
	hello[1][1] = 10;
	cout << hello.at(1).at(1);	

}
