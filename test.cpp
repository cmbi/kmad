#include <iostream>
#include <vector>
int main(){
	std::vector<std::string> myvector ={"a","s","i","a"};
	std::reverse(myvector.begin(),myvector.end());
	for (int i = 0; i < myvector.size(); i++) std::cout << myvector.at(i) << "\n";

}
