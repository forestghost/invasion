#include <algorithm>
#include <ctime>
#include <map>
#include <cmath>
#include <random>
#include "nr3.h"

using std::string;
using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::map;
using std::time;
using std::random_shuffle;
using std::ifstream;
using std::ofstream;
using std::istringstream;

int main() {
	vector<Int> results;
	default_random_engine e(time(0));
	vector<Int> testpattern = {1,2,3,4,100};
	discrete_distribution<Int> u(testpattern.begin(), testpattern.end());
	//discrete_distribution<Int> u{0.5, 0.25, 0.1, 0.};
	for (Int i=0; i<1000; i++) 
		results.push_back(u(e));
	for (auto iter = results.begin(); iter!=results.end(); iter++)
		cout << *iter << endl;

	return 0;
}