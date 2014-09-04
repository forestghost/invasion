#include <algorithm>
#include <ctime>
#include <map>
#include <cmath>
#include <random>
#include <bitset>
#include <string>
#include <iterator>
#include <R.h>
#include "nr3.h"
#include "erf.h"
#include "exponential.h"
#include "landscape.h"
#include "matrixmult.h"
#include "walker.h"
#include "population2.h"
#include "math.h"
#include "cauchy.h"

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
using std::stringstream;
using std::max_element;
using std::isdigit;

extern "C" void getlcps() {   // pass it the river matrix  and the lcps_to_find matrix

	default_random_engine e(time(0));
	default_random_engine &ernie = e;

	Int n = 720;   //// NOTE, these will change if a different preloaded landscape is used!!!!!!
	Int m = 666;

	vector<vector<Int> > M1 = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};  // Moore neighborhood (radius 1)

	// read in waterways
	string riverfile("river_raster");  
	MatInt rivers = read_matrix_int(riverfile, n, m);

	// read in the positions of sampled plant cells
	vector<vector<Int>> source_cells;

	ifstream data;
	string fn("find_these_lcps");
	data.open(fn.c_str());

	Int a, b;
	while (data >> a >> b) 
		source_cells.push_back({a,b});
	data.close();

	Int sizey = source_cells.size();
	cout << sizey << " plants to account for" << endl;	
	MatInt lcps = create_emptymatrix_int(sizey, sizey);

	for (Int i = 0; i < sizey-1; i++) {
		for (Int j = i+1; j < sizey; j++) {

			cout << i << " " << j << endl;
			if (i == j) 
				continue;
			vector<Int> destination = source_cells[j];

			// find distances ////
			map<vector<Int>, vector<Int>> exploration;
			exploration[source_cells[i]] = {0,0}; // start with source cell; first value is distance from source, second value is whether its been explored (1) or not (0)
			bool destination_notyet = true;
			bool nonew_riverpixel = false;
			while (destination_notyet) {
				map<vector<Int>, vector<Int>> tempexp;
				for (auto a : exploration) {
					if ( a.second[1]  == 0) { // then not used yet
						for (auto g : M1) {  // determine distances for adjacent (by first order Moore neighborhood) cells
							Int query1 = a.first[0] + g[0];
							Int query2 = a.first[1] + g[1];
							vector<Int> query;
							query.push_back(query1);
							query.push_back(query2);
							if (exploration.count(query) == 0 && rivers[query1][query2] == 0) {  // then not in exploration yet AND river pixel: add to tempexp map, and, eventually, exploration map
								tempexp[query] = {a.second[0] + 1, 0};  // one further than current pixel
							}
							if (query == destination) {
								lcps[i][j] = a.second[0]+1;
								lcps[j][i] = a.second[0]+1;
								destination_notyet = false;
								exploration.insert(tempexp.begin(), tempexp.end());
								break;
							}	
						}
						if (destination_notyet == false)
							break;
						a.second[1] = 1; // mark as used
					}	
				}
				if (tempexp.size() == 0) {
					nonew_riverpixel = true;
					break;
				}
				exploration.insert(tempexp.begin(), tempexp.end());
			}	
/*
			if (nonew_riverpixel == true) {  // then isolated waterway, determine new, truncated destination
				vector<vector<Int>> terminal_nodes;
				for (auto a : exploration) {
					bool terminal = true;
					for (auto g : M1) {  // determine distances for adjacent (by first order Moore neighborhood) cells
						Int query1 = a.first[0] + g[0];
						Int query2 = a.first[1] + g[1];
						vector<Int> query = { query1, query2};
						if (rivers[query1][query2] == 0 && exploration.count(query) > 0 && exploration[query][0] > a.second[0] )
							terminal = false;
					}
					if (terminal == true)
						terminal_nodes.push_back(a.first);
				}
				Doub min_distance = 1e06;
				vector<Int> truncated_destination;
				for (auto b : terminal_nodes) {
					Doub disty = sqrt ( pow(b[0] - destination[0], 2) + pow(b[1]-destination[1], 2)   );
					if (disty < min_distance) {
						min_distance = disty;
						truncated_destination = b;
					}
				}		
				destination = truncated_destination;	
			}
			// trace path backwards  ///// 
			vector<Int> current = destination;
	//cout <<  current[0] << " " << current[1] << endl;
			Int distance = exploration[current][0];
			bool notsource = true;
			if (exploration.size() == 1) { // then source and destination are the same;
				outfile << source[0] << " " << source[1] << " " << source[0] << " " << source[1] << endl;
				flow_directions[source] = source;
			} else {
				while (notsource) {
					for (auto g : M1) {  // determine distances for adjacent (by first order Moore neighborhood) cells
						Int query1 = current[0] + g[0];
						Int query2 = current[1] + g[1];
						vector<Int> query;
						query.push_back(query1);
						query.push_back(query2);
						if (exploration.count(query)>0) { // then key is in exploration map
							if (exploration[query][0] == distance -1 ) {
								if (flow_directions.count(query) == 0)  // then cell's flow-to cell has not been determined and should print to file
									outfile << query[0] << " " << query[1] << " " << current[0] << " " << current[1] << endl;
								flow_directions[query] = current;
								current = query;
								distance--;
								break;
							}
						}
					}		
					if (current == source)
						notsource = false;	
				}
			}*/

		}	
	} // end source_cells loop
	string outfilename("lcps");
	print_matrix_tofile_Int (outfilename, lcps, sizey, sizey);

	//return lcps;

}
