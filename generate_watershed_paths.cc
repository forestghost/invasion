#include <algorithm>
#include <ctime>
#include <map>
#include <cmath>
#include <random>
#include <bitset>
#include <string>
#include <iterator>
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

int main() {

	/*map<vector<Int>, vector<Int>> riverflows;
	string fdfn("flow_directions");
	ifstream flowfile;
	flowfile.open(fdfn.c_str());
	Int q,r,s,t;
	while (flowfile >> q >> r >> s >> t) 
		riverflows[{q,r}] = {s,t}; */


	string outfilename("flow_directions_watershed16");
	ofstream outfile;
	outfile.open(outfilename.c_str());

	default_random_engine e(time(0));
	default_random_engine &ernie = e;

	Int n = 720;   //// NOTE, these will change if a different preloaded landscape is used!!!!!!
	Int m = 666;

	vector<vector<Int> > M1 = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};  // Moore neighborhood (radius 1)

	// filenames of river and watershed rasters
	string riverfile("river_raster");  
	string watershedfile("watersheds_by_number");

	// holds the determined next cell for each riverine cell
	map<vector<Int>, vector<Int> > flow_directions;

	// holds the outflow points of each watershed
/*	map<Int, vector<Int>> outflows;
	outflows[21] = {113,205}; // Cannon
	outflows[20] = {151,282}; // Rush-Vermillion
	outflows[19] = {173,301}; // Zumbro
	outflows[18] = {290,412}; // Root
	outflows[17] = {353,414}; // Upper Iowa
	outflows[16] = {517,455}; // Turkey
	outflows[15] = {629,578}; // Maquoketa
	outflows[14] = {87,304}; // Red Cedar
	outflows[13] = {151,282}; // Little Chippewa
	outflows[12] = {256,400}; // Buffalo-Whitewater
	outflows[11] = {241,383}; // Trempealeau
	outflows[10] = {252,395}; // Black
	outflows[9] = {285,409}; // LaCrosse-Pine
	outflows[8] = {459,428}; // Coon-Yellow
	outflows[7] = {437,467}; // Kickapoo
	outflows[6]  = {459,428}; // Lower Wisconsin
	outflows[5] = {556,513}; // Grant-Little Maquoketa
	outflows[4] = {690,595}; // Apple-Plum
	outflows[3] = {342,696}; // Baraboo
	outflows[2] = {562,762}; // Pecatonica
	outflows[1] = {564,755}; // Sugar
originial above */ 

	map<Int, vector<Int>> outflows;
	outflows[21] = {113,106}; // Cannon
	outflows[20] = {151,183}; // Rush-Vermillion
	outflows[19] = {173,202}; // Zumbro
	outflows[18] = {290,313}; // Root
	outflows[17] = {353,315}; // Upper Iowa
	outflows[16] = {517,356}; // Turkey
	outflows[15] = {629,479}; // Maquoketa
	outflows[14] = {87,205}; // Red Cedar
	outflows[13] = {151,183}; // Little Chippewa
	outflows[12] = {256,301}; // Buffalo-Whitewater
	outflows[11] = {241,284}; // Trempealeau
	outflows[10] = {252,296}; // Black
	outflows[9] = {285,310}; // LaCrosse-Pine
	outflows[8] = {459,329}; // Coon-Yellow
	outflows[7] = {437,368}; // Kickapoo
	outflows[6]  = {459,329}; // Lower Wisconsin
	outflows[5] = {556,414}; // Grant-Little Maquoketa
	outflows[4] = {690,496}; // Apple-Plum
	outflows[3] = {342,597}; // Baraboo
	outflows[2] = {561, 662}; // Pecatonica
	outflows[1] = {564,656}; // Sugar


	MatInt rivers = read_matrix_int(riverfile, n, m);
	MatInt watersheds = read_matrix_int(watershedfile, n, m);

	// populate source_cells with river pixels
	vector<vector<Int>> source_cells;
	/*source_cells.push_back({345,590});
	source_cells.push_back({403,512}); 
	source_cells.push_back({399,503}); 
	source_cells.push_back({530,429}); 	
	source_cells.push_back({531,429}); 	
	source_cells.push_back({530,428}); 	
	source_cells.push_back({530,427}); 
	source_cells.push_back({530,426}); 
	source_cells.push_back({530,425}); 
	source_cells.push_back({529,427}); 
	source_cells.push_back({529,426}); 
	source_cells.push_back({529,425});
	source_cells.push_back({529,424}); 
	source_cells.push_back({528,423}); 
	source_cells.push_back({527,423}); 			 
*/
	for (Int a=0; a < n; a++) {
		for (Int b=0; b<m; b++) {
			if (rivers[a][b] == 0) 
				source_cells.push_back({a,b});
		}	
	}



		
	string fdfn("flow_directions_NEW");
	ifstream flowfile;
	flowfile.open(fdfn.c_str());
	Int q,r,s,t;
	while (flowfile >> q >> r >> s >> t) 
		flow_directions[{q,r}] = {s,t};
	flowfile.close();

	cout << source_cells.size() << " river pixels" << endl;	

	Int runnumber = 0;

	for (auto source : source_cells) {
		runnumber++;
		cout << "RUN NUMBER " << runnumber << " with " << flow_directions.size() << " cells determined." << endl;

		if (flow_directions.count(source) > 0) // already determined, don't need to do again
			continue;	

		Int watershed = watersheds[source[0]][source[1]];
		if (watershed != 16) // then off the map, so don't consider   // SET to ==0 for total run
			continue;

		vector<Int> destination = outflows[watershed];

		cout << "Working on cell " << source[0] << " " << source[1]  << " in watershed " << watersheds[source[0]][source[1]] << endl;

		// find distances ////
		map<vector<Int>, vector<Int>> exploration;
		exploration[source] = {0,0}; // start with source cell; first value is distance from source, second value is whether its been explored (1) or not (0)
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
		}	

	} // end source_cells loop

	outfile.close();

	return 0;

}