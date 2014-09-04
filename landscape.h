#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include "matrixmult.h"

class Landscape {

private:
	const Int n; // number of rows in landscape matrix
	const Int m; // number of cols in landscape matrix
	string fn; // filename of input landscape raster
	string wfn; // watershed by number filename
	string fdfn; // flow directions filename
	MatInt	land; // the landscape matrix
	MatInt watersheds;  // the matrix of watersheds (0 equals "off map")
	map<vector<Int>, vector<Int>> riverflows; // next cell downsteream in river (value) from current cell (key)
	MatInt badland_check; // used in conjunction with badland_prob checks. 
	//default_random_engine &e;


public:

	void print_features_count()
	{
		map<Int, Int> mappy;
		for (Int i=0; i<n ; i++) 
			for (Int j=0; j<m; j++) 
				mappy[land[i][j]]++;
		cout << "LANDSCAPE FEATURES DISTRIBUTION ********" << endl;
		for (auto feature : mappy)    // new C++11 range loop
			cout << feature.first << ": " << feature.second << endl;
	}

	Doub least_cost_path(vector<Int> source, vector<Int> destination)  // coordinates of starting and ending positions
																	// utilizes a form of Dijsktra's algorithm
													// see http://www.geeksforgeeks.org/greedy-algorithms-set-6-dijkstras-shortest-path-algorithm/
	{
		MatDoub dist(n,m,INFINITY);
		map<vector<Int>, vector<Int>> previous; 
		MatInt sptSET(n,m,0);  // holds whether least cost path to cell has been determined or not (0 for not determined)
		//sptSET[source[0]][source[1]] = 1;
		dist[source[0]][source[1]] = 0.; // zero distance from source
		vector<Int> current;
		bool destination_notyet = true;
		
		while (destination_notyet) {

			vector<Int> next = {-1,-1};
			Doub lowdist = INFINITY;
			for (Int i =0; i<n; i++) 
				for (Int j=0; j<m; j++) 
					if (dist[i][j] < lowdist && sptSET[i][j] == 0) {
						lowdist = dist[i][j];
						next = {i,j};
					}
			current = {next[0], next[1]};
			sptSET[current[0]][current[1]] = 1;
			for (auto g : M1) {  // determine istances for adjacent (by first order Moore neighborhood) cells
				Int query1 = current[0] + g[0];
				Int query2 = current[1] + g[1];
				vector<Int> query;
				query.push_back(query1);
				query.push_back(query2);
				//vector<Int> query = toroidal_correction2(query1, query2);
				//dist[query[0]][query[1]] = dist[current[0]][current[1]] + move_costs[ land[query[0]][query[1]] ];
				
				if (dist[query[0]][query[1]] > dist[current[0]][current[1]] + move_costs[ land[query[0]][query[1]] ] ) {
					dist[query[0]][query[1]] = dist[current[0]][current[1]] + move_costs[ land[query[0]][query[1]] ];
					previous[ query ] = {current[0], current[1]};
				}
			}

			if (current == destination)
				destination_notyet = false;
		}
		vector<Int> q = {destination[0], destination[1]};
		while (q != source) {
			q = previous[q];
			cout << q[0] << ", " << q[1] << endl;
		}
		return dist[current[0]][current[1]];
	}

	inline MatInt get_land() {return land;}
	inline Int get_cell(Int x, Int y) {return land[x][y];}
	inline void set_cell(Int x, Int y, Int val) { land[x][y] = val; }
	inline void print_land(string filey) {print_matrix_tofile_Int(filey, land, n, m);}
	inline Int get_badland_check(Int x, Int y) { return badland_check[x][y]; }
	inline void set_badland_check(Int x, Int y) { badland_check[x][y] = 1;}
	inline vector<Int> get_flow(vector<Int> sourcecell) {return riverflows[sourcecell];}
	inline Int get_flow_size() {return riverflows.size();}
	inline Int check_flow(vector<Int> sourcecell) { return riverflows.count(sourcecell);}
	inline Int get_watershed(Int i, Int j) { return watersheds[i][j];}

	Landscape (const Int nn, const Int mm, string ffn, string wwfn, string ffdfn): n(nn), m(mm), fn(ffn), wfn(wwfn), fdfn(ffdfn) {  // constructor
		land = read_matrix_int(fn, n, m);   // will be 0 for river, 1 for good land, and 2 for bad land
		watersheds = read_matrix_int(wfn, n, m);
		ifstream flowfile;
		flowfile.open(fdfn.c_str());
		Int q,r,s,t;
		while (flowfile >> q >> r >> s >> t) 
			riverflows[{q,r}] = {s,t};
		flowfile.close();
		string filey("landscape");
		print_matrix_tofile_Int(filey, land, n, m);
		badland_check = create_emptymatrix_int(n,m);  // starts as all zeroes, will be set to 1 if tested for badland
	}


	static vector<Doub> move_costs;
	static vector<vector<Int> > M1;  // Moore neighborhood (radius 1)


};


#endif