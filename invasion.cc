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

	default_random_engine e(time(0));
	default_random_engine &ernie = e;

	// READ parameters  /////////////
	Int numrow, numcol, numloci, numyears, reps;
	Doub goodland_cost, badland_cost, badland_prob, germination_probability;
	stringstream ss(stringstream::in | stringstream::out);
	cout << "reading parmameter values ... " << endl;
	map<string,string> parameters;
	string line;
	string infile("parameters");
	ifstream paramfile;
	paramfile.open(infile.c_str());
	while (paramfile >> line) {
		string par = line;
		string value;
		if (par == "POLLEN_DISPERSAL") {
			bool good = 1;
			while (good) {
				paramfile >> line;
				if (line == "//") 
					good  = 0;
				if (good) 
					Population::pollen_dispersal.push_back(atof(line.c_str()));
			}	
		} else if (par == "SEED_GERMINATION_PROBS") {
			bool good = 1;
			while (good) {
				paramfile >> line;
				if (line == "//")
					good = 0;
				if (good) 
					Population::seed_germination_probs.push_back(atof(line.c_str()));
			}
		} else if (par == "SEED_DISPERSAL") {
			bool good = 1;
			while (good) {
				paramfile >> line;
				if (line == "//") 
					good  = 0;
				if (good) 
					Population::seed_dispersal.push_back(atof(line.c_str()));
			}	
		} 	else if (par == "POLLEN_NUMBER_PER_DAY") {
			bool good = 1;
			while (good) {
				paramfile >> line;
				if (line == "//") 
					good  = 0;
				if (good)
					Population::pollen_number_per_day.push_back(atof(line.c_str()));
			}	
		}	else if (par == "FLOWER_NUMBER") {
			bool good = 1;
			while (good) {
				paramfile >> line;
				if (line == "//") 
					good  = 0;
				if (good) 
					Population::flower_number.push_back(atof(line.c_str()));
			}	
		} else 	{
			paramfile >> line;
			parameters[par] = line;
		}
	}
	paramfile.close();

	string landscape_filename =  parameters["LANDSCAPE_RASTER_FILE"];
	string beginning_seeds_filename = parameters["BEGINNING_SEEDS_FILE"];
	string watersheds_filename = parameters["WATERSHEDS_RASTER_FILE"];
	string flowdirections_filename = parameters["FLOW_DIRECTIONS_FILE"];
	numrow = 720;   //// NOTE, these will change if a different preloaded landscape is used!!!!!!
	numcol = 666;
	cout << "landscape_file: " << landscape_filename << endl;
	cout << "landscape size: " << numrow << "x" << numcol << endl;
	numyears = atoi(parameters["NUM_YEARS"].c_str());
	cout << "years to be simulated: " << numyears << endl;
	numloci  = atoi(parameters["NUM_LOCI"].c_str());
	cout << "loci number: " << numloci << endl;
	cout << "beginning seeds file: " << beginning_seeds_filename << endl;
	goodland_cost = atof(parameters["GOODLAND_COST"].c_str());
	badland_cost = atof(parameters["BADLAND_COST"].c_str());
	reps = atoi(parameters["REPLICATES"].c_str());
	Population::badland_prob = atof(parameters["BADLAND_PROB"].c_str());
	Population::flow_parameter = atof(parameters["FLOW_PARAMETER"].c_str());
	Population::flood_year = atoi(parameters["FLOOD_YEAR"].c_str());
	Population::flood_flow_parameter = atof(parameters["FLOOD_FLOW_PARAMETER"].c_str());
	Walker::costs.push_back(goodland_cost);
	Walker::costs.push_back(badland_cost);
	Landscape::move_costs.push_back(goodland_cost);
	Landscape::move_costs.push_back(badland_cost);
	cout << "goodland cost: " << goodland_cost << endl;
	cout << "badland cost: " << badland_cost << endl;
	Population::carrying_capacity = atoi(parameters["CARRYING_CAPACITY"].c_str());
	cout << "carrying capacity: " << Population::carrying_capacity << endl;


	cout << "seed germination probabilities (years 0,1,2): ";
	copy(Population::seed_germination_probs.begin(), Population::seed_germination_probs.end(), ostream_iterator<Doub>(cout, "\t")); 
	cout << endl;
	cout << "pollen dispersal parameters (mean, sd): ";
	copy(Population::pollen_dispersal.begin(), Population::pollen_dispersal.end(), ostream_iterator<Doub>(cout, "\t")); 
	cout << endl;
	cout << "seed dispersal parameters (mean, variance): ";
	copy(Population::seed_dispersal.begin(), Population::seed_dispersal.end(), ostream_iterator<Doub>(cout, "\t"));
	cout << endl;
	cout << "pollen number parameters (mean, variance): ";
	copy(Population::pollen_number_per_day.begin(), Population::pollen_number_per_day.end(), ostream_iterator<Doub>(cout, "\t"));
	cout << endl;
	cout << "flower number parameters (mean, variance): ";
	copy(Population::flower_number.begin(), Population::flower_number.end(), ostream_iterator<Doub>(cout, "\t"));
	cout << endl;


	for (Int r = 0; r<reps; r++) {  // replicate loop
		// create landscape
		Landscape rivers(numrow, numcol, landscape_filename, watersheds_filename, flowdirections_filename);
		Landscape &landie = rivers;

		stringstream dd;
		dd << r;
		string reppy = dd.str();

		// create population from starting seeds
		Population pop (numloci, numrow, numcol, e, landie, beginning_seeds_filename, reppy);
		pop.testbigauss();
		pop.testbicauchy();
		for (Int i = 0; i <= numyears; i++) {
			cout << i << ": " << pop.get_numwalkers() << endl;
			pop.generation(i);
			string ofie("occupancy_");
			string underscore("_");
			stringstream qqq;
			qqq << i;
			ofie += qqq.str();	
			ofie += underscore;
			ofie += reppy;
			pop.print_occupancy_matrix(ofie);
			landie.print_features_count();
		}
	}	
	return 0;	
}

// initialize static Landscape class variables 
vector<Doub> Landscape::move_costs;
vector<vector<Int> > Landscape::M1 = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};  // Moore neighborhood (radius 1)


// initialize static Population class variables
Int Population::carrying_capacity;
Doub Population::badland_prob;
Doub Population::flow_parameter;
Int Population::flood_year;
Doub Population::flood_flow_parameter;
vector<Doub> Population::pollen_dispersal;
vector<Doub> Population::seed_dispersal;
vector<Doub> Population::pollen_number_per_day;
vector<Doub> Population::flower_number;
vector<Doub> Population::seed_germination_probs;
vector<vector<Int> > Population::M1 = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};  // Moore neighborhood (radius 1)

// initialize static Walker class variables
vector<Doub> Walker::costs;


