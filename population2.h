#ifndef POPULATION_H
#define POPULATION_H

#include "cauchy.h"

class Population {

private:
Int n;
Int m;
string rep;
string fn;
Int numwalkers;
Int numloci;
vector<Walker> walkers;
default_random_engine &e;
Landscape &l;
string filename;
vector<Int> theliving; // indices of currently living individuals in vector called "walkers"
vector<Int> thedead;
MatInt occupancy_matrix;
uniform_real_distribution<Doub> d;
Normaldist flowernum_gauss;
Normaldist pollennum_gauss;
Cauchydist seed_bicauchy;
Normaldist pollen_bigauss;
static vector<vector<Int> > M1;  // Moore neighborhood (radius 1)
uniform_real_distribution<Doub> ang;
string flowyname;
ofstream flowy;

void mate (Walker &w1, Walker &w2, Int year)
{
	Doub flow;
	if (year == flood_year)
		flow = flood_flow_parameter;
	else
		flow = flow_parameter;
	w1.reduce_flowernum();  // whether or not seed turns out to be viable
	vector<Int> boo = w1.get_position();
	// determine seed dispersal and whether or not it will land somewhere viable
	bool move_forward = false;
	vector<Int> seed_landing;
	// first find distance from parent plant using cauchy distribution
	Doub dist = seed_bicauchy.invcdf(d(e));
	// second find the angle from origin
	Doub angle = ang(e);
	seed_landing.push_back(boo[0] + dist*(cos(angle)));
	seed_landing.push_back(boo[1] + dist*(sin(angle)));
	// off the map?? 
	if (seed_landing[0] < 0 || seed_landing[0] >= n)
		return;
	if (seed_landing[1] < 0 || seed_landing[1] >= m)
		return;
	if (l.get_cell(seed_landing[0], seed_landing[1]) == 2) { //previously determined to be adjacent to water and FALLOW
		return;
	} else if (l.get_cell(seed_landing[0], seed_landing[1]) == 1) { // on land; adjacent to water?
		for (auto g : M1) {
			Int p1 = seed_landing[0]+g[0]; Int p2 = seed_landing[1]+g[1] ;
			if ( p1 <0 || p1 >= n || p2 < 0 || p2 >= m )   // is the adjacent cell on map????
				continue;
			if  (l.get_cell(p1, p2) == 0) {  // is it next to water???
				move_forward = true;
				break;
			}	
		}	
		if (move_forward && (! l.get_badland_check(seed_landing[0], seed_landing[1]) ) )  {
			l.set_badland_check(seed_landing[0], seed_landing[1]);
			if (d(e) <= badland_prob) {
				l.set_cell(seed_landing[0], seed_landing[1], 2);
				return;
			}
		}
	}	else if (l.get_cell(seed_landing[0], seed_landing[1]) == 0 ) {  // in water, float it to some adjacent land
//flowy.open(flowyname.c_str(), ios::out | ios::app);	
//flowy << seed_landing[0] << " " << seed_landing[1];

		if ( l.check_flow(seed_landing) == 0 )   // in a non-watershed stream
				return;			
		bool stillwater = true;
//Int cnt = 0;		
		while (stillwater) {

//cnt++;
//if (cnt > 20000) {
//	cout << seed_landing[0] << " " << seed_landing[1] << endl;
//}
			// new method 		
			if ( l.check_flow(seed_landing) == 0	)   // in a non-watershed stream
				return;	
			vector<Int> next_cell = l.get_flow(seed_landing);
vector<Int> orig_next = next_cell;			
			Int diff1 = next_cell[0] - seed_landing[0];
			Int diff2 = next_cell[1] - seed_landing[1];
			Doub D1 = d(e);		
			if (D1 <= flow) {
 				;   // next_cell is realized
			} else {
				vector<Int> alt_cell1;
				vector<Int> alt_cell2;
				if (diff1 !=0 && diff2 != 0) { // different on both axes
					alt_cell1 = {seed_landing[0]+diff1, seed_landing[1]};
					alt_cell2 = {seed_landing[0], seed_landing[1] + diff2};
				} else if (diff1  != 0 && diff2 == 0) { 
					alt_cell1 =  {next_cell[0], next_cell[1] + 1  } ;
					alt_cell2 =  {next_cell[0],  next_cell[1] -1  }  ;
				} else {  // second axis changes, but not first.
					alt_cell1 = {next_cell[0] +1 , next_cell[1]};
					alt_cell2 = {next_cell[0] -1, next_cell[1]};
				}
				// determine if either alt_cell is offmap and/or upstream and therefore not possible

				vector<Doub> choose_probs= {0., 0.};
				
				if (l.get_watershed(alt_cell1[0], alt_cell1[1]) > 0 ) {	// then in a watershed
					if ( l.check_flow(alt_cell1) > 0 )  {  // then alt cell is in water and must check if upstream
						if (l.get_flow(alt_cell1) == seed_landing) 
							; // then upstream, can't use
						else 
							choose_probs[0] = 0.5;
					} else  // a land pixel, so possible  
						choose_probs[0] = 0.5;	
				}		

				if (l.get_watershed(alt_cell2[0], alt_cell2[1]) > 0) { // then in a watershed
					if (l.check_flow(alt_cell1) > 0 )  {  // then alt cell is watery and must check if upstream
						if (l.get_flow(alt_cell2) == seed_landing)
							; // then upstream, can't use
						else 
							choose_probs[1] = 0.5;
					} else // a land pixel, so possible
						choose_probs[1] = 0.5;
				}

				Doub summy = choose_probs[0] + choose_probs[1];
				if (summy != 0.) {
					choose_probs[0] /= summy;
					choose_probs[1] /= summy;
					choose_probs[1] = choose_probs[0] + choose_probs[1];
				}

				Doub D2 = d(e);
				if (D2 <= choose_probs[0]) 
					next_cell = alt_cell1;
				else if ( D2 <= choose_probs[1] )
					next_cell =  alt_cell2;
				// otherwise remains to-flow cel

			}
			seed_landing = next_cell;
//flowy << " " << seed_landing[0] << " " << seed_landing[1];			

			// off the map???
			if (seed_landing[0] < 0 || seed_landing[0] > n)
				return;
			if (seed_landing[1] < 0 || seed_landing[1] > m)
				return;
			if (l.get_cell(seed_landing[0], seed_landing[1]) == 2 )  {  // bad land
//flowy << endl;
				return;
			}
			if (l.get_cell(seed_landing[0], seed_landing[1] ) == 1) {   // good land
				stillwater = false;
				move_forward = true;
				if (! l.get_badland_check(seed_landing[0], seed_landing[1]) ) {  // haven't checked before
					l.set_badland_check(seed_landing[0], seed_landing[1]);
					if (d(e) <= badland_prob) {
						l.set_cell(seed_landing[0], seed_landing[1], 2);
//flowy << endl;						
						return;
					}   // otherwise badland check passed and allow offspring plant production
				}	
			}
		}				
//flowy << endl;
//flowy.close();
	}
	if (move_forward && occupancy_matrix[seed_landing[0]][seed_landing[1]] < carrying_capacity ) { // if above 
		// check if seed will germinate, conditioned on landing in a proper area
		Doub germination_age = -1;
		for (Int y=0; y<3; y++) {
			Doub r = d(e);
			if (r < seed_germination_probs[y]) {
				germination_age = y;
				break;
			}
		}
		if (germination_age >= 0) {
			Walker baby(seed_landing, w1, w2, numloci, numwalkers, germination_age, pollennum_gauss.invcdf(d(e)), flowernum_gauss.invcdf(d(e))); // reproduction constructor
			occupancy_matrix[seed_landing[0]][seed_landing[1]]++;
			walkers.push_back(baby);  // all are living
			numwalkers++;
		}	
	}
}

public:
void testbicauchy() 
{
	ofstream dude;
	string ofile("bicauchy");
	dude.open(ofile.c_str());
	for(Int q = 0; q < 20000; q++)  {
		// first find distance from parent plant using cauchy distribution
		Doub dist = seed_bicauchy.invcdf(d(e));
		// second find the angle from origin
		uniform_real_distribution<Doub> a(0,360);
		Doub angle = a(e);
		dude << dist*(cos(angle)) << "\t" << dist*(sin(angle)) << endl;;
	}	
	dude.close();
}

void testbigauss()
{
	ofstream dudette;
	string ofile("bigauss");
	dudette.open(ofile.c_str());
	for (Int q = 0; q < 20000; q++) {
		dudette << pollen_bigauss.invcdf(d(e)) << " " << pollen_bigauss.invcdf(d(e)) << endl;
	}
	dudette.close();
}

void generation (Int year) 
{
	// reset occupancy_matrix
	MatInt o(n,m,0);
	occupancy_matrix = o;
	// germination
	ofstream yearfile;
	string filename("yearfile_");
	stringstream ss;
	ss << year;
	filename += ss.str();	
	filename += "_";
	filename += rep;
	yearfile.open(filename.c_str());
	vector<Int> female_plant_indices;
	vector<Int> male_plant_indices;
	map<vector<Int>, vector<Int> > female_plant_positions; // key is position, value is list of plant indices at that position
	for (Int i=0; i<numwalkers; i++) {
		if (walkers[i].check_germination()) {  // then seed has become plant  
			vector<Int> plantpos = walkers[i].get_position();
			vector<string> haps = walkers[i].get_haplotypes();
			walkers[i].set_death();  /// so it can be killed at the end of the reporoductive season
			if (walkers[i].get_gender() == 0  )
				male_plant_indices.push_back(i);
			else {
				female_plant_indices.push_back(i);
				female_plant_positions[walkers[i].get_position()].push_back(i);
			}
			yearfile << plantpos[0] << " " << plantpos[1] << " " << haps[0] << " " << haps[1] << endl;
		}		
	}
	yearfile.close();
	// reproduction and seed dispersal
	for (Int day = 0; day < 60 ; day++) {   // by day
		for (auto iter=male_plant_indices.begin(); iter != male_plant_indices.end(); iter++) {  // by male plant
			vector<Int> cur = walkers[*iter].get_position();
			for (Int pgrain = 0; pgrain < walkers[*iter].get_pollennum(); pgrain++) {   // by pollen grain
				vector<Int> pollenlanding;
				pollenlanding.push_back(cur[0] + pollen_bigauss.invcdf(d(e)));
				pollenlanding.push_back(cur[1] + pollen_bigauss.invcdf(d(e)));
				if (female_plant_positions.count(pollenlanding)) {// check if the pollen's landing position holds female plant
					// find eligible mate(s)
					vector<Int> eligible;
					for (auto iter2 = female_plant_positions[pollenlanding].begin(); iter2 != female_plant_positions[pollenlanding].end(); iter2++) {
						if (walkers[*iter2].get_flowernum() > 0 )
							eligible.push_back(*iter2); // opposite gender and still has unmated flowers
						else
							(iter2 = female_plant_positions[pollenlanding].erase(iter2))--;
						if (female_plant_positions[pollenlanding].size() == 0) { // then get rid of this female plant position
							female_plant_positions.erase(pollenlanding);
							break;
						}	
					}
					Int mate_index;					
					Doub nummy = static_cast<Doub>(eligible.size());
					if (nummy > 1) {	
						Doub ran = d(e);
						Doub cumulative = 0.;
						for (Int q : eligible) {
							cumulative += 1./nummy;
							if (ran <= cumulative) {
								mate_index = q;
								break;
							}
						}
					} else if (nummy == 1) {
						mate_index = eligible[0];
					} else {
						continue;   // no mate so move on
					}
					// mating
					mate(walkers[mate_index], walkers[*iter], year);  // so first walker is always the female plant
				}
			}	
		}
		Int viable_flowers = 0;
		for (auto iter=female_plant_indices.begin(); iter!= female_plant_indices.end(); iter++) {
			viable_flowers += walkers[*iter].get_flowernum();
		}
		cout << "year: " << year << ", day: " << day << ", numwalkers: " << numwalkers << ", viable flowers: " << viable_flowers << endl;
		if (viable_flowers == 0)
			break;		
	}
	// winter death and age advancement of seeds
	 		// kill off plants that grew this year // add year to non-germinated seeds
	cout << "beginning death sequence ... " << endl;
	cout << "before: " << walkers.size();	
	bool stop = false;
	vector<Walker>::iterator finaliter;
	Int count = 0;
	for (auto iter=walkers.begin(); iter != walkers.end(); iter++) {
		if (!stop) {
			count++;
			if (!(*iter).get_death()) {
				finaliter = iter;
				stop = true;
			}	
		}
		(*iter).add_year();
	}	
	// eliminate the block of walkers past their prime
	walkers.erase(walkers.begin(), finaliter);
	cout << ", after: " << walkers.size() << endl;
	cout << "elimination count: " << count << endl;
	if  (walkers.size() == 0) 
		throw std::runtime_error("EXTINCTION of populuation!\n");
	numwalkers = walkers.size();
}

void print_occupancy_matrix (string occupy)
{
	ofstream ostt;
	ostt.open(occupy.c_str());
	for (Int i = 0; i<n; i++) {
		for (Int j=0; j<m-1; j++)
			ostt << occupancy_matrix[i][j] << " ";
		ostt << occupancy_matrix[i][m-1] << endl;
	}
	ostt.close();
}

vector<Doub> get_frequencies()  //frequencies of the derived allele are returned
{
	vector<Doub> freqs(numloci, 0.);
	for(Int j : theliving) {
		for (Int k = 0; k<2; k++) {
			for (Int l = 0; l<numloci; l++)
				freqs[l] += static_cast<Doub>(walkers[j].get_allele(l,k));
		}
	}
	for (Int l=0; l<numloci; l++)
		freqs[l] /= (2*theliving.size());
	return freqs;
}

inline Int get_numwalkers() { return numwalkers; }

Population (Int nloci, Int numrow, Int numcol, default_random_engine &engine, Landscape &landy, string filename, string rr): numloci(nloci), n(numrow), m(numcol), e(engine), l(landy), fn(filename), rep(rr) {


	//starting allele frequencies for each locus
	uniform_real_distribution<Doub> d(0,1);

	MatInt o(n,m,0);
	occupancy_matrix = o;

//flowyname = "flowlines";

	// create probabilitiy distributions
	Normaldist dummy(pollen_dispersal[0], pollen_dispersal[1]);
	pollen_bigauss = dummy;
	Cauchydist dummy2(seed_dispersal[0], seed_dispersal[1]);
	seed_bicauchy = dummy2;
	Normaldist dummy3(pollen_number_per_day[0], pollen_number_per_day[1]);
	pollennum_gauss = dummy3;
	Normaldist dummy4(flower_number[0], flower_number[1]);
	flowernum_gauss = dummy4;
	uniform_real_distribution<Doub> a(0,360);
	ang = a;

	// create walkers
	walkers.reserve(100000);
	map<vector<Int>, Int> seed_positions;
	map<vector<Int>, vector<Int>> means;
	map<vector<Int>, vector<Doub>> sds;
	   // create starting seeds
	ifstream seedfile;
	seedfile.open(fn.c_str());
	Int flipper = 0;
	vector<Int> coords;
	string line;

	Int q,r,s, m1, m2, m3, m4, m5, cp, mt;  // where the m's are the mean allele sizes of the three SSR loci
	Doub sd1, sd2, sd3, sd4, sd5;   // where these are the std devs on the three SSR loci
	while (seedfile >> q >> r >> s >> mt >> cp >> m1 >> sd1 >> m2 >> sd2 >> m3 >> sd3 >> m4 >> sd4 >> m5 >> sd5) {
		seed_positions[{q,r}] = s;
		means[{q,r}] = {mt, cp, m1, m2, m3, m4, m5};
		sds[{q,r}]  = {sd1, sd2, sd3, sd4, sd5};
		l.set_badland_check(q,r); // plants already growing there, so declare permanent goodland
	}	
	// create initial seeds;
	numwalkers = 0;
	for (auto iter = seed_positions.begin(); iter != seed_positions.end() ; iter++) { 
		normal_distribution<Doub> ssr1(means[iter->first][2], sds[iter->first][0]);
		normal_distribution<Doub> ssr2(means[iter->first][3], sds[iter->first][1]);
		normal_distribution<Doub> ssr3(means[iter->first][4], sds[iter->first][2]);
		normal_distribution<Doub> ssr4(means[iter->first][5], sds[iter->first][3]);
		normal_distribution<Doub> ssr5(means[iter->first][6], sds[iter->first][4]);
		for (Int y=0; y<iter->second; y++) {
			//determine germination age if any and create walker if so
			Doub germination_age = -1;
			for (Int y=0; y<3; y++) {
				Doub r = d(e);
				if (r < seed_germination_probs[y]) {
					germination_age = y;
					break;
				}
			}
			vector<Int> genetics1 = { means[iter->first][0] , means[iter->first][1] };
			vector<Int> genetics2 = { means[iter->first][0] , means[iter->first][1] };

			// determine STR alleles
			Doub str1 = ssr1(e);
			Doub str2 = ssr2(e);
			Doub str3 = ssr3(e);
			Doub str4 = ssr4(e);
			Doub str5 = ssr5(e);
			genetics1.push_back(str1-floor(str1)>0.5?ceil(str1):floor(str1));
			genetics1.push_back(str2-floor(str2)>0.5?ceil(str2):floor(str2));
			genetics1.push_back(str3-floor(str3)>0.5?ceil(str3):floor(str3));
			genetics1.push_back(str4-floor(str4)>0.5?ceil(str4):floor(str4));
			genetics1.push_back(str5-floor(str5)>0.5?ceil(str5):floor(str5));

			str1 = ssr1(e);
			str2 = ssr2(e);
			str3 = ssr3(e);
			str4 = ssr4(e);
			str5 = ssr5(e);
			genetics2.push_back(str1-floor(str1)>0.5?ceil(str1):floor(str1));
			genetics2.push_back(str2-floor(str2)>0.5?ceil(str2):floor(str2));
			genetics2.push_back(str3-floor(str3)>0.5?ceil(str3):floor(str3));
			genetics2.push_back(str4-floor(str4)>0.5?ceil(str4):floor(str4));
			genetics2.push_back(str5-floor(str5)>0.5?ceil(str5):floor(str5));

			map<Int, vector<Int>> genetics;
			genetics[0] = genetics1;
			genetics[1] = genetics2;

			if (germination_age >= 0) {
				Walker w(iter->first, genetics, numwalkers, germination_age, pollennum_gauss.invcdf(d(e)), flowernum_gauss.invcdf(d(e)));
				walkers.push_back(w); 
				vector<Int> p = walkers[numwalkers].get_position();
				occupancy_matrix[p[0]][p[1]]++;
				numwalkers++;
			}
		}

	}
	ofstream outtie;
	string ou("seedy0");
	outtie.open(ou.c_str());
	for (Int p=0; p<numwalkers;p++) {
		vector<Int> c = walkers[p].get_position();
		outtie << static_cast<Doub>(c[0])/numrow << "\t" << static_cast<Doub>(c[1])/numcol << endl;
	}
	outtie.close();


}

static vector<Doub> pollen_dispersal;
static vector<Doub> seed_dispersal;
static vector<Doub> pollen_number_per_day;
static vector<Doub> flower_number;
static vector<Doub> seed_germination_probs;
static Int carrying_capacity;
static Doub badland_prob;
static Doub flow_parameter;
static Int flood_year;
static Doub flood_flow_parameter;

};

#endif
