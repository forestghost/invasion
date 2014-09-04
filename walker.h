#ifndef WALKER_H
#define WALKER_H

#include "matrixmult.h"

class Walker {

private:
default_random_engine e;
vector<Int> pos;  // pos is current position of the walker (in cell number)
vector<Int> genotype;
map<Int, vector<Int>> haplotypes;
Int age; // age in years
Int germination_age; // the age at which the seed will germinate
bool dead; // true if dead
Int numloci;
Int flowernum;
Int pollennum;
uniform_real_distribution<Doub> u;
Int gender; // 0 for male, 1 for female
Int s;


public:
static vector<Doub> costs;

vector<string> get_haplotypes()
{
	vector<string> haps;

	// first haplotype
	stringstream ss;

	for(vector<Int>::iterator iter = (haplotypes[0]).begin(); iter != (haplotypes[0]).end()-1; iter++)
		ss << *iter << " ";
	ss << haplotypes[0].back();

	string hap = ss.str();

	haps.push_back(hap);

	ss.str(string()); // clear stringstream variable


	// second haplotype
	for(vector<Int>::iterator iter = (haplotypes[1]).begin();  iter != (haplotypes[1]).end()-1; iter++)
		ss << *iter << " ";
	ss << haplotypes[1].back();
	hap = ss.str();

	haps.push_back(hap);

	return (haps);
}

bool check_germination() {
	if (age == germination_age)
		return true;
	else 
		return false;
}

inline vector<Int> get_position( ) {return pos;} // returns cell on the landscape (row, col)
inline Int get_allele(Int locnum, Int hapnum) {return haplotypes[hapnum][locnum];}
inline Int get_pollennum() { return pollennum; }
inline Int get_gender() { return gender; }
inline Int get_flowernum() { return flowernum; }
inline Int reduce_flowernum() { flowernum -= 1; }
inline Int set_death() { dead = true ; }
inline Int add_year() { age += 1; }
inline Int get_death() { return dead; }

Walker () = default;
// constructor for initialization of the population
Walker (vector<Int> position, map<Int, vector<Int>> ggenetics, Int seedaugment, Int germage, Doub polnum, Doub flownum):pos(position), haplotypes(ggenetics), s(seedaugment), germination_age(germage), pollennum(polnum), flowernum(flownum) { // constructor for random starting position
	default_random_engine e(time(0)+s);
	uniform_real_distribution<Doub> u(0,1);
	age = 0;
	dead = false;
	if (u(e) < 0.5) {
		gender = 0;  // "male" seed
		flowernum = 0;
	} 	else  {
		gender = 1;  // "female" seed
		pollennum = 0;
	}	
			
}

// constructor for intrasimulation walker birth
Walker (vector<Int> ppos, Walker &w1, Walker &w2, Int nloci, Int seedaugment, Int germage, Doub polnum, Doub flownum): pos(ppos), numloci(nloci), s(seedaugment), germination_age(germage), pollennum(polnum), flowernum(flownum) {  
	default_random_engine e(time(0) + s);
	uniform_real_distribution<Doub> u(0,1); 
	age = -1;
	dead = false;
		if (u(e) < 0.5) {
		gender = 0;  // "male" seed
		flowernum = 0;
	} 	else  {
		gender = 1;  // "female" seed
		pollennum = 0;
	}
	// determine haplotypes of offspring
	vector<Int> b1;
	vector<Int> b2;
	//mtDNA and cpDNA
	if (w1.get_gender() == 1) {  // then w1 is the mother
		b1.push_back(w1.haplotypes[0][0]);
		b2.push_back(w1.haplotypes[0][0]);
		b1.push_back(w1.haplotypes[0][1]);
		b2.push_back(w1.haplotypes[0][1]);
	} else {
		cout << "ERROR: male listed as first parent" << endl;
	}
	//SSRs   // for now simple SMM for mutation
	for (Int i=2; i<7; i++) {
		// copy from first parent
		Int allele_size;
		if (u(e) < 0.5) {
			allele_size = w1.haplotypes[0][i];
			// check for mutation
			if (u(e) <= 0.0001) {
cout << "MUTEY!" << endl;				
				if (u(e) < 0.5)
					allele_size++;
				else 
					allele_size--;
			}	
		} else { 
			allele_size = w1.haplotypes[1][i];
			// check for mutation
			if (u(e) <= 0.0001) {
cout << "MUTEY!" << endl;				
				if (u(e) < 0.5)
					allele_size++;
				else 
					allele_size--;
			}
		}
		b1.push_back(allele_size);
		// copy from second parent
		if (u(e) < 0.5) { 
			allele_size = w2.haplotypes[0][i];
			// check for mutation
			if (u(e) <= 0.0001) {
cout << "MUTEY!" << endl;				
				if (u(e) < 0.5)
					allele_size++;
				else 
					allele_size--;
			} 
		} else  {
			allele_size = w2.haplotypes[1][i];
			// check for mutation
			if (u(e) <= 0.0001) {
cout << "MUTEY!" << endl;				
				if (u(e) < 0.5)
					allele_size++;
				else 
					allele_size--;
			}
		}	
		b2.push_back(allele_size);
	}

	haplotypes[0] = b1;
	haplotypes[1] = b2;
}

};

#endif