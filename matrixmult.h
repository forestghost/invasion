#ifndef MATRIXMULT_H
#define MATRIXMULT_H



VecDoub matrix_by_vector (MatDoub matty, VecDoub vecky, Int n);
MatDoub scalar_by_matrix (Doub scalar, MatDoub matty, Int n, Int m);
Doub scalar_product (MatDoub matty1, MatDoub matty2, Int dim1, Int dim4); // the Matrices may be VEctors
	// dim1 is row dimension of the first matrix (or vector)
	// dim4 is column dimension of second matrix (or vector)
MatDoub matrix_product (MatDoub matty1, MatDoub matty2, Int dim1, Int dim4);  // MAY WORK for any combination of vector and matrix

void print_matrix (MatDoub matty, Int n, Int m);
void print_vector (VecDoub vecky, Int n);
void print_vector_tofile(string filename, VecDoub vecky, Int n, int gen);
void print_matrix_tofile (string filename, MatDoub matty, Int n, Int m);
void print_matrix_tofile_Int (string filename, MatInt matty, Int n, Int m);
void print_map_tofile (string filename, map<Doub, Doub> mappy);
void clear_file(string filename);

MatDoub create_emptymatrix(Int n, Int m);
MatDoub read_matrix(string filename, Int n, Int m);
//MatDoub_I read_matrix(string filename, Int n, Int m);
VecDoub create_emptyvector(Int n);
VecDoub read_vector(string filename, Int n);

VecDoub generate_marginal(VecDoub freqs, MatDoub geno_fitness, Int n, Int m);
MatDoub generate_gmatrix(VecDoub freqs, Int n, Int m);

MatInt read_matrix_int(string filename, Int n, Int m);
MatInt create_emptymatrix_int(Int n, Int m);
VecInt read_vector_int(string filename, Int n);

Doub calc_mean(VecDoub freqs, VecDoub marginal_fitness, Int n);

VecDoub matrix_by_vector(MatDoub matty, VecDoub vecky, Int n) {  // for square matrices only 
	
	VecDoub return_vecky(n);


	for (Int i=0; i<n; ++i) {
		return_vecky[i] = 0;
		for (Int j=0; j<n; ++j) {
			return_vecky[i] += matty[i][j]*vecky[j];
		}
	}

	return return_vecky;


}

MatInt create_emptymatrix_int(Int n, Int m) {
 
	MatInt return_matty(n,m);
	for (Int i=0; i<n; ++i) 
		for (Int j = 0; j<m; j++) 
			return_matty[i][j] = 0;
	return return_matty;

}

MatDoub create_emptymatrix(Int n, Int m) {
 
	MatDoub return_matty(n,m);
	for (Int i=0; i<n; ++i) 
		for (Int j = 0; j<m; j++) 
			return_matty[i][j] = 0.;
	return return_matty;

}


VecDoub create_emptyvector(Int n) {

	VecDoub return_vecky(n);
	for (Int j=0; j<n; j++) 
		return_vecky[j] = 0.;
	return return_vecky;

}


MatDoub scalar_by_matrix (Doub scalar, MatDoub matty, Int n, Int m) {

	for (Int i = 0; i<n; ++i) 
		for (Int j=0; j<m; ++j) 
			matty[i][j] *= scalar;
	
	return matty;

}

void print_matrix (MatDoub matty, Int n, Int m) {
	
	for (Int i = 0; i<n; ++i) {
		for (Int j=0; j<m; ++j) 
			cout << matty[i][j] << "\t";
		cout << endl;
	}

}

void print_matrix_tofile (string filename, MatDoub matty, Int n, Int m) {
	ofstream outstar;
	outstar.open(filename.c_str());

	for (Int i = 0; i<n; ++i) {
		for (Int j=0; j<m; ++j) 
			outstar << matty[i][j] << "\t";
		outstar << endl;
	}

	outstar.close();
	
}

void print_matrix_tofile_Int (string filename, MatInt matty, Int n, Int m) {
	ofstream outstar;
	outstar.open(filename.c_str());

	for (Int i = 0; i<n; ++i) {
		for (Int j=0; j<m; ++j) 
			outstar << matty[i][j] << "\t";
		outstar << endl;
	}

	outstar.close();
	
}


void print_vector (VecDoub vecky, Int n) {

	for (int i =0; i<n; ++i) 
		cout << vecky[i] << "\t";

	cout << endl;

}

void clear_file(string filename) {
	ofstream outstar;
	outstar.open(filename.c_str());
	outstar.close();
}

void print_vector_tofile(string filename, VecDoub vecky, Int n, int gen) {
	
	ofstream outstar;
	//outstar.open(filename.c_str());
	outstar.open(filename.c_str(), ios::app);
	outstar << gen << "\t";
	for (int d=0; d<n; d++) 
		outstar << vecky[d] << "\t";
	outstar << endl;
	outstar.close();
}	

void print_map_tofile (string filename, map<Doub, Doub> mappy) {
	ofstream outstar;
	outstar.open(filename.c_str());
	for (map<Doub, Doub>::iterator iter = mappy.begin(); iter != mappy.end(); ++iter) {
		outstar << iter->first << "\t" << iter->second << endl;
	}
	outstar.close();
}


MatDoub read_matrix(string filename, Int n, Int m) {
	ifstream infile;
	infile.open(filename.c_str()); 
	Int total = 0;
	Int colcounter = 0;
	Int rowcounter = 0;

	MatDoub matty(n,m);

	while (total < n*m) {
		++total;
		Doub datum;
		string dummy; 
		infile >> dummy;
		istringstream ss(dummy);
		ss >> datum;
		matty[rowcounter][colcounter] = datum;
		++colcounter;
		if (colcounter == m) {
			colcounter = 0;
			++rowcounter;
		}
	}


	infile.close();
	return matty;
}

MatInt read_matrix_int(string filename, Int n, Int m) {
	ifstream infile;
	infile.open(filename.c_str()); 
	Int total = 0;
	Int colcounter = 0;
	Int rowcounter = 0;

	MatInt matty(n,m);

	while (total < n*m) {
		++total;
		Int datum;
		string dummy; 
		infile >> dummy;
		matty[rowcounter][colcounter] = atoi(dummy.c_str());
		++colcounter;
		if (colcounter == m) {
			colcounter = 0;
			++rowcounter;
		}
	}
	infile.close();
	return matty;
}

VecDoub read_vector(string filename, Int n) {

	ifstream infile;
	infile.open(filename.c_str());
	Int total = 0;
	Int counter = 0;
	VecDoub vecky(n);

	while (total < n) {

		++total;
		Doub datum;
		string dummy;
		infile >> dummy;
		istringstream ss(dummy);
		ss >> datum;
		vecky[counter] = datum;
		++counter;
	}

	infile.close();

	return vecky;

}

VecInt read_vector_int(string filename, Int n) {

	ifstream infile;
	infile.open(filename.c_str());
	Int total = 0;
	Int counter = 0;
	VecInt vecky(n);

	while (total < n) {

		++total;
		Int datum;
		string dummy;
		infile >> dummy;
		istringstream ss(dummy);
		ss >> datum;
		vecky[counter] = datum;
		++counter;
	}

	infile.close();

	return vecky;

}

VecDoub generate_marginal(VecDoub freqs, MatDoub geno_fitness, Int n, Int m) {
	VecDoub margey(n);
	for (Int i=0; i<n; ++i) {
		Doub marge = 0;
		for (Int j=0; j<n; ++j) {
			marge += freqs[j]*geno_fitness[i][j];		
		}	
		margey[i] = marge;
	}
	return margey;
}

MatDoub generate_gmatrix(VecDoub freqs, Int n, Int m) {
	MatDoub matcity(n,m);
	for (Int i=0; i<n; ++i) {
		for (Int j=0; j<n; ++j) {
			if (i == j) 
				matcity[i][j] = freqs[i]*(1-freqs[i]);
			else 
				matcity[i][j] = -1*freqs[i]*freqs[j];
		}	
	}	
	return matcity;
}

Doub calc_mean(VecDoub freqs, VecDoub marginal_fitness, Int n) {
	Doub mean = 0;
	for (Int i=0; i<n; ++i) 
		mean += freqs[i]*marginal_fitness[i];
	return mean;
}


MatDoub matrix_product (MatDoub matty1, MatDoub matty2, Int dim1, Int dim4) {
	MatDoub outty = create_emptymatrix(dim1, dim4);
	for (Int i=0; i<dim1; i++) {
		for (Int j=0; j<dim4; j++) {
			for (Int k=0; k<dim4; k++) { 
				outty[i][j] += matty1[i][k]*matty2[k][j];
			}
		}		
	}
	return outty;
}

Doub scalar_product (MatDoub matty1, MatDoub matty2, Int dim1, Int dim4) {
	Doub outty = 0.;
	for (Int i=0; i<dim1; i++) {
		for (Int j=0; j<dim4; j++) 
			outty += matty1[i][j]*matty2[j][i];
	}
	return outty;
}


#endif
