#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

struct Sparsematrix
{
	Int nrows;
	Int ncols;
	Int nvals;
	VecInt col_ptr;
	VecInt row_ind;
	VecDoub val;

	Sparsematrix(Int m, Int n, Int nnvals) : nrows(m), ncols(n), nvals(nnvals), col_ptr(n+1,0), row_ind(nnvals,0), val(nnvals,0.) {}
	

	VecDoub ax(const VecDoub &x) const{
		VecDoub y(nrows,0.);
		for (Int j=0; j<ncols; j++) {
			if (col_ptr[j] != 999999) { // then there's at least one nnz in this col
				Int next_index_that_matters;
				for (Int q=j+1; q<ncols+1; q++) {
					if (col_ptr[q] != 999999) {
						next_index_that_matters = q;
						break;
					}
				}
				for (Int i=col_ptr[j]; i<col_ptr[next_index_that_matters]; i++) 
					y[row_ind[i]] += val[i]*x[j];
			}
		}
		return y;
	} 

};



#endif
