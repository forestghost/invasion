require(Matrix);
d = read.table(file = "river_raster");

build_pathmatrix <- function()
{
	pm = Matrix(0, nrow= 20000, ncol = 20000, sparse = T);
	#pm = matrix(NA, nrow  = 20000, ncol = 20000);
	apply(pm, c(1,2), function (x) if (x == 0 ) x=NA);
	return(pm);
}

