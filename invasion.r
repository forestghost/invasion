require(raster);
require(RColorBrewer);
require(MASS);
dyn.load("getlcps.so");
print ("reading raster file ... ");
ras = raster(nrow = 720, ncol = 666, xmn = 1, xmx = 667, ymn = -721, ymx = -1);
rr = read.table(file = "river_raster");
cmp = read.table(file = "corrected_mississippi_flows", header = F, sep = "");
print ("loading rivers ... ");

#################### onetime function
read_river_points <- function(fn)
{
	a = read.table(file = fn, header = F, sep = "");
	nrow = length(a[,1]);
	ncol = length(a[1,]);
	c = matrix(nrow = 49860, ncol = 2);
	current = 1;
	for (i in 1:nrow) {
		for (j in 1:ncol) {
			if (a[i,j] == 0) {
				c[current,1] = j + 0.5;
				c[current,2] = (i*-1) - 0.5;
				current = current+1;
			}	
		}
	}
	return(c);
}

rp = read_river_points("river_raster");
values(ras) = as.matrix(rr);

convert_points <- function (year, low, high, numloci) #year, lowest replicate number, highest replicate number
{
	elements = high-low+1;
	listy = vector("list", elements);
	for (i in low:high) {
		print (i);	
		filename = paste("yearfile_", year, "_", i, sep = "");
		d = read.table(file = filename, header = F, sep = "");#, colClasses = c("numeric", "numeric", "character", "character") );
		d$V1 = ( d$V1 * -1) - 0.5 -1;  # 0.5 factor to center point, -1 to correct for zero-based indexing in c++
		d$V2 = d$V2 + 0.5 + 1;  # 0.5 factor center point , -1 to correct for zero-based indexing in C++ vs. 1-based indexing in R

		l = length(d$V1);
		z = as.matrix(d);
		listy[[i+1]] = z;
	}
	return(listy);
}

zoomy <- function(a, r)   # a is list of converted plant tables, r is the replicate desired
{
	zoom(ras, col = c("black", "wheat2"));
	points(a[[r]][,2], a[[r]][,1], col = "magenta", pch = ".", cex = 3);

}

read_occupancies <- function (year, low, high) #year, lowest replicate number, highest replicate number
{
	elements = high-low+1;
	listy = vector("list", elements);
	for (i in low:high) {
		print(i);
		filename = paste("occupancy_",year,"_",i, sep = "");
		z = as.matrix(read.table(file = filename, header =F, sep = ""));
		listy[[i+1]] = z;
	}
	return(listy);
}

draw_sample_by_rectangle<-function(a, r, n)   # a is list of converted plant tables,  r is the replicate (i.e., list element) to use, n is number of individuals to be sampled 
{
	q = zoom(ras, col = c("black", "wheat2") );
	points(a[[r]][,2], a[[r]][,1], col = "cyan", pch = ".", cex = 3);
	xminnie = floor(q@xmin);
	xmaxie = floor(q@xmax);
	yminnie = floor(q@ymin);
	ymaxie = floor(q@ymax);
	h = as.table(a[[r]]);
	g = h[h[,2] >= xminnie & h[,2] <= xmaxie & h[,1] >= yminnie & h[,1] <=  ymaxie,];
	len = length(g[,1]);
	print (len);
	if (len < n) {
		print ("too few plants to produce a sample of requested size");
		return();
	} else {
		q = seq(1:len);
		b1 = sample(q, n, replace = F);
		h2 = g[b1,];
print (dim(h2));		
		x = get_allele_frequencies(h2);
		return(x);
	}
}

get_allele_frequencies<-function(s)
{
	q = vector();
	for (i in 3:9) {    # need to change range IF MORE THAN FIVE loci
		count = 0;
		for(j in 1:length(s[,1])) {
			count = count + s[j,i];
		}
		count  = count / length(s[,1]);
		q = c(q,count);
	}
	return(q);
}

get_genotypes<-function(s)   #assumes cp, mt haplotpyes plus 5 diploid STRs
									# use returned matrix to print output to a file
{
	len = length(s[,1]);
	matty = matrix(nrow = len, ncol = 1)
	for (j in 1:len) {
		x = paste(s[j,1], s[j,2], s[j,3], s[j,4], s[j,5], sep = " ");
		x = paste(x, s[j,12], sep  = "/");
		x = paste(x, s[j,6], sep = " ");
		x = paste(x, s[j,13], sep = "/");
		x = paste(x, s[j,7], sep = " ");
		x = paste(x, s[j,14], sep = "/");
		x = paste(x, s[j,8], sep = " ");
		x = paste(x, s[j,15], sep = "/");
		x = paste(x, s[j,9], sep = " ");
		x = paste(x, s[j,16], sep  = "/");
		matty[j,1] = x;
	}
	return(matty);
}

get_calculable_genotypes<-function(s)   #assumes cp, mt haplotpyes plus 5 diploid STRs
										# use returned matrix for distance calculations and other popgen analysis
{
	len = length(s[,1]);
	matty = matrix(nrow = len, ncol = 14)
	for (j in 1:len) {
		matty[j,1:5] = s[j,1:5];
		matty[j,6] = s[j,12];
		matty[j,7] = s[j,6];
		matty[j,8] = s[j,13]
		matty[j,9] = s[j,7];
		matty[j,10] = s[j,14];
		matty[j,11] = s[j,8];
		matty[j,12] = s[j,15];
		matty[j,13] = s[j,9];
		matty[j,14] = s[j, 16];
	}
	return(matty);
}

distance_matrices<-function(m) #m is the matrix of positional and genetical data produced by get_calculable_genotypes
{
	print("euclidean matrix");
	len  = length(m[,1]);
	num.pc = choose(len,2);
	matties=list();
	eucl = matrix(nrow=len, ncol=len);  # matrix of euclidean distances
	for( i in 1:len) {
		for (j in 1:len) {
			if (i == j) {
				eucl[i,j] = 0;
			} else {
				eucl[i,j] = sqrt( (m[i,1] - m[j,1])^2 + (m[i,2] - m[j,2])^2      ); 
			}
		}
	}
	matties[[1]] = eucl;

#### look at ade4 manual for genetic distance calculations
	print ("genetic distance matrix");
	genet = matrix(0,nrow = len, ncol = len); # matrix of genetic distances, based on sqared difference between mean allele sizes of all five STRs
	for (i in 1:(len-1)) {
		for (j in (i+1):len) {
			genet[i,j] = 0;
			if (i != j) {
				genet[i,j] = genet[i,j] + (mean(c(m[i,5],m[i,6]))-mean(c(m[j,5],m[j,6])))^2;
				genet[i,j] = genet[i,j] + (mean(c(m[i,7],m[i,8]))-mean(c(m[j,7],m[j,8])))^2;
				genet[i,j] = genet[i,j] + (mean(c(m[i,9],m[i,10]))-mean(c(m[j,9],m[j,10])))^2;
				genet[i,j] = genet[i,j] + (mean(c(m[i,11],m[i,12]))-mean(c(m[j,11],m[j,12])))^2;
				genet[i,j] = genet[i,j] + (mean(c(m[i,13],m[i,14]))-mean(c(m[j,13],m[j,14])))^2;			
				genet[j,i] = genet[i,j];
			} 
		}
	}
	matties[[2]] = genet;

	print ("least cost path matrix");
	lcp = matrix(nrow = len, ncol = len);
	#convert m positional information to cell information in river raster matrix
	for (i in 1:len) {
		m[i,1] = (m[i,1]+1.5) *-1;
		m[i,2] = m[i,2] - 1.5;   
	}

	# now convert to nearest water cell
	for (i in 1:len) {
		found = 0;
		for (h in seq(-1,1,1)) {
			for (j in seq(-1,1,1)) {
				if (rr[ (m[i,1]+h) , (m[i,2]+j) ] == 0) { # then water
					m[i,1] = m[i,1]+h-1; 
					m[i,2] = m[i,2]+j-1; # final -1 to convert from one-based indexing of R to zero-based indexing of the called C++ function
					found = 1;
					break;
				}
			}
			if (found) {
				break;
			}
		}
	}	

	write.table(m[,1:2], file = "find_these_lcps", row.names = F, col.names = F, sep = " ", quote = F);
	.C("getlcps");
	lcp = as.matrix(read.table(file = "lcps", header = F, sep = ""));
	matties[[3]] = lcp;


	return(matties);
}

perform_mantel <- function (a,b,n) # two distance matrices and number of permutations
{
	require(vegan);
	mantel(a,b, method = "spear", permutations = n);
}

get_obs_het <- function(a, nloci)
{
	het = vector();
	for (i in 1:nloci) {
		if (a[i] == 0. | a[i] == 1.) {
			het[i] = 0.;
		} else {
			het[i] = 2*a[i]*(1-a[i]);
		}
	}
	return(het);
}

get_fst <- function(a1, a2, numloci) 
{
	vecky = vector();
	for (i in 1:numloci) {
		if (a1[i] == a2[i]) {
			vecky[i] = 0.;
		} else {
			q = (a1[i]+a2[i]) / 2;
			p = 1-q;
			ht = 2*p*q;
			hs1 = 2*a1[i]*(1-a1[i]);
			hs2 = 2*a2[i]*(1-a2[i]);
			hs = (hs1+hs2)/2;
			fst = (ht-hs)/ht;
			vecky[i] = fst;
		}
	}
	return(vecky);
}



########################### visualization
plot_bicauchy <- function(ofst)
{
	bc = read.table(file = "bicauchy", header = F, sep = "");
	bc[,2] = bc[,2] +500;
	bc[,1] = bc[,1] - 500;
	outlier = bc[bc[,2] > (500+ofst) | bc[,2] < (500-ofst) | bc[,1] > (-500+ofst) | bc[,1] < (-500-ofst),];
	plot(ras, col = c("snow4", "white"), xlim = c(200, 700), ylim = c(-700, -300));
	points(bc[,2], bc[,1], col = "magenta", cex = 1, pch = ".");
	points(outlier[,2], outlier[,1], col = "magenta", cex = 4, pch = "."); 
	print (length(outlier[,2]));
}

plot_bigauss <- function(ofst, mode)  # mode is 1 for dots and anything else for contour plot
{
	bc = read.table(file = "bigauss", header = F, sep = "");
	bc[,2] = bc[,2] +500;
	bc[,1] = bc[,1] - 500;
	outlier = bc[bc[,2] > (500+ofst) | bc[,2] < (500-ofst) | bc[,1] > (-500+ofst) | bc[,1] < (-500-ofst),];
	if (mode == 1) {
				plot(ras, col = c("snow4", "white"), xlim = c(475, 525), ylim = c(-525, -475));
	#	plot(ras, col = c("snow4", "white"), xlim = c(200, 700), ylim = c(-700, -300));
		points(bc[,2], bc[,1], col = "magenta", cex = 1, pch = ".");
		points(outlier[,2], outlier[,1], col = "magenta", cex = 4, pch = "."); 
		print (length(outlier[,2]));
	} else {
		plot(ras, col = c("snow4", "white"), xlim = c(475, 525), ylim = c(-525, -475));
		bub = kde2d(bc[,2], bc[,1], n =100);
		bub[["z"]] = bub[["z"]]*10000;
		contour(bub, add = T, col = "magenta")
#		contour(bub, add = T, col = "magenta", levels = c(0.0005, 0.001, 0.01, 0.25, 0.5, 0.75, 0.95, 0.99), labels = c("0.9995", "0.999", "0.99","0.75", "0.5", "0.25", "0.05", "0.01"))
	}
}

plot_plants <- function (obj, b, size)  # object is list of converted point tables, b is the replicate to plot, size is the size of plant points on the map
{
	plot(ras, col = c("black", "wheat2"), xlim = c(0,700), ylim = c(-700, 0), font.lab = 2, xlab= "East-West", ylab = "North-South");
	points(obj[[b]][,2], obj[[b]][,1], col = "magenta", pch = ".", cex = size);
}

draw_state_lines <- function (colory)
{
	lines (c(0,305), c(-344, -344), lwd = 2, lty = 2, col = colory); # IA/MN state line
	lines (c(414, 650), c(-554, -554), lwd =  2, lty = 2, col = colory); # WI/IL state line
}

draw_rectangle <-function(xl, xh, yl, yh)
{
	rect(xl, yl, xh, yh, border = "black", lwd = 1);
}

get_occupancy_raster <- function (obj, b)
{
	d2 = ras;
	values(d2) = obj[[b]];
	return(d2);
}

add_rivers <- function()
{
	points(rp[,1], rp[,2], col = "navy", cex = .6, pch = ".");
}

zoom_w_rivers_ByExtent <- function (oj, xmin, xmax, ymin, ymax)  # handy if you want to repeatedly bracket the same space on the map
{
	mypalette<-c(colors()[345], brewer.pal(9,"YlOrRd"));
	nr = ymax -ymin;
	nc = xmax - xmin;
	d2 = raster(nrow = nr, ncol = nc, xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax);
	values(d2) = oj[   (-1*ymax):((-1*ymax)+(nr-1))  , xmin:(xmin+(nc-1))  ];
	plot(d2, col = mypalette);	
	l = length(rp[,1]);
	for (i in 1:l) {
		if (rp[i,1] >= xmin && rp[i,1] <= xmax && rp[i,2] >= ymin && rp[i,2] <= ymax) {
		rect(rp[i,1]-0.5, rp[i,2]-0.5, rp[i,1]+0.5, rp[i,2]+0.5, col = colors()[354], border = F);
		}
	}
}

zoom_w_rivers <-function (oj,mode)  # oj is a replicate element of the list returned by the function read_occupancies -- e.g., o[[3]] for replicate 3
									# mode is 1 for ponits and 2 for squares
{
	mypalette<-c(colors()[348], brewer.pal(9,"YlOrRd"));
	d2 = ras;
	values(d2) = oj;
	plot(d2, col = mypalette);
	q = zoom(d2, col = mypalette);
	if (mode == 1) {
		points(rp[,1], rp[,2], col = "navy", cex = .6, pch = ".");
	}
	if (mode == 2) {
		l = length(rp[,1]);
		for (i in 1:l) {
			if (rp[i,1] >= q@xmin && rp[i,1] <= q@xmax && rp[i,2] >= q@ymin && rp[i,2] <= q@ymax) {
				rect(rp[i,1]-0.5, rp[i,2]-0.5, rp[i,1]+0.5, rp[i,2]+0.5, col = colors()[354], border = F);
			}
		}
	}
}

zoom_and_draw_sample <-function (oj, ind, mode, n)  # oj is a replicate element of the list returned by the function read_occupancies -- e.g., o[[3]] for replicate 3
									#inds is a replicate element of the list returned by the function convert_points
									# mode is 1 for ponits and 2 for squares
									# n is number of samples to draw
{
	mypalette<-c(colors()[348], brewer.pal(9,"YlOrRd"));
	d2 = ras;
	values(d2) = oj;
	plot(d2, col = mypalette);
	#q = zoom(d2, col = mypalette);
	e = drawExtent(col = "purple");
	rascrop = crop(d2, e);
	plot(rascrop, col = mypalette);
	if (mode == 1) {
		points(rp[,1], rp[,2], col = "navy", cex = .6, pch = ".");
	}
	if (mode == 2) {
		l = length(rp[,1]);
		for (i in 1:l) {
			if (rp[i,1] >= e@xmin && rp[i,1] <= e@xmax && rp[i,2] >= e@ymin && rp[i,2] <= e@ymax) {
				rect(rp[i,1]-0.5, rp[i,2]-0.5, rp[i,1]+0.5, rp[i,2]+0.5, col = colors()[354], border = F);
			}
		}
	}
	count = 0;
	dataa = matrix();
	while (x <- readline("Take another sample? (y/n)") == "y") {
		qq = drawExtent(col = "purple");
		xminnie = floor(qq@xmin);
		xmaxie = floor(qq@xmax);
		yminnie = floor(qq@ymin);
		ymaxie = floor(qq@ymax);
		h = as.table(ind);
		g = h[h[,2] >= xminnie & h[,2] <= xmaxie & h[,1] >= yminnie & h[,1] <=  ymaxie,];
		len = length(g[,1]);
		print (len);
		if (len < n) {
			print ("too few plants to produce a sample of requested size");
			return();
		} else {
			q = seq(1:len);
			b1 = sample(q, n, replace = F);
			h2 = g[b1,];
			print (dim(h2));		
			x = get_calculable_genotypes(h2);
			if (count == 0) {
				dataa = x;
			} else { 	
				dataa = rbind(dataa, x);
			}	
			#return(x);
		}
		count = count+1;
	}
	return (dataa);
}

# accessory function to previous function


zoom_w_riversBinary <-function (oj, threshold, mode)  # oj is a replicate element of the list returned by the function read_occupancies -- e.g., o[[3]] for replicate 3
									# threshold is the minimum number of plants in cell, which is the condition on which it's shown
									# mode is 1 for ponits and 2 for squares
{
	mypalette<-c(colors()[348], brewer.pal(9,"Oranges"));
	d2 = ras;
	oj = apply(as.matrix(oj), c(1,2), function(x) {if(x>=threshold) {x=1;} else {b=0}      });
	values(d2) = oj;
	plot(d2, col = mypalette);
	q = zoom(d2, col = mypalette);
	if (mode == 1) {
		points(rp[,1], rp[,2], col = "navy", cex = .6, pch = ".");
	}
	if (mode == 2) {
		l = length(rp[,1]);
		for (i in 1:l) {
			if (rp[i,1] >= q@xmin && rp[i,1] <= q@xmax && rp[i,2] >= q@ymin && rp[i,2] <= q@ymax) {
				rect(rp[i,1]-0.5, rp[i,2]-0.5, rp[i,1]+0.5, rp[i,2]+0.5, col = colors()[354], border = F);
			}
		}
	}
}

zoom_w_riversMean <-function (obj,mode)  # obj is the list returned by the function read_occupancies (not a single replicate as in the previous function)
									# mode is 1 for ponits and 2 for squares
{
	oj = Reduce("+", obj) / length(obj);
	mypalette<-c(colors()[249], brewer.pal(9,"Oranges"));
	d2 = ras;
	values(d2) = oj;
	plot(d2, col = mypalette);
	q = zoom(d2, col = mypalette);
	if (mode == 1) {
		points(rp[,1], rp[,2], col = "navy", cex = .6, pch = ".");
	}
	if (mode == 2) {
		l = length(rp[,1]);
		for (i in 1:l) {
			if (rp[i,1] >= q@xmin && rp[i,1] <= q@xmax && rp[i,2] >= q@ymin && rp[i,2] <= q@ymax) {
				rect(rp[i,1]-0.5, rp[i,2]-0.5, rp[i,1]+0.5, rp[i,2]+0.5, col = colors()[131], border = F);
			}
		}
	}
}

plot_initial_seeds <- function (filename, sizey) 
{
	d = as.matrix(read.table(file = filename, header = F, sep = ""));
	points (d[,2]+1.5, -1.5+(d[,1]*-1), pch = 1,  cex = sizey, col = "magenta");
}

read_flowlines <- function ()  # probably only needed for debugging
{
	fc = file("flowlines");
	l = strsplit(readLines(fc), " ");
	close(fc);
	l = lapply(l, as.numeric);
	print ("length of list");
	print (length(l));
	return(l);
}

draw_flowline <- function (ly, num)  # list with one entry per flowline, num is the entry number you want to draw
								#AGAIN, probably only need this for debuggin
{
	b = ly[[num]];
	q = length(b);
	if (q%%2 != 0 )
		return("not divisible by 2!");
	for (i in seq(1,q,2)) {
		a = b[i+1]+1.5;
		c  = -1.5 + b[i]*-1 ;
		print (a); 
		print (c);
		points(a, c, pch = 16, cex = .55, col = "green");
	}	
}