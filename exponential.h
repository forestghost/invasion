#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

struct Expondist {
	Doub bet;
	Expondist(Doub bbet) : bet(bbet) {
		if (bet <= 0.) throw("bad parameter in Expondist");
	}
	Doub p(Doub x) {
		if (x <0.) throw("bad x in Expondist");
		return bet*exp(-bet*x);
	}
	Doub cdf(Doub x) {
		if (x <0.) throw("bad x in Expondist");
		return 1.  - exp(-bet*x);
	}
	Doub invcdf(Doub p) {
		if (p <0. || p>=1.) throw ("bad p in Expondist");
		return -log(1.-p)/bet;
	}
};
#endif
