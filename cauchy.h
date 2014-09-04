#ifndef CAUCHY_H
#define CAUCHY_H

struct Cauchydist {
	Doub mu, sig;
	Cauchydist(Doub mmu = 0., Doub ssig = 1.) : mu(mmu), sig(ssig) {
		if (sig < 0.) throw ("bad sig in Cauchydist");
	}
	Doub p (Doub x) {
		return 0.318309886183790671/(sig*(1.+SQR((x-mu)/sig)));
	}
	Doub cdf(Doub x) {
		return 0.5+0.318309886183790671*atan2(x-mu,sig);
	}
	Doub invcdf(Doub p) {
		if (p <=0. || p >=1.) throw ("bad p in Cauchydist");
		return mu + sig*tan(3.14159265358979324*(p-0.5));
	}
};

#endif