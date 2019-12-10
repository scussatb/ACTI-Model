#ifndef __PARETO_DISTRIBUTION_HPP__
#define __PARETO_DISTRIBUTION_HPP__

class pareto_distribution {
protected:
	double xStar;
	double invd;
	std::uniform_real_distribution<double> dist;

public:
	pareto_distribution(double xStar, double d):
       		dist(0.0, 1.0)	{
		this->xStar = xStar;
		invd = -1.0/d;
	}

	double operator()(std::default_random_engine& g) {
		return xStar*pow(dist(g), invd);
	}	
};

#endif
