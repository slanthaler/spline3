/*
  random.hpp

  Simple wrapper for random number generation with <random>
*/
#pragma once
#include<random>

namespace Random{
	
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_real_distribution<> dist_uniform(0.0,1.0);
	
	double uniform(double min=0.0, double max=1.0){
		return min + (max-min)*dist_uniform(rng);
	};

	double sign(){
		double rnd = uniform(-1.0,1.0);
		return (rnd >= 0.0) - (rnd < 0.0);
	}
	
	void fill_uniform(std::vector<double>& tofill,
					  double min=0.0, double max=1.0){
		for(int i=0; i<tofill.size(); ++i){
			tofill[i] = uniform(min,max);
		}
	}

}
