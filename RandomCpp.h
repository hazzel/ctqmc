#pragma once
#include <random>

class RandomCpp
{
	public:
		using generator_t = std::mt19937_64;
		using dist_t = std::uniform_real_distribution<double>;
		using seed_t = generator_t::result_type;

		RandomCpp(seed_t seed = generator_t::default_seed)
			: rng(generator_t(seed)), dist(dist_t(0.0, 1.0)), seed(seed)
		{}
		
		RandomCpp(double min, double max, seed_t seed = generator_t::default_seed)
			: rng(generator_t(seed)), dist(dist_t(min, max)), seed(seed)
		{}

		double operator()()
		{
			return dist(rng);
		}

		seed_t Seed()
		{
			return seed;
		}

	private:
		generator_t	rng;
		dist_t dist;
		seed_t seed;
};