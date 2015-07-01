#pragma once
#include "random.h"

class Random
{
	public:
		using generator_t = randomnumbergenerator;
		using seed_t = std::size_t;

		Random()
		{
			rng = new generator_t();
			seed = rng->seed();
		}
		
		Random(seed_t seed)
			: seed(seed)
		{
			rng = new generator_t(seed);
		}
		
		~Random()
		{
			delete rng;
		}

		inline double operator()()
		{
			return rng->d();
		}

		seed_t Seed()
		{
			return seed;
		}
		
		void NewRng()
		{
			delete rng;
			rng = new generator_t();
		}
		void NewRng(seed_t seed)
		{
			delete rng;
			rng = new generator_t(seed);
		}
		generator_t* RngHandle()
		{
			return rng;
		}

	private:
		generator_t* rng;
		seed_t seed;
};
