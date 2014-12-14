#pragma once
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <cstdint>
#include "LookUpTable.h"

template<typename RNG, typename Uint_t = std::uint_fast32_t, typename Int_t = std::int_fast32_t>
class ObcChain1D
{
	public:
		typedef Uint_t index_t;
		typedef Int_t int_t;
		enum SublatticeType {A, B};
		
		ObcChain1D()
		{}
		
		ObcChain1D(index_t l)
		{
			Resize(l);
		}
		
		void Resize(index_t l, RNG& rng)
		{
			L = l;
			nSites = l;
			nBonds = l - 1;
			latticeDirection = {{ -1, 1 }};
			distanceHistogram.resize(nSites, 0);
			GenerateDistanceHistogram();
		}
		
		int_t Distance(index_t s1, index_t s2)
		{
			return (s1 > s2) ? s1 - s2 : s2 - s1;
		}

		int_t DistanceHistogram(int_t distance)
		{
			return distanceHistogram[distance];
		}
		
		bool IsNeighbor(index_t s1, index_t s2)
		{
			return Distance(s1, s2) == 1;
		}

		SublatticeType Sublattice(index_t site)
		{
			return SublatticeType::A;
		}

		index_t Sites()
		{
			return nSites;
		}

		index_t Bonds()
		{
			return nBonds;
		}
		
		index_t MaxDistance()
		{
			index_t i = distanceHistogram.size() - 1;
			while (i >= 0 && distanceHistogram[i] == 0)
				--i;
			return i;
		}
		
		index_t ShiftSite(index_t site, int_t direction, int_t distance = 1)
		{
			index_t newSite = site;
			for (int_t i = 0; i < distance; ++i)
			{
				if (newSite != 0 && newSite != nSites - 1)
					newSite = newSite + latticeDirection[direction];
			}
			return newSite;
		}
		
		index_t RandomWalk(index_t site, int_t distance, RNG& rng)
		{
			index_t newSite = site;
			int_t lastDir = latticeDirection.size();
			for (int j = 0; j < distance; ++j)
			{
				int_t newDir = static_cast<int_t>(rng() * latticeDirection.size());
				while(lastDir == newDir)
					newDir = static_cast<int_t>(rng() * latticeDirection.size());
				newSite = ShiftSite(newSite, newDir);
				lastDir = newDir;
			}
			return newSite;
		}
		
		index_t FromNeighborhood(index_t site, int_t distance, RNG& rng)
		{
			double r = rng();
			if (r < 1.0 / 3.0)
				return site;
			else if (r < 2.0 / 3.0)
				return ShiftSite(site, 0);
			else if (r < 3.0 / 3.0)
				return ShiftSite(site, 1);
		}
				
		index_t RandomSite(RNG& rng)
		{
			return static_cast<index_t>(rng() * nSites);
		}
	private:
		void GenerateDistanceHistogram()
		{
			for (index_t i = 0; i < nSites; ++i)
				for (index_t j = 0; j < nSites; ++j)
					distanceHistogram[Distance(i, j)] += 1;
			distanceHistogram[0] = nSites;
		}
	private:
		using array_t = std::array < int_t, 3 > ;
		using histogram_t = std::vector < int_t > ;
		
		index_t L;
		index_t nSites;
		index_t nBonds;
		array_t latticeDirection;
		histogram_t distanceHistogram;
};