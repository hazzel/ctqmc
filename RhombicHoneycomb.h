#pragma once
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <cstdint>
#include "LookUpTable.h"

template<typename RNG, typename Uint_t = std::uint_fast32_t, typename Int_t = std::int_fast32_t>
class RhombicHoneycomb
{
	public:
		typedef Uint_t index_t;
		typedef Int_t int_t;
		enum SublatticeType {A, B};
		
		RhombicHoneycomb()
		{}
		
		void Resize(index_t l, RNG& rng)
		{
			L = l;
			nSites = 2 * L * L;
			nBonds = 3 * nSites / 2;
			latticeDirection = {{ -1, 1, 2 * static_cast<int>(L)-1 }};
			distanceMap.AllocateTable(nSites, nSites);
			distanceHistogram.resize(nSites, 0);
			BuildLookUpTable(rng);
			GenerateDistanceHistogram();
		}
		
		int_t Distance(index_t s1, index_t s2)
		{
			return distanceMap[s1][s2];
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
			return site % 2 == 0 ? SublatticeType::A : SublatticeType::B;
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
				if (direction != 2)
					newSite = static_cast<index_t>(newSite + latticeDirection[direction] + nSites) % nSites;
				else
					newSite = static_cast<index_t>(newSite + latticeDirection[direction] * std::pow(-1.0, newSite) + nSites) % nSites;
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
			return static_cast<index_t>(rng() * nSites);
		}

		index_t RandomSite(RNG& rng)
		{
			return static_cast<index_t>(rng() * nSites);
		}
	private:
		void BuildLookUpTable(RNG& rng)
		{
			for (index_t i = 0; i < nSites; ++i)
			{
				for (index_t j = 0; j < i; ++j)
				{
					distanceMap[i][j] = SimulateDistance(i, j, rng);
					distanceMap[j][i] = distanceMap[i][j];
				}
				distanceMap[i][i] = 0;
			}
		}

		int_t SimulateDistance(index_t i, index_t j, RNG& rng)
		{
			int_t shortestPath = nSites;
			int_t nRuns = 2500 * nSites;
			for (int_t n = 0; n < nRuns; ++n)
			{
				int_t path = 0;
				index_t pos = i;
				while (path < shortestPath)
				{
					pos = RandomWalk(pos, 1, rng);
					++path;
					if (pos == j)
					{
						shortestPath = path;
						break;
					}
				}
			}
			return shortestPath;
		}

		void GenerateDistanceHistogram()
		{
			for (index_t i = 0; i < nSites; ++i)
				for (index_t j = 0; j < nSites; ++j)
					distanceHistogram[Distance(i, j)] += 1;
			distanceHistogram[0] = nSites;
		}
	private:
		using array_t = std::array < int_t, 3 > ;
		using lookup_t = LookUpTable < int_t, index_t, 2 > ;
		using histogram_t = std::vector < int_t > ;
		
		index_t L;
		index_t nSites;
		index_t nBonds;
		array_t latticeDirection;
		lookup_t distanceMap;
		histogram_t distanceHistogram;
};