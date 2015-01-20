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
		using vector_t = std::vector< index_t >;
		enum SublatticeType {A, B};
		
		RhombicHoneycomb()
		{}
		
		void Resize(index_t l, RNG& rng)
		{
			L = l;
			nSites = 2 * L * L;
			nBonds = 3 * nSites / 2;
			nDirections = (L == 1 ? 1 : 3);
			DeallocateNeighborList();
			AllocateNeighborList();
			distanceMap.AllocateTable(nSites, nSites);
			distanceHistogram.resize(nSites, 0);
			BuildLookUpTable(rng);
			GenerateDistanceHistogram();
			numNeighborhood.resize(maxDistance + 1, 0);
			CountNeighborhood();
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
			return ((site % 2) == 0 ? SublatticeType::A : SublatticeType::B);
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
			return maxDistance;
		}
		
		index_t ShiftSite(index_t siteIndex, int_t direction, int_t distance = 1)
		{
			index_t newSite = siteIndex;
			for (int_t i = 0; i < distance; ++i)
				newSite = neighborList[newSite][direction];
			return newSite;
		}
		
		index_t RandomWalk(index_t site, int_t distance, RNG& rng)
		{
			index_t newSite = site;
			for (int j = 0; j < distance; ++j)
			{
				int_t newDir = static_cast<int_t>(rng() * nDirections);
				newSite = ShiftSite(newSite, newDir);
			}
			return newSite;
		}

		index_t FromNeighborhood(index_t site, int_t distance, RNG& rng)
		{
			index_t s = RandomSite(rng);
			while (Distance(s, site) > distance)
				s = RandomSite(rng);
			return s;
		}
		
		index_t NeighborhoodCount(int_t distance)
		{
			return numNeighborhood[distance];
		}

		index_t RandomSite(RNG& rng)
		{
			return static_cast<index_t>(rng() * nSites);
		}
	private:
		void AllocateNeighborList()
		{
			neighborList = new index_t*[nSites];
			for (index_t i = 0; i < nSites; ++i)
			{
				neighborList[i] = new index_t[nDirections + 1];
			}
			for (index_t i = 0; i < nSites; ++i)
				for (index_t j = 0; j < nDirections + 1; ++j)
					neighborList[i][j] = 0;
		}

		void DeallocateNeighborList()
		{
			if (neighborList == 0)
				return;
			for (index_t i = 0; i < nSites; ++i)
				delete[] neighborList[i];
			delete[] neighborList;
			neighborList = 0;
		}
		
		index_t ShiftSiteHardCode(index_t site, int_t direction, int_t distance = 1)
		{
			index_t newSite = site;
			for (int_t i = 0; i < distance; ++i)
			{
				if (Sublattice(newSite) == SublatticeType::B)
				{
					if ((newSite + 1) % (2*L) == 0)
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<index_t>(newSite - 4 * L + 1 + nSites) % nSites;
								break;
							case 1:
								newSite = static_cast<index_t>(newSite - 1 + nSites) % nSites;
								break;
							case 2:
								newSite = static_cast<index_t>(newSite - 2 * L + 1 + nSites) % nSites;
								break;
						}
					}
					else
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<index_t>(newSite - 2 * L + 1 + nSites) % nSites;
								break;
							case 1:
								newSite = static_cast<index_t>(newSite + 1 + nSites) % nSites;
								break;
							case 2:
								newSite = static_cast<index_t>(newSite - 1 + nSites) % nSites;
								break;
						}
					}
				}
				else
				{
					if (newSite % (2*L) == 0)
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<index_t>(newSite + 1 + nSites) % nSites;
								break;
							case 1:
								newSite = static_cast<index_t>(newSite + 2 * L - 1 + nSites) % nSites;
								break;
							case 2:
								newSite = static_cast<index_t>(newSite + 4 * L - 1 + nSites) % nSites;
								break;
						}
					}
					else
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<index_t>(newSite + 1 + nSites) % nSites;
								break;
							case 1:
								newSite = static_cast<index_t>(newSite - 1 + nSites) % nSites;
								break;
							case 2:
								newSite = static_cast<index_t>(newSite + 2 * L - 1 + nSites) % nSites;
								break;
						}
					}
				}
			}
			return newSite;
		}
		
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
			int_t nRuns = 1000 * nSites;
			for (int_t n = 0; n < nRuns; ++n)
			{
				int_t path = 0;
				index_t pos = i;
				while (path < shortestPath)
				{
					pos = ShiftSiteHardCode(pos, static_cast<int_t>(rng() * nDirections));
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
			{
				for (index_t j = 0; j < nSites; ++j)
				{
					distanceHistogram[Distance(i, j)] += 1;
					if (Distance(i, j) == 1)
					{
						neighborList[i][neighborList[i][nDirections]] = j;
						++neighborList[i][nDirections];
					}
				}
			}
				
			index_t i = distanceHistogram.size() - 1;
			while (i >= 0 && distanceHistogram[i] == 0)
				--i;
			maxDistance = i;
		}
		
		void CountNeighborhood()
		{
			index_t i = 0;
			for (index_t j = 0; j < nSites; ++j)
				numNeighborhood[Distance(i, j)] += 1;
			for (index_t j = 1; j <= maxDistance; ++j)
				numNeighborhood[j] += numNeighborhood[j-1];
		}
	private:
		using lookup_t = LookUpTable < int_t, index_t, 2 > ;
		using histogram_t = std::vector < int_t > ;
		
		index_t L;
		index_t nSites;
		index_t nBonds;
		int_t nDirections;
		index_t maxDistance;
		lookup_t distanceMap;
		histogram_t distanceHistogram;
		index_t** neighborList;
		vector_t numNeighborhood;
};