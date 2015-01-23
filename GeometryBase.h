#pragma once
#include <algorithm>
#include <cmath>
#include <array>
#include <map>
#include <tuple>
#include <iostream>
#include <vector>
#include <cstdint>
#include "LookUpTable.h"

template<typename RNG, typename Uint_t = std::uint_fast32_t, typename Int_t = std::int_fast32_t>
class GeometryBase
{
	public:
		typedef Uint_t index_t;
		typedef Int_t int_t;
		using site_t = std::tuple < int_t, int_t, int_t >;
		using lookup_t = LookUpTable < int_t, index_t, 2 >;
		using histogram_t = std::vector < index_t >;
		using vector_t = std::vector< index_t >;
		enum SublatticeType {A, B};
		
	public:
		GeometryBase()
			: neighborList(0)
		{}
		
		virtual ~GeometryBase()
		{
			DeallocateNeighborList();
		}
		
		virtual void Resize(index_t l, RNG& rng) = 0;
		virtual SublatticeType Sublattice(index_t site) = 0;
		
		int_t Distance(index_t s1, index_t s2)
		{
			return distanceMap[s1][s2];
		}

		index_t DistanceHistogram(int_t distance)
		{
			return distanceHistogram[distance];
		}

		bool IsNeighbor(index_t s1, index_t s2)
		{
			return Distance(s1, s2) == 1;
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
				index_t newDir = static_cast<index_t>(rng() * nDirections);
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
	protected:
		void AllocateNeighborList()
		{
			neighborList = new index_t*[nSites];
			for (index_t i = 0; i < nSites; ++i)
			{
				neighborList[i] = new index_t[nDirections+1];
			}
			for (index_t i = 0; i < nSites; ++i)
				for (index_t j = 0; j < nDirections+1; ++j)
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
		
	protected:
		int_t L;
		index_t nSites;
		index_t nBonds;
		int_t nDirections;
		index_t maxDistance;
		lookup_t distanceMap;
		histogram_t distanceHistogram;
		index_t** neighborList;
		vector_t numNeighborhood;
};