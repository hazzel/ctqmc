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
class HexagonalHoneycomb
{
	public:
		typedef Uint_t index_t;
		typedef Int_t int_t;
		using site_t = std::tuple < int_t, int_t, int_t >;
		using lookup_t = LookUpTable < int_t, index_t, 2 >;
		using histogram_t = std::vector < index_t >;
		using index_map_t = std::map < index_t, site_t >;
		using reverse_map_t = std::map < site_t, index_t >;
		using vector_t = std::vector< index_t >;
		enum SublatticeType {A, B};
		
	public:
		HexagonalHoneycomb()
			: neighborList(0)
		{}
		
		~HexagonalHoneycomb()
		{
			DeallocateNeighborList();
		}
		
		void Resize(index_t l, RNG& rng)
		{
			L = l;
			nSites = 6 * L * L;
			nBonds = 9 * L * L;
			nDirections = 3;
			DeallocateNeighborList();
			AllocateNeighborList();
			BuildIndexMap();
			distanceMap.AllocateTable(nSites, nSites);
			distanceHistogram.resize(nSites, 0);
			BuildLookUpTable();
			GenerateDistanceHistogram();
			numNeighborhood.resize(maxDistance + 1, 0);
			CountNeighborhood();
		}
		
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

		SublatticeType Sublattice(index_t site)
		{
			site_t s = indexMap[site];
			if (std::get<0>(s) + std::get<1>(s) + std::get<2>(s) == 1)
				return SublatticeType::A;
			else
				return SublatticeType::B;
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
	private:
		void AllocateNeighborList()
		{
			neighborList = new index_t*[nSites];
			for (index_t i = 0; i < nSites; ++i)
			{
				neighborList[i] = new index_t[4];
			}
			for (index_t i = 0; i < nSites; ++i)
				for (index_t j = 0; j < 4; ++j)
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
		
		void BuildIndexMap()
		{
			index_t cnt = 0;
			for (int_t i = -L + 1; i <= L; ++i)
			{
				for (int_t j = -L + 1; j <= L; ++j)
				{
					for (int_t k = -L + 1; k <= L; ++k)
					{
						int_t sum = i + j + k;
						if (sum == 1 || sum == 2)
						{
							indexMap.insert(std::make_pair(cnt, std::make_tuple(i, j, k)));
							reverseMap.insert(std::make_pair(std::make_tuple(i, j, k), cnt));
							++cnt;
						}
					}
				}
			}
		}

		void BuildLookUpTable()
		{
			for (index_t i = 0; i < nSites; ++i)
			{
				for (index_t j = 0; j < i; ++j)
				{
					site_t s1 = indexMap[i];
					site_t s2 = indexMap[j];
					int_t u1 = std::get<0>(s1);
					int_t v1 = std::get<1>(s1);
					int_t w1 = std::get<2>(s1);
					int_t u2 = std::get<0>(s2);
					int_t v2 = std::get<1>(s2);
					int_t w2 = std::get<2>(s2);
					int_t d0 = std::abs(u2 - u1) + std::abs(v2 - v1) + std::abs(w2 - w1);
					int_t d1 = std::abs(u2 - u1 + 2 * L) + std::abs(v2 - v1 - L) + std::abs(w2 - w1 - L);
					int_t d2 = std::abs(u2 - u1 - 2 * L) + std::abs(v2 - v1 + L) + std::abs(w2 - w1 + L);
					int_t d3 = std::abs(u2 - u1 - L) + std::abs(v2 - v1 + 2 * L) + std::abs(w2 - w1 - L);
					int_t d4 = std::abs(u2 - u1 + L) + std::abs(v2 - v1 - 2 * L) + std::abs(w2 - w1 + L);
					int_t d5 = std::abs(u2 - u1 - L) + std::abs(v2 - v1 - L) + std::abs(w2 - w1 + 2 * L);
					int_t d6 = std::abs(u2 - u1 + L) + std::abs(v2 - v1 + L) + std::abs(w2 - w1 - 2 * L);
					distanceMap[i][j] = std::min({ d0, d1, d2, d3, d4, d5, d6 });
					distanceMap[j][i] = distanceMap[i][j];
				}
				distanceMap[i][i] = 0;
			}
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
						neighborList[i][neighborList[i][3]] = j;
						++neighborList[i][3];
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
		int_t L;
		index_t nSites;
		index_t nBonds;
		int_t nDirections;
		index_t maxDistance;
		lookup_t distanceMap;
		histogram_t distanceHistogram;
		index_map_t indexMap;
		reverse_map_t reverseMap;
		index_t** neighborList;
		vector_t numNeighborhood;
};