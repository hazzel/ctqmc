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
		using site_t = std::tuple < int_t, int_t, int_t > ;
		using lookup_t = LookUpTable < int_t, index_t, 2 >;
		using histogram_t = std::vector < index_t >;
		using index_map_t = std::map < index_t, site_t >;
		using reverse_map_t = std::map < site_t, index_t >;
		enum SublatticeType {A, B};
		
	public:
		HexagonalHoneycomb()
		{}
		
		void Resize(index_t l, RNG& rng)
		{
			L = l;
			nSites = 6 * L * L;
			nBonds = 9 * L * L;
			nDirections = 3;
			BuildIndexMap();
			distanceMap.AllocateTable(nSites, nSites);
			distanceHistogram.resize(nSites, 0);
			BuildLookUpTable();
			GenerateDistanceHistogram();
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
			index_t i = distanceHistogram.size() - 1;
			while (i >= 0 && distanceHistogram[i] == 0)
				--i;
			return i;
		}

		index_t ShiftSite(index_t siteIndex, int_t direction, int_t distance = 1)
		{
			site_t newSite = indexMap[siteIndex];
			for (int_t i = 0; i < distance; ++i)
			{
				int_t& u = std::get<0>(newSite);
				int_t& v = std::get<1>(newSite);
				int_t& w = std::get<2>(newSite);
				
				//(-t + 1, t + 1 - p, p) and (t, -p + 1, p - t)
				//(p, -t + 1, t + 1 - p) and (p - t, t, -p + 1)
				//(t + 1 - p, p, -t + 1) and (-p + 1, p - t, t)
				if ((direction == 0) && (u == 1 - L) && (v == 1 + L - w) && (w >= 1))
				{
					u = L;
					v -= L;
					w -= L;
				}
				else if ((direction == 0) && (u == L) && (v <= 0) && (w == -v + 1 - L))
				{
					u = 1 - L;
					v += L;
					w += L;
				}
				else if ((direction == 1) && (u >= 1) && (v == 1 - L) && (w == 1 + L - u))
				{
					u -= L;
					v = L;
					w -= L;
				}
				else if ((direction == 1) && (u == -w + 1 - L) && (v == L) && (w <= 0))
				{
					u += L;
					v = 1 - L;
					w += L;
				}
				else if ((direction == 2) && (u == 1 + L - v) && (v >= 1) && (w == 1 - L))
				{
					u -= L;
					v -= L;
					w = L;
				}
				else if ((direction == 2) && (u <= 0) && (v == -u + 1 - L) && (w == L))
				{
					u += L;
					v += L;
					w = 1 - L;
				}
				else if(direction == 0)
				{
					std::get<0>(newSite) += (u + v + w == 1 ? 1 : -1);
				}
				else if(direction == 1)
				{
					std::get<1>(newSite) += (u + v + w == 1 ? 1 : -1);
				}
				else if(direction == 2)
				{
					std::get<2>(newSite) += (u + v + w == 1 ? 1 : -1);
				}
			}
			return reverseMap[newSite];
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
			return static_cast<index_t>(rng() * nSites);
		}

		index_t RandomSite(RNG& rng)
		{
			return static_cast<index_t>(rng() * nSites);
		}
	private:
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
				for (index_t j = 0; j < nSites; ++j)
					distanceHistogram[Distance(i, j)] += 1;
		}
		
	private:
		int_t L;
		index_t nSites;
		index_t nBonds;
		int_t nDirections;
		lookup_t distanceMap;
		histogram_t distanceHistogram;
		index_map_t indexMap;
		reverse_map_t reverseMap;
};