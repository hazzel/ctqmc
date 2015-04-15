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
#include "GeometryBase.h"

template<typename RNG, typename Int_t = std::int_fast32_t>
class HexagonalHoneycomb : public GeometryBase<RNG, Int_t>
{
	public:
		typedef Int_t int_t;
		using site_t = std::tuple < int_t, int_t, int_t >;
		using index_map_t = std::map < int_t, site_t >;
		using reverse_map_t = std::map < site_t, int_t >;
		typedef typename GeometryBase<RNG, Int_t>::SublatticeType SublatticeType;
		
	public:
		HexagonalHoneycomb() {}
		~HexagonalHoneycomb() {}
		
		void Resize(int_t l, RNG& rng)
		{
			this->L = l;
			this->nSites = 6 * this->L * this->L;
			this->nBonds = 9 * this->L * this->L;
			this->nDirections = 3;
			this->DeallocateNeighborList();
			this->AllocateNeighborList();
			BuildIndexMap();
			this->distanceMap.AllocateTable(this->nSites, this->nSites);
			this->distanceHistogram.resize(this->nSites, 0);
			this->BuildLookUpTable();
			this->GenerateDistanceHistogram();
			this->numNeighborhood.resize(this->maxDistance + 1, 0);
			this->CountNeighborhood();
			this->sublatVector.resize(this->nSites);
			for (int_t i = 0; i < this->nSites; ++i)
				this->sublatVector[i] = GetSublattice(i);
		}
		
		double Parity(int_t site)
		{
			return (Sublattice(site) == SublatticeType::A ? 1.0 : -1.0);
		}
	private:
		SublatticeType GetSublattice(int_t site)
		{
			site_t s = this->indexMap[site];
			if (std::get<0>(s) + std::get<1>(s) + std::get<2>(s) == 1)
				return SublatticeType::A;
			else
				return SublatticeType::B;
		}
		
		void BuildIndexMap()
		{
			int_t cnt = 0;
			for (int_t i = -this->L + 1; i <= this->L; ++i)
			{
				for (int_t j = -this->L + 1; j <= this->L; ++j)
				{
					for (int_t k = -this->L + 1; k <= this->L; ++k)
					{
						int_t sum = i + j + k;
						if (sum == 1 || sum == 2)
						{
							this->indexMap.insert(std::make_pair(cnt, std::make_tuple(i, j, k)));
							this->reverseMap.insert(std::make_pair(std::make_tuple(i, j, k), cnt));
							++cnt;
						}
					}
				}
			}
		}

		void BuildLookUpTable()
		{
			for (int_t i = 0; i < this->nSites; ++i)
			{
				for (int_t j = 0; j < i; ++j)
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
					int_t d1 = std::abs(u2 - u1 + 2 * this->L) + std::abs(v2 - v1 - this->L) + std::abs(w2 - w1 - this->L);
					int_t d2 = std::abs(u2 - u1 - 2 * this->L) + std::abs(v2 - v1 + this->L) + std::abs(w2 - w1 + this->L);
					int_t d3 = std::abs(u2 - u1 - this->L) + std::abs(v2 - v1 + 2 * this->L) + std::abs(w2 - w1 - this->L);
					int_t d4 = std::abs(u2 - u1 + this->L) + std::abs(v2 - v1 - 2 * this->L) + std::abs(w2 - w1 + this->L);
					int_t d5 = std::abs(u2 - u1 - this->L) + std::abs(v2 - v1 - this->L) + std::abs(w2 - w1 + 2 * this->L);
					int_t d6 = std::abs(u2 - u1 + this->L) + std::abs(v2 - v1 + this->L) + std::abs(w2 - w1 - 2 * this->L);
					this->distanceMap[i][j] = std::min({ d0, d1, d2, d3, d4, d5, d6 });
					this->distanceMap[j][i] = this->distanceMap[i][j];
				}
				this->distanceMap[i][i] = 0;
			}
		}	
	private:
		index_map_t indexMap;
		reverse_map_t reverseMap;
};