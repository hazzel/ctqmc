#pragma once
#include <iostream>
#include <algorithm>
#include <array>
#include <vector>
#include <cstdint>
#include <string>
#include "LookUpTable.h"
#include "GeometryBase.h"

template<typename RNG, typename Int_t = std::int_fast32_t>
class RhombicHoneycomb : public GeometryBase<RNG, Int_t>
{
	public:
		typedef Int_t int_t;
		typedef typename GeometryBase<RNG, Int_t>::SublatticeType SublatticeType;
		
		RhombicHoneycomb(const std::string& filename)
			: filename(filename)
		{}
		~RhombicHoneycomb() {}

		void Resize(int_t l, RNG& rng)
		{
			this->L = l;
			this->nSites = 2 * this->L * this->L;
			this->nBonds = 3 * this->nSites / 2;
			this->nDirections = (this->L == 1 ? 1 : 3);
			this->DeallocateNeighborList();
			this->AllocateNeighborList();
			this->distanceMap.AllocateTable(this->nSites, this->nSites);
			this->distanceHistogram.resize(this->nSites, 0);
			if (!FileExists(filename))
			{
				this->BuildLookUpTable(rng);
				this->GenerateDistanceHistogram();
				this->numNeighborhood.resize(this->maxDistance + 1, 0);
				this->CountNeighborhood();
				this->SaveToFile(filename);
			}
			else
			{
				this->ReadFromFile(filename);
			}
		}

		SublatticeType Sublattice(int_t site)
		{
			return ((site % 2) == 0 ? SublatticeType::A : SublatticeType::B);
		}
	private:		
		int_t ShiftSiteHardCode(int_t site, int_t direction, int_t distance = 1)
		{
			int_t newSite = site;
			for (int_t i = 0; i < distance; ++i)
			{
				if (Sublattice(newSite) == SublatticeType::B)
				{
					if ((newSite + 1) % (2*this->L) == 0)
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<int_t>(newSite - 4 * this->L + 1 + this->nSites) % this->nSites;
								break;
							case 1:
								newSite = static_cast<int_t>(newSite - 1 + this->nSites) % this->nSites;
								break;
							case 2:
								newSite = static_cast<int_t>(newSite - 2 * this->L + 1 + this->nSites) % this->nSites;
								break;
						}
					}
					else
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<int_t>(newSite - 2 * this->L + 1 + this->nSites) % this->nSites;
								break;
							case 1:
								newSite = static_cast<int_t>(newSite + 1 + this->nSites) % this->nSites;
								break;
							case 2:
								newSite = static_cast<int_t>(newSite - 1 + this->nSites) % this->nSites;
								break;
						}
					}
				}
				else
				{
					if (newSite % (2*this->L) == 0)
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<int_t>(newSite + 1 + this->nSites) % this->nSites;
								break;
							case 1:
								newSite = static_cast<int_t>(newSite + 2 * this->L - 1 + this->nSites) % this->nSites;
								break;
							case 2:
								newSite = static_cast<int_t>(newSite + 4 * this->L - 1 + this->nSites) % this->nSites;
								break;
						}
					}
					else
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<int_t>(newSite + 1 + this->nSites) % this->nSites;
								break;
							case 1:
								newSite = static_cast<int_t>(newSite - 1 + this->nSites) % this->nSites;
								break;
							case 2:
								newSite = static_cast<int_t>(newSite + 2 * this->L - 1 + this->nSites) % this->nSites;
								break;
						}
					}
				}
			}
			return newSite;
		}
		
		void BuildLookUpTable(RNG& rng)
		{
			#pragma omp parallel for
			for (int_t i = 0; i < this->nSites; ++i)
			{
				for (int_t j = 0; j < i; ++j)
				{
					this->distanceMap[i][j] = this->SimulateDistance(i, j, rng);
					this->distanceMap[j][i] = this->distanceMap[i][j];
				}
				this->distanceMap[i][i] = 0;
			}
		}

		int_t SimulateDistance(int_t i, int_t j, RNG& rng)
		{
			if (i == j)
				return 0;
			int_t shortestPath = this->nSites;
			int_t nRuns = 100000 * this->nSites;
			for (int_t n = 0; n < nRuns; ++n)
			{
				int_t path = 0;
				int_t pos = i;
				while (path < shortestPath)
				{
					pos = ShiftSiteHardCode(pos, static_cast<int_t>(rng() * this->nDirections));
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
	private:
		std::string filename;
};