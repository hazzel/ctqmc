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
		
		RhombicHoneycomb(const std::string& filename, bool fileIO)
			: filename(filename)
		{
			this->fileIO = fileIO;
		}
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
			if (this->fileIO && FileExists(filename))
			{
				this->ReadFromFile(filename);
				//std::cout << "Distance Histogram:" << std::endl;
				//for (int_t i = 0; i <= this->maxDistance; ++i)
				//	std::cout << i << " : " << this->distanceHistogram[i] << std::endl;
				//std::cout << "Neighborhood Histogram:" << std::endl;
				//for (int_t i = 0; i <= this->maxDistance; ++i)
				//	std::cout << i << " : " << this->numNeighborhood[i] << std::endl;
				//std::cout << "Dist(0, N/2) = " << this->distanceMap[0][this->nSites / 2] << std::endl;
			}
			else
			{
				this->BuildLookUpTable(rng);
				this->GenerateDistanceHistogram();
				this->numNeighborhood.resize(this->maxDistance + 1, 0);
				this->CountNeighborhood();
				if (this->fileIO)
					this->SaveToFile(filename);
				std::cout << "Distance Histogram:" << std::endl;
				for (int_t i = 0; i <= this->maxDistance; ++i)
					std::cout << this->distanceHistogram[i] << std::endl;
				std::cout << "Dist(0, N/2) = " << this->distanceMap(0, this->nSites / 2) << std::endl;
			}
			this->sublatVector.resize(this->nSites);
			for (int_t i = 0; i < this->nSites; ++i)
				this->sublatVector[i] = GetSublattice(i);
		}
		
		double Parity(int_t site)
		{
			return ((site % 2) == 0 ? 1.0 : -1.0);
		}
	private:
		SublatticeType GetSublattice(int_t site)
		{
			return ((site % 2) == 0 ? SublatticeType::A : SublatticeType::B);
		}
		
		int_t ShiftSiteHardCode(int_t site, int_t direction, int_t distance = 1)
		{
			int_t newSite = site;
			for (int_t i = 0; i < distance; ++i)
			{
				if (GetSublattice(newSite) == SublatticeType::B)
				{
					if ((newSite + 1) % (2*this->L) == 0)
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<int_t>(newSite - 4 * this->L + 1 + this->nSites) % this->nSites;
								break;
							case 1:
								newSite = static_cast<int_t>(newSite - 2 * this->L + 1 + this->nSites) % this->nSites;
								break;
							case 2:
								newSite = static_cast<int_t>(newSite - 1 + this->nSites) % this->nSites;
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
								newSite = static_cast<int_t>(newSite + 4 * this->L - 1 + this->nSites) % this->nSites;
								break;
							case 1:
								newSite = static_cast<int_t>(newSite + 2 * this->L - 1 + this->nSites) % this->nSites;
								break;
							case 2:
								newSite = static_cast<int_t>(newSite + 1 + this->nSites) % this->nSites;
								break;
						}
					}
					else
					{
						switch (direction)
						{
							case 0:
								newSite = static_cast<int_t>(newSite + 2 * this->L - 1 + this->nSites) % this->nSites;
								break;
							case 1:
								newSite = static_cast<int_t>(newSite - 1 + this->nSites) % this->nSites;
								break;
							case 2:
								newSite = static_cast<int_t>(newSite + 1 + this->nSites) % this->nSites;
								break;
						}
					}
				}
			}
			return newSite;
		}
		
		
		void BuildLookUpTable(RNG& rng)
		{
			for (int_t i = 0; i < this->nSites; ++i)
				for (int_t j = 0; j < this->nSites; ++j)
					this->distanceMap(i, j) = -1;
			int_t cnt = 0;

			for (int_t i = 0; i < 2; ++i)
			{
				for (int_t j = 0; j < this->nSites; ++j)
				{
					this->distanceMap(i, j) = this->SimulateDistance(i, j, rng);
					this->distanceMap(j, i) = this->distanceMap(i, j);
					cnt += 2;
				}
			}
			for (int_t k = 0; k < 10000000; ++k)
			{
				int_t dir = this->RandomDirection(rng);
				int_t dist = 1;
				int_t i = this->RandomSite(rng);
				int_t j = this->RandomSite(rng);
				while(this->distanceMap(i, j) < 0 || (GetSublattice(i) != GetSublattice(j)))
				{
					i = this->RandomSite(rng);
					j = this->RandomSite(rng);
				}
				int_t u = ShiftSiteHardCode(i, dir, dist);
				int_t v = ShiftSiteHardCode(j, dir, dist);
				if (this->distanceMap(u, v) < 0)
				{
					this->distanceMap(u, v) = this->distanceMap(i, j);
					this->distanceMap(v, u) = this->distanceMap(i, j);
					cnt += 2;
				}
			}
			
			for (int_t i = 0; i < this->nSites; ++i)
			{
				for (int_t j = 0; j < i; ++j)
				{
					if (this->distanceMap(i, j) < 0)
					{
						this->distanceMap(i, j) = this->SimulateDistance(i, j, rng);
						this->distanceMap(j, i) = this->distanceMap(i, j);
						cnt += 2;
						//std::cout << cnt << std::endl;
					}
				}
				if (this->distanceMap(i, i) < 0)
				{
					this->distanceMap(i, i) = 0;
					++cnt;
				}
			}
		}

		int_t SimulateDistance(int_t i, int_t j, RNG& rng)
		{
			if (i == j)
				return 0;
			int_t shortestPath = 2 * std::sqrt(this->nSites);
			int_t nRuns = 10000 * this->nSites;
			for (int_t n = 0; n < nRuns; ++n)
			{
				int_t path = 0;
				int_t pos = i;
				int_t prevPos = i;
				while (path < shortestPath)
				{
					int_t newPos = prevPos;
					while (newPos == prevPos)
						newPos = ShiftSiteHardCode(pos, static_cast<int_t>(rng() * this->nDirections));
					prevPos = pos;
					pos = newPos;
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
