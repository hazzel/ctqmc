#pragma once
#include <algorithm>
#include <cmath>
#include <array>
#include <map>
#include <iostream>
#include <vector>
#include <cstdint>
#include <fstream>
#include "LookUpTable.h"
#include <sys/stat.h>

inline bool FileExists(const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0); 
}

inline bool DirExists(const std::string& name)
{
	struct stat buffer;

	if(stat (name.c_str(), &buffer) != 0)
		return 0;
	else if(buffer.st_mode & S_IFDIR)
		return 1;
	else
		return 0;
}

template<typename RNG, typename Int_t = int32_t>
class GeometryBase
{
	public:
		typedef Int_t int_t;
		typedef LookUpTable < int_t, 2 > lookup_t;
		typedef std::vector< int_t > vector_t;
		enum SublatticeType {A, B};
		typedef std::vector< SublatticeType > sub_lat_vector_t;

		//using site_t = std::tuple < int_t, int_t, int_t >;
		//using lookup_t = LookUpTable < int_t, 2 >;
		//using vector_t = std::vector< int_t >;
		//enum SublatticeType {A, B};
		//using sub_lat_vector_t = std::vector< SublatticeType >;
		
	public:
		GeometryBase()
			: neighborList(0)
		{}
		
		virtual ~GeometryBase()
		{
			DeallocateNeighborList();
		}
		
		virtual void Resize(int_t l, RNG& rng) = 0;
		
		inline SublatticeType Sublattice(int_t site)
		{
			return sublatVector[site];
		}
		
		inline int_t Distance(int_t s1, int_t s2)
		{
			return distanceMap(s1, s2);
		}

		inline int_t DistanceHistogram(int_t distance)
		{
			return distanceHistogram[distance];
		}

		inline bool IsNeighbor(int_t s1, int_t s2)
		{
			return Distance(s1, s2) == 1;
		}

		inline int_t Sites()
		{
			return nSites;
		}

		inline int_t Bonds()
		{
			return nBonds;
		}
		
		inline int_t MaxDistance()
		{
			return maxDistance;
		}

		inline int_t RandomDirection(RNG& rng)
		{
			return static_cast<int_t>(rng() * nDirections);
		}

		int_t ShiftSite(int_t siteIndex, int_t direction, int_t distance = 1)
		{
			int_t newSite = siteIndex;
			for (int_t i = 0; i < distance; ++i)
				newSite = neighborList[newSite][direction];
			return newSite;
		}

		int_t RandomWalk(int_t site, int_t distance, RNG& rng)
		{
			int_t newSite = site;
			for (int j = 0; j < distance; ++j)
			{
				int_t newDir = static_cast<int_t>(rng() * nDirections);
				newSite = ShiftSite(newSite, newDir);
			}
			return newSite;
		}
		
		int_t FromNeighborhood(int_t site, int_t distance, RNG& rng)
		{
			int_t s = RandomSite(rng);
			while (Distance(s, site) > distance)
				s = RandomSite(rng);
			return s;
		}

		int_t FromDistance(int_t site, int_t distance, RNG& rng)
		{
			int_t s = RandomSite(rng);
			while (Distance(s, site) != distance)
				s = RandomSite(rng);
			return s;
		}
		
		inline int_t NeighborhoodCount(int_t distance)
		{
			return numNeighborhood[distance];
		}

		inline int_t DistanceCount(int_t distance)
		{
			return numDistance[distance];
		}

		inline int_t RandomSite(RNG& rng)
		{
			return static_cast<int_t>(rng() * nSites);
		}

		void SaveToFile(const std::string& filename)
		{
			if (FileExists(filename))
				return;
			std::ofstream os(filename.c_str(), std::ofstream::binary);
			os.write((char*)&maxDistance, sizeof(maxDistance));
			os.write((char*)&nSites, sizeof(nSites));
			for (int_t i = 0; i < nSites; ++i)
			{
				for (int_t j = 0; j < nSites; ++j)
				{
					os.write((char*)&distanceMap(i, j), sizeof(distanceMap(i, j)));
				}
				for (int_t j = 0; j < nDirections; ++j)
				{
					os.write((char*)&neighborList[i][j], sizeof(neighborList[i][j]));
				}
				os.write((char*)&distanceHistogram[i], sizeof(distanceHistogram[i]));
			}
			for (int_t i = 0; i <= maxDistance; ++i)
			{
				os.write((char*)&numNeighborhood[i], sizeof(numNeighborhood[i]));
			}
			for (int_t i = 0; i <= maxDistance; ++i)
			{
				os.write((char*)&numDistance[i], sizeof(numDistance[i]));
			}
			os.close();
		}

		void ReadFromFile(const std::string& filename)
		{
			std::ifstream is(filename.c_str(), std::ifstream::binary);
			if (is.is_open())
			{
				while(is.good())
				{
					is.read((char*)&maxDistance, sizeof(maxDistance));
					is.read((char*)&nSites, sizeof(nSites));
					this->numNeighborhood.resize(this->maxDistance + 1, 0);
					this->numDistance.resize(this->maxDistance + 1, 0);
					
					for (int_t i = 0; i < nSites; ++i)
					{
						for (int_t j = 0; j < nSites; ++j)
						{
							is.read((char*)&distanceMap(i, j), sizeof(distanceMap(i, j)));
						}
						for (int_t j = 0; j < nDirections; ++j)
						{
							is.read((char*)&neighborList[i][j], sizeof(neighborList[i][j]));
						}
						is.read((char*)&distanceHistogram[i], sizeof(distanceHistogram[i]));
					}
					for (int_t i = 0; i <= maxDistance; ++i)
					{
						is.read((char*)&numNeighborhood[i], sizeof(numNeighborhood[i]));
					}
					for (int_t i = 0; i <= maxDistance; ++i)
					{
						is.read((char*)&numDistance[i], sizeof(numDistance[i]));
					}
					is.close();
					CountNeighborhood();
				}
			}
		}
	protected:
		void AllocateNeighborList()
		{
			neighborList = new int_t*[nSites];
			for (int_t i = 0; i < nSites; ++i)
			{
				neighborList[i] = new int_t[nDirections+1];
			}
			for (int_t i = 0; i < nSites; ++i)
				for (int_t j = 0; j < nDirections+1; ++j)
					neighborList[i][j] = 0;
		}

		void DeallocateNeighborList()
		{
			if (neighborList == 0)
				return;
			for (int_t i = 0; i < nSites; ++i)
				delete[] neighborList[i];
			delete[] neighborList;
			neighborList = 0;
		}

		void GenerateDistanceHistogram()
		{
			for (int_t i = 0; i < nSites; ++i)
			{
				for (int_t j = 0; j < nSites; ++j)
				{
					distanceHistogram[Distance(i, j)] += 1;
					if (Distance(i, j) == 1)
					{
						neighborList[i][neighborList[i][nDirections]] = j;
						++neighborList[i][nDirections];
					}
				}
			}
			int_t i = distanceHistogram.size() - 1;
			while (i >= 0 && distanceHistogram[i] == 0)
				--i;
			maxDistance = i;
		}
		
		void CountNeighborhood()
		{
			int_t i = 0;
			for (int_t j = 0; j <= maxDistance; ++j)
			{
				numNeighborhood[j] = 0;
				numDistance[j] = 0;
			}
			for (int_t j = 0; j < nSites; ++j)
			{
				numNeighborhood[Distance(i, j)] += 1;
				numDistance[Distance(i, j)] += 1;
			}
			for (int_t j = 1; j <= maxDistance; ++j)
				numNeighborhood[j] += numNeighborhood[j-1];
		}
		
	protected:
		int_t L;
		int_t nSites;
		int_t nBonds;
		int_t nDirections;
		int_t maxDistance;
		sub_lat_vector_t sublatVector;
		lookup_t distanceMap;
		vector_t distanceHistogram;
		int_t** neighborList;
		vector_t numNeighborhood;
		vector_t numDistance;
		bool fileIO;
};