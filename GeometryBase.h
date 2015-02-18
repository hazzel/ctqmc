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
#include <sys/stat.h>
#include <mpi.h>

inline bool FileExists(const std::string& name)
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

template<typename RNG, typename Int_t = std::int_fast32_t>
class GeometryBase
{
	public:
		typedef Int_t int_t;
		using site_t = std::tuple < int_t, int_t, int_t >;
		using lookup_t = LookUpTable < int_t, int_t, 2 >;
		using vector_t = std::vector< int_t >;
		enum SublatticeType {A, B};
		
	public:
		GeometryBase()
			: neighborList(0)
		{}
		
		virtual ~GeometryBase()
		{
			DeallocateNeighborList();
		}
		
		virtual void Resize(int_t l, RNG& rng) = 0;
		virtual SublatticeType Sublattice(int_t site) = 0;
		
		int_t Distance(int_t s1, int_t s2)
		{
			return distanceMap[s1][s2];
		}

		int_t DistanceHistogram(int_t distance)
		{
			return distanceHistogram[distance];
		}

		bool IsNeighbor(int_t s1, int_t s2)
		{
			return Distance(s1, s2) == 1;
		}

		int_t Sites()
		{
			return nSites;
		}

		int_t Bonds()
		{
			return nBonds;
		}
		
		int_t MaxDistance()
		{
			return maxDistance;
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
		
		int_t NeighborhoodCount(int_t distance)
		{
			return numNeighborhood[distance];
		}

		int_t RandomSite(RNG& rng)
		{
			return static_cast<int_t>(rng() * nSites);
		}

		void SaveToFile(const std::string& filename)
		{
			if (FileExists(filename))
				return;
			std::ofstream os(filename, std::ofstream::binary);
			os.write((char*)&maxDistance, sizeof(maxDistance));
			os.write((char*)&nSites, sizeof(nSites));
			for (int_t i = 0; i < nSites; ++i)
			{
				for (int_t j = 0; j < nSites; ++j)
				{
					os.write((char*)&distanceMap[i][j], sizeof(distanceMap[i][j]));
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
			os.close();
		}

		void ReadFromFile(const std::string& filename)
		{
			std::ifstream is(filename, std::ofstream::binary);
			if (is.is_open())
			{
				while(is.good())
				{
					is.read((char*)&maxDistance, sizeof(maxDistance));
					is.read((char*)&nSites, sizeof(nSites));
					this->numNeighborhood.resize(this->maxDistance + 1, 0);
					
					for (int_t i = 0; i < nSites; ++i)
					{
						for (int_t j = 0; j < nSites; ++j)
						{
							is.read((char*)&distanceMap[i][j], sizeof(distanceMap[i][j]));
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
					is.close();
				}
			}
		}

		void SyncMPI(const std::string& filename, const std::string& type)
		{
			int file_free = 0;
			int np;
			int proc_id;
			MPI_Comm_size(MPI_COMM_WORLD, &np);
			MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
			MPI_Status status;

			if (proc_id == 1)
				file_free = 1;
			else
				MPI_Recv(&file_free, 1, MPI_INT, proc_id-1, 1, MPI_COMM_WORLD, &status);
			if (file_free == 1)
			{
				if (type == "save")
					SaveToFile(filename);
				else if (type == "read")
					ReadFromFile(filename);
			}
			if (proc_id != np-1)
				MPI_Send(&file_free, 1, MPI_INT, proc_id+1, 1, MPI_COMM_WORLD);
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
			for (int_t j = 0; j < nSites; ++j)
				numNeighborhood[Distance(i, j)] += 1;
			for (int_t j = 1; j <= maxDistance; ++j)
				numNeighborhood[j] += numNeighborhood[j-1];
		}
		
	protected:
		int_t L;
		int_t nSites;
		int_t nBonds;
		int_t nDirections;
		int_t maxDistance;
		lookup_t distanceMap;
		vector_t distanceHistogram;
		int_t** neighborList;
		vector_t numNeighborhood;
		bool fileIO = true;
};