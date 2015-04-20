#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <list>
#include <utility>
#include <cmath>
#include <numeric>
#include <cstdint>
#include "LookUpTable.h"
#include "measurements.h"
#include "random.h"
#include "parser.h"
#include "types.h"
#include "ConfigSpace.h"
#include "UpdateHandler.h"
#include "Eigen/Dense"

template<typename Index_t, typename Value_t>
struct Node
{
	Node()
		: Site(0), Tau(0.0), Worm(false)
	{}
	Node(Index_t site, Value_t tau, bool worm)
		: Site(site), Tau(tau), Worm(worm)
	{}
	
	Index_t Site;
	Value_t Tau;
	bool Worm;
};

template<typename ConfigSpace_t>
class VertexHandler
{
	public:
		typedef typename ConfigSpace_t::uint_t uint_t;
		typedef typename ConfigSpace_t::int_t int_t;
		typedef typename ConfigSpace_t::value_t value_t;
		typedef Node<uint_t, value_t> node_t;
		
		struct Configuration
		{
			std::vector<node_t> nodes;
			std::vector<node_t> wormNodes;
		};
		
		VertexHandler(ConfigSpace_t& configSpace)
			: configSpace(configSpace)
		{
			std::size_t maxBufferSize = 100;
			nodeBuffer.resize(maxBufferSize);
			nodeBufferEnd = nodeBuffer.end();
			indexBuffer.resize(maxBufferSize);
			indexBufferEnd = indexBuffer.end();
			nodes.resize(500);
		}
		
		void PrintVertices()
		{
			value_t minTauDiff = configSpace.beta;
			value_t minWormTauDiff = configSpace.beta;
			for (uint_t i = 0; i < nodeNumber; i+=2)
			{
				std::cout << "(" << nodes[i].Site << " , " << nodes[i+1].Site << ", " << nodes[i].Tau << ", " << nodes[i].Worm << ") ";
			}
			std::cout << std::endl;
			/*
			for (uint_t i = 2; i < nodeNumber; i+=2)
			{
				value_t diff = std::abs(nodes[i-2].Tau - nodes[i].Tau);
				if (diff < minTauDiff)
					minTauDiff = diff;
			}
			std::cout << "MinTauDiff: " << minTauDiff / configSpace.dtau << std::endl;
			*/
		}

		void PrintWormVertices()
		{
			for (auto node : wormNodes)
				std::cout << node << " , ";
			std::cout << " -> ";
			for (auto node : wormNodes)
				std::cout << "(" << nodes[node].Site << " , " << nodes[node].Tau << ") ";
			std::cout << std::endl;
		}
		
		void PrintVertexBuffer()
		{
			for (auto node : nodeBuffer)
				std::cout << "(" << node.Site << " , " << node.Tau << ") ";
			std::cout << std::endl;
		}
		
		void PrintIndexBuffer()
		{
			for (auto index : indexBuffer)
				std::cout << index << " ";
			std::cout << std::endl;
		}

		template<typename Matrix_t>
		void PrintMatrix(const Matrix_t& M)
		{
			using namespace Eigen;
			IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
			IOFormat CleanFmt(1, 0, ", ", "\n", "[", "]");
			IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
			IOFormat HeavyFmt(FullPrecision, 0, ", ", ",\n", "{", "}", "{", "}");
			//std::cout << M.format(CommaInitFmt) << std::endl;
			std::cout << M.format(CleanFmt) << std::endl;
			//std::cout << M.format(OctaveFmt) << std::endl;
			//std::cout << M.format(HeavyFmt) << std::endl;
		}
		
		void AddBufferedVertices(bool isWorm)
		{
			if (nodeNumber == nodes.size())
				nodes.resize(nodeNumber + 500);
			for (auto it = nodeBuffer.begin(); it != nodeBufferEnd; ++it)
			{
				it->Worm = isWorm;
				nodes[nodeNumber] = *it;
				++nodeNumber;
			}
			if (isWorm)
			{
				uint_t n = std::distance(nodeBuffer.begin(), nodeBufferEnd);
				for (uint_t i = 0; i < n; ++i)
					wormNodes.push_back(nodeNumber - n + i);
				std::sort(wormNodes.begin(), wormNodes.end());
			}
		}

		void OpenUpdate()
		{
			for (auto it = indexBuffer.begin(); it != indexBufferEnd; ++it)
			{
				nodes[*it].Worm = true;
				wormNodes.push_back(*it);
			}
		}

		void CloseUpdate()
		{
			for (auto it = indexBuffer.begin(); it != indexBufferEnd; ++it)
				nodes[wormNodes[*it]].Worm = false;
			//wormNodes only contains one worm : W2 -> Z
			wormNodes.clear();
		}
		
		value_t VertexBufferParity()
		{
			value_t parity = 1.0;
			for (auto it = nodeBuffer.begin(); it != nodeBufferEnd; ++it)
				parity *= (configSpace.lattice->Sublattice(it->Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			return parity;
		}

		value_t WormParity()
		{
			value_t parity = 1.0;
			for (auto node : wormNodes)
				parity *= (configSpace.lattice->Sublattice(nodes[node].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			return parity;
		}
	
		value_t WormIndexBufferParity()
		{
			value_t parity = 1.0;
			for (auto it = indexBuffer.begin(); it != indexBufferEnd; ++it)
				parity *= (configSpace.lattice->Sublattice(nodes[wormNodes[*it]].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			return parity;
		}
		
		value_t WormShiftParity()
		{
			return wormShiftParity;
		}
		
		uint_t WormDistance()
		{
			return configSpace.lattice->Distance(nodes[wormNodes[0]].Site, nodes[wormNodes[1]].Site);
		}

		template<int_t N>
		bool WormIndexBufferDistance()
		{
			if (N == 1)
				return configSpace.lattice->Distance(nodes[wormNodes[indexBuffer[0]]].Site, nodes[wormNodes[indexBuffer[1]]].Site) <= configSpace.nhoodDist;
			else if(N == 2)
			{			
				uint_t dist[4];
				uint_t r = configSpace.rng() * 2;
				for (uint_t i = 0; i < 4; ++i)
					dist[i] = configSpace.lattice->Distance(nodes[wormNodes[indexBuffer[2*r]]].Site, nodes[wormNodes[indexBuffer[i]]].Site);
				return (*std::max_element(dist, dist+4)) <= configSpace.nhoodDist;
			}
		}
		
		template<int_t N>
		void AddRandomVerticesToBuffer()
		{
			for (uint_t i = 0; i < N; ++i)
			{
				uint_t site = configSpace.lattice->RandomSite(configSpace.rng);
				value_t tau = configSpace.rng() * configSpace.beta;
				nodeBuffer[2*i] = node_t(site, tau, false);
				nodeBuffer[2*i+1] = node_t(configSpace.lattice->RandomWalk(site, 1, configSpace.rng), tau, false);
			}
			nodeBufferEnd = nodeBuffer.begin() + 2 * N;
		}
		
		template<int_t N>
		void AddRandomWormsToBuffer(uint_t nhoodDist)
		{
			uint_t site = configSpace.lattice->RandomSite(configSpace.rng);
			value_t tau;
			if (wormNodes.size() > 0)
				tau = nodes[wormNodes[0]].Tau;
			else
				tau = configSpace.rng() * configSpace.beta;
			nodeBuffer[0] = node_t(site, tau, true);
			for (uint_t i = 1; i < 2 * N; ++i)
			{
				uint_t nsite = configSpace.lattice->FromNeighborhood(site, nhoodDist, configSpace.rng);
				nodeBuffer[i] = node_t(nsite, tau, true);
			}
			nodeBufferEnd = nodeBuffer.begin() + 2 * N;
		}
		
		void RemoveBufferedVertices(bool isWorm)
		{
			uint_t n = std::distance(indexBuffer.begin(), indexBufferEnd);
			for (uint_t i = 0; i < n; ++i)
			{
				if (isWorm)
				{
					if (nodes[nodeNumber - i - 1].Worm)
					{
						auto wit = std::find(wormNodes.begin(), wormNodes.end(), nodeNumber - i - 1);
						*wit = wormNodes[indexBuffer[n - i - 1]];
					}
					std::swap(nodes[wormNodes[indexBuffer[n - i - 1]]], nodes[nodeNumber - i - 1]);
				}
				else
				{
					if (nodes[nodeNumber - i - 1].Worm)
					{
						auto wit = std::find(wormNodes.begin(), wormNodes.end(), nodeNumber - i - 1);
						*wit = indexBuffer[n - i - 1];
					}
					std::swap(nodes[indexBuffer[n - i - 1]], nodes[nodeNumber - i - 1]);
				}
			}
			nodeNumber -= n;
			if (isWorm)
			{
				for (auto it = indexBufferEnd; it != indexBuffer.begin(); --it)
					wormNodes.erase(wormNodes.begin() + *(it-1));
			}
			std::sort(wormNodes.begin(), wormNodes.end());
		}
		
		template<int_t N>
		void AddRandomIndicesToBuffer()
		{
			for (uint_t i = 0; i < 2 * N; ++i)
				indexBuffer[i] = nodeNumber;
			for (uint_t i = 0; i < N;)
			{
				uint_t r = configSpace.rng() * nodeNumber / 2;
				if ((!nodes[2 * r].Worm) && std::find(indexBuffer.begin(), indexBuffer.begin() + 2*N, 2 * r) == indexBuffer.begin() + 2*N)
				{
					indexBuffer[2*i] = 2 * r;
					indexBuffer[2*i + 1] = 2 * r + 1;
					++i;
				}
			}
			indexBufferEnd = indexBuffer.begin() + 2 * N;
			std::sort(indexBuffer.begin(), indexBufferEnd);
		}
		
		template<int_t N>
		void AddRandomWormIndicesToBuffer()
		{
			if (2 * N == wormNodes.size())
			{
				for (uint_t i = 0; i < 2 * N; ++i)
					indexBuffer[i] = i;
			}
			//must be N=1, W=2 -> chose one worm at random
			else
			{
				uint_t r = configSpace.rng() * wormNodes.size() / 2;
				indexBuffer[0] = 2 * r;
				indexBuffer[1] = 2 * r + 1;
			}
			indexBufferEnd = indexBuffer.begin() + 2 * N;
		}

		void ShiftWormToBuffer()
		{
			uint_t l = wormNodes.size();
			for (uint_t i = 0; i < l; ++i)
			{
				nodeBuffer[i] = nodes[wormNodes[i]];
				nodeBuffer[i+l] = nodes[wormNodes[i]];
				indexBuffer[i] = i;
			}
			nodeBufferEnd = nodeBuffer.begin() + l;
			indexBufferEnd = indexBuffer.begin() + l;

			uint_t r = static_cast<uint_t>(configSpace.rng() * l);
			wormShiftParity = 1.0;
			wormShiftParity *= (configSpace.lattice->Sublattice(nodeBuffer[r].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			
			nodeBuffer[r].Site = configSpace.lattice->RandomWalk(nodeBuffer[r].Site, 1, configSpace.rng);
			wormShiftParity *= (configSpace.lattice->Sublattice(nodeBuffer[r].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			nodeBuffer[r].Tau += -0.05 * configSpace.beta + configSpace.rng() * 0.1 * configSpace.beta;
			if (nodeBuffer[r].Tau > configSpace.beta)
				nodeBuffer[r].Tau -= configSpace.beta;
			else if (nodeBuffer[r].Tau < 0.0)
				nodeBuffer[r].Tau += configSpace.beta;
			
			for (uint_t i = 0; i < l; ++i)
				nodeBuffer[i].Tau = nodeBuffer[r].Tau;
		}

		template<int_t W>
		void RestoreAfterShift()
		{
			uint_t l = 2 * W;
			for (uint_t i = 0; i < l; ++i)
				nodeBuffer[i] = nodeBuffer[i+l];
			nodeBufferEnd = nodeBuffer.begin() + l;
		}
		
		void ApplyWormShift()
		{
			for (uint_t i = 0; i < wormNodes.size(); ++i)
				nodes[wormNodes[i]] = nodeBuffer[i];
		}
		
		std::size_t Vertices()
		{
			return (nodeNumber - wormNodes.size()) / 2;
		}
		
		std::size_t Worms()
		{
			return wormNodes.size() / 2;
		}

		void Clear()
		{
			nodes.clear();
			wormNodes.clear();
		}
		
		template<typename Matrix>
		void PropagatorMatrix(Matrix& G)
		{
			for (uint_t i = 0; i < nodeNumber; ++i)
			{
				for (uint_t j = 0; j < i; ++j)
				{
					G(j, i) = configSpace.LookUpG0(nodes[j].Site, nodes[i].Site, nodes[j].Tau - nodes[i].Tau + configSpace.infinTau);
					value_t sign = (configSpace.lattice->Sublattice(nodes[i].Site) == configSpace.lattice->Sublattice(nodes[j].Site) ? -1.0 : 1.0);
					G(i, j) = G(j, i) * sign;
				}
				G(i, i) = 0.0;
			}
		}
		
		template<typename U, typename V, typename A>
		void WoodburyAddVertices(U& u, V& v, A& a)
		{
			uint_t k = nodeNumber;
			uint_t n = a.n_cols;
			for (uint_t i = 0; i < n; ++i)
			{
				for (uint_t j = 0; j < k; ++j)
				{
					u(j, i) = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, nodes[j].Tau - nodeBuffer[i].Tau + configSpace.infinTau);
					value_t sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodes[j].Site) ? -1.0 : 1.0);
					v(i, j) = u(j, i) * sign;
				}
				for (uint_t j = 0; j < n; ++j)
				{
					if (i < j)
					{
						a(i, j) = configSpace.LookUpG0(nodeBuffer[i].Site, nodeBuffer[j].Site, nodeBuffer[i].Tau - nodeBuffer[j].Tau + configSpace.infinTau);
						value_t sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodeBuffer[j].Site) ? -1.0 : 1.0);
						a(j, i) = a(i, j) * sign;
					}
				}
				a(i, i) = 0.0;
			}
		}

		template<typename U, typename V, typename A>
		void WoodburyWorm(U& u, V& v, A& a)
		{
			uint_t k = nodeNumber;
			uint_t l = wormNodes.size();
			for (uint_t i = 0; i < l; ++i)
			{
				uint_t n = 0;
				for (uint_t j = 0; j < k; ++j)
				{
					if (!nodes[j].Worm)
					{
						u(n, i) = configSpace.LookUpG0(nodes[j].Site, nodes[wormNodes[i]].Site, nodes[j].Tau - nodes[wormNodes[i]].Tau + configSpace.infinTau);
						v(i, n) = configSpace.LookUpG0(nodes[wormNodes[i]].Site, nodes[j].Site, nodes[wormNodes[i]].Tau - nodes[j].Tau - configSpace.infinTau);
						++n;
					}
				}
				for (uint_t j = 0; j < l; ++j)
				{
					if (i < j)
					{
						a(i, j) = configSpace.LookUpG0(nodes[wormNodes[i]].Site, nodes[wormNodes[j]].Site, nodes[wormNodes[i]].Tau - nodes[wormNodes[j]].Tau + configSpace.infinTau);
						a(j, i) = configSpace.LookUpG0(nodes[wormNodes[j]].Site, nodes[wormNodes[i]].Site, nodes[wormNodes[j]].Tau - nodes[wormNodes[i]].Tau - configSpace.infinTau);
					}
				}
				a(i, i) = 0.0;
			}
		}

		template<typename U, typename V, typename A>
		void WoodburyShiftWorm(U& u, V& v, A& a)
		{
			uint_t k = nodeNumber;
			uint_t l = wormNodes.size();
			for (uint_t i = 0; i < l; ++i)
			{
				uint_t n = 0;
				for (uint_t j = 0; j < k; ++j)
				{
					if (!nodes[j].Worm)
					{
						u(n, i) = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, nodes[j].Tau - nodeBuffer[i].Tau + configSpace.infinTau);
						v(i, n) = configSpace.LookUpG0(nodeBuffer[i].Site, nodes[j].Site, nodeBuffer[i].Tau - nodes[j].Tau - configSpace.infinTau);
						++n;
					}
				}
				for (uint_t j = 0; j < l; ++j)
				{
					if (i < j)
					{
						a(i, j) = configSpace.LookUpG0(nodeBuffer[i].Site, nodeBuffer[j].Site, nodeBuffer[i].Tau - nodeBuffer[j].Tau + configSpace.infinTau);
						a(j, i) = configSpace.LookUpG0(nodeBuffer[j].Site, nodeBuffer[i].Site, nodeBuffer[j].Tau - nodeBuffer[i].Tau - configSpace.infinTau);
					}
				}
				a(i, i) = 0.0;
			}
		}
		
		template<typename P>
		void PermutationMatrix(P& perm, bool isWorm)
		{
			for (uint_t i = 0; i < perm.size(); ++i)
				perm[i] = i;
			uint_t n = std::distance(indexBuffer.begin(), indexBufferEnd);
			for (uint_t i = 0; i < n; ++i)
			{
				if (isWorm)
				{	
					perm[wormNodes[indexBuffer[n - i - 1]]] = perm[nodeNumber - i - 1];
					perm[nodeNumber - i - 1] = wormNodes[indexBuffer[n - i - 1]];
				}
				else
				{
					perm[indexBuffer[n - i - 1]] = perm[nodeNumber - i - 1];
					perm[nodeNumber - i - 1] = indexBuffer[n - i - 1];
				}
			}
		}
		/*
		template<typename P>
		void PermutationMatrix(P& perm, bool isWorm)
		{
			for (uint_t i = 0; i < perm.size(); ++i)
				perm[i] = i;
			uint_t n = std::distance(indexBuffer.begin(), indexBufferEnd);
			for (uint_t i = 0; i < n; i+=2)
			{
				if (isWorm)
				{
					perm[wormNodes[indexBuffer[n - i - 2]]] = perm[perm.size() - i - 2];
					perm[perm.size() - i - 2] = wormNodes[indexBuffer[n - i - 2]];
				
					perm[wormNodes[indexBuffer[n - i - 1]]] = perm[perm.size() - i - 1];
					perm[perm.size() - i - 1] = wormNodes[indexBuffer[n - i - 1]];
				}
				else
				{
					perm[indexBuffer[n - i - 2]] = perm[perm.size() - i - 2];
					perm[perm.size() - i - 2] = indexBuffer[n - i - 2];
					
					perm[indexBuffer[n - i - 1]] = perm[perm.size() - i - 1];
					perm[perm.size() - i - 1] = indexBuffer[n - i - 1];
				}
			}
		}
		*/
			
		void Serialize(odump& d)
		{
			d.write(nodeNumber);
			d.write(nodes);
			d.write(wormNodes);
		}
		
		void Serialize(idump& d)
		{
			d.read(nodeNumber);
			d.read(nodes);
			d.read(wormNodes);
		}
	private:
		ConfigSpace_t& configSpace;
		std::vector<node_t> nodes;
		std::vector<std::size_t> wormNodes;
		std::size_t nodeNumber = 0;
		std::size_t wormNumber = 0;
		std::vector<node_t> nodeBuffer;
		typename std::vector<node_t>::iterator nodeBufferEnd;
		std::vector<std::size_t> indexBuffer;
		typename std::vector<std::size_t>::iterator indexBufferEnd;
		value_t wormShiftParity;
};