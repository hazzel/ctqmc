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
			std::size_t maxBufferSize = 4;
			nodeBuffer.resize(maxBufferSize);
			nodeBufferEnd = nodeBuffer.end();
			indexBuffer.resize(maxBufferSize);
			indexBufferEnd = indexBuffer.end();
		}
		
		void PrintVertices()
		{
			value_t minTauDiff = configSpace.beta;
			value_t minWormTauDiff = configSpace.beta;
			for (uint_t i = 0; i < nodes.size(); i+=2)
			{
				std::cout << "(" << nodes[i].Site << " , " << nodes[i+1].Site << ", " << nodes[i].Tau << ", " << nodes[i].Worm << ") ";
			}
			for (uint_t i = 2; i < nodes.size(); i+=2)
			{
				value_t diff = std::abs(nodes[i-2].Tau - nodes[i].Tau);
				if (diff < minTauDiff)
					minTauDiff = diff;
			}
			std::cout << std::endl;
			std::cout << "MinTauDiff: " << minTauDiff / configSpace.dtau << std::endl;
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
		
		void AddBufferedVertices(bool isWorm)
		{
			nodes.insert(nodes.end(), nodeBuffer.begin(), nodeBufferEnd);
			if (isWorm)
			{
				uint_t n = std::distance(nodeBuffer.begin(), nodeBufferEnd);
				for (uint_t i = 0; i < n; ++i)
					wormNodes.push_back(nodes.size() - n + i);
			}
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
		void AddRandomWormsToBuffer()
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
				uint_t nsite = configSpace.lattice->FromNeighborhood(site, configSpace.nhoodDist, configSpace.rng);
				nodeBuffer[i] = node_t(nsite, tau, true);
			}
			nodeBufferEnd = nodeBuffer.begin() + 2 * N;
		}
		
		void RemoveBufferedVertices(bool isWorm)
		{
			if (isWorm)
			{
				for (auto it = indexBufferEnd; it != indexBuffer.begin(); --it)
				{
					nodes.erase(nodes.begin() + wormNodes[*(it-1)]);
				}
				for (auto it = indexBufferEnd; it != indexBuffer.begin(); --it)
				{
					for (uint_t i = 0; i < wormNodes.size(); ++i)
					{
						if (wormNodes[i] > *(it-1))
							--wormNodes[i];
					}
				}
				for (auto it = indexBufferEnd; it != indexBuffer.begin(); --it)
					wormNodes.erase(wormNodes.begin() + *(it-1));
			}
			else
			{
				for (auto it = indexBufferEnd; it != indexBuffer.begin(); --it)
				{
					for (uint_t i = 0; i < wormNodes.size(); ++i)
					{
						if (wormNodes[i] > *(it-1))
							--wormNodes[i];
					}
					nodes.erase(nodes.begin() + *(it-1));
				}
			}
		}
		
		template<int_t N>
		void AddRandomIndicesToBuffer()
		{
			for (uint_t i = 0; i < 2 * N; ++i)
				indexBuffer[i] = nodes.size();
			for (uint_t i = 0; i < N;)
			{
				uint_t r = configSpace.rng() * nodes.size() / 2;
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
		
		//TODO: simpplify
		template<int_t N>
		void AddRandomWormIndicesToBuffer()
		{
			for (uint_t i = 0; i < 2 * N; ++i)
				indexBuffer[i] = wormNodes.size();
			for (uint_t i = 0; i < N;)
			{
				uint_t r = configSpace.rng() * wormNodes.size() / 2;
				if (std::find(indexBuffer.begin(), indexBuffer.begin() + 2*N, 2 * r) == indexBuffer.begin() + 2*N)
				{
					indexBuffer[2*i] = 2 * r;
					indexBuffer[2*i + 1] = 2 * r + 1;
					++i;
				}
			}
			indexBufferEnd = indexBuffer.begin() + 2 * N;
			std::sort(indexBuffer.begin(), indexBufferEnd);
		}

		void ShiftWorm()
		{
			uint_t l = wormNodes.size();
			uint_t r = static_cast<uint_t>(configSpace.rng() * l);
			nodeBuffer[0] = wormNodes[r];
			nodeBufferEnd = nodeBuffer.begin() + 1;
			indexBuffer[0] = r;
			indexBufferEnd = indexBuffer.begin() + 1;
			wormShiftParity = 1.0;
			wormShiftParity *= (configSpace.lattice->Sublattice(wormNodes[r].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			
			wormNodes[r].Site = configSpace.lattice->RandomWalk(wormNodes[r].Site, 1, configSpace.rng);
			wormShiftParity *= (configSpace.lattice->Sublattice(wormNodes[r].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			wormNodes[r].Tau += -0.05 * configSpace.beta + configSpace.rng() * 0.1 * configSpace.beta;
			if (wormNodes[r].Tau > configSpace.beta)
				wormNodes[r].Tau -= configSpace.beta;
			else if (wormNodes[r].Tau < 0.0)
				wormNodes[r].Tau += configSpace.beta;
			
			for (uint_t i = 0; i < wormNodes.size(); ++i)
				wormNodes[i].Tau = wormNodes[r].Tau;
		}
		
		void UndoWormShift()
		{
			/*
			wormNodes[indexBuffer[0]].Site = nodeBuffer[0].Site;
			for (uint_t i = 0; i < wormNodes.size(); ++i)
				wormNodes[i].Tau = nodeBuffer[0].Tau;
			*/
		}
		
		std::size_t Vertices()
		{
			return (nodes.size() - wormNodes.size()) / 2;
		}
		
		std::size_t Worms()
		{
			return wormNodes.size() / 2;
		}
		
		template<typename Matrix>
		void PropagatorMatrix(Matrix& G)
		{
			for (uint_t i = 0; i < nodes.size(); ++i)
			{
				for (uint_t j = 0; j < i; ++j)
				{
					G(j, i) = configSpace.LookUpG0(nodes[j].Site, nodes[i].Site, nodes[j].Tau - nodes[i].Tau + configSpace.infinTau);
					G(i, j) = configSpace.LookUpG0(nodes[i].Site, nodes[j].Site, nodes[i].Tau - nodes[j].Tau - configSpace.infinTau);
				}
				G(i, i) = 0.0;
			}
		}
		
		template<typename U, typename V, typename A>
		void WoodburyAddVertices(U& u, V& v, A& a)
		{
			uint_t k = nodes.size();
			uint_t n = a.cols();
			for (uint_t i = 0; i < n; ++i)
			{
				for (uint_t j = 0; j < k; ++j)
				{
					u(j, i) = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, nodes[j].Tau - nodeBuffer[i].Tau + configSpace.infinTau);
					v(i, j) = configSpace.LookUpG0(nodeBuffer[i].Site, nodes[j].Site, nodeBuffer[i].Tau - nodes[j].Tau - configSpace.infinTau);
				}
				for (uint_t j = 0; j < n; ++j)
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
		
		template<typename U, typename V, typename A, typename P>
		void WoodburyRemoveVertices(U& u, V& v, A& a, const P& perm)
		{
			uint_t k = nodes.size();
			uint_t l = wormNodes.size();
			uint_t n = a.cols();
			for (uint_t i = 0; i < n; ++i)
			{
				for (uint_t j = 0; j < k - n; ++j)
				{
					u(j, i) = configSpace.LookUpG0(nodes[perm[j]].Site, nodes[indexBuffer[i]].Site, nodes[perm[j]].Tau - nodes[indexBuffer[i]].Tau + configSpace.infinTau);
					v(i, j) = configSpace.LookUpG0(nodes[indexBuffer[i]].Site, nodes[perm[j]].Site, nodes[indexBuffer[i]].Tau - nodes[perm[j]].Tau - configSpace.infinTau);
				}
				for (uint_t j = 0; j < n; ++j)
				{
					if (i < j)
					{
						a(i, j) = configSpace.LookUpG0(nodes[indexBuffer[i]].Site, nodes[indexBuffer[j]].Site, nodes[indexBuffer[i]].Tau - nodes[indexBuffer[j]].Tau + configSpace.infinTau);
						a(j, i) = configSpace.LookUpG0(nodes[indexBuffer[j]].Site, nodes[indexBuffer[i]].Site, nodes[indexBuffer[j]].Tau - nodes[indexBuffer[i]].Tau - configSpace.infinTau);
					}
				}
				a(i, i) = 0.0;
			}
		}

		template<typename U, typename V>
		void WoodburyShiftWorm(U& u, V& v)
		{
			uint_t k = nodes.size();
			for (uint_t i = 0; i < k - 1; ++i)
			{
				u(i, 0) = configSpace.LookUpG0(nodes[i].Site, nodeBuffer[0].Site, nodes[i].Tau - nodeBuffer[0].Tau + configSpace.infinTau);
				v(0, i) = configSpace.LookUpG0(nodeBuffer[0].Site, nodes[i].Site, nodeBuffer[0].Tau - nodes[i].Tau - configSpace.infinTau);
			}
		}
		
		template<typename P>
		void PermutationMatrix(P& perm, bool isWorm)
		{
			uint_t cnt = 0;
			if (isWorm)
			{
				std::vector<std::size_t> buf;
				for (auto it = indexBuffer.begin(); it != indexBufferEnd; ++it)
					buf.push_back(wormNodes[*it]);

				for (uint_t i = 0; i < perm.size(); ++i)
				{
					if (find(buf.begin(), buf.end(), i) == buf.end())
					{
						perm[cnt] = i;
						++cnt;
					}
				}
				int i = 0;
				for (auto it = buf.begin(); it != buf.end(); ++it)
				{
					perm[cnt + i] = *it;
					++i;
				}
			}
			else
			{
				for (uint_t i = 0; i < perm.size(); ++i)
				{
					if (find(indexBuffer.begin(), indexBufferEnd, i) == indexBufferEnd)
					{
						perm[cnt] = i;
						++cnt;
					}
				}
				int i = 0;
				for (auto it = indexBuffer.begin(); it != indexBufferEnd; ++it)
				{
					perm[cnt + i] = *it;
					++i;
				}
			}
		}
			
		void Serialize(odump& d)
		{
			d.write(nodes);
			d.write(wormNodes);
		}
		
		void Serialize(idump& d)
		{
			d.read(nodes);
			d.read(wormNodes);
		}
	private:
		ConfigSpace_t& configSpace;
		std::vector<node_t> nodes;
		std::vector<std::size_t> wormNodes;
		std::vector<node_t> nodeBuffer;
		typename std::vector<node_t>::iterator nodeBufferEnd;
		std::vector<std::size_t> indexBuffer;
		typename std::vector<std::size_t>::iterator indexBufferEnd;
		value_t wormShiftParity;
};