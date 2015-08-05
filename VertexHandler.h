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
		
		VertexHandler(ConfigSpace_t& configSpace)
			: configSpace(configSpace)
		{
			std::size_t maxBufferSize = 20;
			nodeBuffer.resize(40);
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
			for (uint_t i = 0; i < nodeBuffer.size(); i+=2)
				std::cout << "(" << nodeBuffer[i].Site << " , " << nodeBuffer[i+1].Site << ", " << nodeBuffer[i].Tau << ", D = " << configSpace.lattice->Distance(nodeBuffer[i].Site, nodeBuffer[i+1].Site) << ") ";
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
			for (auto it = nodeBuffer.begin(); it != nodeBufferEnd; ++it)
				it->Worm = isWorm;
			nodes.insert(nodes.end(), nodeBuffer.begin(), nodeBufferEnd);
			if (isWorm)
			{
				uint_t n = std::distance(nodeBuffer.begin(), nodeBufferEnd);
				for (uint_t i = 0; i < n; ++i)
					wormNodes.push_back(nodes.size() - n + i);
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
			wormNodes.clear();
			/*
			for (auto it = indexBufferEnd; it != indexBuffer.begin();)
			{
				--it;
				wormNodes.erase(it);
			}
			*/
		}
		
		template<uint_t N>
		bool ReplaceWorm()
		{
			AddRandomIndicesToBuffer<N>();
			if (WormDistance() == 1)
			{
				nodes[wormNodes[0]].Worm = false;
				nodes[wormNodes[1]].Worm = false;
				nodes[indexBuffer[0]].Worm = true;
				nodes[indexBuffer[1]].Worm = true;
				wormNodes[0] = indexBuffer[0];
				wormNodes[1] = indexBuffer[1];
				return true;
			}
			else
				return false;
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
		int_t WormIndexBufferDistance()
		{
			if (N == 1)
				return configSpace.lattice->Distance(nodes[wormNodes[indexBuffer[0]]].Site, nodes[wormNodes[indexBuffer[1]]].Site);
			else if(N == 2)
			{			
				uint_t dist[4];
				uint_t r = configSpace.rng() * 2;
				for (uint_t i = 0; i < 4; ++i)
					dist[i] = configSpace.lattice->Distance(nodes[wormNodes[indexBuffer[2*r]]].Site, nodes[wormNodes[indexBuffer[i]]].Site);
				return *std::max_element(dist, dist+4);
			}
		}
		
		template<int_t N>
		bool WormIndexBufferDistance(uint_t distance)
		{
			if (N == 1)
			{
				return configSpace.lattice->Distance(nodes[wormNodes[indexBuffer[0]]].Site, nodes[wormNodes[indexBuffer[1]]].Site) == distance;
			}
			else if(N == 2)
			{			
				return (configSpace.lattice->Distance(nodes[wormNodes[indexBuffer[0]]].Site, nodes[wormNodes[indexBuffer[1]]].Site) == distance) && (configSpace.lattice->Distance(nodes[wormNodes[indexBuffer[2]]].Site, nodes[wormNodes[indexBuffer[3]]].Site) == distance);
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
			int_t r = configSpace.rng() * (2 * N - 1) + 1;
			for (uint_t i = 1; i < 2 * N; ++i)
			{
				uint_t nsite;
				if (N == 1)
				{
					nsite = configSpace.lattice->FromDistance(site, nhoodDist, configSpace.rng);
				}
				else if(N == 2)
				{
					if (i == r)
						nsite = configSpace.lattice->FromDistance(site, nhoodDist, configSpace.rng);
					else
						nsite = configSpace.lattice->FromNeighborhood(site, nhoodDist, configSpace.rng);
				}
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
						if (wormNodes[i] > wormNodes[*(it-1)])
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
		
		void RemoveBufferedVertices2(bool isWorm)
		{
			uint_t k = nodes.size();
			uint_t n = std::distance(indexBuffer.begin(), indexBufferEnd);
			for (uint_t i = 0; i < n; ++i)
			{
				if (isWorm)
				{
					if (nodes[k - i - 1].Worm)
					{
						auto wit = std::find(wormNodes.begin(), wormNodes.end(), k - i - 1);
						*wit = wormNodes[indexBuffer[n - i - 1]];
					}
					std::swap(nodes[wormNodes[indexBuffer[n - i - 1]]], nodes[k - i - 1]);
				}
				else
				{
					if (nodes[k - i - 1].Worm)
					{
						auto wit = std::find(wormNodes.begin(), wormNodes.end(), k - i - 1);
						*wit = indexBuffer[n - i - 1];
					}
					std::swap(nodes[indexBuffer[n - i - 1]], nodes[k - i - 1]);
				}
			}
			for (uint_t i = 0; i < n; ++i)
				nodes.erase(nodes.end() - 1);
			for (auto it = indexBufferEnd; it != indexBuffer.begin(); --it)
			{
				if (isWorm)
					wormNodes.erase(wormNodes.begin() + *(it-1));
			}
			std::sort(wormNodes.begin(), wormNodes.end());
		}
		
		void PermuteVertices(bool isWorm)
		{
			uint_t k = nodes.size();
			uint_t n = std::distance(indexBuffer.begin(), indexBufferEnd);
			for (uint_t i = 0; i < n; ++i)
			{
				if (isWorm)
				{
					if (nodes[k - i - 1].Worm)
					{
						auto wit = std::find(wormNodes.begin(), wormNodes.end(), k - i - 1);
						*wit = wormNodes[indexBuffer[n - i - 1]];
					}
					std::swap(nodes[wormNodes[indexBuffer[n - i - 1]]], nodes[k - i - 1]);
					wormNodes[indexBuffer[n - i - 1]] = k - i - 1;
				}
				else
				{
					if (nodes[k - i - 1].Worm)
					{
						auto wit = std::find(wormNodes.begin(), wormNodes.end(), k - i - 1);
						*wit = indexBuffer[n - i - 1];
					}
					std::swap(nodes[indexBuffer[n - i - 1]], nodes[k - i - 1]);
					if (nodes[k - i - 1].Worm)
					{
						auto wit = std::find(wormNodes.begin(), wormNodes.end(), indexBuffer[n - i - 1]);
						*wit = k - i - 1;
					}
				}
			}
			std::sort(wormNodes.begin(), wormNodes.end());
		}
		
		void RandomSwap()
		{
			uint_t r = configSpace.rng() * (nodes.size() / 2);
			while (!nodes[2*r].Worm)
			{
				r = configSpace.rng() * (nodes.size() / 2);
			}
			indexBuffer[0] = 2*r;
			indexBuffer[1] = 2*r+1;
			indexBufferEnd = indexBuffer.begin() + 2;
		}
		
		void RandomSwapWorm()
		{
			indexBuffer[0] = 0;
			indexBuffer[1] = 1;
			indexBufferEnd = indexBuffer.begin() + 2;
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

			/*
			for (uint_t i = 0; i < l; i+=2)
			{
				uint_t r = static_cast<uint_t>(configSpace.rng() * 2);
				wormShiftParity = 1.0;
				wormShiftParity *= (configSpace.lattice->Sublattice(nodeBuffer[i+r].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			
				nodeBuffer[i+r].Site = configSpace.lattice->RandomWalk(nodeBuffer[i+r].Site, 1, configSpace.rng);
				wormShiftParity *= (configSpace.lattice->Sublattice(nodeBuffer[i+r].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			}
			*/
			
			uint_t r = static_cast<uint_t>(configSpace.rng() * l);
			wormShiftParity = 1.0;
			wormShiftParity *= (configSpace.lattice->Sublattice(nodeBuffer[r].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			
			//nodeBuffer[r].Site = configSpace.lattice->RandomWalk(nodeBuffer[r].Site, 1, configSpace.rng);
			nodeBuffer[r].Site = configSpace.lattice->FromNeighborhood(nodeBuffer[r].Site, configSpace.nhoodDist, configSpace.rng);
			wormShiftParity *= (configSpace.lattice->Sublattice(nodeBuffer[r].Site) == ConfigSpace_t::Geometry_t::SublatticeType::A ? 1.0 : -1.0);
			
			
			value_t tau = nodeBuffer[0].Tau - 0.05 * configSpace.beta + configSpace.rng() * 0.1 * configSpace.beta;
			if (tau > configSpace.beta)
				tau -= configSpace.beta;
			else if (tau < 0.0)
				tau += configSpace.beta;
			for (uint_t i = 0; i < l; ++i)
				nodeBuffer[i].Tau = tau;
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
			{
				nodes[wormNodes[i]] = nodeBuffer[i];
			}
		}
		
		template<typename P>
		void ApplyWormShift(P& perm)
		{
			uint_t k = nodes.size();
			std::vector<node_t> nodesTmp(nodes);
			for (uint_t i = 0; i < k; ++i)
				nodes[i] = nodesTmp[perm[i]];
			
			for (uint_t i = 0; i < wormNodes.size(); ++i)
			{
				wormNodes[i] = k - wormNodes.size() + i;
				nodes[wormNodes[i]] = nodeBuffer[i];
			}
		}
		
		std::size_t Vertices()
		{
			return (nodes.size() - wormNodes.size()) / 2;
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
			for (uint_t i = 0; i < nodes.size(); ++i)
			{
				for (uint_t j = 0; j < i; ++j)
				{
					value_t dt = nodes[j].Tau - nodes[i].Tau;
					if (std::abs(dt) < configSpace.infinTau)
						dt += configSpace.infinTau;
					G(j, i) = configSpace.LookUpG0(nodes[j].Site, nodes[i].Site, dt);
					value_t sign = (configSpace.lattice->Sublattice(nodes[i].Site) == configSpace.lattice->Sublattice(nodes[j].Site) ? -1.0 : 1.0);
					G(i, j) = G(j, i) * sign;
					//G(i, j) = configSpace.LookUpG0(nodes[i].Site, nodes[j].Site, nodes[i].Tau - nodes[j].Tau - configSpace.infinTau);
				}
				G(i, i) = 0.0;
			}
		}
		
		template<typename U, typename V, typename A>
		void WoodburyAddVertices(U& u, V& v, A& a)
		{
			uint_t k = nodes.size();
			uint_t n = std::distance(nodeBuffer.begin(), nodeBufferEnd);
			for (uint_t i = 0; i < n; ++i)
			{
				for (uint_t j = 0; j < k; j+=2)
				{
					value_t dt = nodes[j].Tau - nodeBuffer[i].Tau;
					if (std::abs(dt) < configSpace.infinTau)
						dt += configSpace.infinTau;
					u(j, i) = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, dt);
					value_t sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodes[j].Site) ? -1.0 : 1.0);
					v(i, j) = u(j, i) * sign;
					
					dt = nodes[j + 1].Tau - nodeBuffer[i].Tau;
					if (std::abs(dt) < configSpace.infinTau)
						dt += configSpace.infinTau;
					u(j + 1, i) = configSpace.LookUpG0(nodes[j + 1].Site, nodeBuffer[i].Site, dt);
					sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodes[j + 1].Site) ? -1.0 : 1.0);
					v(i, j + 1) = u(j + 1, i) * sign;
				}
				for (uint_t j = 0; j < n; ++j)
				{
					if (i < j)
					{
						value_t dt = nodeBuffer[i].Tau - nodeBuffer[j].Tau;
						if (std::abs(dt) < configSpace.infinTau)
							dt += configSpace.infinTau;
						a(i, j) = configSpace.LookUpG0(nodeBuffer[i].Site, nodeBuffer[j].Site, dt);
						value_t sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodeBuffer[j].Site) ? -1.0 : 1.0);
						a(j, i) = a(i, j) * sign;
					}
				}
				a(i, i) = 0.0;
			}
		}
		
		template<typename U, typename V, typename A>
		void WoodburyAddShiftedVertices(U& u, V& v, A& a)
		{
			uint_t k = nodes.size();
			uint_t n = std::distance(nodeBuffer.begin(), nodeBufferEnd);
			for (uint_t i = 0; i < n; ++i)
			{
				uint_t c = 0;
				for (uint_t j = 0; j < k; j+=2)
				{
					if (nodes[j].Worm)
						continue;
					else
						c+=2;
					value_t dt = nodes[j].Tau - nodeBuffer[i].Tau;
					if (std::abs(dt) < configSpace.infinTau)
						dt += configSpace.infinTau;
					u(c, i) = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, dt);
					value_t sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodes[j].Site) ? -1.0 : 1.0);
					v(i, c) = u(c, i) * sign;
					
					dt = nodes[j + 1].Tau - nodeBuffer[i].Tau;
					if (std::abs(dt) < configSpace.infinTau)
						dt += configSpace.infinTau;
					u(c + 1, i) = configSpace.LookUpG0(nodes[j + 1].Site, nodeBuffer[i].Site, dt);
					sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodes[j + 1].Site) ? -1.0 : 1.0);
					v(i, c + 1) = u(c + 1, i) * sign;
				}
				for (uint_t j = 0; j < n; ++j)
				{
					if (i < j)
					{
						value_t dt = nodeBuffer[i].Tau - nodeBuffer[i].Tau;
						if (std::abs(dt) < configSpace.infinTau)
							dt += configSpace.infinTau;
					
						a(i, j) = configSpace.LookUpG0(nodeBuffer[i].Site, nodeBuffer[j].Site, dt);
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
			uint_t k = nodes.size();
			uint_t l = wormNodes.size();
			for (uint_t i = 0; i < l; ++i)
			{
				uint_t n = 0;
				for (uint_t j = 0; j < k; ++j)
				{
					if (!nodes[j].Worm)
					{
						value_t dt = nodes[j].Tau - nodes[wormNodes[i]].Tau;
						if (std::abs(dt) < configSpace.infinTau)
							dt += configSpace.infinTau;
						u(n, i) = configSpace.LookUpG0(nodes[j].Site, nodes[wormNodes[i]].Site, dt);
						value_t sign = (configSpace.lattice->Sublattice(nodes[wormNodes[i]].Site) == configSpace.lattice->Sublattice(nodes[j].Site) ? -1.0 : 1.0);
						v(i, n) = u(n, i) * sign;
						++n;
					}
				}
				for (uint_t j = 0; j < l; ++j)
				{
					if (i < j)
					{
						value_t dt = nodes[wormNodes[i]].Tau - nodes[wormNodes[j]].Tau;
						if (std::abs(dt) < configSpace.infinTau)
							dt += configSpace.infinTau;
						a(i, j) = configSpace.LookUpG0(nodes[wormNodes[i]].Site, nodes[wormNodes[j]].Site, dt);
						value_t sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodeBuffer[j].Site) ? -1.0 : 1.0);
						a(j, i) = a(i, j) * sign;
					}
				}
				a(i, i) = 0.0;
			}
		}

		template<typename U, typename V, typename A>
		void WoodburyShiftWorm(U& u, V& v, A& a)
		{
			uint_t k = nodes.size();
			uint_t l = wormNodes.size();
			for (uint_t i = 0; i < l; ++i)
			{
				uint_t n = 0;
				for (uint_t j = 0; j < k; ++j)
				{
					if (!nodes[j].Worm)
					{
						value_t dt = nodes[j].Tau - nodeBuffer[i].Tau;
						if (std::abs(dt) < configSpace.infinTau)
							dt += configSpace.infinTau;
						u(n, i) = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, dt);
						value_t sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodes[j].Site) ? -1.0 : 1.0);
						v(i, n) = u(n, i) * sign;
						++n;
					}
				}
				for (uint_t j = 0; j < l; ++j)
				{
					if (i < j)
					{
						value_t dt = nodeBuffer[i].Tau - nodeBuffer[j].Tau;
						if (std::abs(dt) < configSpace.infinTau)
							dt += configSpace.infinTau;
						a(i, j) = configSpace.LookUpG0(nodeBuffer[i].Site, nodeBuffer[j].Site, dt);
						value_t sign = (configSpace.lattice->Sublattice(nodeBuffer[i].Site) == configSpace.lattice->Sublattice(nodeBuffer[j].Site) ? -1.0 : 1.0);
						a(j, i) = a(i, j) * sign;
					}
				}
				a(i, i) = 0.0;
			}
		}
		
		template<typename U, typename V>
		void WoodburyShiftRowsCols(U& u, V& v)
		{
			uint_t k = nodes.size();
			uint_t l = wormNodes.size();
			value_t shiftedU, U;
			for (uint_t i = 0; i < l; ++i)
			{
				uint_t n = 0;
				for (uint_t j = 0; j < k; ++j)
				{
					if (j == wormNodes[i])
					{
						u(j, i) = 0.0;
						v(i, j) = 0.0;
						++n;
					}
					else if (!nodes[j].Worm)
					{
						if (j < wormNodes[i])
						{
							value_t dt = nodes[j].Tau - nodeBuffer[i].Tau;
							if (std::abs(dt) < configSpace.infinTau)
								dt += configSpace.infinTau;
							shiftedU = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, dt);
							
							dt = nodes[j].Tau - nodes[wormNodes[i]].Tau;
							if (std::abs(dt) < configSpace.infinTau)
								dt += configSpace.infinTau;
							U = configSpace.LookUpG0(nodes[j].Site, nodes[wormNodes[i]].Site, dt);
							u(j, i) = shiftedU - U;
						}
						else
						{
							value_t dt = nodes[j].Tau - nodeBuffer[i].Tau;
							if (std::abs(dt) < configSpace.infinTau)
								dt -= configSpace.infinTau;
							shiftedU = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, dt);
							
							dt = nodes[j].Tau - nodes[wormNodes[i]].Tau;
							if (std::abs(dt) < configSpace.infinTau)
								dt -= configSpace.infinTau;
							U = configSpace.LookUpG0(nodes[j].Site, nodes[wormNodes[i]].Site, dt);
							u(j, i) = shiftedU - U;
						}
						value_t shiftedV = shiftedU * (configSpace.lattice->Sublattice(nodes[j].Site) == configSpace.lattice->Sublattice(nodeBuffer[i].Site) ? -1.0 : 1.0);
						value_t V = U * (configSpace.lattice->Sublattice(nodes[j].Site) == configSpace.lattice->Sublattice(nodes[wormNodes[i]].Site) ? -1.0 : 1.0);
						v(i, j) = shiftedV - V;
					}
					else
					{
						if (j < wormNodes[i])
						{
							value_t dt = nodeBuffer[n].Tau - nodeBuffer[i].Tau;
							if (std::abs(dt) < configSpace.infinTau)
								dt += configSpace.infinTau;
							shiftedU = configSpace.LookUpG0(nodeBuffer[n].Site, nodeBuffer[i].Site, dt);
							
							dt = nodes[j].Tau - nodes[wormNodes[i]].Tau;
							if (std::abs(dt) < configSpace.infinTau)
								dt += configSpace.infinTau;
							U = configSpace.LookUpG0(nodes[j].Site, nodes[wormNodes[i]].Site, dt);
							u(j, i) = shiftedU - U;
						}
						else
						{
							value_t dt = nodeBuffer[n].Tau - nodeBuffer[i].Tau;
							if (std::abs(dt) < configSpace.infinTau)
								dt -= configSpace.infinTau;
							shiftedU = configSpace.LookUpG0(nodeBuffer[n].Site, nodeBuffer[i].Site, dt);
							
							dt = nodes[j].Tau - nodes[wormNodes[i]].Tau;
							if (std::abs(dt) < configSpace.infinTau)
								dt -= configSpace.infinTau;
							U = configSpace.LookUpG0(nodes[j].Site, nodes[wormNodes[i]].Site, dt);
							u(j, i) = shiftedU - U;
						}
						v(i, j) = 0.0;
						++n;
					}
				}
			}
		}
		
		void SwapVertexPosition(uint_t pos1, uint_t pos2)
		{
			if (nodes[pos1].Worm)
			{
				uint_t wormpos = std::find(wormNodes.begin(), wormNodes.end(), pos1) - wormNodes.begin();
				wormNodes[wormpos] = pos2;
				wormNodes[wormpos+1] = pos2+1;
			}
			if (nodes[pos2].Worm)
			{
				uint_t wormpos = std::find(wormNodes.begin(), wormNodes.end(), pos2) - wormNodes.begin();
				wormNodes[wormpos] = pos1;
				wormNodes[wormpos+1] = pos1+1;
			}
			std::swap(nodes[pos1], nodes[pos2]);
			std::swap(nodes[pos1+1], nodes[pos2+1]);
		}
		
		std::vector<std::size_t>& WormPositions()
		{
			return wormNodes;
		}
		
		std::vector<std::size_t>& IndexBuffer()
		{
			return indexBuffer;
		}
		
		template<typename M>
		void SwapRowsCols(M& m, uint_t i, uint_t j)
		{
			uint_t k = nodes.size();
			M cols = m.submat(0, i, k - 1, i + 1);
			m.submat(0, i, k - 1, i + 1) = m.submat(0, j, k - 1, j + 1);
			m.submat(0, j, k - 1, j + 1) = cols;
			M rows = m.submat(i, 0, i + 1, k - 1);
			m.submat(i, 0, i + 1, k - 1) = m.submat(j, 0, j + 1, k - 1);
			m.submat(j, 0, j + 1, k - 1) = rows;
		}
		
		template<typename M>
		void PermuteProgagatorMatrix(M& m, bool isWorm)
		{
			uint_t k = nodes.size();
			uint_t n = std::distance(indexBuffer.begin(), indexBufferEnd);
			for (uint_t i = 0; i < n; i+=2)
			{
				if (isWorm)
				{
					SwapRowsCols(m, wormNodes[indexBuffer[n - i - 2]], k - i - 2);
				}
				else
				{
					SwapRowsCols(m, indexBuffer[n - i - 2], k - i - 2);
				}
			}
		}
		
		template<typename M>
		void PermuteBackProgagatorMatrix(M& m, bool isWorm)
		{
			uint_t k = nodes.size();
			uint_t n = std::distance(indexBuffer.begin(), indexBufferEnd);
			for (uint_t i = 0; i < n; i+=2)
			{
				if (isWorm)
				{
					SwapRowsCols(m, wormNodes[indexBuffer[i]], k - n + i);
				}
				else
				{
					SwapRowsCols(m, indexBuffer[i], k - n + i);
				}
			}
		}
		
		template<typename M, typename N>
		void FillSMatrix(M& S, N& invG, bool isWorm)
		{
			uint_t n = std::distance(indexBuffer.begin(), indexBufferEnd);
			for (uint_t i = 0; i < n; i+=2)
			{
				for (uint_t j = 0; j < n; j+=2)
				{
					if (isWorm)
						//S.template block<2, 2>(i, j) = invG.template block<2, 2>(wormNodes[indexBuffer[i]], wormNodes[indexBuffer[j]]);
						S.submat(i, j, i+1, j+1) = invG.submat(wormNodes[indexBuffer[i]], wormNodes[indexBuffer[j]], wormNodes[indexBuffer[i]]+1, wormNodes[indexBuffer[j]]+1);
					else
						//S.template block<2, 2>(i, j) = invG.template block<2, 2>(indexBuffer[i], indexBuffer[j]);
						S.submat(i, j, i+1, j+1) = invG.submat(indexBuffer[i], indexBuffer[j], indexBuffer[i]+1, indexBuffer[j]+1);
				}
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
		uint_t sum = 0;
};
