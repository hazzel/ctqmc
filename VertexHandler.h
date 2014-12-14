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
//#define EIGEN_USE_MKL_ALL
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"
#include <flens/flens.cxx>
#include "Permutation.h"
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
		: Site(0), Tau(0.0)
	{}
	Node(Index_t site, Value_t tau)
		: Site(site), Tau(tau)
	{}
	
	Index_t Site;
	Value_t Tau;
};

template<typename ConfigSpace_t>
class VertexHandler
{
	public:
		typedef typename ConfigSpace_t::uint_t uint_t;
		typedef typename ConfigSpace_t::int_t int_t;
		typedef typename ConfigSpace_t::value_t value_t;
		typedef Node<uint_t, value_t> node_t;
		
		typedef flens::GeMatrix< flens::FullStorage<value_t, flens::ColMajor> > GeMatrix;
		typedef typename GeMatrix::IndexType IndexType;
		typedef flens::DenseVector< flens::Array<IndexType> > IndexVector;
		
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
			std::cout << std::endl;
			for (uint_t i = 0; i < nodes.size(); i+=2)
				std::cout << "(" << nodes[i].Site << " , " << nodes[i+1].Site << ", " << nodes[i].Tau << ") ";
			std::cout << std::endl;
		}
		
		void PrintVertexBuffer()
		{
			std::cout << std::endl;
			for (auto node : nodeBuffer)
				std::cout << "(" << node.Site << " , " << node.Tau << ") ";
			std::cout << std::endl;
		}
		
		void PrintIndexBuffer()
		{
			std::cout << std::endl;
			for (auto index : indexBuffer)
				std::cout << index << " ";
			std::cout << std::endl;
		}
		
		void AddBufferedVertices()
		{
			nodes.insert(nodes.end(), nodeBuffer.begin(), nodeBufferEnd);
		}
		
		template<int_t N>
		void AddRandomVerticesToBuffer()
		{
			for (uint_t i = 0; i < N; ++i)
			{
				uint_t site = configSpace.lattice.RandomSite(configSpace.rng);
				value_t tau = configSpace.rng() * configSpace.beta;
				nodeBuffer[2*i] = node_t(site, tau);
				nodeBuffer[2*i+1] = node_t(configSpace.lattice.RandomWalk(site, 1, configSpace.rng), tau);
			}
			nodeBufferEnd = nodeBuffer.begin() + 2 * N;
		}
		
		void RemoveBufferedVertices()
		{
			for (auto it = indexBufferEnd; it != indexBuffer.begin(); --it)
				nodes.erase(nodes.begin() + *(it-1));
		}
		
		template<int_t N>
		void AddRandomIndicesToBuffer()
		{
			for (uint_t i = 0; i < 2 * N; ++i)
				indexBuffer[i] = nodes.size();
			for (uint_t i = 0; i < N;)
			{
				uint_t r = configSpace.rng() * nodes.size() / 2;
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
		
		std::size_t Vertices()
		{
			return nodes.size() / 2;
		}
		
		template<int_t N>
		void ComputeWoodburyMatrices(GeMatrix& u, GeMatrix& v, GeMatrix& a)
		{
			for (uint_t i = 0; i < 2 * N; ++i)
			{
				for (uint_t j = 0; j < nodes.size(); ++j)
				{
					u(1+j, 1+i) = configSpace.LookUpG0(nodes[j].Site, nodeBuffer[i].Site, nodes[j].Tau - nodeBuffer[i].Tau + configSpace.infinTau);
					v(1+i, 1+j) = configSpace.LookUpG0(nodeBuffer[i].Site, nodes[j].Site, nodeBuffer[i].Tau - nodes[j].Tau - configSpace.infinTau);
				}
				for (uint_t j = 0; j < 2 * N; ++j)
				{
					if (i < j)
					{
						a(1+i, 1+j) = configSpace.LookUpG0(nodeBuffer[i].Site, nodeBuffer[j].Site, nodeBuffer[i].Tau - nodeBuffer[j].Tau + configSpace.infinTau);
						a(1+j, 1+i) = configSpace.LookUpG0(nodeBuffer[j].Site, nodeBuffer[i].Site, nodeBuffer[j].Tau - nodeBuffer[i].Tau - configSpace.infinTau);
					}
				}
				a(1+i, 1+i) = 0.0;
			}
		}
		
		void PermutationMatrix(GeMatrix& perm)
		{
			perm = 0.0;
			int cnt = perm.firstRow();
			for (int i = perm.firstRow(); i <= perm.lastRow(); ++i)
			{
				if (find(indexBuffer.begin(), indexBufferEnd, i-1) == indexBufferEnd)
				{
					perm(i, cnt) = 1.0;
					++cnt;
				}
			}
			int i = 0;
			for (auto it = indexBuffer.begin(); it != indexBufferEnd; ++it)
			{
				perm(*it+1, cnt + i) = 1.0;
				++i;
			}
		}
	private:
		ConfigSpace_t& configSpace;
		std::vector<node_t> nodes;
		std::vector<node_t> nodeBuffer;
		typename std::vector<node_t>::iterator nodeBufferEnd;
		std::vector<std::size_t> indexBuffer;
		typename std::vector<std::size_t>::iterator indexBufferEnd;
};