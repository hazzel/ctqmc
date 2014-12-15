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
//#define EIGEN_USE_MKL_ALL
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"

template<typename Value_t, int N, int M>
class MemoryPool
{
	using matrix_t = Eigen::Matrix<value_t, N, M>;
	public:
		MemoryPool() {}
	private:
		std::vector<matrix_t>3
};