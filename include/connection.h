#pragma once

#include <Eigen/Sparse>

// Connection define links between agents
// Не стоит делать кучу всего ради экономии двух строчек кода
// особенно, когда эти связи нужно писать руками

class PermanentConnection {
public:
	Eigen::SparseMatrix<double> adjacency_matrix;

public:
	PermanentConnection(size_t n): adjacency_matrix(n, n){};
	void add_connection(size_t i, size_t j, double w)
	{
		adjacency_matrix.coeffRef(i, j) = w;
	}
	void remove_connection(size_t i, size_t j)
	{
		adjacency_matrix.coeffRef(i, j) = 0.;
		adjacency_matrix.prune(0.);
	}
	Eigen::SparseMatrix<double>::InnerIterator neighbours(size_t i)
	{
		Eigen::SparseMatrix<double>::InnerIterator it(adjacency_matrix, i);
		return it;
	}
};
