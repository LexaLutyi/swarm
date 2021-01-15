#pragma once

#include <Eigen/Sparse>

// Connection define links between agents

class Connection {
protected:
	/* data */
public:
	Connection(/* args */);
};

class ConstantConnection {
	// adjacency(i, j) <> 0 <=> exists link i -> j
	// adjacency(i, j) = 0 <=> there is no link from i to j

private:
	Eigen::SparseMatrix<double> adjacency;

public:
	ConstantConnection() {}
	void add_link(size_t i, size_t j, double weight = 1)
	{
		// if (weight == 0.) raise exception
		adjacency.coeffRef(i, j) = weight;
	}
	void remove_link(size_t i, size_t j)
	{
		adjacency.coeffRef(i, j) = 0.;
		adjacency.prune(0.0);
	}
};
