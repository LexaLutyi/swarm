#pragma once

#include <vector>
#include <set>

#include <Eigen/Sparse>

// #include "agent.h"

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
	Eigen::SparseMatrix<double>::InnerIterator neighbours(size_t i) const
	{
		Eigen::SparseMatrix<double>::InnerIterator it(adjacency_matrix, i);
		return it;
	}
};

class MetricConnection {
private:
	double x_start;
	double x_end;
	// n_start = 0
	size_t n_end;

	double x_delta;

public:
	// v_i = {j_1, ..., j_m}
	// if position(j_1) < x_i
	// x0, x1, ..., x_n -> n + 1 points -> n + 2 sets

	std::vector<std::set<size_t>> adjacency_sets;

	// id -> v_i
	std::map<size_t, size_t> id_map;

	// must be declared outside
	// double position(Agent &a) { return a.get_state()[1]; }

	MetricConnection(double x0, double x1, size_t np)
	  : x_start(x0), x_end(x1), n_end(np)
	{
		x_delta = (x_end - x_start) / n_end;
		adjacency_sets.resize(n_end + 2);
	}

	void update(size_t id, double pos)
	{
		auto i = (int)floor((pos - x_start) / x_delta);
		i = i < 0 ? 0 : i;
		i = i > n_end + 1 ? n_end + 1 : i;
		auto it = id_map.find(id);
		if (it == id_map.end()) {
			id_map[id] = i;
			adjacency_sets[i].insert(id);
		}
		else if (it->second != i) {
			adjacency_sets[it->second].erase(it->first);
			it->second = i;
			adjacency_sets[i].insert(it->first);
		}
	}

	void del(size_t id)
	{
		auto it = id_map.find(id);
		if (it != id_map.end()) {
			adjacency_sets[it->second].erase(it->first);
			id_map.erase(it);
		}
	}

	double delta() const { return x_delta; }
};

// for fixed i and radius
// return every possible neighbouring agent
class MetricIterator {
	// agent_i identifire
	const size_t id_i;
	// radius of interest
	const double r;
	// position data
	const MetricConnection &mc;
	// radius in cells
	size_t r_cells;
	// agent_i's cell
	size_t cell_i;
	// cell to start
	size_t cell_start;
	// cell to stop, works like end(): real end + 1
	size_t cell_stop;

	// current cell
	size_t cell_j;
	// current set iterator
	std::set<size_t>::const_iterator it;
	// are there any neighbours
	bool is_not_empty;

public:
	MetricIterator(MetricConnection &metric_connection, size_t id, double radius)
	  : mc(metric_connection), id_i(id), r(radius)
	{
		auto it_i = mc.id_map.find(id_i);
		if (it_i == mc.id_map.cend()) {
			is_not_empty = false;
			return;
		}
		cell_i = it_i->second;
		double delta = mc.delta();
		r_cells = (size_t)(radius / delta) + 1;
		ptrdiff_t q = cell_i - r_cells;
		cell_start = q < 0 ? 0 : q;
		q = cell_i + r_cells;
		cell_stop =
		  q + 1 < mc.adjacency_sets.size() ? q + 1 : mc.adjacency_sets.size();
		for (auto j = cell_start; j < cell_stop; j++)
			if (!mc.adjacency_sets[j].empty()) {
				is_not_empty = true;
				it = mc.adjacency_sets[j].cbegin();
				cell_j = j;
				return;
			}
		is_not_empty = false;
	}

	explicit operator bool() const { return is_not_empty; }

	MetricIterator &operator++()
	{
		it++;
		if (it == mc.adjacency_sets[cell_j].end()) {
			for (auto j = cell_j + 1; j < cell_stop; j++) {
				if (!mc.adjacency_sets[j].empty()) {
					it = mc.adjacency_sets[j].cbegin();
					cell_j = j;
					return *this;
				}
			}
			is_not_empty = false;
		}
		return *this;
	}

	size_t j() const { return *it; }
};