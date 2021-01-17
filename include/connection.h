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
	size_t idx_i;
	// cell to start
	size_t coords_start;
	// cell to stop, works like end(): real end + 1
	size_t coords_stop;

	// current cell
	size_t coords_cur;
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
		idx_i = it_i->second;
		double delta = mc.delta();
		r_cells = (size_t)(radius / delta) + 1;
		ptrdiff_t q = idx_i - r_cells;
		coords_start = q < 0 ? 0 : q;
		q = idx_i + r_cells;
		coords_stop =
		  q + 1 < mc.adjacency_sets.size() ? q + 1 : mc.adjacency_sets.size();
		for (auto j = coords_start; j < coords_stop; j++)
			if (!mc.adjacency_sets[j].empty()) {
				is_not_empty = true;
				it = mc.adjacency_sets[j].cbegin();
				coords_cur = j;
				return;
			}
		is_not_empty = false;
	}

	explicit operator bool() const { return is_not_empty; }

	MetricIterator &operator++()
	{
		it++;
		if (it == mc.adjacency_sets[coords_cur].end()) {
			for (auto j = coords_cur + 1; j < coords_stop; j++) {
				if (!mc.adjacency_sets[j].empty()) {
					it = mc.adjacency_sets[j].cbegin();
					coords_cur = j;
					return *this;
				}
			}
			is_not_empty = false;
		}
		return *this;
	}

	size_t j() const { return *it; }
};

class MetricConnectionXD {
private:
	// dimension number
	size_t n_dim;
	// start points
	std::vector<double> x_start;
	// end points
	std::vector<double> x_end;
	// number of intervals (size)
	// depends on is_cycle
	std::vector<size_t> x_size;
	// size of intervals
	std::vector<double> x_delta;

	// boundaries type
	std::vector<bool> is_cycle;

public:
	// xyz_i = {j_1, ..., j_m}
	// if position(j_1) < x_i
	// x0, x1, ..., x_n -> n + 1 points -> n + 2 sets

	std::map<size_t, std::set<size_t>> adjacency_sets;
	// std::multimap<size_t, size_t> adjacency_sets;

	// id -> xyz_i
	std::map<size_t, size_t> id_map;

	// must be declared outside
	// double position(Agent &a) { return a.get_state()[1]; }

	MetricConnectionXD(size_t dim, double x_0, double x_n, size_t n,
	                   bool cycle = false)
	  : n_dim(dim)
	{
		auto x_del = (x_n - x_0) / n;
		for (auto dim = 0; dim < n_dim; dim++) {
			x_start.push_back(x_0);
			x_end.push_back(x_n);
			x_delta.push_back(x_del);
			is_cycle.push_back(cycle);
			x_size.push_back(cycle ? n : n + 2);
		}
	}

	// index from good coords
	size_t get_index(const std::vector<size_t> &coords) const
	{
		size_t idx = 0;
		size_t power = 1;
		for (auto dim = n_dim; dim > 0;) {
			dim--;
			idx += coords[dim] * power;
			power *= x_size[dim];
		}
		return idx;
	}

	// coords from index
	void get_coords(size_t idx, std::vector<size_t> &coords) const
	{
		for (size_t dim = n_dim; dim > 0;) {
			dim--;
			size_t n = x_size[dim];
			size_t i = idx % n;

			coords[dim] = i;
			idx = (idx - i) / n;
		}
	}

	// return out of range coords for cycle connection
	// good coords otherwise
	void sum_coords(const std::vector<size_t> &left,
	                const std::vector<ptrdiff_t> &right,
	                std::vector<ptrdiff_t> &result) const
	{
		for (size_t dim = 0; dim < n_dim; dim++) {
			ptrdiff_t res = left[dim] + right[dim];
			result[dim] = is_cycle[dim] ? res : no_cycle(res, x_size[dim]);
		}
	}

	// index from bad coords
	size_t get_index(const std::vector<ptrdiff_t> &coords) const
	{
		size_t idx = 0;
		size_t power = 1;
		for (auto dim = n_dim; dim > 0;) {
			dim--;
			size_t good_coord = is_cycle[dim] ? cycle(coords[dim], x_size[dim]) :
			                                    no_cycle(coords[dim], x_size[dim]);
			idx += good_coord * power;
			power *= x_size[dim];
		}
		return idx;
	}

	// i = 0 ... n-1
	size_t cycle(ptrdiff_t i, size_t n) const
	{
		if (i < 0) {
			return n - (-i) % n;
		}
		else {
			return i % n;
		}
	}
	// i = 0 ... n+1
	size_t no_cycle(ptrdiff_t i, size_t n) const
	{
		i = i < 0 ? 0 : i;
		i = i > n - 1 ? n - 1 : i;
		return i;
	}

	void update(size_t id, const std::vector<double> &pos)
	{
		std::vector<ptrdiff_t> coords(n_dim);
		for (auto dim = 0; dim < n_dim; dim++) {
			double q = (pos[dim] - x_start[dim]) / x_delta[dim];
			if (q >= 0) {
				coords[dim] = (ptrdiff_t)floor(q);
			}
			else {
				coords[dim] = (ptrdiff_t)floor(q) - 1;
			}
		}
		auto i = get_index(coords);
		auto it = id_map.find(id);
		if (it == id_map.end()) {
			id_map[id] = i;
			adjacency_sets[i].insert(id);
		}
		else if (it->second != i) {
			adjacency_sets[it->second].erase(id);
			it->second = i;
			adjacency_sets[i].insert(id);
		}
	}

	void del(size_t id)
	{
		auto it = id_map.find(id);
		if (it != id_map.end()) {
			adjacency_sets[it->second].erase(id);
			id_map.erase(it);
		}
	}

	const std::vector<double> &delta() const { return x_delta; }
};

// for fixed i and radius
// return every possible neighbouring agent
class MetricIteratorXD {
	// agent_i identifire
	const size_t id_i;
	// radius of interest
	const double r;
	// position data
	const MetricConnectionXD &mc;
	// agent_i's cell index
	size_t idx_i;
	// dimension number
	size_t n_dim;
	// cell to start
	std::vector<ptrdiff_t> coords_start;
	// cell to stop, works like end(): real end + 1
	std::vector<ptrdiff_t> coords_stop;

	// current cell
	std::vector<ptrdiff_t> coords_cur;
	// current cell iterator
	std::map<size_t, std::set<size_t>>::const_iterator cell_it;
	// current set iterator
	std::set<size_t>::const_iterator it;
	// are there any neighbours
	bool is_not_empty;

	void next_cell()
	{
		// for initialization call
		if (!is_not_empty) {
			coords_cur = coords_start;
			size_t j = mc.get_index(coords_cur);
			cell_it = mc.adjacency_sets.find(j);
			if (cell_it != mc.adjacency_sets.cend()) {
				// found new cell
				is_not_empty = true;
				return;
			}
		}

		// for iteration calls

		for (auto dim = n_dim; dim > 0;) {
			dim--;

			coords_cur[dim]++;
			for (; coords_cur[dim] != coords_stop[dim]; coords_cur[dim]++) {
				if (dim + 1 < n_dim) {
					// we always iterate on last dimension
					dim++;
				}
				size_t j = mc.get_index(coords_cur);
				cell_it = mc.adjacency_sets.find(j);
				if (cell_it != mc.adjacency_sets.cend()) {
					// found new cell
					is_not_empty = true;
					return;
				}
			}
			coords_cur[dim] = coords_start[dim];
		}
		// no cells left
		is_not_empty = false;
		return;
	}

public:
	MetricIteratorXD(const MetricConnectionXD &metric_connection, size_t id,
	                 double radius)
	  : mc(metric_connection), id_i(id), r(radius)
	{
		auto it_i = mc.id_map.find(id_i);
		if (it_i == mc.id_map.cend()) {
			is_not_empty = false;
			return;
		}
		idx_i = it_i->second;

		std::vector<double> delta = mc.delta();
		n_dim = delta.size();

		std::vector<ptrdiff_t> from(n_dim);
		std::vector<ptrdiff_t> to(n_dim);

		for (size_t dim = 0; dim < n_dim; dim++) {
			ptrdiff_t r = (ptrdiff_t)((radius / delta[dim]) + 1);
			from[dim] = -r;
			to[dim] = r + 1;
		};

		std::vector<size_t> cur(n_dim);
		mc.get_coords(idx_i, cur);
		coords_start.resize(n_dim);
		mc.sum_coords(cur, from, coords_start);
		coords_stop.resize(n_dim);
		mc.sum_coords(cur, to, coords_stop);

		is_not_empty = false;
		next_cell();
		for (; is_not_empty; next_cell()) {
			it = cell_it->second.cbegin();
			if (it != cell_it->second.cend())
				return;
		}
	}

	explicit operator bool() const { return is_not_empty; }

	MetricIteratorXD &operator++()
	{
		it++;
		if (it != cell_it->second.cend())
			return *this;

		next_cell();
		for (; is_not_empty; next_cell()) {
			it = cell_it->second.cbegin();
			if (it != cell_it->second.cend()) {
				return *this;
			}
		}

		return *this;
	}

	size_t j() const { return *it; }
};