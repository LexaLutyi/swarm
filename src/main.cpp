#include <iostream>
#include <vector>
#include <random>
#include <memory>
#include <chrono>
#include <algorithm>

#include <Eigen/Dense>

#include "agent.h"
#include "connection.h"
#include "model.h"

int main(int argc, char *argv[])
{
	auto start = std::chrono::high_resolution_clock::now();
	// initial state distribution
	std::default_random_engine gen(200);
	std::uniform_real_distribution<> nd{0., 20.};

	// vector of agents
	std::vector<std::unique_ptr<Agent>> agents;
	size_t agents_number = 1000;

	// initialization of model
	for (auto i = 0; i < agents_number; i++) {
		SimpleAgent sa(nd(gen), 0);
		agents.push_back(move(std::make_unique<SimpleAgent>(sa)));
	}

	// disturbance is zero
	Eigen::Matrix<double, 1, 1> v;
	v << 0;

	// time of integration and step definition
	double t0 = 0., t1 = 10.;
	size_t n_steps = 201;
	double dt = (t1 - t0) / (n_steps - 1);

	// Permanent neighbours
	PermanentConnection pc(agents_number);
	for (auto i = 0; i < agents_number; i++) {
		// if (i != 0)
		// 	pc.add_connection(0, i, 1.);
		// if (i + 1 != agents_number)
		// 	pc.add_connection(agents_number - 1, i, 1.);
		// if (i + 10 < agents_number)
		// 	pc.add_connection(i + 9, i, 1.);
	}
	MetricConnectionXD mc(2, 0, 20, 10, true);
	for (auto i = 0; i < agents_number; i++)
		if (agents[i]->has_position)
			mc.update(i, {agents[i]->position()[0], 19});

	//
	// std::cout << "********"
	//           << "\n";
	// for (auto &s : mc.adjacency_sets)
	// 	std::cout << s.size() << "\n";
	// std::cout << "********"
	//           << "\n";
	// for (auto &a : agents)
	// 	std::cout << a->state()[0] << "\n";
	// std::cout << "********"
	//           << "\n";
	//

	auto start_int = std::chrono::high_resolution_clock::now();
	// integration
	for (size_t nt = 0; nt < n_steps; nt++) {
		auto t = t0 + nt * dt;

		// computing dx for each agent
		// in PermanentConnection
		for (auto i = 0; i < agents_number; i++) {
			auto &x_i = agents[i]->state();
			Eigen::VectorXd u_i = Eigen::VectorXd::Zero(1);
			// computing control
			// ui = sum(w_ij * (xj - xi))
			// for (auto it = pc.neighbours(i); it; ++it) {
			// 	auto j = it.row();
			// 	auto w_ij = it.value();
			// 	auto &xj = agents[j]->state();
			// 	ui += w_ij * (xj - xi);
			// }

			for (auto it_j = MetricIteratorXD(mc, i, 1.); it_j; ++it_j) {
				auto j = it_j.j();
				auto &x_j = agents[j]->state();
				u_i += x_j - x_i;
			}
			if (u_i.norm() > 1)
				u_i.normalize();
			agents[i]->update_agent(u_i, v);
		}

		// moving agents towards dx
		for (auto i = 0; i < agents_number; i++) {
			agents[i]->move_agent(dt);
			if (agents[i]->has_position)
				mc.update(i, {agents[i]->position()[0], 19});
		}
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration1 =
	  std::chrono::duration_cast<std::chrono::microseconds>(start_int - start);
	auto duration2 =
	  std::chrono::duration_cast<std::chrono::microseconds>(stop - start_int);
	std::cout << duration1.count() << "\n" << duration2.count() << "\n";
	std::cout << pc.adjacency_matrix.nonZeros() << "\n";

	// for (auto &a : agents)
	// 	std::cout << a->state()[0] << "\n";
	std::cout << "***********"
	          << "\n";
	for (auto &s : mc.adjacency_sets)
		for (auto i : s.second) {
			std::vector<size_t> ccc(2);
			mc.get_coords(s.first, ccc);
			std::cout << s.first << "\t:\t" << i << "\t" << agents[i]->state()
			          << "\ti:\t" << ccc[0] << "\tj:\t" << ccc[1] << "\n";
			;
		}
	// std::cout << "***********"
	//           << "\n";
}
