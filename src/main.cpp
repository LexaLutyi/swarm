#include <iostream>
#include <vector>
#include <random>
#include <memory>
#include <chrono>

#include <Eigen/Dense>

#include "agent.h"
#include "connection.h"
#include "model.h"

int main(int argc, char *argv[])
{
	auto start = std::chrono::high_resolution_clock::now();
	// initial state distribution
	std::default_random_engine gen(1);
	std::normal_distribution<> nd{10, 10};

	// vector of agents
	std::vector<std::unique_ptr<Agent>> agents;
	size_t agents_number = 100000;

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
		if (i != 0)
			pc.add_connection(0, i, 1.);
		if (i + 1 != agents_number)
			pc.add_connection(agents_number - 1, i, 1.);
		if (i + 10 < agents_number)
			pc.add_connection(i + 9, i, 1.);
	}
	auto start_int = std::chrono::high_resolution_clock::now();
	// integration
	for (size_t nt = 0; nt < n_steps; nt++) {
		auto t = t0 + nt * dt;

		// computing dx for each agent
		for (auto i = 0; i < agents_number; i++) {
			auto &xi = agents[i]->get_state();
			Eigen::VectorXd ui = Eigen::VectorXd::Zero(1);
			// ui = sum(w_ij * (xj - xi))
			for (auto it = pc.neighbours(i); it; ++it) {
				auto j = it.row();
				auto w_ij = it.value();
				auto &xj = agents[j]->get_state();
				ui += w_ij * (xj - xi);
			}

			agents[i]->update_agent(ui, v);
		}

		// moving agents towards dx
		for (auto i = 0; i < agents_number; i++) {
			agents[i]->move_agent(dt);
			// if (round(t) == t) {
			// 	if (i == 0)
			// 		std::cout << "t = " << t << "\n";
			// 	std::cout << agents[i]->get_state() << "\n";
			// }
		}
		// std::cout << t << "\n";
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration1 =
	  std::chrono::duration_cast<std::chrono::microseconds>(start_int - start);
	auto duration2 =
	  std::chrono::duration_cast<std::chrono::microseconds>(stop - start_int);
	std::cout << duration1.count() << "\n" << duration2.count() << "\n";
	std::cout << pc.adjacency_matrix.nonZeros() << "\n";
}
