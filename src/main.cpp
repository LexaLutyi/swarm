#include <iostream>
#include <vector>
#include <random>
#include <memory>

#include <Eigen/Dense>

#include "agent.h"

int main(int argc, char *argv[])
{
	std::default_random_engine gen(1);
	std::normal_distribution<> nd{10, 1};

	// vector of agents
	std::vector<std::unique_ptr<Agent>> agents;

	for (int i = 0; i < 10; i++) {
		SimpleAgent sa(nd(gen), 0);
		agents.push_back(move(std::make_unique<SimpleAgent>(sa)));
	}
	for (auto &a : agents) {
		std::cout << *a << "\n";
	}
	std::cout << "------------------------------"
	          << "\n";
	Eigen::Matrix<double, 1, 1> v;
	v << 0;
	double dt = 0.05;

	for (double t = 0; t < 10; t += dt) {
		std::cout << "t = " << t << "\n";
		for (auto &a : agents) {
			Eigen::VectorXd xi = a->get_state();
			Eigen::VectorXd ui(1);
			ui << 0;

			for (auto &aj : agents) {
				if (aj == a)
					continue;
				ui += aj->get_state();
				ui -= 2 * xi;
			}
			if (a == agents.back())
				ui << sin(t);
			a->update_agent(ui, v);
			std::cout << a->get_state() << "\n";
		}
	}
}
