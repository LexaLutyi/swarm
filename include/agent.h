#pragma once

#include <Eigen/Dense>
#include <iostream>

class Agent {
protected:
	size_t id;
	Eigen::VectorXd x;
	Eigen::VectorXd p;
	virtual void update_dx() = 0;

	Eigen::VectorXd u;
	Eigen::VectorXd v;
	Eigen::VectorXd dx;
	double t;

public:
	Agent(size_t identifier, double time = 0): id(identifier), t(time) {}

	// update dx before move
	void update_agent(const Eigen::VectorXd &control,
	                  const Eigen::VectorXd &disturbance)
	{
		u = control;
		v = disturbance;
		update_dx();
	}
	// move the agent after updating all agents
	// standart Euler scheme
	void move_agent(double dt)
	{
		x += dx * dt;
		t += dt;
	}
	const Eigen::VectorXd &get_state() { return x; }
	friend std::ostream &operator<<(std::ostream &out, const Agent &a);
};

// dx/dt = u
class SimpleAgent : public Agent {
private:
	void update_dx() { dx = u; }

public:
	SimpleAgent(double state, double time): Agent(time)
	{
		x.resize(1);
		dx.resize(1);
		u.resize(1);
		p.resize(0);
		v.resize(0);

		x(0, 0) = state;
		dx(0, 0) = 0;
	}
};
