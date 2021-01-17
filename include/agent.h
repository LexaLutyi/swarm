#pragma once

#include <Eigen/Dense>
#include <iostream>

class Agent {
protected:
	Eigen::VectorXd x;
	Eigen::VectorXd pos;
	Eigen::VectorXd p;
	virtual void update_dx() = 0;
	virtual void update_pos() = 0;

	Eigen::VectorXd u;
	Eigen::VectorXd v;
	Eigen::VectorXd dx;
	double t;

public:
	Agent(double time = 0, bool positionable = false)
	  : t(time), has_position(positionable)
	{
	}

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
		update_pos();
		t += dt;
	}
	const Eigen::VectorXd &state() { return x; }

	const bool has_position;

	const Eigen::VectorXd &position() { return pos; }

	friend std::ostream &operator<<(std::ostream &out, const Agent &a);
};

// dx/dt = u
class SimpleAgent : public Agent {
protected:
	void update_dx() { dx = u; }
	void update_pos() { pos = x; }

public:
	SimpleAgent(double state, double time): Agent(time, true)
	{
		x.resize(1);
		pos.resize(1);
		dx.resize(1);
		u.resize(1);

		p.resize(0);
		v.resize(0);

		x << state;
		pos << state;
		dx(0, 0) = 0;
	}
};

class SimpleAgent2 : public Agent {
protected:
	void update_pos() { pos << -x[0], -x[0]; }
	void update_dx() { dx = u; }

public:
	SimpleAgent2(double state, double time): Agent(time, true)
	{
		x.resize(1);
		pos.resize(2);
		dx.resize(1);
		u.resize(1);
		p.resize(0);
		v.resize(0);

		x << state;
		pos << -state, -state;
		dx << 0;
	}
};