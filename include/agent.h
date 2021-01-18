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

// position = [p_1, p_2] \in rectangle:(0, d_1)x(0, d_2)
// cyclic boundaries
// dx = [v sin(phi), v cos(phi)]
class Agent2D : public Agent {
protected:
	void update_pos()
	{
		// validate periodic variables x1, x2, phi
		x[0] = periodic(x[0], p[0]);
		x[1] = periodic(x[1], p[1]);
		x[3] = periodic(x[3], 8 * atan(1));
		// update position
		pos << x[0], x[1];
	}
	void update_dx()
	{
		// dx_1 = v sin(phi)
		dx[0] = x[2] * sin(x[3]);
		// dx_2 = v cos(phi)
		dx[1] = x[2] * cos(x[3]);
		// dv = u_1
		dx[2] = u[1];
		// dphi = u_2
		dx[3] = u[2];
	}

	// return r >= 0:
	// x = period * n + r
	double periodic(double x, double period) const
	{
		if (x >= 0)
			return fmod(x, period);
		else
			return period - fmod(-x, period);
	}
	// return r2: -T/2 <= r2 < T/2:
	// x + T/2 = T * n + r
	// r2 = r - T/2
	double periodic2(double x, double period) const
	{
		double semiperiod = 0.5 * period;
		double x2 = x + semiperiod;
		if (x2 >= 0)
			return fmod(x2, period) - semiperiod;
		else
			return semiperiod - fmod(-x2, period);
	}

public:
	Agent2D(const Eigen::Vector4d &state, const Eigen::Vector2d &params,
	        double time = 0)
	  : Agent(time, true)
	{
		x = state;
		pos.resize(2);
		dx.resize(4);
		u.resize(2);
		// rectangle parameters;
		p = params;
		v.resize(0);

		update_pos();
	}

	// compute periodic distanse based on agent's parameters
	Eigen::Vector4d diff_state(const Eigen::Vector4d &x) const
	{
		Eigen::Vector4d res = x;
		res[0] = periodic(res[0], p[0]);
		res[1] = periodic(res[1], p[1]);
		res[3] = periodic(res[3], 8 * atan(4));
		return res;
	}
};