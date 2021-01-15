#pragma once

#include <vector>
#include <numeric>
#include <cmath>

class Environment {
private:
	std::vector<double> mass;
	double x_start;
	double x_end;

public:
	Environment(double x0, double x1, int n, double m);
	void action(double dt);
	void consume(double place, double amt);
	double get(double place) const;
	double get_all() const;
};

Environment::Environment(double x0, double x1, int n, double m)
{
	mass.resize(n, 0);
	x_start = x0;
	x_end = x1;
	for (double &x : mass)
		x = m / n;
	// mass[4] = 1;
}

void Environment::action(double dt)
{
	size_t n = mass.size();
	double tmp = mass[0];
	mass[0] += (2 * mass[0] + mass[1]) * (1 - mass[0]) * dt;
	for (auto i = 1; i < n - 1; i++) {
		double tmp2 = mass[i];
		mass[i] += (tmp + 2 * mass[i] + mass[i + 1]) * (1 - mass[i]) * dt;
		tmp = tmp2;
	}
	mass[n - 1] += (tmp + 2 * mass[n - 1]) * (1 - mass[n - 1]) * dt;
}

double Environment::get_all() const
{
	return std::accumulate(mass.begin(), mass.end(), 0.) * (x_end - x_start) /
	       (mass.size() - 1);
}

double Environment::get(double place) const
{
	size_t n = mass.size();
	double x_step = (x_end - x_start) / (n - 1);
	if (place < x_start || place > x_end)
		return -1;
	if (place == x_end)
		return mass[n - 1];
	double i;
	double alpha = modf((place - x_start) / x_step, &i);

	return (1 - alpha) * mass[i] + alpha * mass[i + 1];
}

void Environment::consume(double place, double percent)
{
	double left = 1 - percent;
	size_t n = mass.size();
	double x_step = (x_end - x_start) / (n - 1);
	if (place < x_start || place > x_end)
		return;
	if (place == x_end) {
		mass[n - 1] *= left;
		return;
	}
	double i;
	double alpha = modf((place - x_start) / x_step, &i);

	// double q = 0.8 * left / (1 - alpha);

	mass[i] *= left;
	mass[i + 1] *= left;

	// (1 - a) * x *a *p + a *y *(1 - a) *q = 0.5 * (1 - a) * x + 0.5 * a * y;
	// p = (0.5 * (1 - a) * x + 0.5 * a * y - a *y *(1 - a) *q) / ((1 - a) * x *a)
}
