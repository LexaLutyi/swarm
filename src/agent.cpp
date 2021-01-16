#include <iostream>
// #include <vector>

#include "agent.h"

std::ostream &operator<<(std::ostream &out, const Agent &a)
{
	// Поскольку operator<< является другом класса Point, то мы имеем прямой
	// доступ к членам Point
	out << "t: " << a.t << "\n";
	out << "x:\n" << a.x.transpose() << "\n";
	out << "dx:\n" << a.dx.transpose() << "\n";

	return out;
}
