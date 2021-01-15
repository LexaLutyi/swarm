#pragma once

#include <map>
#include "agent.h"

class Model {
private:
	std::map<size_t id, std::unique_ptr<Agent>> agents;

public:
};