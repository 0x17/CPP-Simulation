#pragma once

#include <vector>

class AbstractSimulation;

void solveWithGurobi(AbstractSimulation &sim, std::vector<std::vector<int>> &scenarios);

