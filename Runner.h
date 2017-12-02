//
// Created by André Schnabel on 23.05.17.
//

#pragma once

#include <string>

class Runner {
public:
	static void commandLine(const std::list<std::string> &args);
	static void benchmark(const std::string &dir);
	static void runOptimizers();
};

