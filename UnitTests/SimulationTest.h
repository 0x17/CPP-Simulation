//
// Created by André Schnabel on 27.11.17.
//

#pragma once

#include <gtest/gtest.h>

class MultiClassSimulation;

class MultiClassSimulationTest : public testing::Test {
protected:
    std::unique_ptr<MultiClassSimulation> mcs;

    void SetUp() override;
};

