#pragma once

#include <gtest/gtest.h>
#include "TestHelpers.h"
#include "../Matrix3D.h"
#include <memory>

class Matrix3DTest : public testing::Test {
protected:
	std::unique_ptr<Matrix3D<int>> m;

	void SetUp() override {
		m = std::unique_ptr<Matrix3D<int>>(new Matrix3D<int>({{
			{1,2,3,4},
			{5,6,7,8},
			{9,10,11,12}
		}}));
	}
};
