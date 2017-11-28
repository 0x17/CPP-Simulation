//
// Created by André Schnabel on 27.11.17.
//

#include <gtest/gtest.h>
#include "../Helpers.h"

TEST(HelpersTest, testVecAverageSubRange) {
    std::vector<double> nums = {1,2,3,4,5,6};
    ASSERT_FLOAT_EQ(2.0, Helpers::vecAverageSubRange(nums, 0, 3));
}