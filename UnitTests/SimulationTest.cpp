//
// Created by André Schnabel on 27.11.17.
//

#include "SimulationTest.h"
#include "../Simulation.h"
#include "../Helpers.h"

#include <boost/filesystem.hpp>
#include <cmath>

const static std::string    TMP_DATA_FILENAME = "temp_data.json",
                            TMP_DATA_CONTENTS = "{\n"
        "  \"capacity\": 20,\n"
        "  \"clients\": [\n"
        "    {\n"
        "      \"name\": \"A\",\n"
        "      \"expD\": 2.0,\n"
        "      \"devD\": 1.0,\n"
        "      \"description\": \"kurzfristig, 2 Wochen vorher, trifft spät ein\",\n"
        "      \"consumptionPerReq\": 3.0,\n"
        "      \"revenuePerReq\": 6.0\n"
        "    },\n"
        "    {\n"
        "      \"name\": \"B\",\n"
        "      \"expD\": 3.0,\n"
        "      \"devD\": 1.0,\n"
        "      \"description\": \"langfristig, trifft früh ein\",\n"
        "      \"consumptionPerReq\": 2.0,\n"
        "      \"revenuePerReq\": 3.0\n"
        "    },\n"
        "    {\n"
        "      \"name\": \"C\",\n"
        "      \"expD\": 19.0,\n"
        "      \"devD\": 1.0,\n"
        "      \"description\": \"sehr langfristig, tritt als erstes ein\",\n"
        "      \"consumptionPerReq\": 1.0,\n"
        "      \"revenuePerReq\": 1.0\n"
        "    }\n"
        "  ]\n"
        "}";

void MultiClassSimulationTest::SetUp() {
    Helpers::spit(TMP_DATA_CONTENTS, TMP_DATA_FILENAME);
    Toggles toggles;
    mcs = std::make_unique<MultiClassSimulation>(TMP_DATA_FILENAME, toggles);
    boost::filesystem::remove(TMP_DATA_FILENAME);
}

TEST_F(MultiClassSimulationTest, testPickDemands) {
    double eps = 2.0;
    DemandScenario scenario = mcs->pickDemands();
    for(int cix = 0; cix < mcs->getNumClasses(); cix++) {
        Customer c = mcs->getCustomer(cix);
        ASSERT_TRUE(std::fabs((double)c.expD - scenario[cix]) < eps);
    }
}

TEST_F(MultiClassSimulationTest, testPickDemandsDescriptive) {
    double eps = 2.0;

    std::vector<DistParameters> distParams = mcs->demandDistributionParametersForCustomers();
    auto lutList = generateLookupTableList(mcs->getNumClasses(), 1000, distParams);
    DemandScenario scenario = mcs->pickDemandsDescriptive(lutList);
    for(int cix = 0; cix < mcs->getNumClasses(); cix++) {
        Customer c = mcs->getCustomer(cix);
        ASSERT_TRUE(std::fabs((double)c.expD - scenario[cix]) <= eps);
    }
}