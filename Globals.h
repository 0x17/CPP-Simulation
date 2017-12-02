//
// Created by André Schnabel on 28.05.17.
//

#pragma once

namespace globals {
	const double EPSILON = 1.0; // 0.00000001;
	const double EPSILON2 = 0.00001;//std::numeric_limits<double>::min();
	
	const double TIME_LIMIT = 2.0;
	const bool TRACING_ENABLED = false;
	
	const bool ECONOMY_OF_SCALE_ENABLED = false;
	const bool CONDITIONAL_VALUE_AT_RISK_ENABLED = false;
	const bool STOCHASTIC_CONSUMPTIONS_ENABLED = false;

	const double PROFIT_WEIGHT = 0.5;
	const double ALPHA = 0.7;
}
