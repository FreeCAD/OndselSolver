#include "PosICDragLimitNewtonRaphson.h"
#include "SystemSolver.h"
#include "SimulationStoppingError.h"
#include "Part.h"
#include "Constraint.h"
#include <algorithm>

using namespace MbD;

std::shared_ptr<PosICDragLimitNewtonRaphson> MbD::PosICDragLimitNewtonRaphson::With()
{
	auto newtonRaphson = std::make_shared<PosICDragLimitNewtonRaphson>();
	newtonRaphson->initialize();
	return newtonRaphson;
}

void MbD::PosICDragLimitNewtonRaphson::preRun()
{
	std::string str("MbD: Assembling system with limits. ");
	system->logString(str);
	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->prePosIC(); });

}

void MbD::PosICDragLimitNewtonRaphson::run()
{
	preRun();
	system->deactivateLimits();
	if (system->limitsSatisfied()) {
		std::string str("MbD: No limits reached. ");
		system->logString(str);
		return;
	}
	auto limits = system->limits();
	for (size_t i = 0; i < limits->size(); i++) {
		bool activatedLimit = false;
		for (auto& limit : *limits) {
			// Keep earlier violated limits active while adding newly violated ones.
			if (!limit->active && !limit->satisfied()) {
				limit->activate();
				activatedLimit = true;
			}
		}
		if (!activatedLimit) {
			break;
		}

		preRun();
		initializeLocally();
		initializeGlobally();
		iterate();
		postRun();
		if (system->limitsSatisfied()) {
			system->deactivateLimits();
			return;
		}
	}
	system->deactivateLimits();
	throw SimulationStoppingError("Limits cannot be satisfied.");
}
