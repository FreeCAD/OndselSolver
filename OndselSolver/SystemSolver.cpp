/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include <vector>
#include <set>
#include <algorithm>

#include "SystemSolver.h"
#include "DynIntegrator.h"
#include "LimitIJ.h"
#include "System.h"
#include "NewtonRaphson.h"
#include "PosICNewtonRaphson.h"
#include "CREATE.h"
#include "RedundantConstraint.h"
#include "NotKinematicError.h"
#include "ICKineIntegrator.h"
#include "KineIntegrator.h"
#include "DiscontinuityError.h"
#include "PosICKineNewtonRaphson.h"
#include "PosKineNewtonRaphson.h"
#include "VelICSolver.h"
#include "AccICNewtonRaphson.h"
#include "VelKineSolver.h"
#include "AccKineNewtonRaphson.h"
#include "VelICKineSolver.h"
#include "AccICKineNewtonRaphson.h"
#include "PosICDragNewtonRaphson.h"
#include "PosICDragLimitNewtonRaphson.h"

using namespace MbD;

//class PosICNewtonRaphson;

void SystemSolver::setSystem(Solver*)
{
	//Do not use
	throw SimulationStoppingError("To be implemented.");
}

void SystemSolver::initialize()
{
	tstartPasts = std::make_shared<std::vector<double>>();
}

void SystemSolver::initializeLocally()
{
	setsOfRedundantConstraints = std::make_shared<std::vector<std::shared_ptr<std::set<std::string>>>>();
	direction = (tstart < tend) ? 1.0 : -1.0;
	toutFirst = tstart + (direction * hout);
}

void SystemSolver::initializeGlobally()
{
}

void SystemSolver::runAllIC()
{
	while (true)
	{
		initializeLocally();
		initializeGlobally();
		runPosIC();
		while (needToRedoPosIC())
		{
			runPosIC();
		}
		runVelIC();
		runAccIC();
		auto discontinuities = system->discontinuitiesAtIC();
		if (discontinuities->size() == 0) break;
		if (std::find(discontinuities->begin(), discontinuities->end(), "REBOUND") != discontinuities->end())
		{
			preCollision();
			runCollisionDerivativeIC();
			runBasicCollision();
		}
	}
}

void SystemSolver::runPosIC()
{
	icTypeSolver = CREATE<PosICNewtonRaphson>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::runVelIC()
{
	icTypeSolver = CREATE<VelICSolver>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::runAccIC()
{
	icTypeSolver = CREATE<AccICNewtonRaphson>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::releaseSeparatingLimits()
{
    // Removing a tensile reaction redistributes the others. Re-solve after
    // each batch; at least one active limit is removed per iteration.
    while (true) {
        bool released = false;
        for (const auto& limit : *system->limits) {
            if (limit->hasTensileReaction()) {
                limit->deactivate();
                released = true;
            }
        }
        if (!released) return;
        runAccIC();
    }
}

bool SystemSolver::needToRedoPosIC()
{
	auto allRedunCons = this->allRedundantConstraints();
	auto newSet = std::make_shared<std::set<std::string>>();
	for (auto& con : *allRedunCons) {
		auto aaa = std::static_pointer_cast<RedundantConstraint>(con);
		auto& bbb = aaa->constraint->name;
		newSet->insert(bbb);
	}
	//std::transform(allRedunCons->begin(), allRedunCons->end(), newSet->begin(), [](auto con) {
	//	return std::static_pointer_cast<RedundantConstraint>(con)->constraint->name;
	//	});
	if (newSet->empty()) return false;
	auto itr = std::find_if(setsOfRedundantConstraints->begin(), setsOfRedundantConstraints->end(), [&](auto& set) {
		for (auto& name : *set) {
			if (newSet->find(name) == newSet->end()) return false;
		}
		return true;
		});
	if (itr != setsOfRedundantConstraints->end()) {
		//"Same set of redundant constraints found."
		setsOfRedundantConstraints->push_back(newSet);
		return false;
	}
	if (setsOfRedundantConstraints->size() >= 2) {
		auto it = std::find_if(setsOfRedundantConstraints->begin(), setsOfRedundantConstraints->end(), [&](auto set) {
			return set->size() == newSet->size();
			});
		if (it != setsOfRedundantConstraints->end()) {
			//"Equal number of redundant constraints found."
			setsOfRedundantConstraints->push_back(newSet);
			return false;
		}
	}
	setsOfRedundantConstraints->push_back(newSet);
	this->partsJointsMotionsDo([](auto item) { item->reactivateRedundantConstraints(); });
	return true;
}

void SystemSolver::preCollision()
{
}

void SystemSolver::runCollisionDerivativeIC()
{
}

void SystemSolver::runBasicCollision()
{
}

void SystemSolver::runBasicDynamic()
{
    while (direction * tstart < direction * tend) {
        try {
            basicIntegrator = DynIntegrator::With();
            basicIntegrator->setSystem(this);
            basicIntegrator->run();
            break;
        }
        catch (const DiscontinuityError& error) {
            const auto& types = error.types();
            const bool event = types && std::find(types->begin(), types->end(), EVENT) != types->end();
            if (!types
                || (!event && std::find(types->begin(), types->end(), TOUCHDOWN) == types->end()
                    && std::find(types->begin(), types->end(), LIFTOFF) == types->end())) {
                throw;
            }
            // A joint stop is a perfectly inelastic impact.  Project the
            // interpolated velocity onto the newly active constraint, then
            // restart the DAE with its new equation count.
            const auto restartStage = [this](const char* stage, const auto& operation) {
                try {
                    operation();
                }
                catch (const std::exception& restartError) {
                    throw SimulationStoppingError(
                        std::string("Dynamic discontinuity ") + stage + " failed: " + restartError.what()
                    );
                }
            };
            // Interpolation locates the impact, and position IC removes the
            // remaining interpolation error before the new constraint is
            // used by the restarted DAE.
            do {
                restartStage("position projection", [this] { runPosIC(); });
                restartStage("velocity projection", [this] { runVelIC(); });
                restartStage("acceleration restart", [this] { runAccIC(); });
                restartStage("limit release", [this] { releaseSeparatingLimits(); });
            } while (event && system->dynamicEvents->settle(tstart));
            if (event) output();
        }
    }
}

void SystemSolver::runBasicKinematic()
{
	if (tstart == tend) return;
	try {
		basicIntegrator = CREATE<KineIntegrator>::With();
		basicIntegrator->setSystem(this);
		basicIntegrator->run();
	}
	catch (const NotKinematicError& ex) {
		this->runQuasiKinematic();
	}
}

void SystemSolver::runPreDrag()
{
	initializeLocally();
	initializeGlobally();
	//Redundant constraints are removed here.
	runPosIC();
}

void MbD::SystemSolver::runDragStep(std::shared_ptr<std::vector<std::shared_ptr<Part>>> dragParts)
{
	//Assume no redundant constraints
	runPosICDrag(dragParts);
	runPosICDragLimit(dragParts);
}

void SystemSolver::runQuasiKinematic()
{
	try {
		basicIntegrator = CREATE<ICKineIntegrator>::With();
		basicIntegrator->setSystem(this);
		basicIntegrator->run();
	}
	catch (const DiscontinuityError& ex) {
		this->discontinuityBlock();
	}
}

void SystemSolver::runPosKine()
{
	icTypeSolver = CREATE<PosKineNewtonRaphson>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::runVelKine()
{
	icTypeSolver = CREATE<VelKineSolver>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::runAccKine()
{
	icTypeSolver = CREATE<AccKineNewtonRaphson>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void MbD::SystemSolver::runPosICDrag(std::shared_ptr<std::vector<std::shared_ptr<Part>>> dragParts)
{
	//Assume no redundant constraints
	auto newtonRaphson = PosICDragNewtonRaphson::With();
	newtonRaphson->setdragParts(dragParts);
	icTypeSolver = newtonRaphson;
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void MbD::SystemSolver::runPosICDragLimit(
	std::shared_ptr<std::vector<std::shared_ptr<Part>>> dragParts
)
{
	//Assume no redundant constraints
	auto newtonRaphson = PosICDragLimitNewtonRaphson::With();
	newtonRaphson->setdragParts(dragParts);
	icTypeSolver = newtonRaphson;
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::runPosICKine()
{
	icTypeSolver = CREATE<PosICKineNewtonRaphson>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::runVelICKine()
{
	icTypeSolver = CREATE<VelICKineSolver>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::runAccICKine()
{
	icTypeSolver = CREATE<AccICKineNewtonRaphson>::With();
	icTypeSolver->setSystem(this);
	icTypeSolver->run();
}

void SystemSolver::partsJointsMotionsDo(const std::function<void(std::shared_ptr<Item>)>& f)
{
	if (system->runMode == System::RunMode::Dynamic) system->partsJointsMotionsLimitsDo(f);
	else system->partsJointsMotionsDo(f);
}

void SystemSolver::logString(const std::string& str)
{
	system->logString(str);
}

std::shared_ptr<std::vector<std::shared_ptr<Part>>> SystemSolver::parts()
{
	return system->parts;
}

std::shared_ptr<std::vector<std::shared_ptr<LimitIJ>>> MbD::SystemSolver::limits()
{
	return system->limits;
}

std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> SystemSolver::essentialConstraints()
{
	return system->essentialConstraints();
}

std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> SystemSolver::displacementConstraints()
{
	return system->displacementConstraints();
}

std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> SystemSolver::perpendicularConstraints()
{
	return system->perpendicularConstraints();
}

std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> SystemSolver::allRedundantConstraints()
{
	return system->allRedundantConstraints();
}

std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> SystemSolver::allConstraints()
{
	return system->allConstraints();
}

std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> MbD::SystemSolver::allConstraintsLimits()
{
	return system->allConstraintsLimits();
}

void SystemSolver::postNewtonRaphson()
{
	throw SimulationStoppingError("To be implemented.");
}

void SystemSolver::partsJointsMotionsForcesTorquesDo(const std::function<void(std::shared_ptr<Item>)>& f)
{
	if (system->runMode == System::RunMode::Dynamic) system->partsJointsMotionsLimitsForcesTorquesDo(f);
	else system->partsJointsMotionsForcesTorquesDo(f);
}

void MbD::SystemSolver::partsJointsMotionsLimitsDo(const std::function<void(std::shared_ptr<Item>)>& f)
{
	system->partsJointsMotionsLimitsDo(f);
}

void MbD::SystemSolver::partsJointsMotionsLimitsForcesTorquesDo(const std::function<void(std::shared_ptr<Item>)>& f)
{
	system->partsJointsMotionsLimitsForcesTorquesDo(f);
}

void SystemSolver::discontinuityBlock()
{
	throw SimulationStoppingError("To be implemented.");
}

double SystemSolver::startTime()
{
	return tstart;
}

double SystemSolver::outputStepSize()
{
	return hout;
}

double SystemSolver::maxStepSize()
{
	return hmax;
}

double SystemSolver::minStepSize()
{
	return hmin;
}

double SystemSolver::firstOutputTime()
{
	return toutFirst;
}

double SystemSolver::endTime()
{
	return tend;
}

void SystemSolver::settime(double tnew)
{
	system->mbdTimeValue(tnew);
}

void SystemSolver::tstartPastsAddFirst(double tstartPast)
{
	tstartPasts->insert(tstartPasts->begin(), tstartPast);
}

void SystemSolver::output()
{
	system->outputFor(DYNAMIC);
}

void SystemSolver::time(double t)
{
	system->mbdTimeValue(t);
}

bool MbD::SystemSolver::limitsSatisfied()
{
	return 	system->limitsSatisfied();
}

void MbD::SystemSolver::deactivateLimits()
{
	system->deactivateLimits();
}
