/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include <string>
#include <cassert>

#include "IntegratorInterface.h"
#include "SystemSolver.h"
#include "BasicQuasiIntegrator.h"
#include "NormalBasicDAEIntegrator.h"
#include "StartingBasicDAEIntegrator.h"

using namespace MbD;

std::shared_ptr<IntegratorInterface> IntegratorInterface::With()
{
    //Should not create abstract class.
    throw SimulationStoppingError("To be implemented.");
    return std::shared_ptr<IntegratorInterface>();
}

void IntegratorInterface::initializeGlobally()
{
    tstart = system->startTime();
    hout = system->outputStepSize();
    hmax = system->maxStepSize();
    hmin = system->minStepSize();
    tout = system->firstOutputTime();
    tend = system->endTime();
    direction = (tstart < tend) ? 1.0 : -1.0;
}

void IntegratorInterface::preRun()
{
    //Subclasses must implement.
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::checkForDiscontinuity()
{
    //Subclasses must implement.
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::setSystem(Solver* sys)
{
    system = static_cast<SystemSolver*>(sys);
}

void IntegratorInterface::logString(const std::string& str)
{
    system->logString(str);
}

size_t IntegratorInterface::orderMax() const
{
    return system->orderMax;
}

void IntegratorInterface::incrementTime(double tnew)
{
    system->settime(tnew);
}

void IntegratorInterface::postFirstStep()
{
    throw SimulationStoppingError("To be implemented.");
}

double IntegratorInterface::suggestSmallerOrAcceptFirstStepSize(double hnew)
{
    //Subclasses must implement.
    throw SimulationStoppingError("To be implemented.");
    return 0.0;
}

double IntegratorInterface::suggestSmallerOrAcceptStepSize(double hnew)
{
    //Subclasses must implement.
    throw SimulationStoppingError("To be implemented.");
    return 0.0;
}

void IntegratorInterface::checkForOutputThrough(double t)
{
    //Subclasses must implement.
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::interpolateAt(double tArg)
{
    //"Interpolate for system state at tArg and leave system in that state."
    throw SimulationStoppingError("To be implemented.");
    //auto yout = integrator->yDerivat(0, tArg);
    //auto ydotout = integrator->yDerivat(1, tArg);
    //auto yddotout = integrator->yDerivat(2, tArg);
    //system->time(tArg);
    //system->y(yout);
    //system->ydot(ydotout);
    //system->yddot(yddotout);
    //system->simUpdateAll();
}

void IntegratorInterface::fillF(FColDsptr vecF)
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::fillpFpy(SpMatDsptr mat)
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::fillpFpydot(SpMatDsptr mat)
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::changeTime(double t)
{
    system->settime(t);
}

void IntegratorInterface::y(FColDsptr col)
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::ydot(FColDsptr col)
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::updateForDAECorrector()
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::useTrialStepStats(std::shared_ptr<SolverStatistics> stats)
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::useDAEStepStats(std::shared_ptr<SolverStatistics> stats)
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::useQuasiStepStats(std::shared_ptr<SolverStatistics> stats)
{
    throw SimulationStoppingError("To be implemented.");
}

void IntegratorInterface::run()
{
	this->preRun();
	this->initializeLocally();
	this->initializeGlobally();
	if (hout > (4 * std::numeric_limits<double>::epsilon()) && (direction * tout < (direction * (tend + (0.1 * direction * hout))))) {
		integrator->run();
	}
	this->finalize();
	this->reportStats();
	this->postRun();
}
