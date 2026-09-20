/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <memory>

#include "DAEIntegrator.h"
#include "Item.h"
#include "SystemSolver.h"
#include "BasicQuasiIntegrator.h"
#include "SingularMatrixError.h"
#include "SimulationStoppingError.h"
#include "TooSmallStepSizeError.h"
#include "TooManyTriesError.h"
#include "SingularMatrixError.h"
#include "DiscontinuityError.h"
#include "StartingBasicDAEIntegrator.h"
#include "NormalBasicDAEIntegrator.h"

using namespace MbD;

std::shared_ptr<DAEIntegrator> DAEIntegrator::With()
{
    auto inst = std::make_shared<DAEIntegrator>();
    inst->initialize();
    return inst;
}

void DAEIntegrator::initialize()
{
    Solver::initialize();
    integrator = StartingBasicDAEIntegrator::With();
    integrator->setSystem(this);
}

void DAEIntegrator::initializeGlobally()
{
    IntegratorInterface::initializeGlobally();
    tout = std::min(tout, tend);
    assignEquationNumbers();
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->useEquationNumbers(); });
}

void DAEIntegrator::preFirstStep()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->preDynFirstStep(); });
}

void DAEIntegrator::checkForOutputThrough(double t)
{
    throw SimulationStoppingError("To be implemented.");
}

void DAEIntegrator::preRun()
{
    throw SimulationStoppingError("To be implemented.");
}

void DAEIntegrator::preStep()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->preDynStep(); });
}

double DAEIntegrator::suggestSmallerOrAcceptStepSize(double hnew)
{
    throw SimulationStoppingError("To be implemented.");    //May not be used?
    return 0.0;
    //auto hnew2 = system->suggestSmallerOrAcceptStepSize(hnew);
    //if (hnew2 > hmax) {
    //    hnew2 = hmax;
    //    system->logString("MbD: Step size is at user specified maximum.");
    //}
    //if (hnew2 < hmin) {
    //    std::stringstream ss;
    //    ss << "MbD: Step size " << hnew2 << " < " << hmin << "  user specified minimum.";
    //    auto str = ss.str();
    //    system->logString(str);
    //    throw TooSmallStepSizeError("");
    //}
    //return hnew2;
}

void DAEIntegrator::postRun()
{
    throw SimulationStoppingError("To be implemented.");
}

void DAEIntegrator::runInitialConditionTypeSolution()
{
    throw SimulationStoppingError("To be implemented.");
}

void DAEIntegrator::iStep(size_t i)
{
    throw SimulationStoppingError("To be implemented.");
}

void DAEIntegrator::selectOrder()
{
    throw SimulationStoppingError("To be implemented.");
}

void DAEIntegrator::checkForDiscontinuity()
{
    throw SimulationStoppingError("To be implemented.");
}

double DAEIntegrator::suggestSmallerOrAcceptFirstStepSize(double hnew)
{
    throw SimulationStoppingError("To be implemented.");
    return 0.0;
}

FColDsptr DAEIntegrator::integrationRelativeTolerance()
{
    //"Answer column of tolerances used by the integration error estimator."
    //"Algebraic variables are not included in the error estimator."

    auto relTol = system->integrationRelativeTolerance();
    auto col = std::make_shared<FullColumn<double>>(neqn);
    for (size_t i = 0; i < neqn - ncon; i++)
    {
        col->atiput(i, relTol);
    }
    for (size_t i = neqn - ncon; i < neqn; i++)
    {
        col->atiput(i, std::numeric_limits<double>::min());
    }
    return col;
}

FColDsptr DAEIntegrator::integrationAbsoluteTolerance()
{
    auto absTol = system->integrationAbsoluteTolerance();
    auto col = std::make_shared<FullColumn<double>>(neqn);
    for (size_t i = 0; i < neqn - ncon; i++)
    {
        col->atiput(i, absTol);
    }
    for (size_t i = neqn - ncon; i < neqn; i++)
    {
        col->atiput(i, std::numeric_limits<double>::min());
    }
    return col;
}

FColDsptr DAEIntegrator::correctorRelativeTolerance()
{
    auto corRelTol = system->correctorRelativeTolerance();
    auto col = std::make_shared<FullColumn<double>>(neqn);
    for (size_t i = 0; i < neqn; i++)
    {
        col->atiput(i, corRelTol);
    }
    return col;
}

FColDsptr DAEIntegrator::correctorAbsoluteTolerance()
{
    auto corAbsTol = system->correctorAbsoluteTolerance();
    auto col = std::make_shared<FullColumn<double>>(neqn);
    for (size_t i = 0; i < neqn; i++)
    {
        col->atiput(i, corAbsTol);
    }
    return col;
}

void DAEIntegrator::y(FColDsptr col)
{
    throw SimulationStoppingError("To be implemented.");
    //system->y(col);
}

void DAEIntegrator::ydot(FColDsptr col)
{
    throw SimulationStoppingError("To be implemented.");
    //system->ydot(col);
}

void DAEIntegrator::preDAEOutput()
{
    throw SimulationStoppingError("To be implemented.");
}

void DAEIntegrator::postDAEOutput()
{
    throw SimulationStoppingError("To be implemented.");
}

void DAEIntegrator::useTrialStepStats(std::shared_ptr<SolverStatistics> stats)
{
    throw SimulationStoppingError("To be implemented.");
    //system->useDAETrialStepStats(stats);
}

void DAEIntegrator::useDAEStepStats(std::shared_ptr<SolverStatistics> stats)
{
    system->useDAEStepStats(stats);
}

void DAEIntegrator::run()
{
    preRun();
    initializeLocally();
    initializeGlobally();
    if (tout <= tend) {
        integrator->run();
        auto startingintegrator = std::dynamic_pointer_cast<StartingBasicDAEIntegrator>(integrator);
        auto normalIntegrator = std::make_shared<NormalBasicDAEIntegrator>(startingintegrator);
        integrator = normalIntegrator;
        integrator->run();
    }
    finalize();
    reportStats();
    postRun();
}
