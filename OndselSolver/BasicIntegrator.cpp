/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include "BasicIntegrator.h"
#include "StableBackwardDifference.h"
#include "IntegratorInterface.h"
#include "SimulationStoppingError.h"

#include "CREATE.h"

using namespace MbD;

std::shared_ptr<BasicIntegrator> BasicIntegrator::With()
{
    auto inst = std::make_shared<BasicIntegrator>();
    inst->initialize();
    return inst;
}

void BasicIntegrator::initialize()
{
    Solver::initialize();
    //statistics = IdentityDictionary new.
    tpast = std::make_shared<std::vector<double>>();
    opBDF = CREATE<StableBackwardDifference>::With();
    opBDF->timeNodes = tpast;
}

void BasicIntegrator::initializeLocally()
{
    _continue = true;
}

void BasicIntegrator::iStep(size_t integer)
{
    istep = integer;
    opBDF->setiStep(integer);
}

void BasicIntegrator::postFirstStep()
{
    t = tnew;
    system->postFirstStep();
}

void BasicIntegrator::postRun()
{
    //Do nothing.
}

void BasicIntegrator::postStep()
{
    t = tnew;
    system->postStep();
}

void BasicIntegrator::initializeGlobally()
{
    //"Get info from system and prepare for start of simulation."
    //"Integrator asks system for info. Not system setting integrator."

    sett(system->tstart);
    direction = system->direction;
    orderMax = system->orderMax();
}

void BasicIntegrator::setSystem(Solver* sys)
{
    system = static_cast<IntegratorInterface*>(sys);
}

void BasicIntegrator::calcOperatorMatrix()
{
    opBDF->calcOperatorMatrix();
}

void BasicIntegrator::incrementTime()
{
    tpast->insert(tpast->begin(), t);

    if (tpast->size() > (orderMax + 1)) { tpast->pop_back(); }
    auto istepNew = istep + 1;
    iStep(istepNew);
    setorder(orderNew);
    h = hnew;
    settnew(t + (direction * h));
    calcOperatorMatrix();
    system->incrementTime(tnew);
}

void BasicIntegrator::incrementTry()
{
    iTry++;
}

void BasicIntegrator::logString(const std::string& str)
{
    system->logString(str);
}

void BasicIntegrator::run()
{
    preRun();
    initializeLocally();
    initializeGlobally();
    firstStep();
    subsequentSteps();
    finalize();
    reportStats();
    postRun();
}

void BasicIntegrator::selectOrder()
{
    //"Increase order consecutively with step."
    if (iTry == 1) orderNew = std::min(istep + 1, orderMax);
}

void BasicIntegrator::preFirstStep()
{
    system->preFirstStep();
}

void BasicIntegrator::preRun()
{
    //Do nothing.
}

void BasicIntegrator::preStep()
{
    system->preStep();
}

void BasicIntegrator::reportStats()
{
    //Do nothing.
}

void BasicIntegrator::setorder(size_t o)
{
    order = o;
    opBDF->setorder(o);
}

void BasicIntegrator::settnew(double t)
{
    tnew = t;
    settime(t);
}

void BasicIntegrator::sett(double tt)
{
    t = tt;
    opBDF->settime(tt);
}

void BasicIntegrator::settime(double tt)
{
    opBDF->settime(tt);
}

double BasicIntegrator::tprevious() const
{
    return tpast->at(0);
}

FColDsptr BasicIntegrator::yDerivat(size_t n, double time)
{
    throw SimulationStoppingError("To be implemented.");
    return FColDsptr();
}

void BasicIntegrator::subsequentSteps()
{
    while (_continue) { nextStep(); }
}
