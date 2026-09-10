/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "Integrator.h"
#include <cassert>
#include "SimulationStoppingError.h"

using namespace MbD;

std::shared_ptr<Integrator> Integrator::With()
{
    auto inst = std::make_shared<Integrator>();
    inst->initialize();
    return inst;
}

void Integrator::setSystem(Solver*)
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::run()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::firstStep()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::preFirstStep()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::postFirstStep()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::subsequentSteps()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::nextStep()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::preStep()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::postStep()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::runInitialConditionTypeSolution()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::iStep(size_t i)
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::selectOrder()
{
    throw SimulationStoppingError("To be implemented.");
}

void Integrator::selectStepSize()
{
    throw SimulationStoppingError("To be implemented.");
}

size_t Integrator::iterMax()
{
    throw SimulationStoppingError("To be implemented.");
    return 0;
}
