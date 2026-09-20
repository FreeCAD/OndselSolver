/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#pragma once

#include <memory>
#include "Solver.h"

namespace MbD {
    class SystemSolver;

    class Integrator : public Solver
    {
        //system direction 
    public:
        static std::shared_ptr<Integrator> With();

        void setSystem(Solver* sys) override;
        void run() override;
        virtual void firstStep();
        virtual void preFirstStep();
        virtual void postFirstStep();
        virtual void subsequentSteps();
        virtual void nextStep();
        virtual void preStep();
        virtual void postStep();
        virtual void runInitialConditionTypeSolution();
        virtual void iStep(size_t i);
        virtual void selectOrder();
        virtual void selectStepSize();
        virtual size_t iterMax();

        double direction = 1;
    };
}

