/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "BasicDAEIntegrator.h"
#include "StableStartingBDF.h"

namespace MbD {
    class StartingBasicDAEIntegrator : public BasicDAEIntegrator
    {
        //
    public:
        static std::shared_ptr<StartingBasicDAEIntegrator> With();
        void initialize() override;

        void initializeLocally() override;
        FColDsptr yDerivat(size_t _order, double tout) override;
        std::shared_ptr<LinearMultiStepMethod> correctorBDF() override;
        void calcOperatorMatrix() override;
        void setorder(size_t o) override;
        void settime(double t) override;
        void iStep(size_t i) override;
        FColDsptr yDeriv(size_t order);
        FColDsptr dyOrderPlusOnedt() override;
        void run() override;
        void firstSteps();
        void earlySteps();

        std::shared_ptr<StableStartingBDF> startingBDF;
    };
}
