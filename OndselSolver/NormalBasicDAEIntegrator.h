/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "BasicDAEIntegrator.h"

namespace MbD {
    class StableBackwardDifference;
    class StartingBasicDAEIntegrator;

    class NormalBasicDAEIntegrator : public BasicDAEIntegrator
    {
        //
    public:
        NormalBasicDAEIntegrator() {}
        NormalBasicDAEIntegrator(std::shared_ptr<StartingBasicDAEIntegrator> startingBasicDAEIntegrator);
        static std::shared_ptr<NormalBasicDAEIntegrator> With();
        void initialize() override;

        void initializeLocally() override;
        FColDsptr yDerivat(size_t _order, double tout) override;
        std::shared_ptr<LinearMultiStepMethod> correctorBDF() override;
        void calcOperatorMatrix() override;
        void settime(double t) override;
        void iStep(size_t i) override;
        void setorder(size_t o) override;
        FColDsptr yDeriv(size_t order);
        FColDsptr dyOrderPlusOnedt() override;
        void run() override;
        void selectOrder() override;
        void selectOrderNormal();
        void selectStepSize() override;

        std::shared_ptr<StableBackwardDifference> opBDFhigher;
    };
}
