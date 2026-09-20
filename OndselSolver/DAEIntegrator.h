/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

//#include <vector>

#include "IntegratorInterface.h"
#include "enum.h"
#include "FullColumn.h"

namespace MbD {
    class DAEIntegrator : public IntegratorInterface
    {
        //neqn ncon
    public:
        static std::shared_ptr<DAEIntegrator> With();
        void initialize() override;

        void initializeGlobally() override;
        void preFirstStep() override;
        void checkForOutputThrough(double t) override;
        void preRun() override;
        void preStep() override;
        double suggestSmallerOrAcceptStepSize(double hnew) override;
        //void incrementTime(double tnew) override;
        //void throwDiscontinuityError(const std::string& str, std::shared_ptr<std::vector<DiscontinuityType>> discontinuityTypes);
        //void interpolateAt(double t) override;
        void postRun() override;
        void runInitialConditionTypeSolution() override;
        void iStep(size_t i) override;
        void selectOrder() override;
        void checkForDiscontinuity() override;
        double suggestSmallerOrAcceptFirstStepSize(double hnew) override;
        FColDsptr integrationRelativeTolerance();
        FColDsptr integrationAbsoluteTolerance();
        FColDsptr correctorRelativeTolerance();
        FColDsptr correctorAbsoluteTolerance();
        void y(FColDsptr col) override;
        void ydot(FColDsptr col) override;
        virtual void preDAEOutput();
        virtual void postDAEOutput();
        void useTrialStepStats(std::shared_ptr<SolverStatistics> stats) override;
        void useDAEStepStats(std::shared_ptr<SolverStatistics> stats) override;
        void run() override;

        size_t neqn = SIZE_MAX, ncon = SIZE_MAX;

    };
}
