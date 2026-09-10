/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

//#include <vector>

#include "DAEIntegrator.h"
#include "FullMatrix.h"
#include "enum.h"

namespace MbD {
    class DynIntegrator : public DAEIntegrator
    {
        //
    public:
        static std::shared_ptr<DynIntegrator> With();

        void assignEquationNumbers() override;
        void checkForDiscontinuity() override;
        void checkForOutputThrough(double t) override;
        void fillF(FColDsptr vecF) override;
        void fillpFpy(SpMatDsptr mat) override;
        void fillpFpydot(SpMatDsptr mat) override;
        void fillY(FColDsptr y);
        void fillYdot(FColDsptr ydot);
        void incrementTime(double t) override;
        void interpolateAt(double t) override;
        size_t iterMax() override;
        void postDAECorrector();
        bool acceptTrial() const;
        void postDAECorrectorIteration();
        void postDAEFirstStep();
        void postDAEOutput() override;
        void postDAEPredictor();
        void postDAEStep();
        void postRun() override;
        void preDAECorrector();
        void preDAECorrectorIteration();
        void preDAEFirstStep();
        void preDAEOutput() override;
        void preDAEPredictor();
        void preDAEStep();
        void preRun() override;
        void run() override;
        double suggestSmallerOrAcceptFirstStepSize(double hnew) override;
        double suggestSmallerOrAcceptStepSize(double hnew) override;
        void updateForDAECorrector();
        void y(FColDsptr col) override;
        void ydot(FColDsptr col) override;
        void throwDiscontinuityError(const std::string& str, std::shared_ptr<std::vector<DiscontinuityType>> discontinuityTypes);
        void useTrialStepStats(std::shared_ptr<SolverStatistics> stats);
        void useDAEStepStats(std::shared_ptr<SolverStatistics> stats);
        void reportStats() override;

    };
}
