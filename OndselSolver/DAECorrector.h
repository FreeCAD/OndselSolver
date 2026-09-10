/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "VectorNewtonRaphson.h"
#include "SparseMatrix.h"

namespace MbD {
    class BasicDAEIntegrator;

    class DAECorrector : public VectorNewtonRaphson
    {
        //
    public:
        static std::shared_ptr<DAECorrector> With();

        void iterate();
        void fillPyPx() override;
        void fillY() override;
        void passRootToSystem() override;
        void basicSolveEquations() override;
        void initializeGlobally() override;
        void calcdxNorm() override;
        void handleSingularMatrix() override;
        void run() override;
        void preRun() override;
        void askSystemToUpdate() override;
        bool isConverged() override;
        void postRun() override;
        void setSystem(Solver* sys) override;
        void reportStats() override;
        void outputSpreadsheet();

        SpMatDsptr pypx;
        BasicDAEIntegrator* daeSystem;
    };
}
