/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#pragma once

#include "Integrator.h"
#include "FullColumn.h"
#include "SparseMatrix.h"

namespace MbD {
    class BasicIntegrator;
    class StartingBasicDAEIntegrator;
    class NormalBasicDAEIntegrator;

    class IntegratorInterface : public Integrator
    {
        //tout hout hmin hmax tstart tend integrator 
    public:
        static std::shared_ptr<IntegratorInterface> With();

        void initializeGlobally() override;
        void run() override;
        virtual void preRun() override;
        virtual void checkForDiscontinuity();
        
        void setSystem(Solver* sys) override;
        void logString(const std::string& str) override;
        size_t orderMax() const;
        virtual void incrementTime(double tnew);

        void postFirstStep() override;
        virtual double suggestSmallerOrAcceptFirstStepSize(double hnew);
        virtual double suggestSmallerOrAcceptStepSize(double hnew);
        virtual void checkForOutputThrough(double t);
        virtual void interpolateAt(double t);
        virtual void fillF(FColDsptr vecF);
        virtual void fillpFpy(SpMatDsptr mat);
        virtual void fillpFpydot(SpMatDsptr mat);
        void changeTime(double t);
        virtual void y(FColDsptr col);
        virtual void ydot(FColDsptr col);
        virtual void updateForDAECorrector();
        virtual void useTrialStepStats(std::shared_ptr<SolverStatistics> stats);
        virtual void useDAEStepStats(std::shared_ptr<SolverStatistics> stats);
        virtual void useQuasiStepStats(std::shared_ptr<SolverStatistics> stats);

        SystemSolver* system = nullptr;
        double tout = 0.0, hout = 0.0, hmin = 0.0, hmax = 0.0, tstart = 0.0, tend = 0.0;
        std::shared_ptr<BasicIntegrator> integrator;
    };
}

