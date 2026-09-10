/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "BasicIntegrator.h"
#include "FullColumn.h"
#include "SparseMatrix.h"
#include "Extrapolator.h"
#include "DAECorrector.h"

namespace MbD {
    class BasicDAEIntegrator : public BasicIntegrator
    {
        //y ydot dy ypast ydotpast aF pFpy pFpydot alp aG extrapolator newtonRaphson corAbsTol corRelTol corOK integAbsTol integRelTol truncError
    public:
        static std::shared_ptr<BasicDAEIntegrator> With();
        void initialize() override;

        void initializeGlobally() override;
        void firstStep() override;
        bool isRedoingFirstStep();
        bool isRedoingStep() const;
        void nextStep() override;
        void reportStepStats();
        void reportTrialStepStats();
        void runInitialConditionTypeSolution() override;
        void selectFirstStepSize();
        void selectStepSize() override;
        void predictFirstStep();
        void correctFirstStep();
        void changeTime();
        double selectBasicStepSize();
        double integErrorNormFromwrt(FColDsptr err, FColDsptr ref);
        double rmswrtrelativeTolabsoluteTol(FColDsptr vector, FColDsptr baseVector, FColDsptr relativeTol, FColDsptr absoluteTol);
        void predictValuesAtFirstStep();
        void predictValuesAtNextStep();
        void incrementTime() override;
        void correctValuesAtFirstStep() const;
        void correctValuesAtNextStep() const;
        void preDAECorrector();
        void postDAECorrector();
        FColDsptr fillF();
        size_t iterMax() override;
        SpMatDsptr calcG();
        virtual std::shared_ptr<LinearMultiStepMethod> correctorBDF();
        void calcOperatorMatrix() override;
        void setorder(size_t o) override;
        void settime(double t) override;
        void iStep(size_t i) override;
        double corErrorNormFromwrt(FColDsptr error, FColDsptr ref);
        void updateForDAECorrector();
        virtual FColDsptr yDeriv(size_t order);
        void calcTruncError();
        virtual FColDsptr dyOrderPlusOnedt();
        bool isConvergedForand(size_t iterNo, std::shared_ptr<std::vector<double>> dyNorms) const;
        void postFirstStep();
        void postStep();
        void predict();
        void correct();
        void selectStepSizeNormal();
        void useDAECorrectorStats(std::shared_ptr<SolverStatistics> stats) const;

        FColDsptr y, ydot, dy, aF, corAbsTol, corRelTol, integAbsTol, integRelTol;
        std::shared_ptr<std::vector<FColDsptr>> ypast, ydotpast;
        SpMatDsptr pFpy, pFpydot, matG;
        double alp = 0.0, truncError = 0.0;
        bool corOK = false;
        bool firstStepGrowthAllowed = true;
        std::shared_ptr<Extrapolator> extrapolator;
        std::shared_ptr<DAECorrector> newtonRaphson;
    };
}
