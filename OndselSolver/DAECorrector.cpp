/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <fstream>

#include "DAECorrector.h"
#include "BasicDAEIntegrator.h"
#include "GESpMatParPvMarkoFast.h"
#include "GESpMatParPvPrecise.h"
#include "SystemSolver.h"
#include "SimulationStoppingError.h"

#include "CREATE.h"

using namespace MbD;

std::shared_ptr<DAECorrector> DAECorrector::With()
{
    auto inst = std::make_shared<DAECorrector>();
    inst->initialize();
    return inst;
}

void DAECorrector::iterate()
{
    //VectorNewtonRaphson::iterate();    //Inlined to help debugging
    iterNo = SIZE_MAX;
    fillY();
    calcyNorm();
    yNorms->push_back(yNorm);

    while (true) {
        incrementIterNo();
        fillPyPx();
        //outputSpreadsheet();
        solveEquations();
        calcDXNormImproveRootCalcYNorm();
        if (isConverged()) {
            //std::cout << "iterNo = " << iterNo << std::endl;
            break;
        }
    }
}

void DAECorrector::fillPyPx()
{
    pypx = daeSystem->calcG();
}

void DAECorrector::fillY()
{
    y = daeSystem->fillF();
}

void DAECorrector::passRootToSystem()
{
    daeSystem->y = x;
}

void DAECorrector::calcdxNorm()
{
    dxNorm = daeSystem->corErrorNormFromwrt(dx, x);
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<double>::max_digits10);
    ss << "MbD: Convergence = " << dxNorm;
    auto str = ss.str();
    daeSystem->logString(str);
}

void DAECorrector::basicSolveEquations()
{
    dx = matrixSolver->solvewithsaveOriginal(pypx, y->negated(), false);
}

void DAECorrector::handleSingularMatrix()
{
    if (std::dynamic_pointer_cast<GESpMatParPvMarkoFast>(matrixSolver)) {
        matrixSolver = CREATE<GESpMatParPvPrecise>::With();
        solveEquations();
    }
    else {
        matrixSolver->throwSingularMatrixError("Singular forward-dynamics corrector matrix");
    }
}

void DAECorrector::initializeGlobally()
{
    iterMax = daeSystem->iterMax();
    x = daeSystem->y;
    matrixSolver = CREATE<GESpMatParPvMarkoFast>::With();
}

void DAECorrector::run()
{
    preRun();
    initializeLocally();
    initializeGlobally();
    iterate();
    finalize();
    reportStats();
    postRun();
}

void DAECorrector::preRun()
{
    //auto basicDAEIntegrator = static_cast<BasicDAEIntegrator*>(system);
    daeSystem->preDAECorrector();
}

void DAECorrector::askSystemToUpdate()
{
    daeSystem->updateForDAECorrector();
}

bool DAECorrector::isConverged()
{
    return daeSystem->isConvergedForand(iterNo, dxNorms);
}

void DAECorrector::postRun()
{
    daeSystem->postDAECorrector();
}

void DAECorrector::setSystem(Solver* sys)
{
    daeSystem = static_cast<BasicDAEIntegrator*>(sys);
}

void DAECorrector::reportStats()
{
    statistics->iterNo = iterNo;
    daeSystem->useDAECorrectorStats(statistics);
}

void DAECorrector::outputSpreadsheet()
{
    std::ofstream os("../testapp/spreadsheetcpp.csv");
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (size_t i = 0; i < pypx->nrow(); i++)
    {
        auto rowi = pypx->at(i);
        for (size_t j = 0; j < pypx->ncol(); j++)
        {
            if (rowi->find(j) == rowi->end()) {
                os << 0.0;
            }
            else {
                os << rowi->at(j);
            }
            os << '\t';
        }
        os << "\t" << y->at(i) << std::endl;
    }
}
