/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <memory>

#include "NormalBasicDAEIntegrator.h"
#include "LinearMultiStepMethod.h"
#include "StableBackwardDifference.h"
#include "StartingBasicDAEIntegrator.h"
#include "SimulationStoppingError.h"

#include "CREATE.h"

using namespace MbD;

NormalBasicDAEIntegrator::NormalBasicDAEIntegrator(std::shared_ptr<StartingBasicDAEIntegrator> startingBasicDAEIntegrator)
{
    //system direction
    //istep iTry maxTry tpast t tnew h hnew order orderNew orderMax opBDF continue
    //y ydot dy ypast ydotpast aF pFpy pFpydot alp aG extrapolator newtonRaphson corAbsTol corRelTol corOK integAbsTol integRelTol truncError
    system = startingBasicDAEIntegrator->system;
    direction = startingBasicDAEIntegrator->direction;
    istep = startingBasicDAEIntegrator->istep;
    iTry = startingBasicDAEIntegrator->iTry;
    maxTry = startingBasicDAEIntegrator->maxTry;
    tpast = startingBasicDAEIntegrator->tpast;
    t = startingBasicDAEIntegrator->t;
    tnew = startingBasicDAEIntegrator->tnew;
    h = startingBasicDAEIntegrator->h;
    hnew = startingBasicDAEIntegrator->hnew;
    order = startingBasicDAEIntegrator->order;
    orderNew = startingBasicDAEIntegrator->orderNew;
    orderMax = startingBasicDAEIntegrator->orderMax;
    opBDF = startingBasicDAEIntegrator->opBDF;
    _continue = startingBasicDAEIntegrator->_continue;
    y = startingBasicDAEIntegrator->y;
    ydot = startingBasicDAEIntegrator->ydot;
    dy = startingBasicDAEIntegrator->dy;
    ypast = startingBasicDAEIntegrator->ypast;
    ydotpast = startingBasicDAEIntegrator->ydotpast;
    aF = startingBasicDAEIntegrator->aF;
    pFpy = startingBasicDAEIntegrator->pFpy;
    pFpydot = startingBasicDAEIntegrator->pFpydot;
    alp = startingBasicDAEIntegrator->alp;
    matG = startingBasicDAEIntegrator->matG;
    extrapolator = startingBasicDAEIntegrator->extrapolator;
    newtonRaphson = startingBasicDAEIntegrator->newtonRaphson;
    newtonRaphson->setSystem(this);
    corAbsTol = startingBasicDAEIntegrator->corAbsTol;
    corRelTol = startingBasicDAEIntegrator->corRelTol;
    corOK = startingBasicDAEIntegrator->corOK;
    integAbsTol = startingBasicDAEIntegrator->integAbsTol;
    integRelTol = startingBasicDAEIntegrator->integRelTol;
    truncError = startingBasicDAEIntegrator->truncError;
    opBDFhigher = CREATE<StableBackwardDifference>::With();
    opBDFhigher->timeNodes = tpast;
    opBDFhigher->time = t;
    opBDFhigher->iStep = istep;
    opBDFhigher->order = order + 1;
    //calcOperatorMatrix();
}

std::shared_ptr<NormalBasicDAEIntegrator> NormalBasicDAEIntegrator::With()
{
    auto inst = std::make_shared<NormalBasicDAEIntegrator>();
    inst->initialize();
    return inst;
}

void NormalBasicDAEIntegrator::initialize()
{
    BasicDAEIntegrator::initialize();
    opBDFhigher = CREATE<StableBackwardDifference>::With();
    opBDFhigher->timeNodes = tpast;
}

void NormalBasicDAEIntegrator::initializeLocally()
{
    //"NormalBasicDAEIntegrator is not used for starting an integration."
    //"Change to StartingBasicDAEIntegrator."
    throw SimulationStoppingError("To be implemented.");
}

FColDsptr NormalBasicDAEIntegrator::yDerivat(size_t n, double time)
{
    // An order-one position polynomial has no second derivative, but the DAE
    // state also carries its first derivative. Interpolate that history to
    // retain valid acceleration output when orderMax == 1.
    if (n > order) {
        return opBDF->derivativeatpresentpast(n - 1, time, ydot, ydotpast);
    }
    return opBDF->derivativeatpresentpast(n, time, y, ypast);
}

std::shared_ptr<LinearMultiStepMethod> NormalBasicDAEIntegrator::correctorBDF()
{
    return opBDF;
}

void NormalBasicDAEIntegrator::calcOperatorMatrix()
{
    BasicDAEIntegrator::calcOperatorMatrix();
    opBDFhigher->calcOperatorMatrix();
}

void NormalBasicDAEIntegrator::settime(double t)
{
    BasicDAEIntegrator::settime(t);
    opBDFhigher->settime(t);
}

void NormalBasicDAEIntegrator::iStep(size_t i)
{
    BasicDAEIntegrator::iStep(i);
    opBDFhigher->setiStep(i);
}

void NormalBasicDAEIntegrator::setorder(size_t o)
{
    BasicDAEIntegrator::setorder(o);
    opBDFhigher->setorder(o + 1);
}

FColDsptr NormalBasicDAEIntegrator::yDeriv(size_t deriv)
{
    return opBDF->derivativepresentpast(deriv, y, ypast);
}

FColDsptr NormalBasicDAEIntegrator::dyOrderPlusOnedt()
{
    return opBDFhigher->derivativepresentpast(order + 1, y, ypast);
}

void NormalBasicDAEIntegrator::run()
{
    subsequentSteps();
    finalize();
    reportStats();
    postRun();
}

void NormalBasicDAEIntegrator::selectOrder()
{
    if (tpast->size() < order + 1) return;    //Needed to transition from Starting to Normal
    selectOrderNormal();
}

void NormalBasicDAEIntegrator::selectStepSize()
{
    if (tpast->size() < order + 1) return;    //Needed to transition from Starting to Normal
    BasicDAEIntegrator::selectStepSize();
}

void NormalBasicDAEIntegrator::selectOrderNormal()
{
    //"Brenan's book pp. 126-7"
    //"Check last nterm of Taylor series plus first term of remainder."

    auto errorTrunc = std::make_shared<FullVector<double>>();
    auto nterm = 3;
    size_t istart;
    FColDsptr yndot;
    double yndotNorm, hpower;
    if (order > nterm) {
        istart = order - nterm + 1;
    }
    else {
        istart = 1;
    }
    for (size_t i = istart; i < order; i++)
    {
        yndot = yDeriv(i);
        yndotNorm = integErrorNormFromwrt(yndot, y);
        hpower = std::pow(h, i);
        errorTrunc->push_back(yndotNorm * hpower);
    }


    yndot = dyOrderPlusOnedt();
    yndotNorm = integErrorNormFromwrt(yndot, y);
    hpower = std::pow(h, order + 1);
    errorTrunc->push_back(yndotNorm * hpower);
    orderNew = order;
    if (errorTrunc->isIncreasing()) {
        if (order > 1) orderNew = order - 1;
    }
    if (errorTrunc->isDecreasingIfExceptionsAreLessThan(0.01)) {
        if (order < orderMax) orderNew = order + 1;
    }
}
