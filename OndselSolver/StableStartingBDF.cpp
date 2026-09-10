/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include "StableStartingBDF.h"
#include "FullColumn.h"

using namespace MbD;

std::shared_ptr<StableStartingBDF> StableStartingBDF::With()
{
    auto inst = std::make_shared<StableStartingBDF>();
    inst->initialize();
    return inst;
}

void StableStartingBDF::initialize()
{
    StableBackwardDifference::initialize();
    iStep = 0;
    order = iStep + 1;
}

void StableStartingBDF::initializeLocally()
{
    initialize();
}

double StableStartingBDF::pvdotpv()
{
    //"pvdotpv = operatorMatrix timesColumn: #(-1.0d ... -1.0d, 0.0d)."

    auto& coeffs = operatorMatrix->at(0);
    auto sum = 0.0;
    for (size_t i = 0; i < order - 1; i++)
    {
        sum -= coeffs->at(i);
    }
    return sum;
}

void StableStartingBDF::formTaylorMatrix()
{
    //"
    //This form is numerically more stable and is prefered over the full Taylor Matrix.
    //For method order 3:
    //|    (t1 - t)    (t1 - t)^2/2!        (t1 - t)^3/3!    |    |qd(t)    |        |q(t1) - q(t)    |
    //|    (t2 - t)    (t2 - t)^2/2!        (t2 - t)^3/3!    |    |qdd(t)    |        |q(t2) - q(t)    |
    //|        1        (t2 - t)            (t2 - t)^2/2!    |    |qddd(t)|        |qd(t2)         |
    //"

    instantiateTaylorMatrix();
    for (size_t i = 0; i < order - 1; i++)
    {
        formTaylorRowwithTimeNodederivative(i, i, 0);
    }
    formTaylorRowwithTimeNodederivative(order - 1, order - 2, 1);
}

void StableStartingBDF::setorder(size_t o)
{
    //"order is controlled by iStep."
    if ((order != o) && (order != o + 1)) throw std::runtime_error("iStep and order must be consistent.");
}

void StableStartingBDF::setiStep(size_t i)
{
    //"iStep is the current step of interest."
    //"iStep must increase consecutively."

    auto iStepNew = iStep + 1;
    if (iStepNew == i) {
        iStep = iStepNew;
        order = iStep + 1;
    }
    else {
        throw std::runtime_error("Not appropriate iStep");
    }
}

FColDsptr StableStartingBDF::derivativepresentpastpresentDerivativepastDerivative(size_t deriv,
    FColDsptr y, std::shared_ptr<std::vector<FColDsptr>> ypast,
    FColDsptr ydot, std::shared_ptr<std::vector<FColDsptr>> ydotpast)
{
    if (deriv == 0) return y->copy();
    auto series = std::make_shared<std::vector<FColDsptr>>(order);
    for (size_t j = 0; j < order - 1; j++)
    {
        series->at(j) = ypast->at(j)->minusFullColumn(y);
    }
    series->at(order - 1) = ydotpast->at(order - 2);
    auto& coeffs = operatorMatrix->at(deriv - 1);
    auto answer = coeffs->dot(series);
    return std::static_pointer_cast<FullColumn<double>>(answer);
}
