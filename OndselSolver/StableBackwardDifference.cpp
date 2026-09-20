/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include "StableBackwardDifference.h"
#include "FullColumn.h"

using namespace MbD;

void StableBackwardDifference::formTaylorMatrix()
{
	//This form is numerically more stable and is preferred over the full Taylor Matrix.
	//For method order 3:
	//| (t1 - t)	(t1 - t) ^ 2 / 2!	(t1 - t) ^ 3 / 3!|	|qd(t)  | = | q(t1) - q(t)	|
	//|	(t2 - t)	(t2 - t) ^ 2 / 2!	(t2 - t) ^ 3 / 3!|	|qdd(t)	|	|q(t2) - q(t)	|
	//|	(t3 - t)	(t3 - t) ^ 2 / 2!	(t3 - t) ^ 3 / 3!|	|qddd(t)|	|q(t3) - q(t)	|

	this->instantiateTaylorMatrix();
	for (size_t i = 0; i < order; i++)
	{
		this->formTaylorRowwithTimeNodederivative(i, i, 0);
	}
}

double MbD::StableBackwardDifference::pvdotpv()
{
	//"pvdotpv = operatorMatrix timesColumn: #(-1.0d ... -1.0d)."

	auto& coeffs = operatorMatrix->at(0);
	auto sum = 0.0;
	for (size_t i = 0; i < order; i++)
	{
		sum -= coeffs->at(i);
	}
	return sum;
}

FColDsptr MbD::StableBackwardDifference::derivativepresentpastpresentDerivativepastDerivative(
        size_t, FColDsptr, std::shared_ptr<std::vector<FColDsptr>>, FColDsptr,
        std::shared_ptr<std::vector<FColDsptr>>)
{
	throw SimulationStoppingError("To be implemented.");
	return FColDsptr();
}

void StableBackwardDifference::instantiateTaylorMatrix()
{
	if (taylorMatrix == nullptr || (taylorMatrix->nrow() != (order))) {
		taylorMatrix = std::make_shared<FullMatrix<double>>(order, order);
	}
}

void StableBackwardDifference::formTaylorRowwithTimeNodederivative(size_t i, size_t ii, size_t k)
{
	//| rowi hi hipower aij |
	auto& rowi = taylorMatrix->at(i);
	if (k > 0) {
		for (ssize_t j = 0; j < (ssize_t)k - 2; j++)
		{
			rowi->at(j) = 0.0;
		}
		rowi->at(k - 1) = 1.0;
	}

	auto hi = timeNodes->at(ii) - time;
	auto hipower = 1.0;
	for (size_t j = k; j < order; j++)
	{
		hipower *= hi;
		auto aij = hipower * OneOverFactorials->at(j - k + 1);
		rowi->at(j) = aij;
	}
}

FColDsptr MbD::StableBackwardDifference::derivativepresentpast(size_t deriv, FColDsptr y, std::shared_ptr<std::vector<FColDsptr>> ypast)
{
	//"Answer ith derivative given present value and past values."

	if (deriv == 0) {
		return y->copy();
	}
	else {
		if (deriv <= order) {
			auto series = std::make_shared<std::vector<FColDsptr>>(order);
			for (size_t i = 0; i < order; i++)
			{
				series->at(i) = ypast->at(i)->minusFullColumn(y);
			}
			auto& coeffs = operatorMatrix->at(deriv - 1);
			auto answer = coeffs->dot(series);
			return std::static_pointer_cast<FullColumn<double>>(answer);
		}
		else {
            auto ySize = y->size();
			return std::make_shared<FullColumn<double>>(ySize, 0.0);
		}
	}
}

FColDsptr StableBackwardDifference::derivativeatpresentpastpresentDerivativepastDerivative(size_t n, double t, FColDsptr y, std::shared_ptr<std::vector<FColDsptr>> ypast, FColDsptr ydot, std::shared_ptr<std::vector<FColDsptr>> ydotpast)
{
    //"Interpolate or extrapolate."
    //"dfdt(t) = df0dt + d2f0dt2*(t - t0) + d3f0dt3*(t - t0)^2 / 2! + ..."

    auto answer = derivativepresentpastpresentDerivativepastDerivative(n, y, ypast, ydot, ydotpast);
    if (t != time) {
        auto dt = t - time;
        auto dtpower = 1.0;
        for (size_t i = n + 1; i <= order; i++)
        {
            auto diydti = derivativepresentpastpresentDerivativepastDerivative(i, y, ypast, ydot, ydotpast);
            dtpower = dtpower * dt;
            answer->equalSelfPlusFullColumntimes(diydti, dtpower * OneOverFactorials->at(i - n));
        }
    }
    return answer;
}
