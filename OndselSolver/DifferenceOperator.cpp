/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <cmath>

#include "DifferenceOperator.h"
#include "CREATE.h"
#include "SingularMatrixError.h"
#include "LDUFullMatParPv.h"
#include "FullRow.h"

using namespace MbD;

void DifferenceOperator::formDegenerateTaylorRow(size_t i) const
{
    auto row = taylorMatrix->at(i);
    row->zeroSelf();
    row->at(0) = 1.0;
}

FColDsptr DifferenceOperator::valueWith(std::shared_ptr<std::vector<FColDsptr>> series)
{
    return derivativewith(0, series);
}

FColDsptr DifferenceOperator::derivativewith(size_t deriv, std::shared_ptr<std::vector<FColDsptr>> series) const
{
    const auto coefficients = operatorMatrix->at(deriv);
    auto result = series->at(0)->times(coefficients->at(0));
    for (size_t i = 1; i < coefficients->size(); ++i)
        result->equalSelfPlusFullVectortimes(series->at(i), coefficients->at(i));
    return result;
}

FRowDsptr DifferenceOperator::OneOverFactorials = []() {
	auto oneOverFactorials = std::make_shared<FullRow<double>>(10);
	for (size_t i = 0; i < oneOverFactorials->size(); i++)
	{
		oneOverFactorials->at(i) = 1.0 / std::tgamma(i + 1);
	}
	return oneOverFactorials;
}();

void DifferenceOperator::calcOperatorMatrix()
{
	//Compute operatorMatrix such that 
	//value(time) : = (operatorMatrix at : 1) timesColumn : series.
	//valuedot(time) : = (operatorMatrix at : 2) timesColumn : series.
	//valueddot(time) : = (operatorMatrix at : 3) timesColumn : series.

    formTaylorMatrix();
    // Taylor columns scale as powers of the time step. Equilibrate them before
    // pivoting instead of ignoring a singular-matrix exception (which can leave
    // a stale operator of the previous order). Rescale inverse rows afterwards.
    const auto n = taylorMatrix->nrow();
    auto scaled = std::make_shared<FullMatrix<double>>(n, n);
    std::vector<double> scales(n, 0.0);
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i)
            scales[j] = std::max(scales[j], std::abs(taylorMatrix->at(i)->at(j)));
        if (!(scales[j] > 0) || !std::isfinite(scales[j]))
            throw SingularMatrixError("Invalid time nodes in the integration operator");
        for (size_t i = 0; i < n; ++i)
            scaled->at(i)->at(j) = taylorMatrix->at(i)->at(j) / scales[j];
    }
    auto inverse = CREATE<LDUFullMatParPv>::With()->inversesaveOriginal(scaled, false);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            inverse->at(i)->at(j) /= scales[i];
    operatorMatrix = inverse;
}

void DifferenceOperator::initialize()
{
	//Do nothing
}

void MbD::DifferenceOperator::initializeLocally()
{
	throw SimulationStoppingError("To be implemented.");
}

void DifferenceOperator::setiStep(size_t i)
{
	iStep = i;
}

void DifferenceOperator::setorder(size_t o)
{
	order = o;
}

void DifferenceOperator::instantiateTaylorMatrix()
{
	if (taylorMatrix == nullptr || (taylorMatrix->nrow() != (order + 1))) {
		taylorMatrix = std::make_shared<FullMatrix<double>>(order + 1, order + 1);
	}
}

void DifferenceOperator::formTaylorRowwithTimeNodederivative(size_t i, size_t ii, size_t k)
{
	//| rowi hi hipower aij |
	auto& rowi = taylorMatrix->at(i);
	for (size_t j = 0; j < k; j++)
	{
		rowi->at(j) = 0.0;
	}
	rowi->at(k) = 1.0;
	auto hi = timeNodes->at(ii) - time;
	auto hipower = 1.0;
	for (size_t j = k + 1; j < order + 1; j++)
	{
		hipower = hipower * hi;
		auto aij = hipower * OneOverFactorials->at(j - k);
		rowi->atiput(j, aij);
	}
}

void DifferenceOperator::settime(double t)
{
	time = t;
}
