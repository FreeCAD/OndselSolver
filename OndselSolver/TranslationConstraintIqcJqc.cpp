/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "TranslationConstraintIqcJqc.h"
#include "DispCompIeqcJeqcKeqc.h"
#include "EndFrameqc.h"
#include "CREATE.h"

using namespace MbD;

namespace {
void addOuterProduct(
	SpMatDsptr mat,
	size_t rowStart,
	const FRowDsptr& row,
	size_t columnStart,
	const FRowDsptr& column,
	double factor
)
{
	for (size_t i = 0; i < row->size(); ++i) {
		for (size_t j = 0; j < column->size(); ++j) {
			mat->atijplusNumber(rowStart + i, columnStart + j, factor * row->at(i) * column->at(j));
		}
	}
}
}

TranslationConstraintIqcJqc::TranslationConstraintIqcJqc(EndFrmsptr frmi, EndFrmsptr frmj, size_t axisi) :
	TranslationConstraintIqcJc(frmi, frmj, axisi)
{
}

void TranslationConstraintIqcJqc::initriIeJeIe()
{
	riIeJeIe = CREATE<DispCompIeqcJeqcKeqc>::With(frmI, frmJ, frmI, axisI);
}

void TranslationConstraintIqcJqc::calcPostDynCorrectorIteration()
{
	TranslationConstraintIqcJc::calcPostDynCorrectorIteration();
	pGpXJ = riIeJeIe->pvaluepXJ();
	pGpEJ = riIeJeIe->pvaluepEJ();
	ppGpEIpXJ = riIeJeIe->ppvaluepXJpEK()->transpose();
	ppGpEIpEJ = riIeJeIe->ppvaluepEJpEK()->transpose();
	ppGpEJpEJ = riIeJeIe->ppvaluepEJpEJ();
}

void TranslationConstraintIqcJqc::useEquationNumbers()
{
	TranslationConstraintIqcJc::useEquationNumbers();
	auto frmJeqc = std::static_pointer_cast<EndFrameqc>(frmJ);
	iqXJ = frmJeqc->iqX();
	iqEJ = frmJeqc->iqE();
}

std::string MbD::TranslationConstraintIqcJqc::constraintSpec()
{
	return "TranslationConstraintI" + MbDMath::xyzFromInt(axisI) + "J";
}

void TranslationConstraintIqcJqc::fillPosICError(FColDsptr col)
{
	TranslationConstraintIqcJc::fillPosICError(col);
	col->atiplusFullVectortimes(iqXJ, pGpXJ, lam);
	col->atiplusFullVectortimes(iqEJ, pGpEJ, lam);
}

void TranslationConstraintIqcJqc::fillPosICJacob(SpMatDsptr mat)
{
	TranslationConstraintIqcJc::fillPosICJacob(mat);
	mat->atijplusFullRow(iG, iqXJ, pGpXJ);
	mat->atijplusFullColumn(iqXJ, iG, pGpXJ->transpose());
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
	mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
	auto ppGpEIpXJlam = ppGpEIpXJ->times(lam);
	mat->atijplusFullMatrix(iqEI, iqXJ, ppGpEIpXJlam);
	mat->atijplusTransposeFullMatrix(iqXJ, iqEI, ppGpEIpXJlam);
	auto ppGpEIpEJlam = ppGpEIpEJ->times(lam);
	mat->atijplusFullMatrix(iqEI, iqEJ, ppGpEIpEJlam);
	mat->atijplusTransposeFullMatrix(iqEJ, iqEI, ppGpEIpEJlam);
	mat->atijplusFullMatrixtimes(iqEJ, iqEJ, ppGpEJpEJ, lam);
}

void TranslationConstraintIqcJqc::fillPosKineJacob(SpMatDsptr mat)
{
	TranslationConstraintIqcJc::fillPosKineJacob(mat);
	mat->atijplusFullRow(iG, iqXJ, pGpXJ);
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
}

void TranslationConstraintIqcJqc::fillVelICJacob(SpMatDsptr mat)
{
	TranslationConstraintIqcJc::fillVelICJacob(mat);
	mat->atijplusFullRow(iG, iqXJ, pGpXJ);
	mat->atijplusFullColumn(iqXJ, iG, pGpXJ->transpose());
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
	mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
}

void TranslationConstraintIqcJqc::fillAccICIterError(FColDsptr col)
{
	TranslationConstraintIqcJc::fillAccICIterError(col);
	col->atiplusFullVectortimes(iqXJ, pGpXJ, lam);
	col->atiplusFullVectortimes(iqEJ, pGpEJ, lam);
	auto efrmIqc = std::static_pointer_cast<EndFrameqc>(frmI);
	auto efrmJqc = std::static_pointer_cast<EndFrameqc>(frmJ);
	auto qEdotI = efrmIqc->qEdot();
	auto qXdotJ = efrmJqc->qXdot();
	auto qEdotJ = efrmJqc->qEdot();
	double sum = pGpXJ->timesFullColumn(efrmJqc->qXddot());
	sum += pGpEJ->timesFullColumn(efrmJqc->qEddot());
	sum += 2.0 * (qEdotI->transposeTimesFullColumn(ppGpEIpXJ->timesFullColumn(qXdotJ)));
	sum += 2.0 * (qEdotI->transposeTimesFullColumn(ppGpEIpEJ->timesFullColumn(qEdotJ)));
	sum += qEdotJ->transposeTimesFullColumn(ppGpEJpEJ->timesFullColumn(qEdotJ));
	col->atiplusNumber(iG, sum);
}

void TranslationConstraintIqcJqc::fillpFpy(SpMatDsptr mat)
{
    TranslationConstraintIqcJc::fillpFpy(mat);
    mat->atijplusFullRow(iG, iqXJ, pGpXJ);
    mat->atijplusFullRow(iG, iqEJ, pGpEJ);
    auto ppGpEIpXJlam = ppGpEIpXJ->times(lam);
    mat->atijplusFullMatrix(iqEI, iqXJ, ppGpEIpXJlam);
    mat->atijplusTransposeFullMatrix(iqXJ, iqEI, ppGpEIpXJlam);
    auto ppGpEIpEJlam = ppGpEIpEJ->times(lam);
    mat->atijplusFullMatrix(iqEI, iqEJ, ppGpEIpEJlam);
    mat->atijplusTransposeFullMatrix(iqEJ, iqEI, ppGpEIpEJlam);
    mat->atijplusFullMatrixtimes(iqEJ, iqEJ, ppGpEJpEJ, lam);
}

void TranslationConstraintIqcJqc::fillpFpydot(SpMatDsptr mat)
{
    TranslationConstraintIqcJc::fillpFpydot(mat);
    mat->atijplusFullColumn(iqXJ, iG, pGpXJ->transpose());
    mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
}

double TranslationConstraintIqcJqc::constraintVelocity() const
{
	auto frameJ = std::static_pointer_cast<EndFrameqc>(frmJ);
	return TranslationConstraintIqcJc::constraintVelocity()
		+ pGpXJ->timesFullColumn(frameJ->qXdot())
		+ pGpEJ->timesFullColumn(frameJ->qEdot());
}

void TranslationConstraintIqcJqc::fillGeneralizedForce(FColDsptr col, double multiplier)
{
	TranslationConstraintIqcJc::fillGeneralizedForce(col, multiplier);
	col->atiplusFullVectortimes(iqXJ, pGpXJ, multiplier);
	col->atiplusFullVectortimes(iqEJ, pGpEJ, multiplier);
}

void TranslationConstraintIqcJqc::fillGeneralizedForcePositionJacobian(
	SpMatDsptr mat,
	double multiplier,
	double derivative
)
{
	TranslationConstraintIqcJc::fillGeneralizedForcePositionJacobian(mat, multiplier, derivative);
	auto eixjHessian = ppGpEIpXJ->times(multiplier);
	mat->atijplusFullMatrix(iqEI, iqXJ, eixjHessian);
	mat->atijplusTransposeFullMatrix(iqXJ, iqEI, eixjHessian);
	auto eiejHessian = ppGpEIpEJ->times(multiplier);
	mat->atijplusFullMatrix(iqEI, iqEJ, eiejHessian);
	mat->atijplusTransposeFullMatrix(iqEJ, iqEI, eiejHessian);
	mat->atijplusFullMatrixtimes(iqEJ, iqEJ, ppGpEJpEJ, multiplier);

	const std::pair<size_t, FRowDsptr> iSegments[] = {{iqXI, pGpXI}, {iqEI, pGpEI}};
	const std::pair<size_t, FRowDsptr> jSegments[] = {{iqXJ, pGpXJ}, {iqEJ, pGpEJ}};
	for (const auto& i : iSegments) {
		for (const auto& j : jSegments) {
			addOuterProduct(mat, i.first, i.second, j.first, j.second, derivative);
			addOuterProduct(mat, j.first, j.second, i.first, i.second, derivative);
		}
	}
	for (const auto& row : jSegments) {
		for (const auto& column : jSegments) {
			addOuterProduct(mat, row.first, row.second, column.first, column.second, derivative);
		}
	}
}

void TranslationConstraintIqcJqc::fillGeneralizedForceVelocityJacobian(
	SpMatDsptr mat,
	double derivative
)
{
	TranslationConstraintIqcJc::fillGeneralizedForceVelocityJacobian(mat, derivative);
	const std::pair<size_t, FRowDsptr> iSegments[] = {{iqXI, pGpXI}, {iqEI, pGpEI}};
	const std::pair<size_t, FRowDsptr> jSegments[] = {{iqXJ, pGpXJ}, {iqEJ, pGpEJ}};
	for (const auto& i : iSegments) {
		for (const auto& j : jSegments) {
			addOuterProduct(mat, i.first, i.second, j.first, j.second, derivative);
			addOuterProduct(mat, j.first, j.second, i.first, i.second, derivative);
		}
	}
	for (const auto& row : jSegments) {
		for (const auto& column : jSegments) {
			addOuterProduct(mat, row.first, row.second, column.first, column.second, derivative);
		}
	}
}
