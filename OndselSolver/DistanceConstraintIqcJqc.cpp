/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "DistanceConstraintIqcJqc.h"
#include "EndFrameqc.h"
#include "CREATE.h"
#include "DistIeqcJeqc.h"

using namespace MbD;

namespace
{
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
            mat->atijplusNumber(
                rowStart + i,
                columnStart + j,
                factor * row->at(i) * column->at(j)
            );
        }
    }
}
}  // namespace

MbD::DistanceConstraintIqcJqc::DistanceConstraintIqcJqc(EndFrmsptr frmi, EndFrmsptr frmj) : DistanceConstraintIqcJc(frmi, frmj)
{
}

void MbD::DistanceConstraintIqcJqc::calcPostDynCorrectorIteration()
{
	DistanceConstraintIqcJc::calcPostDynCorrectorIteration();
	pGpXJ = distIeJe->pvaluepXJ();
	pGpEJ = distIeJe->pvaluepEJ();
	ppGpXIpXJ = distIeJe->ppvaluepXIpXJ();
	ppGpEIpXJ = distIeJe->ppvaluepEIpXJ();
	ppGpXJpXJ = distIeJe->ppvaluepXJpXJ();
	ppGpXIpEJ = distIeJe->ppvaluepXIpEJ();
	ppGpEIpEJ = distIeJe->ppvaluepEIpEJ();
	ppGpXJpEJ = distIeJe->ppvaluepXJpEJ();
	ppGpEJpEJ = distIeJe->ppvaluepEJpEJ();
}

void MbD::DistanceConstraintIqcJqc::fillAccICIterError(FColDsptr col)
{
	DistanceConstraintIqcJc::fillAccICIterError(col);
	col->atiplusFullVectortimes(iqXJ, pGpXJ, lam);
	col->atiplusFullVectortimes(iqEJ, pGpEJ, lam);
	auto frmIeqc = std::static_pointer_cast<EndFrameqc>(frmI);
	auto frmJeqc = std::static_pointer_cast<EndFrameqc>(frmJ);
	auto qXdotI = frmIeqc->qXdot();
	auto qEdotI = frmIeqc->qEdot();
	auto qXdotJ = frmJeqc->qXdot();
	auto qEdotJ = frmJeqc->qEdot();
	double sum = 0.0;
	sum += pGpXJ->timesFullColumn(frmJeqc->qXddot());
	sum += pGpEJ->timesFullColumn(frmJeqc->qEddot());
	sum += 2.0 * (qXdotI->transposeTimesFullColumn(ppGpXIpXJ->timesFullColumn(qXdotJ)));
	sum += 2.0 * (qEdotI->transposeTimesFullColumn(ppGpEIpXJ->timesFullColumn(qXdotJ)));
	sum += qXdotJ->transposeTimesFullColumn(ppGpXJpXJ->timesFullColumn(qXdotJ));
	sum += 2.0 * (qXdotI->transposeTimesFullColumn(ppGpXIpEJ->timesFullColumn(qEdotJ)));
	sum += 2.0 * (qEdotI->transposeTimesFullColumn(ppGpEIpEJ->timesFullColumn(qEdotJ)));
	sum += 2.0 * (qXdotJ->transposeTimesFullColumn(ppGpXJpEJ->timesFullColumn(qEdotJ)));
	sum += qEdotJ->transposeTimesFullColumn(ppGpEJpEJ->timesFullColumn(qEdotJ));
	col->atiplusNumber(iG, sum);
}

void MbD::DistanceConstraintIqcJqc::fillPosICError(FColDsptr col)
{
	DistanceConstraintIqcJc::fillPosICError(col);
	col->atiplusFullVectortimes(iqXJ, pGpXJ, lam);
	col->atiplusFullVectortimes(iqEJ, pGpEJ, lam);
}

void MbD::DistanceConstraintIqcJqc::fillPosICJacob(SpMatDsptr mat)
{
	DistanceConstraintIqcJc::fillPosICJacob(mat);
	mat->atijplusFullRow(iG, iqXJ, pGpXJ);
	mat->atijplusFullColumn(iqXJ, iG, pGpXJ->transpose());
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
	mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
	auto ppGpXIpXJlam = ppGpXIpXJ->times(lam);
	mat->atijplusFullMatrix(iqXI, iqXJ, ppGpXIpXJlam);
	mat->atijplusTransposeFullMatrix(iqXJ, iqXI, ppGpXIpXJlam);
	auto ppGpEIpXJlam = ppGpEIpXJ->times(lam);
	mat->atijplusFullMatrix(iqEI, iqXJ, ppGpEIpXJlam);
	mat->atijplusTransposeFullMatrix(iqXJ, iqEI, ppGpEIpXJlam);
	mat->atijplusFullMatrixtimes(iqXJ, iqXJ, ppGpXJpXJ, lam);
	auto ppGpXIpEJlam = ppGpXIpEJ->times(lam);
	mat->atijplusFullMatrix(iqXI, iqEJ, ppGpXIpEJlam);
	mat->atijplusTransposeFullMatrix(iqEJ, iqXI, ppGpXIpEJlam);
	auto ppGpEIpEJlam = ppGpEIpEJ->times(lam);
	mat->atijplusFullMatrix(iqEI, iqEJ, ppGpEIpEJlam);
	mat->atijplusTransposeFullMatrix(iqEJ, iqEI, ppGpEIpEJlam);
	auto ppGpXJpEJlam = ppGpXJpEJ->times(lam);
	mat->atijplusFullMatrix(iqXJ, iqEJ, ppGpXJpEJlam);
	mat->atijplusTransposeFullMatrix(iqEJ, iqXJ, ppGpXJpEJlam);
	mat->atijplusFullMatrixtimes(iqEJ, iqEJ, ppGpEJpEJ, lam);
}

void MbD::DistanceConstraintIqcJqc::fillPosKineJacob(SpMatDsptr mat)
{
	DistanceConstraintIqcJc::fillPosKineJacob(mat);
	mat->atijplusFullRow(iG, iqXJ, pGpXJ);
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
}

void MbD::DistanceConstraintIqcJqc::fillVelICJacob(SpMatDsptr mat)
{
	DistanceConstraintIqcJc::fillVelICJacob(mat);
	mat->atijplusFullRow(iG, iqXJ, pGpXJ);
	mat->atijplusFullColumn(iqXJ, iG, pGpXJ->transpose());
	mat->atijplusFullRow(iG, iqEJ, pGpEJ);
	mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
}

void MbD::DistanceConstraintIqcJqc::init_distIeJe()
{
	distIeJe = CREATE<DistIeqcJeqc>::With(frmI, frmJ);
}

void MbD::DistanceConstraintIqcJqc::useEquationNumbers()
{
	DistanceConstraintIqcJc::useEquationNumbers();
	auto frmJeqc = std::static_pointer_cast<EndFrameqc>(frmJ);
	iqXJ = frmJeqc->iqX();
	iqEJ = frmJeqc->iqE();
}

std::string MbD::DistanceConstraintIqcJqc::constraintSpec()
{
	return "DistanceConstraintIJ";
}

void DistanceConstraintIqcJqc::fillpFpy(SpMatDsptr mat)
{
    DistanceConstraintIqcJc::fillpFpy(mat);
    mat->atijplusFullRow(iG, iqXJ, pGpXJ);
    mat->atijplusFullRow(iG, iqEJ, pGpEJ);
    auto ppGpXIpXJlam = ppGpXIpXJ->times(lam);
    mat->atijplusFullMatrix(iqXI, iqXJ, ppGpXIpXJlam);
    mat->atijplusTransposeFullMatrix(iqXJ, iqXI, ppGpXIpXJlam);
    auto ppGpEIpXJlam = ppGpEIpXJ->times(lam);
    mat->atijplusFullMatrix(iqEI, iqXJ, ppGpEIpXJlam);
    mat->atijplusTransposeFullMatrix(iqXJ, iqEI, ppGpEIpXJlam);
    mat->atijplusFullMatrixtimes(iqXJ, iqXJ, ppGpXJpXJ, lam);
    auto ppGpXIpEJlam = ppGpXIpEJ->times(lam);
    mat->atijplusFullMatrix(iqXI, iqEJ, ppGpXIpEJlam);
    mat->atijplusTransposeFullMatrix(iqEJ, iqXI, ppGpXIpEJlam);
    auto ppGpEIpEJlam = ppGpEIpEJ->times(lam);
    mat->atijplusFullMatrix(iqEI, iqEJ, ppGpEIpEJlam);
    mat->atijplusTransposeFullMatrix(iqEJ, iqEI, ppGpEIpEJlam);
    auto ppGpXJpEJlam = ppGpXJpEJ->times(lam);
    mat->atijplusFullMatrix(iqXJ, iqEJ, ppGpXJpEJlam);
    mat->atijplusTransposeFullMatrix(iqEJ, iqXJ, ppGpXJpEJlam);
    mat->atijplusFullMatrixtimes(iqEJ, iqEJ, ppGpEJpEJ, lam);
}

void DistanceConstraintIqcJqc::fillpFpydot(SpMatDsptr mat)
{
    DistanceConstraintIqcJc::fillpFpydot(mat);
    mat->atijplusFullColumn(iqXJ, iG, pGpXJ->transpose());
    mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
}

double DistanceConstraintIqcJqc::constraintVelocity() const
{
    auto frameJ = std::static_pointer_cast<EndFrameqc>(frmJ);
    return DistanceConstraintIqcJc::constraintVelocity()
        + pGpXJ->timesFullColumn(frameJ->qXdot())
        + pGpEJ->timesFullColumn(frameJ->qEdot());
}

void DistanceConstraintIqcJqc::fillGeneralizedForce(FColDsptr col, double multiplier)
{
    DistanceConstraintIqcJc::fillGeneralizedForce(col, multiplier);
    col->atiplusFullVectortimes(iqXJ, pGpXJ, multiplier);
    col->atiplusFullVectortimes(iqEJ, pGpEJ, multiplier);
}

void DistanceConstraintIqcJqc::fillGeneralizedForcePositionJacobian(
    SpMatDsptr mat,
    double multiplier,
    double derivative
)
{
    DistanceConstraintIqcJc::fillGeneralizedForcePositionJacobian(
        mat, multiplier, derivative
    );
    const struct Hessian
    {
        size_t row;
        size_t column;
        FMatDsptr value;
    } crossTerms[] = {
        {iqXI, iqXJ, ppGpXIpXJ},
        {iqEI, iqXJ, ppGpEIpXJ},
        {iqXI, iqEJ, ppGpXIpEJ},
        {iqEI, iqEJ, ppGpEIpEJ},
        {iqXJ, iqEJ, ppGpXJpEJ},
    };
    for (const auto& term : crossTerms) {
        auto hessian = term.value->times(multiplier);
        mat->atijplusFullMatrix(term.row, term.column, hessian);
        mat->atijplusTransposeFullMatrix(term.column, term.row, hessian);
    }
    mat->atijplusFullMatrixtimes(iqXJ, iqXJ, ppGpXJpXJ, multiplier);
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
            addOuterProduct(
                mat,
                row.first,
                row.second,
                column.first,
                column.second,
                derivative
            );
        }
    }
}

void DistanceConstraintIqcJqc::fillGeneralizedForceVelocityJacobian(
    SpMatDsptr mat,
    double derivative
)
{
    DistanceConstraintIqcJc::fillGeneralizedForceVelocityJacobian(mat, derivative);
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
            addOuterProduct(
                mat,
                row.first,
                row.second,
                column.first,
                column.second,
                derivative
            );
        }
    }
}
