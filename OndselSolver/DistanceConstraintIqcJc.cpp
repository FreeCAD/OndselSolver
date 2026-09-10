/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "DistanceConstraintIqcJc.h"
#include "EndFrameqc.h"
#include "CREATE.h"
#include "DistIeqcJec.h"

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

MbD::DistanceConstraintIqcJc::DistanceConstraintIqcJc(EndFrmsptr frmi, EndFrmsptr frmj) : DistanceConstraintIJ(frmi, frmj)
{
}

void MbD::DistanceConstraintIqcJc::addToJointForceI(FColDsptr col)
{
	col->equalSelfPlusFullVectortimes(pGpXI, lam);
}

void MbD::DistanceConstraintIqcJc::addToJointTorqueI(FColDsptr jointTorque)
{
	auto cForceT = pGpXI->times(lam);
	auto frmIeqc = std::static_pointer_cast<EndFrameqc>(frmI);
	auto rIpIeIp = frmIeqc->rpep();
	auto pAOIppEI = frmIeqc->pAOppE();
	auto aBOIp = frmIeqc->aBOp();
	auto fpAOIppEIrIpIeIp = std::make_shared<FullColumn<double>>(4, 0.0);
	for (size_t i = 0; i < 4; i++)
	{
		auto dum = cForceT->timesFullColumn(pAOIppEI->at(i)->timesFullColumn(rIpIeIp));
		fpAOIppEIrIpIeIp->atiput(i, dum);
	}
	auto lampGpE = pGpEI->transpose()->times(lam);
	auto c2Torque = aBOIp->timesFullColumn(lampGpE->minusFullColumn(fpAOIppEIrIpIeIp));
	jointTorque->equalSelfPlusFullColumntimes(c2Torque, 0.5);
}

void MbD::DistanceConstraintIqcJc::calcPostDynCorrectorIteration()
{
	DistanceConstraintIJ::calcPostDynCorrectorIteration();
	pGpXI = distIeJe->pvaluepXI();
	pGpEI = distIeJe->pvaluepEI();
	ppGpXIpXI = distIeJe->ppvaluepXIpXI();
	ppGpXIpEI = distIeJe->ppvaluepXIpEI();
	ppGpEIpEI = distIeJe->ppvaluepEIpEI();
}

void MbD::DistanceConstraintIqcJc::fillAccICIterError(FColDsptr col)
{
	col->atiplusFullVectortimes(iqXI, pGpXI, lam);
	col->atiplusFullVectortimes(iqEI, pGpEI, lam);
	auto efrmIqc = std::static_pointer_cast<EndFrameqc>(frmI);
	auto qXdotI = efrmIqc->qXdot();
	auto qEdotI = efrmIqc->qEdot();
	auto sum = pGpXI->timesFullColumn(efrmIqc->qXddot());
	sum += pGpEI->timesFullColumn(efrmIqc->qEddot());
	sum += qXdotI->transposeTimesFullColumn(ppGpXIpXI->timesFullColumn(qXdotI));
	sum += 2.0 * (qXdotI->transposeTimesFullColumn(ppGpXIpEI->timesFullColumn(qEdotI)));
	sum += qEdotI->transposeTimesFullColumn(ppGpEIpEI->timesFullColumn(qEdotI));
	col->atiplusNumber(iG, sum);
}

void MbD::DistanceConstraintIqcJc::fillPosICError(FColDsptr col)
{
	DistanceConstraintIJ::fillPosICError(col);
	col->atiplusFullVectortimes(iqXI, pGpXI, lam);
	col->atiplusFullVectortimes(iqEI, pGpEI, lam);
}

void MbD::DistanceConstraintIqcJc::fillPosICJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqXI, pGpXI);
	mat->atijplusFullColumn(iqXI, iG, pGpXI->transpose());
	mat->atijplusFullRow(iG, iqEI, pGpEI);
	mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
	mat->atijplusFullMatrixtimes(iqXI, iqXI, ppGpXIpXI, lam);
	auto ppGpXIpEIlam = ppGpXIpEI->times(lam);
	mat->atijplusFullMatrix(iqXI, iqEI, ppGpXIpEIlam);
	mat->atijplusTransposeFullMatrix(iqEI, iqXI, ppGpXIpEIlam);
	mat->atijplusFullMatrixtimes(iqEI, iqEI, ppGpEIpEI, lam);
}

void MbD::DistanceConstraintIqcJc::fillPosKineJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqXI, pGpXI);
	mat->atijplusFullRow(iG, iqEI, pGpEI);
}

void MbD::DistanceConstraintIqcJc::fillVelICJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqXI, pGpXI);
	mat->atijplusFullColumn(iqXI, iG, pGpXI->transpose());
	mat->atijplusFullRow(iG, iqEI, pGpEI);
	mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
}

void MbD::DistanceConstraintIqcJc::init_distIeJe()
{
	distIeJe = CREATE<DistIeqcJec>::With(frmI, frmJ);
}

void MbD::DistanceConstraintIqcJc::useEquationNumbers()
{
	auto frmIeqc = std::static_pointer_cast<EndFrameqc>(frmI);
	iqXI = frmIeqc->iqX();
	iqEI = frmIeqc->iqE();
}

void DistanceConstraintIqcJc::fillpFpy(SpMatDsptr mat)
{
    mat->atijplusFullRow(iG, iqXI, pGpXI);
    mat->atijplusFullRow(iG, iqEI, pGpEI);
    mat->atijplusFullMatrixtimes(iqXI, iqXI, ppGpXIpXI, lam);
    auto ppGpXIpEIlam = ppGpXIpEI->times(lam);
    mat->atijplusFullMatrix(iqXI, iqEI, ppGpXIpEIlam);
    mat->atijplusTransposeFullMatrix(iqEI, iqXI, ppGpXIpEIlam);
    mat->atijplusFullMatrixtimes(iqEI, iqEI, ppGpEIpEI, lam);
}

void DistanceConstraintIqcJc::fillpFpydot(SpMatDsptr mat)
{
    mat->atijplusFullColumn(iqXI, iG, pGpXI->transpose());
    mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
}

double DistanceConstraintIqcJc::constraintVelocity() const
{
    auto frameI = std::static_pointer_cast<EndFrameqc>(frmI);
    return pGpXI->timesFullColumn(frameI->qXdot())
        + pGpEI->timesFullColumn(frameI->qEdot());
}

void DistanceConstraintIqcJc::fillGeneralizedForce(FColDsptr col, double multiplier)
{
    col->atiplusFullVectortimes(iqXI, pGpXI, multiplier);
    col->atiplusFullVectortimes(iqEI, pGpEI, multiplier);
}

void DistanceConstraintIqcJc::fillGeneralizedForcePositionJacobian(
    SpMatDsptr mat,
    double multiplier,
    double derivative
)
{
    mat->atijplusFullMatrixtimes(iqXI, iqXI, ppGpXIpXI, multiplier);
    auto crossHessian = ppGpXIpEI->times(multiplier);
    mat->atijplusFullMatrix(iqXI, iqEI, crossHessian);
    mat->atijplusTransposeFullMatrix(iqEI, iqXI, crossHessian);
    mat->atijplusFullMatrixtimes(iqEI, iqEI, ppGpEIpEI, multiplier);

    const std::pair<size_t, FRowDsptr> segments[] = {{iqXI, pGpXI}, {iqEI, pGpEI}};
    for (const auto& row : segments) {
        for (const auto& column : segments) {
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

void DistanceConstraintIqcJc::fillGeneralizedForceVelocityJacobian(
    SpMatDsptr mat,
    double derivative
)
{
    const std::pair<size_t, FRowDsptr> segments[] = {{iqXI, pGpXI}, {iqEI, pGpEI}};
    for (const auto& row : segments) {
        for (const auto& column : segments) {
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
