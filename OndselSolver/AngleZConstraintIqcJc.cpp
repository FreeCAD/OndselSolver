#include "AngleZConstraintIqcJc.h"
#include "AngleZIeqcJec.h"
#include "EndFrameqc.h"

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

MbD::AngleZConstraintIqcJc::AngleZConstraintIqcJc(EndFrmsptr frmi, EndFrmsptr frmj) : AngleZConstraintIJ(frmi, frmj)
{
	pGpEI = std::make_shared<FullRow<double>>(4);
	ppGpEIpEI = std::make_shared<FullMatrix<double>>(4, 4);
}

void MbD::AngleZConstraintIqcJc::initthezIeJe()
{
	thezIeJe = std::make_shared<AngleZIeqcJec>(frmI, frmJ);
}

void MbD::AngleZConstraintIqcJc::addToJointTorqueI(FColDsptr jointTorque)
{
	auto frmIeqc = std::static_pointer_cast<EndFrameqc>(frmI);
	auto rIpIeIp = frmIeqc->rpep();
	auto pAOIppEI = frmIeqc->pAOppE();
	auto aBOIp = frmIeqc->aBOp();
	auto fpAOIppEIrIpIeIp = std::make_shared<FullColumn<double>>(4, 0.0);
	auto lampGpE = pGpEI->transpose()->times(lam);
	auto c2Torque = aBOIp->timesFullColumn(lampGpE->minusFullColumn(fpAOIppEIrIpIeIp));
	jointTorque->equalSelfPlusFullColumntimes(c2Torque, 0.5);
}

void MbD::AngleZConstraintIqcJc::calc_pGpEI()
{
	pGpEI = thezIeJe->pvaluepEI();
}

void MbD::AngleZConstraintIqcJc::calc_ppGpEIpEI()
{
	ppGpEIpEI = thezIeJe->ppvaluepEIpEI();
}

void MbD::AngleZConstraintIqcJc::calcPostDynCorrectorIteration()
{
	AngleZConstraintIJ::calcPostDynCorrectorIteration();
	this->calc_pGpEI();
	this->calc_ppGpEIpEI();
}

void MbD::AngleZConstraintIqcJc::fillAccICIterError(FColDsptr col)
{
	col->atiplusFullVectortimes(iqEI, pGpEI, lam);
	auto efrmIqc = std::static_pointer_cast<EndFrameqc>(frmI);
	auto qXdotI = efrmIqc->qXdot();
	auto qEdotI = efrmIqc->qEdot();
	auto sum = pGpEI->timesFullColumn(efrmIqc->qEddot());
	sum += qEdotI->transposeTimesFullColumn(ppGpEIpEI->timesFullColumn(qEdotI));
	col->atiplusNumber(iG, sum);
}

void MbD::AngleZConstraintIqcJc::fillPosICError(FColDsptr col)
{
	AngleZConstraintIJ::fillPosICError(col);
	col->atiplusFullVectortimes(iqEI, pGpEI, lam);
}

void MbD::AngleZConstraintIqcJc::fillPosICJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqEI, pGpEI);
	mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
	mat->atijplusFullMatrixtimes(iqEI, iqEI, ppGpEIpEI, lam);
}

void MbD::AngleZConstraintIqcJc::fillPosKineJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqEI, pGpEI);
}

void MbD::AngleZConstraintIqcJc::fillVelICJacob(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqEI, pGpEI);
	mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
}

void MbD::AngleZConstraintIqcJc::fillpFpy(SpMatDsptr mat)
{
	mat->atijplusFullRow(iG, iqEI, pGpEI);
	mat->atijplusFullMatrixtimes(iqEI, iqEI, ppGpEIpEI, lam);
}

void MbD::AngleZConstraintIqcJc::fillpFpydot(SpMatDsptr mat)
{
	mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
}

void MbD::AngleZConstraintIqcJc::useEquationNumbers()
{
	auto frmIeqc = std::static_pointer_cast<EndFrameqc>(frmI);
	iqEI = frmIeqc->iqE();
}

double MbD::AngleZConstraintIqcJc::constraintVelocity() const
{
	auto frameI = std::static_pointer_cast<EndFrameqc>(frmI);
	return pGpEI->timesFullColumn(frameI->qEdot());
}

void MbD::AngleZConstraintIqcJc::fillGeneralizedForce(FColDsptr col, double multiplier)
{
	col->atiplusFullVectortimes(iqEI, pGpEI, multiplier);
}

void MbD::AngleZConstraintIqcJc::fillGeneralizedForcePositionJacobian(
	SpMatDsptr mat,
	double multiplier,
	double derivative
)
{
	mat->atijplusFullMatrixtimes(iqEI, iqEI, ppGpEIpEI, multiplier);
	addOuterProduct(mat, iqEI, pGpEI, iqEI, pGpEI, derivative);
}

void MbD::AngleZConstraintIqcJc::fillGeneralizedForceVelocityJacobian(
	SpMatDsptr mat,
	double derivative
)
{
	addOuterProduct(mat, iqEI, pGpEI, iqEI, pGpEI, derivative);
}
