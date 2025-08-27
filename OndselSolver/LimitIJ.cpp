#include "LimitIJ.h"
#include "Constraint.h"

using namespace MbD;

bool MbD::LimitIJ::satisfied() const
{
	auto& constraint = constraints->front();
	if (type == "=<") {
		return constraint->aG < tol;
	}
	else if (type == "=>") {
		return constraint->aG > -tol;
	}
	assert(false);
	return true;
}

void MbD::LimitIJ::deactivate()
{
	active = false;
}

void MbD::LimitIJ::activate()
{
	active = true;
}

void MbD::LimitIJ::fillConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> allConstraints)
{
	if (active) {
		ConstraintSet::fillConstraints(allConstraints);
	}
}

void MbD::LimitIJ::fillPosICError(FColDsptr col)
{
	if (active) {
		ConstraintSet::fillPosICError(col);
	}
}

void MbD::LimitIJ::fillPosICJacob(SpMatDsptr mat)
{
	if (active) {
		ConstraintSet::fillPosICJacob(mat);
	}
}

void MbD::LimitIJ::fillqsulam(FColDsptr col)
{
	if (active) {
		ConstraintSet::fillqsulam(col);
	}
}

void MbD::LimitIJ::setqsulam(FColDsptr col)
{
	if (active) {
		ConstraintSet::setqsulam(col);
	}
}

void MbD::LimitIJ::useEquationNumbers()
{
	if (active) {
		ConstraintSet::useEquationNumbers();
	}
}
