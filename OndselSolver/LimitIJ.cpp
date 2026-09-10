#include "LimitIJ.h"
#include "Constraint.h"

#include <algorithm>
#include <cmath>

using namespace MbD;

namespace {
struct ContactForce
{
	double multiplier = 0.0;
	double positionDerivative = 0.0;
	double velocityDerivative = 0.0;
};

ContactForce contactForce(const MbD::LimitIJ& limit)
{
	const auto& constraint = limit.constraints->front();
	const double normal = limit.type == "=>" ? 1.0 : limit.type == "=<" ? -1.0 : 0.0;
	if (normal == 0.0) {
		throw MbD::SimulationStoppingError("Unknown joint-limit comparison type");
	}

	const double penetration = -normal * constraint->aG;
	if (penetration <= 0.0) return {};

	const double penetrationVelocity = -normal * constraint->constraintVelocity();
	const double magnitude = limit.stiffness * penetration + limit.damping * penetrationVelocity;
	if (magnitude <= 0.0) return {};

	return {
		normal * magnitude,
		-limit.stiffness,
		-limit.damping,
	};
}
}

bool MbD::LimitIJ::satisfied() const
{
	auto& constraint = constraints->front();
	if (type == "=<") {
		return constraint->aG < tol;
	}
	else if (type == "=>") {
		return constraint->aG > -tol;
	}
	throw SimulationStoppingError("To be implemented.");
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

void MbD::LimitIJ::setCompliant(double newStiffness, double newDamping)
{
	if (!(newStiffness > 0.0) || newDamping < 0.0
		|| !std::isfinite(newStiffness) || !std::isfinite(newDamping)) {
		throw SimulationStoppingError("A compliant joint limit needs positive stiffness and nonnegative damping.");
	}
	compliant = true;
	active = false;
	stiffness = newStiffness;
	damping = newDamping;
}

MbD::LimitIJ::EnergyState MbD::LimitIJ::energyState() const
{
	if (!compliant) return {};
	const auto& constraint = constraints->front();
	const double normal = type == "=>" ? 1.0 : type == "=<" ? -1.0 : 0.0;
	if (normal == 0.0) throw SimulationStoppingError("Unknown joint-limit comparison type");
	const double penetration = -normal * constraint->aG;
	if (penetration <= 0.0) return {};
	const double penetrationVelocity = -normal * constraint->constraintVelocity();
	const double magnitude = stiffness * penetration + damping * penetrationVelocity;
	if (magnitude <= 0.0) return {0.0, 0.5 * stiffness * penetration * penetration, 0.0};
	const double dissipated = damping * penetrationVelocity * penetrationVelocity;
	return {
		-stiffness * penetration * penetrationVelocity - dissipated,
		0.5 * stiffness * penetration * penetration,
		dissipated
	};
}

void MbD::LimitIJ::fillConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> allConstraints)
{
	if (active) {
		ConstraintSet::fillConstraints(allConstraints);
	}
}

void MbD::LimitIJ::fillDispConstraints(
	std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> allConstraints
)
{
	if (active) ConstraintSet::fillDispConstraints(allConstraints);
}

void MbD::LimitIJ::fillEssenConstraints(
	std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> allConstraints
)
{
	if (active) ConstraintSet::fillEssenConstraints(allConstraints);
}

void MbD::LimitIJ::fillPerpenConstraints(
	std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> allConstraints
)
{
	if (active) ConstraintSet::fillPerpenConstraints(allConstraints);
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

void MbD::LimitIJ::fillqsudot(FColDsptr col)
{
	if (active) ConstraintSet::fillqsudot(col);
}

void MbD::LimitIJ::fillqsuddotlam(FColDsptr col)
{
	if (active) ConstraintSet::fillqsuddotlam(col);
}

void MbD::LimitIJ::fillVelICError(FColDsptr col)
{
	if (active) ConstraintSet::fillVelICError(col);
}

void MbD::LimitIJ::fillVelICJacob(SpMatDsptr mat)
{
	if (active) ConstraintSet::fillVelICJacob(mat);
}

void MbD::LimitIJ::fillAccICIterError(FColDsptr col)
{
	if (compliant) {
		auto contact = contactForce(*this);
		constraints->front()->fillGeneralizedForce(col, contact.multiplier);
	}
	else if (active) ConstraintSet::fillAccICIterError(col);
}

void MbD::LimitIJ::fillAccICIterJacob(SpMatDsptr mat)
{
	if (active) ConstraintSet::fillAccICIterJacob(mat);
}

void MbD::LimitIJ::fillpqsumu(FColDsptr col)
{
	if (active) ConstraintSet::fillpqsumu(col);
}

void MbD::LimitIJ::fillpqsumudot(FColDsptr col)
{
	if (active) ConstraintSet::fillpqsumudot(col);
}

void MbD::LimitIJ::fillDynError(FColDsptr col)
{
	if (compliant) {
		auto contact = contactForce(*this);
		constraints->front()->fillGeneralizedForce(col, contact.multiplier);
	}
	else if (active) ConstraintSet::fillDynError(col);
}

void MbD::LimitIJ::fillpFpy(SpMatDsptr mat)
{
	if (compliant) {
		auto contact = contactForce(*this);
		constraints->front()->fillGeneralizedForcePositionJacobian(
			mat,
			contact.multiplier,
			contact.positionDerivative
		);
	}
	else if (active) ConstraintSet::fillpFpy(mat);
}

void MbD::LimitIJ::fillpFpydot(SpMatDsptr mat)
{
	if (compliant) {
		auto contact = contactForce(*this);
		constraints->front()->fillGeneralizedForceVelocityJacobian(
			mat,
			contact.velocityDerivative
		);
	}
	else if (active) ConstraintSet::fillpFpydot(mat);
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

void MbD::LimitIJ::setqsudotlam(FColDsptr col)
{
	if (active) ConstraintSet::setqsudotlam(col);
}

void MbD::LimitIJ::setqsuddotlam(FColDsptr col)
{
	if (active) ConstraintSet::setqsuddotlam(col);
}

void MbD::LimitIJ::setpqsumu(FColDsptr col)
{
	if (active) ConstraintSet::setpqsumu(col);
}

void MbD::LimitIJ::setpqsumudot(FColDsptr col)
{
	if (active) ConstraintSet::setpqsumudot(col);
}

void MbD::LimitIJ::useEquationNumbers()
{
	if (active || compliant) {
		ConstraintSet::useEquationNumbers();
	}
}

namespace {
double clearance(const MbD::LimitIJ& limit)
{
	const double gap = limit.constraints->front()->aG;
	if (limit.type == "=<") return -gap;
	if (limit.type == "=>") return gap;
	throw MbD::SimulationStoppingError("Unknown joint-limit comparison type");
}
}

void MbD::LimitIJ::preDyn()
{
	ConstraintSet::preDyn();
	previousClearance = clearance(*this);
	transitionTime = std::numeric_limits<double>::quiet_NaN();
	previousReaction = normalReaction();
	releasing = false;
}

double MbD::LimitIJ::normalReaction() const
{
	// This solver uses force-minus-inertia residuals: lambda * grad(G)
	// is the physical constraint force. A unilateral stop cannot pull.
	return (type == "=>" ? 1.0 : -1.0) * constraints->front()->lam;
}

bool MbD::LimitIJ::hasTensileReaction() const
{
	return active && !compliant && normalReaction() < -1e-9;
}

void MbD::LimitIJ::preDynStep()
{
	previousClearance = clearance(*this);
	previousReaction = normalReaction();
}

double MbD::LimitIJ::checkForDynDiscontinuityBetweenand(double tprev, double t)
{
	if (compliant || !std::isfinite(previousClearance) || !(t > tprev)) return t;
	if (active) {
		if (!hasTensileReaction()) return t;
		const double current = normalReaction();
		const double denominator = previousReaction - current;
		const double fraction = denominator > 0
			? std::clamp(previousReaction / denominator, 0.0, 1.0) : 0.0;
		transitionTime = tprev + fraction * (t - tprev);
		releasing = true;
		return transitionTime;
	}

	const double currentClearance = clearance(*this);
	if (currentClearance >= -tol || previousClearance < -tol) return t;

	const double denominator = previousClearance - currentClearance;
	const double fraction
		= denominator > 0.0 ? std::clamp(previousClearance / denominator, 0.0, 1.0) : 0.0;
	transitionTime = tprev + fraction * (t - tprev);
	releasing = false;
	return transitionTime;
}

void MbD::LimitIJ::discontinuityAtaddTypeTo(
	double t,
	std::shared_ptr<std::vector<DiscontinuityType>> disconTypes
)
{
	const double scale = std::max({1.0, std::abs(t), std::abs(transitionTime)});
	if ((releasing == active) && std::isfinite(transitionTime)
		&& std::abs(t - transitionTime) <= 32.0 * std::numeric_limits<double>::epsilon() * scale) {
		if (releasing) deactivate();
		else activate();
		transitionTime = std::numeric_limits<double>::quiet_NaN();
		disconTypes->push_back(releasing ? LIFTOFF : TOUCHDOWN);
	}
}
