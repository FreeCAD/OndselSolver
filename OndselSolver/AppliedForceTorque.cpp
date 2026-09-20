// SPDX-License-Identifier: LGPL-2.1-or-later
#include "AppliedForceTorque.h"
#include "EndFrameqc.h"
#include "EndFrameqct.h"
#include "Constraint.h"
#include "Joint.h"
#include "SimulationStoppingError.h"
#include "Symbolic.h"
#include "System.h"
#include "SymTime.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace MbD;
namespace {
using V = AppliedForceTorque::Vector;
V operator+(const V& a, const V& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
V operator-(const V& a, const V& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
V operator*(double s, const V& a) { return {s*a[0], s*a[1], s*a[2]}; }
double dot(const V& a, const V& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
V cross(const V& a, const V& b)
{ return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}; }
V vector(const FColDsptr& a) { return {a->at(0), a->at(1), a->at(2)}; }
bool finite(const V& a) { return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]); }
bool nonnegative(const V& a) { return a[0] >= 0 && a[1] >= 0 && a[2] >= 0; }

// Position and angular-velocity maps for the seven body coordinates. Derivatives
// are analytical, including offset attachment points and quaternion curvature.
struct Endpoint {
    bool moving = false;
    V r{}, v{}, w{};
    std::array<size_t, 7> equation{};
    std::array<V, 7> J{}, W{}, dv{}, dw{};
    std::array<std::array<V, 7>, 7> H{}, dW{};
    explicit Endpoint(const std::shared_ptr<EndFramec>& frame)
    {
        r = vector(frame->rOeO);
        auto qc = std::dynamic_pointer_cast<EndFrameqc>(frame);
        if (!qc) return;
        moving = true;
        v = vector(qc->qXdot());
        for (size_t a = 0; a < 3; ++a) {
            equation[a] = qc->iqX()+a;
            J[a][a] = 1;
        }
        auto B = qc->aBOp();
        for (size_t a = 0; a < 4; ++a) {
            equation[a+3] = qc->iqE()+a;
            for (size_t xyz = 0; xyz < 3; ++xyz) {
                J[a+3][xyz] = qc->prOeOpE->at(xyz)->at(a);
                W[a+3][xyz] = 2*B->at(xyz)->at(a);
            }
            v = v + qc->qEdot()->at(a)*J[a+3];
            for (size_t b = 0; b < 4; ++b) {
                H[a+3][b+3] = vector(qc->pprOeOpEpE->at(a)->at(b));
                dv[a+3] = dv[a+3] + qc->qEdot()->at(b)*H[a+3][b+3];
            }
        }
        for (size_t b = 0; b < 4; ++b) {
            std::array<double, 4> e{};
            e[b] = 2;
            dW[3][b+3] = {e[3], e[2], -e[1]};
            dW[4][b+3] = {-e[2], e[3], e[0]};
            dW[5][b+3] = {e[1], -e[0], e[3]};
            dW[6][b+3] = {-e[0], -e[1], -e[2]};
        }
        for (size_t a = 0; a < 4; ++a) {
            w = w + qc->qEdot()->at(a)*W[a+3];
            for (size_t b = 0; b < 4; ++b)
                dw[b+3] = dw[b+3] + qc->qEdot()->at(a)*dW[a+3][b+3];
        }
    }
};
void validateFrames(const std::shared_ptr<EndFramec>& i, const std::shared_ptr<EndFramec>& j)
{
    if (!i || !j) throw std::invalid_argument("A load requires two attachment markers");
    if (std::dynamic_pointer_cast<EndFrameqct>(i) || std::dynamic_pointer_cast<EndFrameqct>(j))
        throw std::invalid_argument("Loads require markers fixed to their bodies");
}
}

struct AppliedForceTorque::Evaluation {
    Endpoint i, j;
    V separation{}, relativeVelocity{}, relativeAngularVelocity{}, unit{}, f{}, t{}, tj{};
    V xI{}, yI{}, xJ{}, yJ{};
    double length = 0, magnitude = 0, twistSin = 0, twistCos = 0;
    double frictionStaticMagnitude = 0, frictionDynamicMagnitude = 0;
    double storedEnergy = 0, dissipatedPower = 0;
    Evaluation(const std::shared_ptr<EndFramec>& a, const std::shared_ptr<EndFramec>& b) : i(a), j(b) {}
};

AppliedForceTorque::AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j, V f, V t)
    : frameI(i), frameJ(j), force(f), torque(t)
{
    validateFrames(i, j);
    if (!finite(f) || !finite(t)) throw std::invalid_argument("Load components must be finite");
}

AppliedForceTorque::AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j, double k, double c, double l)
    : frameI(i), frameJ(j), spring(true), stiffness(k), damping(c), restLength(l)
{
    validateFrames(i, j);
    if (!std::isfinite(k) || !std::isfinite(c) || !std::isfinite(l) || k < 0 || c < 0 || l < 0)
        throw std::invalid_argument("Spring stiffness, damping and rest length must be finite and nonnegative");
}

AppliedForceTorque::AppliedForceTorque(std::shared_ptr<EndFramec> i,
    std::shared_ptr<EndFramec> j, double k, double c, double angle, bool torsional)
    : frameI(i), frameJ(j), torsionalSpring(torsional), stiffness(k), damping(c), restLength(angle)
{
    validateFrames(i, j);
    if (!torsional || !std::isfinite(k) || !std::isfinite(c) || !std::isfinite(angle)
        || k < 0 || c < 0)
        throw std::invalid_argument("Torsional stiffness and damping must be finite and nonnegative");
}

AppliedForceTorque::AppliedForceTorque(std::shared_ptr<EndFramec> i,
    std::shared_ptr<EndFramec> j, Vector direction,
    std::shared_ptr<Symbolic> magnitude, bool isTorque)
    : frameI(i), frameJ(j), formulaDirection(direction),
      formulaMagnitude(std::move(magnitude)), formulaTorque(isTorque)
{
    validateFrames(i, j);
    if (!finite(direction) || !formulaMagnitude)
        throw std::invalid_argument("Formula load direction and magnitude are required");
}

AppliedForceTorque::AppliedForceTorque(std::shared_ptr<EndFramec> i,
    std::shared_ptr<EndFramec> j, ContactEvaluator evaluator)
    : frameI(i), frameJ(j), contactEvaluator(std::move(evaluator))
{
    validateFrames(i, j);
    if (!contactEvaluator)
        throw std::invalid_argument("A contact load requires an evaluator");
}

void AppliedForceTorque::validateBushingMatrix(const std::array<double, 36>& matrix)
{
    // Normalize by the diagonal before Cholesky: translation and rotation
    // blocks carry different units and may have very different magnitudes.
    double a[6][6] {};
    for (size_t i = 0; i < 6; ++i) {
        if (!std::isfinite(matrix[6*i+i]) || matrix[6*i+i] < 0)
            throw std::invalid_argument("Bushing matrix requires nonnegative finite diagonals");
        for (size_t j = 0; j < 6; ++j) {
            const double value = matrix[6*i+j];
            if (!std::isfinite(value) || value != matrix[6*j+i])
                throw std::invalid_argument("Bushing matrix must be finite and symmetric");
            const double scale = std::sqrt(matrix[6*i+i]) * std::sqrt(matrix[6*j+j]);
            if (scale > 0) a[i][j] = value / scale;
            else if (value != 0)
                throw std::invalid_argument("A zero bushing diagonal requires a zero row and column");
        }
    }
    for (size_t i = 0; i < 6; ++i) {
        if (a[i][i] < -1e-10)
            throw std::invalid_argument("Bushing matrix must be positive semidefinite");
        if (a[i][i] <= 1e-10) {
            for (size_t j = i + 1; j < 6; ++j)
                if (std::abs(a[j][i]) > 1e-10)
                    throw std::invalid_argument("Bushing matrix must be positive semidefinite");
            continue;
        }
        for (size_t j = i + 1; j < 6; ++j)
            for (size_t k = j; k < 6; ++k)
                a[k][j] = a[j][k] = a[k][j] - a[j][i] * a[k][i] / a[i][i];
    }
}

AppliedForceTorque::AppliedForceTorque(std::shared_ptr<EndFramec> i,
    std::shared_ptr<EndFramec> j, BushingParameters parameters)
    : frameI(i), frameJ(j), bushing(true), bushingParameters(parameters)
{
    validateFrames(i, j);
    if (parameters.coupled) {
        validateBushingMatrix(parameters.stiffnessMatrix);
        validateBushingMatrix(parameters.dampingMatrix);
    }
    if (!finite(parameters.linearStiffness) || !finite(parameters.linearDamping)
        || !finite(parameters.angularStiffness) || !finite(parameters.angularDamping)
        || !nonnegative(parameters.linearStiffness) || !nonnegative(parameters.linearDamping)
        || !nonnegative(parameters.angularStiffness) || !nonnegative(parameters.angularDamping))
        throw std::invalid_argument("Bushing stiffness and damping must be finite and nonnegative");
}

AppliedForceTorque::AppliedForceTorque(std::shared_ptr<EndFramec> i,
    std::shared_ptr<EndFramec> j, AxisFrictionParameters parameters)
    : frameI(i), frameJ(j), axisFriction(true), axisFrictionParameters(parameters)
{
    validateFrames(i, j);
    if (!std::isfinite(parameters.staticMagnitude)
        || !std::isfinite(parameters.dynamicMagnitude)
        || !std::isfinite(parameters.transitionVelocity)
        || !std::isfinite(parameters.viscousCoefficient)
        || parameters.staticMagnitude < parameters.dynamicMagnitude
        || parameters.dynamicMagnitude < 0 || !(parameters.transitionVelocity > 0)
        || parameters.viscousCoefficient < 0
        || (parameters.reactionBased
            && (!std::isfinite(parameters.effectiveRadius)
                || !(parameters.effectiveRadius > 0)
                || parameters.reactionJoint.expired())))
        throw std::invalid_argument(
            "Axis friction requires static >= dynamic >= 0, positive transition velocity, "
            "nonnegative viscous damping, and a valid reaction source"
        );
}

AppliedForceTorque::Evaluation AppliedForceTorque::evaluate() const
{
    Evaluation e(frameI, frameJ);
    if (!enabled) return e;
    e.separation = e.j.r-e.i.r;
    e.relativeVelocity = e.j.v-e.i.v;
    e.f = force;
    e.t = torque;
    if (contactEvaluator) {
        const auto state = [](const Endpoint& endpoint, const std::shared_ptr<EndFramec>& frame) {
            FrameState result;
            result.time = frame->root()->time->getValue();
            result.position = endpoint.r;
            result.velocity = endpoint.v;
            result.omega = endpoint.w;
            for (size_t row = 0; row < 3; ++row)
                for (size_t column = 0; column < 3; ++column)
                    result.rotation[3*row+column] = frame->aAOe->at(row)->at(column);
            return result;
        };
        const auto wrench = contactEvaluator(state(e.i, frameI), state(e.j, frameJ));
        if (!finite(wrench.forceOnI) || !finite(wrench.torqueOnI))
            throw SimulationStoppingError("Contact evaluator returned a non-finite wrench");
        e.f = wrench.forceOnI;
        e.t = wrench.torqueOnI;
        e.storedEnergy = wrench.storedEnergy;
        e.dissipatedPower = wrench.dissipatedPower;
    }
    else if (bushing) {
        const auto column = [](const FMatDsptr& matrix, size_t index) {
            return V {matrix->at(0)->at(index), matrix->at(1)->at(index), matrix->at(2)->at(index)};
        };
        std::array<V, 3> axesI, axesJ;
        std::array<double, 6> strains{}, rates{};
        for (size_t axis = 0; axis < 3; ++axis) {
            axesI[axis] = column(frameI->aAOe, axis);
            axesJ[axis] = column(frameJ->aAOe, axis);
            const double displacement = dot(axesI[axis], e.separation);
            // Translational strain is expressed in the rotating I frame.
            const double speed = dot(axesI[axis],
                e.relativeVelocity - cross(e.i.w, e.separation));
            strains[axis] = displacement;
            rates[axis] = speed;
            if (bushingParameters.coupled) continue;
            e.f = e.f + (bushingParameters.linearStiffness[axis] * displacement
                + bushingParameters.linearDamping[axis] * speed) * axesI[axis];
            e.storedEnergy += 0.5 * bushingParameters.linearStiffness[axis]
                * displacement * displacement;
            e.dissipatedPower += bushingParameters.linearDamping[axis] * speed * speed;
        }

        double relative[3][3] {};
        double trace = 0;
        for (size_t row = 0; row < 3; ++row) {
            for (size_t columnIndex = 0; columnIndex < 3; ++columnIndex)
                relative[row][columnIndex] = dot(axesI[row], axesJ[columnIndex]);
            trace += relative[row][row];
        }
        const double cosine = std::clamp(0.5 * (trace - 1), -1.0, 1.0);
        const double angle = std::acos(cosine);
        V rotationVector {
            relative[2][1] - relative[1][2],
            relative[0][2] - relative[2][0],
            relative[1][0] - relative[0][1]
        };
        const double sine = std::sin(angle);
        if (angle < 1e-8) {
            rotationVector = 0.5 * rotationVector;
        }
        else if (std::abs(sine) > 1e-8) {
            rotationVector = (angle / (2 * sine)) * rotationVector;
        }
        else {
            V axis {
                std::sqrt(std::max(0.0, 0.5 * (relative[0][0] + 1))),
                std::sqrt(std::max(0.0, 0.5 * (relative[1][1] + 1))),
                std::sqrt(std::max(0.0, 0.5 * (relative[2][2] + 1)))
            };
            const auto largest = std::distance(axis.begin(), std::max_element(axis.begin(), axis.end()));
            if (axis[largest] > 1e-8) {
                const size_t next = (largest + 1) % 3;
                const size_t last = (largest + 2) % 3;
                axis[next] = (relative[largest][next] + relative[next][largest])
                    / (4 * axis[largest]);
                axis[last] = (relative[largest][last] + relative[last][largest])
                    / (4 * axis[largest]);
            }
            rotationVector = angle * axis;
        }
        e.relativeAngularVelocity = e.j.w - e.i.w;
        // phiDot = J_l(phi)^-1 * R_I^T * (omega_J - omega_I).
        // A torque conjugate to phi is mapped with the transpose Jacobian.
        const double a = angle < 1e-4
            ? 1.0 / 12.0 + angle * angle / 720.0
            : (1.0 - 0.5 * angle / std::tan(0.5 * angle)) / (angle * angle);
        V localOmega {};
        for (size_t axis = 0; axis < 3; ++axis)
            localOmega[axis] = dot(axesI[axis], e.relativeAngularVelocity);
        const V strainRate = localOmega - 0.5 * cross(rotationVector, localOmega)
            + a * cross(rotationVector, cross(rotationVector, localOmega));
        V conjugate {};
        for (size_t axis = 0; axis < 3; ++axis) {
            strains[axis + 3] = rotationVector[axis];
            rates[axis + 3] = strainRate[axis];
            if (bushingParameters.coupled) continue;
            conjugate[axis] = bushingParameters.angularStiffness[axis] * rotationVector[axis]
                + bushingParameters.angularDamping[axis] * strainRate[axis];
            e.storedEnergy += 0.5 * bushingParameters.angularStiffness[axis]
                * rotationVector[axis] * rotationVector[axis];
            e.dissipatedPower += bushingParameters.angularDamping[axis]
                * strainRate[axis] * strainRate[axis];
        }
        if (bushingParameters.coupled) {
            std::array<double, 6> stress{};
            for (size_t row = 0; row < 6; ++row) {
                for (size_t column = 0; column < 6; ++column) {
                    const double elastic = bushingParameters.stiffnessMatrix[6*row+column]
                        * strains[column];
                    const double viscous = bushingParameters.dampingMatrix[6*row+column]
                        * rates[column];
                    stress[row] += elastic + viscous;
                    e.storedEnergy += 0.5 * strains[row] * elastic;
                    e.dissipatedPower += rates[row] * viscous;
                }
            }
            for (size_t axis = 0; axis < 3; ++axis) {
                e.f = e.f + stress[axis] * axesI[axis];
                conjugate[axis] = stress[axis + 3];
            }
        }
        const V localTorque = conjugate + 0.5 * cross(rotationVector, conjugate)
            + a * cross(rotationVector, cross(rotationVector, conjugate));
        // The translational frame-rotation couple belongs to I. Transporting
        // the total opposite wrench below leaves only -rotationalTorque on J.
        e.t = cross(e.separation, e.f);
        for (size_t axis = 0; axis < 3; ++axis)
            e.t = e.t + localTorque[axis] * axesI[axis];
    }
    else if (axisFriction) {
        e.unit = {
            frameI->aAOe->at(0)->at(2),
            frameI->aAOe->at(1)->at(2),
            frameI->aAOe->at(2)->at(2)
        };
        e.relativeAngularVelocity = e.j.w - e.i.w;
        const double speed = dot(
            e.unit,
            axisFrictionParameters.rotational ? e.relativeAngularVelocity : e.relativeVelocity
        );
        double scale = 1;
        if (axisFrictionParameters.reactionBased) {
            const auto joint = axisFrictionParameters.reactionJoint.lock();
            if (!joint) {
                throw SimulationStoppingError("Reaction-based friction lost its source joint");
            }
            const auto reaction = vector(joint->aFX());
            const auto transverse = axisFrictionParameters.axialReaction
                ? dot(reaction, e.unit) * e.unit
                : reaction - dot(reaction, e.unit) * e.unit;
            scale = std::sqrt(dot(transverse, transverse));
            if (axisFrictionParameters.rotational) {
                scale *= axisFrictionParameters.effectiveRadius;
            }
        }
        e.frictionStaticMagnitude = scale * axisFrictionParameters.staticMagnitude;
        e.frictionDynamicMagnitude = scale * axisFrictionParameters.dynamicMagnitude;
        const double ratio = speed / axisFrictionParameters.transitionVelocity;
        const double dryMagnitude = e.frictionDynamicMagnitude
            + (e.frictionStaticMagnitude - e.frictionDynamicMagnitude)
                * std::exp(-(ratio * ratio));
        const double magnitude = dryMagnitude * std::tanh(ratio)
            + axisFrictionParameters.viscousCoefficient * speed;
        e.length = speed;
        e.magnitude = magnitude;
        if (axisFrictionParameters.rotational)
            e.t = magnitude * e.unit;
        else
            e.f = magnitude * e.unit;
        e.dissipatedPower = magnitude * speed;
    }
    else if (formulaMagnitude) {
        const double magnitude = formulaMagnitude->getValue();
        if (!std::isfinite(magnitude))
            throw SimulationStoppingError("Load formula returned a non-finite magnitude");
        if (formulaTorque) e.t = magnitude*formulaDirection;
        else e.f = magnitude*formulaDirection;
    }
    else if (spring) {
        e.length = std::sqrt(dot(e.separation, e.separation));
        if (!(e.length > 0) || !std::isfinite(e.length))
            throw SimulationStoppingError("Spring attachment points must have a finite, nonzero separation");
        e.unit = (1/e.length)*e.separation;
        e.magnitude = stiffness*(e.length-restLength)+damping*dot(e.unit, e.relativeVelocity);
        e.f = e.magnitude*e.unit;
        const double extension = e.length - restLength;
        const double speed = dot(e.unit, e.relativeVelocity);
        e.storedEnergy = 0.5 * stiffness * extension * extension;
        e.dissipatedPower = damping * speed * speed;
    }
    else if (torsionalSpring) {
        const auto column = [](const FMatDsptr& matrix, size_t index) {
            return V {matrix->at(0)->at(index), matrix->at(1)->at(index), matrix->at(2)->at(index)};
        };
        e.xI = column(frameI->aAOe, 0);
        e.yI = column(frameI->aAOe, 1);
        e.xJ = column(frameJ->aAOe, 0);
        e.yJ = column(frameJ->aAOe, 1);
        e.unit = column(frameI->aAOe, 2);
        const double axisLength = std::sqrt(dot(e.unit, e.unit));
        if (!(axisLength > 0) || !std::isfinite(axisLength))
            throw SimulationStoppingError("Torsional spring axis must be finite and nonzero");
        e.unit = (1/axisLength)*e.unit;
        e.twistSin = dot(e.yI, e.xJ)-dot(e.xI, e.yJ);
        e.twistCos = dot(e.xI, e.xJ)+dot(e.yI, e.yJ);
        const double twistProjection = e.twistSin*e.twistSin+e.twistCos*e.twistCos;
        if (!(twistProjection > 1e-20) || !std::isfinite(twistProjection))
            throw SimulationStoppingError("Torsional spring twist is undefined for these attachment frames");
        // Spatial gradient of atan2(sinTwist, cosTwist). For coaxial
        // attachments this is I's Z axis; with swing it is not a unit axis.
        // Its work-conjugate torque preserves the defined elastic energy.
        const V sinGradient = cross(e.xJ, e.yI) - cross(e.yJ, e.xI);
        const V cosGradient = cross(e.xJ, e.xI) + cross(e.yJ, e.yI);
        e.unit = (1 / twistProjection)
            * (e.twistCos * sinGradient - e.twistSin * cosGradient);
        const double angle = continuousTwist(std::atan2(e.twistSin, e.twistCos));
        e.relativeAngularVelocity = e.j.w-e.i.w;
        e.length = angle;
        e.magnitude = stiffness*(angle-restLength)
            + damping*dot(e.unit, e.relativeAngularVelocity);
        e.t = e.magnitude*e.unit;
        const double deflection = angle - restLength;
        const double speed = dot(e.unit, e.relativeAngularVelocity);
        e.storedEnergy = 0.5 * stiffness * deflection * deflection;
        e.dissipatedPower = damping * speed * speed;
    }
    // Transport the opposite wrench to J; the net force AND moment are zero.
    if (follower && !contactEvaluator && !bushing && !axisFriction
        && !spring && !torsionalSpring) {
        const auto world = [&](const V& local) {
            V result{};
            for (size_t row = 0; row < 3; ++row)
                for (size_t col = 0; col < 3; ++col)
                    result[row] += frameI->aAOe->at(row)->at(col) * local[col];
            return result;
        };
        e.f = world(e.f);
        e.t = world(e.t);
    }
    e.tj = cross(e.separation, e.f)-e.t;
    return e;
}

double AppliedForceTorque::continuousTwist(double principal) const
{
    const double time = frameI->root()->time->getValue();
    double reference = restLength;
    for (const auto& sample : acceptedTwists) {
        if (sample.first > time) break;
        reference = sample.second;
    }
    constexpr double twoPi = 6.2831853071795864769;
    return reference + std::remainder(principal - reference, twoPi);
}

void AppliedForceTorque::postDynFirstStep() { postDynStep(); }

void AppliedForceTorque::setContactStepValidator(ContactStepValidator validator)
{
    contactStepValidator = std::move(validator);
}

void AppliedForceTorque::preDynStep()
{
    if (contactStepValidator) {
        previousI = frameState(frameI);
        previousJ = frameState(frameJ);
        previousTime = frameI->root()->time->getValue();
    }
}

bool AppliedForceTorque::acceptDynTrial() const
{
    if (!contactStepValidator) return true;
    const auto currentI = frameState(frameI), currentJ = frameState(frameJ);
    const double dt = std::abs(frameI->root()->time->getValue() - previousTime);
    // Endpoint orientations alone cannot distinguish a full turn from no turn.
    // Keep each contact sweep within a small angular interval, including an
    // accelerating trial. Rejected trials do not overwrite the starting state.
    for (const auto& state : {previousI, previousJ, currentI, currentJ})
        if (dt * std::sqrt(dot(state.omega, state.omega)) > 0.25) return false;
    return contactStepValidator(previousI, currentI, previousJ, currentJ);
}

void AppliedForceTorque::postDynStep()
{
    if (!torsionalSpring) return;
    const double time = frameI->root()->time->getValue();
    const double twist = evaluate().length;
    while (!acceptedTwists.empty() && acceptedTwists.back().first >= time)
        acceptedTwists.pop_back();
    acceptedTwists.emplace_back(time, twist);
    // Output interpolation and discontinuity rollback only need the current
    // accepted interval. Retain an extra predecessor for a restarted first step.
    if (acceptedTwists.size() > 3) acceptedTwists.erase(acceptedTwists.begin());
}

double AppliedForceTorque::suggestSmallerOrAcceptDynFirstStepSize(double h)
{
    return suggestSmallerOrAcceptDynStepSize(h);
}

double AppliedForceTorque::suggestSmallerOrAcceptDynStepSize(double h)
{
    if (!torsionalSpring) return h;
    const Endpoint i(frameI), j(frameJ);
    const V relative = j.w - i.w;
    const double speed = std::sqrt(dot(relative, relative));
    // Sample well before a half turn so a winding cannot be skipped.
    return speed > 0 ? std::min(h, 1.0 / speed) : h;
}

AppliedForceTorque::Vector AppliedForceTorque::forceOnI() const { return evaluate().f; }
AppliedForceTorque::Vector AppliedForceTorque::torqueOnI() const { return evaluate().t; }
AppliedForceTorque::EnergyState AppliedForceTorque::energyState() const
{
    return resultState().energy;
}

AppliedForceTorque::FrameState AppliedForceTorque::frameState(
    const std::shared_ptr<EndFramec>& frame
)
{
    const Endpoint endpoint(frame);
    FrameState result;
    result.time = frame->root()->time->getValue();
    result.position = endpoint.r;
    result.velocity = endpoint.v;
    result.omega = endpoint.w;
    for (size_t row = 0; row < 3; ++row) {
        for (size_t column = 0; column < 3; ++column) {
            result.rotation[3 * row + column] = frame->aAOe->at(row)->at(column);
        }
    }
    // For a time-prescribed end frame, these are deliberately the velocities
    // of the underlying physical body. Including the frame's explicit motion
    // would make actuator work cancel against its own virtual reference side.
    return result;
}

AppliedForceTorque::ResultState AppliedForceTorque::resultState() const
{
    const auto e = evaluate();
    const double power = dot(e.f, e.i.v) + dot(e.t, e.i.w)
        - dot(e.f, e.j.v) + dot(e.tj, e.j.w);
    return {e.f, e.t, {power, e.storedEnergy, e.dissipatedPower}};
}

void AppliedForceTorque::fillAccICIterError(FColDsptr col) { fillDynError(col); }
void AppliedForceTorque::fillDynError(FColDsptr col)
{
    const auto e = evaluate();
    for (size_t a = 0; a < 7; ++a) {
        if (e.i.moving) col->at(e.i.equation[a]) += dot(e.i.J[a], e.f)+dot(e.i.W[a], e.t);
        if (e.j.moving) col->at(e.j.equation[a]) += -dot(e.j.J[a], e.f)+dot(e.j.W[a], e.tj);
    }
}

void AppliedForceTorque::fillpFpy(SpMatDsptr mat) { fillJacobian(mat, false); }
void AppliedForceTorque::fillpFpydot(SpMatDsptr mat) { fillJacobian(mat, true); }
void AppliedForceTorque::fillJacobian(SpMatDsptr mat, bool velocity) const
{
    if (!enabled) return;
    // CAD shape contact is an externally evaluated, nonsmooth penalty force.
    // Re-evaluate it for every residual, but use a lagged tangent so Newton
    // iterations do not multiply expensive collision queries by 14 numerical
    // perturbations. Small integration steps provide the required stability.
    if (contactEvaluator || bushing || follower) return;
    const auto e = evaluate();
    for (size_t side = 0; side < 2; ++side) {
        const auto& source = side == 0 ? e.i : e.j;
        if (!source.moving) continue;
        const double sign = side == 0 ? -1 : 1;
        for (size_t b = 0; b < 7; ++b) {
            const V dD = velocity ? V{} : sign*source.J[b];
            const V dV = sign*(velocity ? source.J[b] : source.dv[b]);
            V dF{}, dT{};
            if (spring) {
                const double dl = dot(e.unit, dD);
                const V du = (1/e.length)*(dD-dl*e.unit);
                dF = (stiffness*dl+damping*(dot(du, e.relativeVelocity)+dot(e.unit, dV)))*e.unit
                    + e.magnitude*du;
            }
            else if (torsionalSpring) {
                V dAxis{};
                double dAngle = 0;
                double dOmega = 0;
                if (velocity) {
                    dOmega = sign*dot(e.unit, source.W[b]);
                }
                else {
                    V dxI{}, dyI{}, dxJ{}, dyJ{};
                    if (side == 0) {
                        dxI = cross(source.W[b], e.xI);
                        dyI = cross(source.W[b], e.yI);
                    }
                    else {
                        dxJ = cross(source.W[b], e.xJ);
                        dyJ = cross(source.W[b], e.yJ);
                    }
                    const double dSin = dot(dyI, e.xJ)+dot(e.yI, dxJ)
                        -dot(dxI, e.yJ)-dot(e.xI, dyJ);
                    const double dCos = dot(dxI, e.xJ)+dot(e.xI, dxJ)
                        +dot(dyI, e.yJ)+dot(e.yI, dyJ);
                    dAngle = (e.twistCos*dSin-e.twistSin*dCos)
                        /(e.twistSin*e.twistSin+e.twistCos*e.twistCos);
                    dOmega = sign*dot(e.unit, source.dw[b]);
                    const V sinGradient = cross(e.xJ, e.yI) - cross(e.yJ, e.xI);
                    const V cosGradient = cross(e.xJ, e.xI) + cross(e.yJ, e.yI);
                    const V dSinGradient = cross(dxJ, e.yI) + cross(e.xJ, dyI)
                        - cross(dyJ, e.xI) - cross(e.yJ, dxI);
                    const V dCosGradient = cross(dxJ, e.xI) + cross(e.xJ, dxI)
                        + cross(dyJ, e.yI) + cross(e.yJ, dyI);
                    const double projection = e.twistSin*e.twistSin + e.twistCos*e.twistCos;
                    dAxis = (1 / projection) * (dCos*sinGradient + e.twistCos*dSinGradient
                        - dSin*cosGradient - e.twistSin*dCosGradient
                        - (2*(e.twistSin*dSin + e.twistCos*dCos))*e.unit);
                    dOmega += dot(dAxis, e.relativeAngularVelocity);
                }
                const double dMagnitude = stiffness*dAngle+damping*dOmega;
                dT = dMagnitude*e.unit+e.magnitude*dAxis;
            }
            else if (axisFriction) {
                // The velocity tangent is the part that controls convergence
                // through the sharp zero-speed transition. Keep the positional
                // orientation dependence lagged; quaternion curvature is still
                // included below in the generalized-force transformation.
                if (velocity) {
                    const V dRelative = sign * (
                        axisFrictionParameters.rotational ? source.W[b] : source.J[b]
                    );
                    const double dSpeed = dot(e.unit, dRelative);
                    const double transition = axisFrictionParameters.transitionVelocity;
                    const double ratio = e.length / transition;
                    const double exponential = std::exp(-(ratio * ratio));
                    const double tangent = std::tanh(ratio);
                    const double dryMagnitude = e.frictionDynamicMagnitude
                        + (e.frictionStaticMagnitude - e.frictionDynamicMagnitude) * exponential;
                    const double drySlope
                        = (e.frictionStaticMagnitude - e.frictionDynamicMagnitude)
                            * exponential * (-2 * ratio / transition);
                    const double slope = drySlope * tangent
                        + dryMagnitude * (1 - tangent*tangent) / transition
                        + axisFrictionParameters.viscousCoefficient;
                    if (axisFrictionParameters.rotational)
                        dT = slope*dSpeed*e.unit;
                    else
                        dF = slope*dSpeed*e.unit;
                }
            }
            const V dTj = cross(dD, e.f)+cross(e.separation, dF);
            for (size_t a = 0; a < 7; ++a) {
                double qi = dot(e.i.J[a], dF)+dot(e.i.W[a], dT);
                double qj = -dot(e.j.J[a], dF)+dot(e.j.W[a], dTj-dT);
                if (!velocity && side == 0)
                    qi += dot(e.i.H[a][b], e.f)+dot(e.i.dW[a][b], e.t);
                if (!velocity && side == 1)
                    qj += -dot(e.j.H[a][b], e.f)+dot(e.j.dW[a][b], e.tj);
                if (e.i.moving && qi != 0) mat->atijplusNumber(e.i.equation[a], source.equation[b], qi);
                if (e.j.moving && qj != 0) mat->atijplusNumber(e.j.equation[a], source.equation[b], qj);
            }
        }
    }
    // Dynamic multipliers are stored in qdot, so this coupling belongs in the
    // velocity Jacobian even though it represents a reaction-force dependency.
    if (velocity && axisFriction && axisFrictionParameters.reactionBased) {
        const auto joint = axisFrictionParameters.reactionJoint.lock();
        const auto reaction = vector(joint->aFX());
        const auto transverse = axisFrictionParameters.axialReaction
            ? dot(reaction, e.unit) * e.unit
            : reaction - dot(reaction, e.unit) * e.unit;
        const double normalLoad = std::sqrt(dot(transverse, transverse));
        if (normalLoad > 1e-14) {
            const double ratio = e.length / axisFrictionParameters.transitionVelocity;
            const double exponential = std::exp(-(ratio * ratio));
            const double coefficient = axisFrictionParameters.dynamicMagnitude
                + (axisFrictionParameters.staticMagnitude
                   - axisFrictionParameters.dynamicMagnitude) * exponential;
            const double speedSign = std::tanh(ratio);
            const double radius = axisFrictionParameters.rotational
                ? axisFrictionParameters.effectiveRadius
                : 1;
            joint->constraintsDo([&](const std::shared_ptr<Constraint>& constraint) {
                if (constraint->iG == SIZE_MAX || std::abs(constraint->lam) < 1e-14) {
                    return;
                }
                auto contribution = std::make_shared<FullColumn<double>>(3, 0.0);
                constraint->addToJointForceI(contribution);
                const V dReaction = (1 / constraint->lam) * vector(contribution);
                const V dTransverse = axisFrictionParameters.axialReaction
                    ? dot(dReaction, e.unit) * e.unit
                    : dReaction - dot(dReaction, e.unit) * e.unit;
                const double dNormal = dot(transverse, dTransverse) / normalLoad;
                const double dMagnitude = radius * coefficient * speedSign * dNormal;
                const V dF = axisFrictionParameters.rotational ? V{} : dMagnitude * e.unit;
                const V dT = axisFrictionParameters.rotational ? dMagnitude * e.unit : V{};
                const V dTj = cross(e.separation, dF) - dT;
                for (size_t a = 0; a < 7; ++a) {
                    const double qi = dot(e.i.J[a], dF) + dot(e.i.W[a], dT);
                    const double qj = -dot(e.j.J[a], dF) + dot(e.j.W[a], dTj);
                    if (e.i.moving && qi != 0) {
                        mat->atijplusNumber(e.i.equation[a], constraint->iG, qi);
                    }
                    if (e.j.moving && qj != 0) {
                        mat->atijplusNumber(e.j.equation[a], constraint->iG, qj);
                    }
                }
            });
        }
    }
}
