// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <array>
#include <functional>
#include <memory>
#include <vector>
#include "ForceTorqueItem.h"
#include "EndFramec.h"

namespace MbD {
class Symbolic;
class Joint;
// A world-resolved wrench on I, with the balancing wrench on J. Linear springs
// act between attachment points; torsional springs act about attachment I's Z axis.
class AppliedForceTorque : public ForceTorqueItem
{
public:
    using Vector = std::array<double, 3>;
    using Matrix = std::array<double, 9>;
    struct FrameState {
        Vector position{};
        Matrix rotation{};
        Vector velocity{};
        Vector omega{};
        double time = 0;
    };
    struct ContactWrench {
        Vector forceOnI{};
        Vector torqueOnI{};
        double storedEnergy = 0;
        double dissipatedPower = 0;
    };
    struct EnergyState {
        // Positive power adds mechanical energy to the connected bodies.
        double power = 0;
        double storedEnergy = 0;
        // Positive dissipated power removes mechanical energy.
        double dissipatedPower = 0;
    };
    struct ResultState {
        Vector forceOnI{};
        Vector torqueOnI{};
        EnergyState energy;
    };
    struct BushingParameters {
        Vector linearStiffness{};
        Vector linearDamping{};
        Vector angularStiffness{};
        Vector angularDamping{};
        // Row-major, work-conjugate [x,y,z,rx,ry,rz] matrices in marker I.
        // Rotations are rotation-vector strains in radians, not Euler angles.
        std::array<double, 36> stiffnessMatrix{};
        std::array<double, 36> dampingMatrix{};
        bool coupled = false;
    };
    static void validateBushingMatrix(const std::array<double, 36>& matrix);
    struct AxisFrictionParameters {
        double staticMagnitude = 0;
        double dynamicMagnitude = 0;
        double transitionVelocity = 0;
        double viscousCoefficient = 0;
        bool rotational = false;
        // When enabled, the two magnitudes are friction coefficients. The
        // normal load is the joint reaction perpendicular to the free axis.
        bool reactionBased = false;
        double effectiveRadius = 1;
        std::weak_ptr<Joint> reactionJoint;
        bool axialReaction = false; // Thrust bearing; otherwise radial bearing load.
    };
    using ContactEvaluator = std::function<ContactWrench(const FrameState&, const FrameState&)>;
    using ContactStepValidator = std::function<bool(const FrameState&, const FrameState&,
                                                    const FrameState&, const FrameState&)>;
    void setContactStepValidator(ContactStepValidator validator);
    void setFollower(bool value) { follower = value; }
    void setEnabled(bool value) { enabled = value; }
    void setRuntimeFormula(Vector direction, std::shared_ptr<Symbolic> magnitude, bool isTorque) {
        formulaDirection = direction;
        formulaMagnitude = std::move(magnitude);
        formulaTorque = isTorque;
    }
    void preDynStep() override;
    bool acceptDynTrial() const override;
    AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j, Vector force, Vector torque);
    AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j, double stiffness,
                       double damping, double restLength);
    AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j, double stiffness,
                       double damping, double restAngle, bool torsional);
    AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j,
                       Vector direction, std::shared_ptr<Symbolic> magnitude, bool torque);
    AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j,
                       ContactEvaluator evaluator);
    AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j,
                       BushingParameters parameters);
    AppliedForceTorque(std::shared_ptr<EndFramec> i, std::shared_ptr<EndFramec> j,
                       AxisFrictionParameters parameters);
    void fillAccICIterError(FColDsptr col) override;
    void fillDynError(FColDsptr col) override;
    void fillpFpy(SpMatDsptr mat) override;
    void fillpFpydot(SpMatDsptr mat) override;
    Vector forceOnI() const;
    Vector torqueOnI() const;
    EnergyState energyState() const;
    ResultState resultState() const;
    void postDynFirstStep() override;
    void postDynStep() override;
    double suggestSmallerOrAcceptDynFirstStepSize(double h) override;
    double suggestSmallerOrAcceptDynStepSize(double h) override;
    static FrameState frameState(const std::shared_ptr<EndFramec>& frame);

private:
    struct Evaluation;
    Evaluation evaluate() const;
    void fillJacobian(SpMatDsptr mat, bool velocity) const;
    std::shared_ptr<EndFramec> frameI, frameJ;
    Vector force{}, torque{};
    Vector formulaDirection{};
    std::shared_ptr<Symbolic> formulaMagnitude;
    ContactEvaluator contactEvaluator;
    ContactStepValidator contactStepValidator;
    FrameState previousI, previousJ;
    double previousTime = 0;
    bool formulaTorque = false;
    bool follower = false;
    bool enabled = true;
    bool spring = false;
    bool torsionalSpring = false;
    bool bushing = false;
    bool axisFriction = false;
    double stiffness = 0, damping = 0, restLength = 0;
    BushingParameters bushingParameters;
    AxisFrictionParameters axisFrictionParameters;
    // Only accepted states are retained. Residual/Jacobian/output evaluation
    // never changes winding, including when a corrector retries a step.
    std::vector<std::pair<double, double>> acceptedTwists;
    double continuousTwist(double principal) const;
};
}
