/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#pragma once

#include "ASMTItemIJ.h"
#include "AppliedForceTorque.h"
#include <array>

namespace MbD {
    class ASMTJoint;
    class ASMTForceTorque : public ASMTItemIJ
    {
        //
    public:
        static std::shared_ptr<ASMTForceTorque> With();
        void initialize() override;
        void clearResults();
        // Constant world-resolved load on I; the balancing wrench acts on J.
        void setForce3D(double x, double y, double z);
        void setTorque3D(double x, double y, double z);
        // When true, constant/formula vector components are marker-I-local.
        void setFollower(bool value) { follower = value; }
        void setForceFormula3D(double x, double y, double z, const std::string& expression);
        void setTorqueFormula3D(double x, double y, double z, const std::string& expression);
        void setSpringDamper(double stiffness, double damping, double restLength);
        void setTorsionalSpringDamper(double stiffness, double damping, double restAngle);
        void setBushing(const std::array<double, 3>& linearStiffness,
                        const std::array<double, 3>& linearDamping,
                        const std::array<double, 3>& angularStiffness,
                        const std::array<double, 3>& angularDamping);
        void setCoupledBushing(const std::array<double, 36>& stiffness,
                               const std::array<double, 36>& damping);
        void setAxisFriction(double staticMagnitude,
                             double dynamicMagnitude,
                             double transitionVelocity,
                             double viscousCoefficient,
                             bool rotational);
        void setReactionAxisFriction(double staticCoefficient,
                                     double dynamicCoefficient,
                                     double transitionVelocity,
                                     double viscousCoefficient,
                                     double effectiveRadius,
                                     bool rotational,
                                     const std::shared_ptr<ASMTJoint>& reactionJoint,
                                     bool axialReaction = false);
        void setContactEvaluator(AppliedForceTorque::ContactEvaluator evaluator);
        void setContactStepValidator(AppliedForceTorque::ContactStepValidator validator);
        void createMbD(std::shared_ptr<System> system, std::shared_ptr<Units> units) override;
        void parseASMT(std::vector<std::string>& lines) override;
        void storeOnLevel(std::ofstream& os, size_t level) override;
        void updateFromMbD() override;
        void compareResults(AnalysisType type) override;
        void outputResults(AnalysisType type) override;
        bool isPassive() const;

        FRowDsptr powers, storedEnergies, dissipatedPowers;


    private:
        std::array<double, 3> force{}, torque{};
        std::array<double, 3> formulaDirection{};
        std::string formula;
        AppliedForceTorque::ContactEvaluator contactEvaluator;
        AppliedForceTorque::ContactStepValidator contactStepValidator;
        bool formulaTorque = false;
        bool follower = false;
        bool spring = false;
        bool torsionalSpring = false;
        bool bushing = false;
        bool axisFriction = false;
        double stiffness = 0, damping = 0, restLength = 0;
        AppliedForceTorque::BushingParameters bushingParameters;
        AppliedForceTorque::AxisFrictionParameters axisFrictionParameters;
        std::weak_ptr<ASMTJoint> reactionJoint;
    };
}

