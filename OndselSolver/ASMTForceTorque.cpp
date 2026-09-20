/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.
 *   This file is part of OndselSolver.
 *   See LICENSE file for details about copyright.
 ***************************************************************************/
#include "ASMTForceTorque.h"
#include "AppliedForceTorque.h"
#include "ASMTAssembly.h"
#include "ASMTMarker.h"
#include "ASMTJoint.h"
#include "System.h"
#include "Units.h"
#include "SymbolicParser.h"
#include "BasicUserFunction.h"
#include "Constant.h"
#include "Joint.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace MbD;

std::shared_ptr<ASMTForceTorque> ASMTForceTorque::With()
{
    auto item = std::make_shared<ASMTForceTorque>();
    item->initialize();
    return item;
}

void ASMTForceTorque::initialize()
{
    ASMTItemIJ::initialize();
    powers = std::make_shared<FullRow<double>>();
    storedEnergies = std::make_shared<FullRow<double>>();
    dissipatedPowers = std::make_shared<FullRow<double>>();
}

void ASMTForceTorque::clearResults()
{
    ASMTItemIJ::clearResults();
    for (auto row : {powers, storedEnergies, dissipatedPowers}) {
        if (row) row->clear();
    }
}

void ASMTForceTorque::setForce3D(double x, double y, double z)
{
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        throw std::invalid_argument("Force components must be finite");
    force = {x, y, z};
    contactEvaluator = {};
    formula.clear();
    spring = false;
    torsionalSpring = false;
    bushing = false;
    axisFriction = false;
}

void ASMTForceTorque::setTorque3D(double x, double y, double z)
{
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        throw std::invalid_argument("Torque components must be finite");
    torque = {x, y, z};
    contactEvaluator = {};
    formula.clear();
    spring = false;
    torsionalSpring = false;
    bushing = false;
    axisFriction = false;
}

namespace {
void setFormula(std::array<double, 3>& direction, std::string& target,
                bool& torque, bool isTorque, double x, double y, double z,
                const std::string& expression)
{
    const double length = std::sqrt(x*x + y*y + z*z);
    if (!std::isfinite(length) || length <= 0 || expression.empty())
        throw std::invalid_argument("A formula load requires a finite direction and expression");
    direction = {x/length, y/length, z/length};
    target = expression;
    torque = isTorque;
}
}

void ASMTForceTorque::setForceFormula3D(double x, double y, double z,
                                        const std::string& expression)
{
    setFormula(formulaDirection, formula, formulaTorque, false, x, y, z, expression);
    contactEvaluator = {};
    spring = false;
    torsionalSpring = false;
    bushing = false;
    axisFriction = false;
}

void ASMTForceTorque::setTorqueFormula3D(double x, double y, double z,
                                         const std::string& expression)
{
    setFormula(formulaDirection, formula, formulaTorque, true, x, y, z, expression);
    contactEvaluator = {};
    spring = false;
    torsionalSpring = false;
    bushing = false;
    axisFriction = false;
}

void ASMTForceTorque::setSpringDamper(double k, double c, double l)
{
    if (!std::isfinite(k) || !std::isfinite(c) || !std::isfinite(l) || k < 0 || c < 0 || l < 0)
        throw std::invalid_argument("Spring stiffness, damping and rest length must be finite and nonnegative");
    stiffness = k;
    damping = c;
    restLength = l;
    contactEvaluator = {};
    formula.clear();
    spring = true;
    torsionalSpring = false;
    bushing = false;
    axisFriction = false;
}

void ASMTForceTorque::setTorsionalSpringDamper(double k, double c, double angle)
{
    if (!std::isfinite(k) || !std::isfinite(c) || !std::isfinite(angle) || k < 0 || c < 0)
        throw std::invalid_argument("Torsional stiffness and damping must be finite and nonnegative");
    stiffness = k;
    damping = c;
    restLength = angle;
    contactEvaluator = {};
    formula.clear();
    spring = false;
    torsionalSpring = true;
    bushing = false;
    axisFriction = false;
}

void ASMTForceTorque::setContactStepValidator(AppliedForceTorque::ContactStepValidator validator)
{
    contactStepValidator = std::move(validator);
}

void ASMTForceTorque::setContactEvaluator(AppliedForceTorque::ContactEvaluator evaluator)
{
    if (!evaluator)
        throw std::invalid_argument("A contact load requires an evaluator");
    contactEvaluator = std::move(evaluator);
    formula.clear();
    spring = false;
    torsionalSpring = false;
    bushing = false;
    axisFriction = false;
}

void ASMTForceTorque::setBushing(const std::array<double, 3>& linearStiffness,
    const std::array<double, 3>& linearDamping,
    const std::array<double, 3>& angularStiffness,
    const std::array<double, 3>& angularDamping)
{
    const auto valid = [](const auto& values) {
        return std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value) && value >= 0;
        });
    };
    if (!valid(linearStiffness) || !valid(linearDamping)
        || !valid(angularStiffness) || !valid(angularDamping))
        throw std::invalid_argument("Bushing stiffness and damping must be finite and nonnegative");
    bushingParameters = {linearStiffness, linearDamping, angularStiffness, angularDamping};
    contactEvaluator = {};
    formula.clear();
    spring = false;
    torsionalSpring = false;
    bushing = true;
    axisFriction = false;
}

void ASMTForceTorque::setCoupledBushing(const std::array<double, 36>& stiffness,
                                      const std::array<double, 36>& damping)
{
    AppliedForceTorque::validateBushingMatrix(stiffness);
    AppliedForceTorque::validateBushingMatrix(damping);
    setBushing({}, {}, {}, {});
    bushingParameters.coupled = true;
    bushingParameters.stiffnessMatrix = stiffness;
    bushingParameters.dampingMatrix = damping;
}

void ASMTForceTorque::setAxisFriction(double staticMagnitude,
    double dynamicMagnitude, double transitionVelocity,
    double viscousCoefficient, bool rotational)
{
    AppliedForceTorque::AxisFrictionParameters parameters {
        staticMagnitude,
        dynamicMagnitude,
        transitionVelocity,
        viscousCoefficient,
        rotational
    };
    // Let the solver object own the common validation contract. Constructing a
    // temporary is not possible before markers exist, so validate scalars here.
    if (!std::isfinite(staticMagnitude) || !std::isfinite(dynamicMagnitude)
        || !std::isfinite(transitionVelocity) || !std::isfinite(viscousCoefficient)
        || staticMagnitude < dynamicMagnitude || dynamicMagnitude < 0
        || !(transitionVelocity > 0) || viscousCoefficient < 0)
        throw std::invalid_argument(
            "Axis friction requires static >= dynamic >= 0, positive transition velocity, "
            "and nonnegative viscous damping"
        );
    axisFrictionParameters = parameters;
    reactionJoint.reset();
    contactEvaluator = {};
    formula.clear();
    spring = false;
    torsionalSpring = false;
    bushing = false;
    axisFriction = true;
}

void ASMTForceTorque::setReactionAxisFriction(double staticCoefficient,
    double dynamicCoefficient, double transitionVelocity,
    double viscousCoefficient, double effectiveRadius, bool rotational,
    const std::shared_ptr<ASMTJoint>& source, bool axialReaction)
{
    if (!source || !std::isfinite(staticCoefficient)
        || !std::isfinite(dynamicCoefficient) || !std::isfinite(transitionVelocity)
        || !std::isfinite(viscousCoefficient) || !std::isfinite(effectiveRadius)
        || staticCoefficient < dynamicCoefficient || dynamicCoefficient < 0
        || !(transitionVelocity > 0) || viscousCoefficient < 0
        || (rotational && !(effectiveRadius > 0))) {
        throw std::invalid_argument(
            "Reaction friction requires static >= dynamic >= 0, positive transition velocity, "
            "nonnegative viscous damping, and a positive rotational effective radius"
        );
    }
    setAxisFriction(staticCoefficient, dynamicCoefficient, transitionVelocity,
        viscousCoefficient, rotational);
    axisFrictionParameters.reactionBased = true;
    if (axialReaction && !rotational)
        throw std::invalid_argument("Thrust bearing friction requires a rotational axis");
    axisFrictionParameters.axialReaction = axialReaction;
    axisFrictionParameters.effectiveRadius = rotational ? effectiveRadius : 1;
    reactionJoint = source;
}

void ASMTForceTorque::createMbD(std::shared_ptr<System> system, std::shared_ptr<Units> units)
{
    auto i = root()->markerAt(markerI);
    auto j = root()->markerAt(markerJ);
    if (!i || !j) throw std::invalid_argument("Load attachment marker was not found");
    auto fi = std::dynamic_pointer_cast<EndFramec>(i->mbdObject);
    auto fj = std::dynamic_pointer_cast<EndFramec>(j->mbdObject);
    std::shared_ptr<AppliedForceTorque> load;
    if (contactEvaluator) {
        load = std::make_shared<AppliedForceTorque>(fi, fj, contactEvaluator);
    }
    else if (bushing) {
        auto parameters = bushingParameters;
        for (size_t row = 0; row < 6; ++row) {
            const double effort = row < 3 ? units->force : units->torque;
            for (size_t col = 0; col < 6; ++col) {
                parameters.stiffnessMatrix[6*row+col]
                    *= (col < 3 ? units->length : units->angle) / effort;
                parameters.dampingMatrix[6*row+col]
                    *= (col < 3 ? units->velocity : units->omega) / effort;
            }
        }
        for (size_t axis = 0; axis < 3; ++axis) {
            parameters.linearStiffness[axis] *= units->length / units->force;
            parameters.linearDamping[axis] *= units->velocity / units->force;
            parameters.angularStiffness[axis] *= units->angle / units->torque;
            parameters.angularDamping[axis] *= units->omega / units->torque;
        }
        load = std::make_shared<AppliedForceTorque>(fi, fj, parameters);
    }
    else if (axisFriction) {
        auto parameters = axisFrictionParameters;
        if (parameters.reactionBased) {
            const auto source = reactionJoint.lock();
            if (!source) {
                throw std::invalid_argument("Reaction friction source joint is no longer available");
            }
            parameters.reactionJoint = std::dynamic_pointer_cast<Joint>(source->mbdObject);
            if (parameters.reactionJoint.expired()) {
                throw std::invalid_argument("Reaction friction source solver joint was not created");
            }
            parameters.effectiveRadius /= units->length;
        }
        if (parameters.rotational) {
            if (!parameters.reactionBased) {
                parameters.staticMagnitude /= units->torque;
                parameters.dynamicMagnitude /= units->torque;
            }
            parameters.transitionVelocity /= units->omega;
            parameters.viscousCoefficient *= units->omega / units->torque;
        }
        else {
            if (!parameters.reactionBased) {
                parameters.staticMagnitude /= units->force;
                parameters.dynamicMagnitude /= units->force;
            }
            parameters.transitionVelocity /= units->velocity;
            parameters.viscousCoefficient *= units->velocity / units->force;
        }
        load = std::make_shared<AppliedForceTorque>(fi, fj, parameters);
    }
    else if (spring) {
        load = std::make_shared<AppliedForceTorque>(fi, fj, stiffness*units->length/units->force,
            damping*units->velocity/units->force, restLength/units->length);
    }
    else if (torsionalSpring) {
        load = std::make_shared<AppliedForceTorque>(fi, fj,
            stiffness*units->angle/units->torque,
            damping*units->omega/units->torque,
            restLength/units->angle, true);
    }
    else if (!formula.empty()) {
        auto parser = std::make_shared<SymbolicParser>();
        parser->owner = this;
        parser->variables->insert(std::make_pair("time", root()->geoTime()));
        auto userFunction = std::make_shared<BasicUserFunction>(formula, 1.0);
        parser->parseUserFunction(userFunction);
        auto magnitude = parser->stack->top();
        const double unit = formulaTorque ? units->torque : units->force;
        magnitude = Symbolic::times(magnitude, sptrConstant(1.0 / unit));
        magnitude->createMbD(system, units);
        magnitude = magnitude->simplified(magnitude);
        load = std::make_shared<AppliedForceTorque>(
            fi, fj, formulaDirection, magnitude, formulaTorque
        );
    }
    else {
        auto f = force, t = torque;
        for (size_t axis = 0; axis < 3; ++axis) {
            f[axis] /= units->force;
            t[axis] /= units->torque;
        }
        load = std::make_shared<AppliedForceTorque>(fi, fj, f, t);
    }
    load->setContactStepValidator(contactStepValidator);
    load->setFollower(follower);
    initialize(); // A rerun replaces, rather than appends to, load histories.
    load->name = fullName("");
    mbdObject = load;
    system->addForceTorque(load);
}

void ASMTForceTorque::updateFromMbD()
{
    const auto load = std::static_pointer_cast<AppliedForceTorque>(mbdObject);
    const auto state = load->resultState();
    const auto& f = state.forceOnI;
    const auto& t = state.torqueOnI;
    const auto& energy = state.energy;
    const auto units = root()->mbdUnits;
    fxs->push_back(f[0]*units->force);
    fys->push_back(f[1]*units->force);
    fzs->push_back(f[2]*units->force);
    txs->push_back(t[0]*units->torque);
    tys->push_back(t[1]*units->torque);
    tzs->push_back(t[2]*units->torque);
    const double energyUnit = units->torque / units->angle;
    const double powerUnit = energyUnit / units->time;
    powers->push_back(energy.power * powerUnit);
    storedEnergies->push_back(energy.storedEnergy * energyUnit);
    dissipatedPowers->push_back(energy.dissipatedPower * powerUnit);
}

void ASMTForceTorque::compareResults(AnalysisType) {}
void ASMTForceTorque::outputResults(AnalysisType) {}

bool ASMTForceTorque::isPassive() const
{
    return spring || torsionalSpring || bushing || axisFriction
        || static_cast<bool>(contactEvaluator);
}

void ASMTForceTorque::storeOnLevel(std::ofstream& os, size_t level)
{
    if (contactEvaluator)
        throw std::invalid_argument("Runtime CAD contact evaluators cannot be serialized");
    if (axisFriction && axisFrictionParameters.reactionBased)
        throw std::invalid_argument("Reaction-based joint friction cannot be serialized as a standalone load");
    storeOnLevelString(os, level, "ForceTorque");
    ASMTItemIJ::storeOnLevel(os, level);
    std::ostringstream data;
    data << std::setprecision(17);
    if (axisFriction) {
        data << (axisFrictionParameters.rotational ? "RotationalAxisFriction "
                                                   : "LinearAxisFriction ")
             << axisFrictionParameters.staticMagnitude << ' '
             << axisFrictionParameters.dynamicMagnitude << ' '
             << axisFrictionParameters.transitionVelocity << ' '
             << axisFrictionParameters.viscousCoefficient;
    }
    else if (bushing && bushingParameters.coupled) {
        data << "CoupledBushing";
        for (const auto& values : {bushingParameters.stiffnessMatrix, bushingParameters.dampingMatrix})
            for (double value : values) data << ' ' << value;
    }
    else if (bushing) {
        data << "Bushing";
        for (const auto& values : {bushingParameters.linearStiffness,
                 bushingParameters.linearDamping, bushingParameters.angularStiffness,
                 bushingParameters.angularDamping})
            for (double value : values) data << ' ' << value;
    }
    else if (spring) data << "LinearSpringDamper " << stiffness << ' ' << damping << ' ' << restLength;
    else if (torsionalSpring) data << "TorsionalSpringDamper " << stiffness << ' ' << damping << ' ' << restLength;
    else if (!formula.empty()) data << (follower
        ? (formulaTorque ? "FollowerTorqueFormula " : "FollowerForceFormula ")
        : (formulaTorque ? "WorldTorqueFormula " : "WorldForceFormula "))
        << formulaDirection[0] << ' ' << formulaDirection[1] << ' ' << formulaDirection[2]
        << ' ' << std::quoted(formula);
    else data << (follower ? "FollowerWrench " : "WorldWrench ") << force[0] << ' ' << force[1] << ' ' << force[2]
        << ' ' << torque[0] << ' ' << torque[1] << ' ' << torque[2];
    storeOnLevelString(os, level+1, data.str());
}

void ASMTForceTorque::parseASMT(std::vector<std::string>& lines)
{
    // This explicit format does not interpret FreeCADMbD load expressions.
    if (lines.size() < 7) throw std::invalid_argument("Incomplete typed load");
    ASMTItemIJ::parseASMT(lines);
    std::istringstream data(lines.front());
    std::string kind, extra;
    double x, y, z, tx, ty, tz;
    data >> kind;
    follower = kind.rfind("Follower", 0) == 0;
    if ((kind == "WorldWrench" || kind == "FollowerWrench")
        && (data >> x >> y >> z >> tx >> ty >> tz) && !(data >> extra)) {
        setForce3D(x, y, z);
        setTorque3D(tx, ty, tz);
    }
    else if (kind == "LinearSpringDamper" && (data >> x >> y >> z) && !(data >> extra)) {
        setSpringDamper(x, y, z);
    }
    else if (kind == "TorsionalSpringDamper" && (data >> x >> y >> z) && !(data >> extra)) {
        setTorsionalSpringDamper(x, y, z);
    }
    else if (kind == "CoupledBushing") {
        std::array<double, 36> k{}, c{};
        for (auto* matrix : {&k, &c})
            for (double& value : *matrix)
                if (!(data >> value)) throw std::invalid_argument("Invalid coupled bushing matrix");
        if (data >> extra) throw std::invalid_argument("Invalid coupled bushing specification");
        setCoupledBushing(k, c);
    }
    else if (kind == "Bushing") {
        std::array<double, 3> linearK, linearC, angularK, angularC;
        const auto read = [&](auto& values) {
            for (double& value : values)
                if (!(data >> value)) throw std::invalid_argument("Invalid bushing specification");
        };
        read(linearK);
        read(linearC);
        read(angularK);
        read(angularC);
        data >> std::ws;
        if (!data.eof()) throw std::invalid_argument("Invalid bushing specification");
        setBushing(linearK, linearC, angularK, angularC);
    }
    else if ((kind == "LinearAxisFriction" || kind == "RotationalAxisFriction")
             && (data >> x >> y >> z >> tx) && !(data >> extra)) {
        setAxisFriction(x, y, z, tx, kind == "RotationalAxisFriction");
    }
    else if (kind == "WorldForceFormula" || kind == "WorldTorqueFormula"
             || kind == "FollowerForceFormula" || kind == "FollowerTorqueFormula") {
        if (!(data >> x >> y >> z >> std::quoted(extra)))
            throw std::invalid_argument("Invalid formula load specification");
        data >> std::ws;
        if (!data.eof())
            throw std::invalid_argument("Invalid formula load specification");
        if (kind == "WorldForceFormula" || kind == "FollowerForceFormula") setForceFormula3D(x, y, z, extra);
        else setTorqueFormula3D(x, y, z, extra);
    }
    else throw std::invalid_argument("Invalid or unsupported typed load specification");
    lines.erase(lines.begin());
}
