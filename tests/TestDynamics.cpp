// SPDX-License-Identifier: LGPL-2.1-or-later
// Independent physical benchmarks. All models use metres, kilograms, seconds,
// radians and inertia about the centre of mass. No CAD or fixture files needed.
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include <ASMTAssembly.h>
#include <ASMTPart.h>
#include <ASMTMarker.h>
#include <ASMTRevoluteJoint.h>
#include <ASMTRotationLimit.h>
#include <ASMTTranslationalJoint.h>
#include <ASMTDistanceLimit.h>
#include <ASMTTranslationLimit.h>
#include <ASMTRotationalMotion.h>
#include <System.h>
#include <ExternalSystem.h>
#include <ASMTPrincipalMassMarker.h>
#include <ASMTForceTorque.h>
#include <AppliedForceTorque.h>
#include <Part.h>
#include <PartFrame.h>
#include <EndFrameqc.h>
#include <fstream>
#include <filesystem>
#include <limits>
#include <StableBackwardDifference.h>
#include <SingularMatrixError.h>
#include <DynamicEvents.h>

TEST(DynamicEvents, EarliestBoundaryAndZeroDelayCascade)
{
    MbD::DynamicEvents engine;
    std::vector<int> order;
    MbD::DynamicEvents::Event later, first, chained;
    later.time = .8; later.action = [&](double) { order.push_back(0); };
    first.time = .3; first.action = [&](double) { order.push_back(1); };
    chained.trigger = MbD::DynamicEvents::Event::After;
    chained.predecessor = 1; chained.action = [&](double) { order.push_back(2); };
    engine.events = {later, first, chained};
    EXPECT_FALSE(engine.initialize(0));
    EXPECT_DOUBLE_EQ(engine.locate(0, 1, [](double) {}), .3);
    EXPECT_TRUE(order.empty()); // Probing must not execute an action.
    engine.apply(.3);
    EXPECT_EQ(order, (std::vector<int>{1,2}));
    EXPECT_DOUBLE_EQ(engine.locate(.3, 1, [](double) {}), .8);
    engine.apply(.8);
    EXPECT_EQ(engine.firings.size(), 3);
    EXPECT_FALSE(engine.initialize(0)); // A new run resets all run-local state.
    EXPECT_TRUE(engine.firings.empty());
}

TEST(DynamicEvents, ThresholdRepeatAndHysteresis)
{
    MbD::DynamicEvents engine;
    double time = 0;
    const auto sample = [&](double t) { time = t; };
    MbD::DynamicEvents::Event event;
    event.trigger = MbD::DynamicEvents::Event::Threshold;
    event.repeat = true;
    event.hysteresis = .2;
    event.measure = [&] { return std::sin(time); };
    event.action = [](double) {};
    event.fireInitially = true;
    engine.events = {event};
    EXPECT_TRUE(engine.initialize(0));
    EXPECT_FALSE(std::isfinite(engine.locate(0, 1, sample)));
    engine.commit(0, 1, sample);
    const double crossing = engine.locate(1, 7, sample);
    EXPECT_NEAR(crossing, 2*std::acos(-1), 1e-8);
    engine.commit(1, crossing, sample);
    engine.apply(crossing);
    EXPECT_EQ(engine.firings.size(), 2);
    EXPECT_FALSE(engine.events[0].armed);
}

TEST(DynamicEvents, FallingCrossingAndInvalidInput)
{
    MbD::DynamicEvents engine;
    double time = 0;
    MbD::DynamicEvents::Event event;
    event.trigger = MbD::DynamicEvents::Event::Threshold;
    event.rising = false;
    event.measure = [&] { return 1-time*time; };
    event.action = [](double) {};
    engine.events = {event};
    EXPECT_FALSE(engine.initialize(0));
    EXPECT_NEAR(engine.locate(0, 2, [&](double t) { time=t; }), 1, 1e-8);
    engine.events[0].hysteresis = -1;
    EXPECT_THROW(engine.initialize(0), std::invalid_argument);
}

TEST(DynamicEvents, DelaysAreNotResetByRepeatedPredecessors)
{
    for (bool repeat : {false, true}) {
        MbD::DynamicEvents engine;
        double time = 0;
        const auto sample = [&](double t) { time = t; };
        MbD::DynamicEvents::Event crossing, delayed;
        crossing.trigger = MbD::DynamicEvents::Event::Threshold;
        crossing.repeat = true;
        crossing.hysteresis = .1;
        crossing.measure = [&] { return std::sin(2*std::acos(-1)*time); };
        crossing.action = [](double) {};
        delayed.trigger = MbD::DynamicEvents::Event::After;
        delayed.predecessor = 0;
        delayed.delay = 1.5;
        delayed.repeat = repeat;
        delayed.action = [](double) {};
        engine.events = {crossing, delayed};
        engine.initialize(0);
        double previous = 0;
        while (previous < 3.6) {
            const double end = std::min(previous+.1, 3.6);
            const double firing = engine.locate(previous, end, sample);
            const double accepted = std::isfinite(firing) ? firing : end;
            engine.commit(previous, accepted, sample);
            if (std::isfinite(firing)) engine.apply(firing);
            previous = accepted;
        }
        std::vector<double> times;
        for (const auto& firing : engine.firings)
            if (firing.event == 1) times.push_back(firing.time);
        ASSERT_EQ(times.size(), repeat ? 2u : 1u);
        EXPECT_NEAR(times[0], 2.5, 1e-7);
        if (repeat) EXPECT_NEAR(times[1], 3.5, 1e-7);
    }
}

namespace {
using namespace MbD;
using Assembly = std::shared_ptr<ASMTAssembly>;
using Body = std::shared_ptr<ASMTPart>;
using Marker = std::shared_ptr<ASMTMarker>;

std::string number(double x)
{
    std::ostringstream s;
    s << std::setprecision(17) << x;
    return s.str();
}

Assembly model(double end = 1.0, double tolerance = 1e-8, double output = 0.02)
{
    auto a = ASMTAssembly::With();
    a->setName("Benchmark");
    a->setVelocity3D(0, 0, 0);
    a->setOmega3D(0, 0, 0);
    auto g = ASMTConstantGravity::With();
    g->setg(0, 0, 0);
    a->setConstantGravity(g);
    auto p = ASMTSimulationParameters::With();
    p->settstart(0);
    p->settend(end);
    p->sethout(output);
    p->sethmax(0.02);
    p->sethmin(1e-12);
    p->seterrorTol(tolerance);
    p->errorTolPosKine = p->errorTolAccKine = tolerance;
    p->corAbsTol = p->corRelTol = p->intAbsTol = p->intRelTol = tolerance;
    p->iterMaxDyn = 25;
    a->setSimulationParameters(p);
    return a;
}

Body body(const Assembly& a, const std::string& name, double mass, double inertia,
          double x = 0, double y = 0, double angle = 0)
{
    auto b = ASMTPart::With();
    b->setName(name);
    b->setPosition3D(x, y, 0);
    const double c = std::cos(angle), s = std::sin(angle);
    b->setRotationMatrix(c, -s, 0, s, c, 0, 0, 0, 1);
    b->setVelocity3D(0, 0, 0);
    b->setOmega3D(0, 0, 0);
    auto mm = ASMTPrincipalMassMarker::With();
    mm->setMass(mass);
    mm->setDensity(1000);
    mm->setMomentOfInertias(inertia, inertia, inertia);
    mm->setPosition3D(0, 0, 0);
    b->setPrincipalMassMarker(mm);
    a->addPart(b);
    return b;
}

Marker marker(const std::shared_ptr<ASMTSpatialContainer>& b, const std::string& name,
              double x = 0, double y = 0, bool alongX = false)
{
    auto m = ASMTMarker::With();
    m->setName(name);
    m->setPosition3D(x, y, 0);
    // Translational joints slide along marker z. Rotate it onto global x.
    if (alongX) m->setRotationMatrix(0, 0, 1, 0, 1, 0, -1, 0, 0);
    b->addMarker(m);
    return m;
}

template<class T>
std::shared_ptr<T> joint(const Assembly& a, const std::string& name,
                         const Marker& i, const Marker& j)
{
    auto q = T::With();
    q->setName(name);
    q->setMarkerI(i->fullName(""));
    q->setMarkerJ(j->fullName(""));
    a->addJoint(q);
    return q;
}

std::shared_ptr<ASMTRotationalMotion> drive(const Assembly& a,
    const std::shared_ptr<ASMTRevoluteJoint>& j, const std::string& expression)
{
    auto m = ASMTRotationalMotion::With();
    m->setName("Driver");
    m->setMotionJoint(j->fullName(""));
    m->setRotationZ(expression);
    a->addMotion(m);
    return m;
}

void run(const Assembly& a, bool dynamic)
{
    if (dynamic) {
        a->runDYNAMIC();
        return;
    }
    // Preserve the legacy kinematic API, but do not swallow its errors here.
    a->mbdSystem = std::make_shared<System>();
    a->mbdSystem->externalSystem->asmtAssembly = a.get();
    a->mbdSystem->runKINEMATIC(a->mbdSystem);
}

TEST(Dynamics, FixedBodyPreservesInputPoseIncludingHalfTurns)
{
    for (const double angle : {0.0, 0.5*std::acos(-1.0), std::acos(-1.0), 2.0*std::acos(-1.0)/3.0}) {
        SCOPED_TRACE(angle);
        auto a = model(0.1);
        a->constantGravity->setg(0, 0, -9.81);
        auto ground = body(a, "Ground", 2, 0.3, 2, -3);
        const double c = std::cos(angle), s = std::sin(angle);
        ground->setPosition3D(2, -3, 4);
        ground->setRotationMatrix(1, 0, 0, 0, c, -s, 0, s, c);
        ground->principalMassMarker->setPosition3D(0.2, -0.1, 0.4);
        ground->isFixed = true;
        auto falling = body(a, "Falling", 1, 0.1);
        a->runDYNAMIC();
        for (size_t k = 0; k < a->times->size(); ++k) {
            EXPECT_NEAR(ground->xs->at(k), 2, 1e-9);
            EXPECT_NEAR(ground->ys->at(k), -3, 1e-9);
            EXPECT_NEAR(ground->zs->at(k), 4, 1e-9);
            auto rotation = ground->getRotationMatrix(k);
            const double expected[3][3] = {{1,0,0}, {0,c,-s}, {0,s,c}};
            for (size_t i=0; i<3; ++i)
                for (size_t j=0; j<3; ++j)
                    EXPECT_NEAR(rotation->at(i)->at(j), expected[i][j], 1e-9);
            const double t = a->times->at(k);
            EXPECT_NEAR(falling->zs->at(k), -4.905*t*t, 1e-6);
        }
    }
}

TEST(Dynamics, SmoothForceStartupDoesNotOscillateFirstStepSize)
{
    auto a = model(0.1, 1e-10, 0.01);
    a->simulationParameters->sethmax(0.005);
    auto b = body(a, "Body", 1, 1);
    auto i = marker(b, "I"), j = marker(a, "J");
    auto load = ASMTForceTorque::With();
    load->setName("Smooth force");
    load->setMarkerI(i->fullName(""));
    load->setMarkerJ(j->fullName(""));
    load->setForceFormula3D(0, 0, 1,
        "piecewise(time,functions(10*(10*(time/0.02)^3-15*(time/0.02)^4+6*(time/0.02)^5),10),transitions(0.02))");
    a->addForceTorque(load);
    a->runDYNAMIC();
    EXPECT_NEAR(b->vzs->back(), 0.9, 1e-5);
}

// Validate the complete output grid, not only the final pose. Backends may emit
// an extra input/initial-condition sample at t=0; it is not a solved sample.
void complete(const Assembly& a, const Body& b)
{
    ASSERT_TRUE(a->times);
    ASSERT_TRUE(b->xs);
    ASSERT_EQ(a->times->size(), b->xs->size());
    const auto& p = *a->simulationParameters;
    const size_t n = static_cast<size_t>(std::llround(p.tend / p.hout));
    ASSERT_GE(a->times->size(), n + 1);
    EXPECT_NEAR(a->times->back(), p.tend, 1e-9);
    size_t next = 1;
    double previous = -1;
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        ASSERT_TRUE(std::isfinite(t));
        EXPECT_GE(t, previous);
        previous = t;
        if (t <= 1e-12) continue;
        EXPECT_NEAR(t, next * p.hout, 1e-9);
        ++next;
        for (const auto& values : {b->xs, b->ys, b->zs, b->vxs, b->vys, b->vzs,
                                  b->axs, b->ays, b->azs, b->bryzs,
                                  b->omexs, b->omeys, b->omezs,
                                  b->alpxs, b->alpys, b->alpzs}) {
            ASSERT_TRUE(values);
            ASSERT_EQ(values->size(), a->times->size());
            EXPECT_TRUE(std::isfinite(values->at(k)));
        }
    }
    EXPECT_EQ(next, n + 1);
}

double reaction(const std::shared_ptr<ASMTItemIJ>& j, size_t k, bool torque, size_t axis)
{
    const auto values = torque ? (axis == 0 ? j->txs : axis == 1 ? j->tys : j->tzs)
                               : (axis == 0 ? j->fxs : axis == 1 ? j->fys : j->fzs);
    if (!values) throw std::runtime_error("Missing reaction output");
    const double value = values->at(k);
    if (!std::isfinite(value)) throw std::runtime_error("Non-finite reaction output");
    return value;
}


TEST(Dynamics, FreeFall)
{
    auto a = model();
    auto b = body(a, "FallingBody", 2, 0.3, 0, 2);
    b->setVelocity3D(0.4, 0.2, 0);
    a->constantGravity->setg(0, -9.81, 0);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        EXPECT_NEAR(b->xs->at(k), 0.4*t, 2e-6);
        EXPECT_NEAR(b->ys->at(k), 2 + 0.2*t - 4.905*t*t, 2e-6);
        EXPECT_NEAR(b->vys->at(k), 0.2 - 9.81*t, 2e-6);
        EXPECT_NEAR(b->ays->at(k), -9.81, 2e-6);
        const double energy = b->vxs->at(k)*b->vxs->at(k)
            + b->vys->at(k)*b->vys->at(k) + 2*9.81*b->ys->at(k);
        EXPECT_NEAR(energy, 0.4*0.4 + 0.2*0.2 + 2*9.81*2, 2e-5);
    }
}

TEST(Dynamics, FirstOrderFreeFallReportsAcceleration)
{
    auto a = model(0.1, 1e-8, 0.01);
    a->simulationParameters->orderMax = 1;
    auto b = body(a, "FallingBody", 2, 0.3, 0, 2);
    a->constantGravity->setg(0, -9.81, 0);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    for (size_t k = 1; k < a->times->size(); ++k) {
        EXPECT_NEAR(b->ays->at(k), -9.81, 2e-6);
    }
}

TEST(Dynamics, TranslationalJointLimitStopsInitialVelocity)
{
    auto a = model(0.1, 1e-8, 0.01);
    auto b = body(a, "Slider", 1, 0.1);
    b->setVelocity3D(-1, 0, 0);
    auto groundMarker = marker(a, "Guide", 0, 0, true);
    auto bodyMarker = marker(b, "Guide", 0, 0, true);
    joint<ASMTTranslationalJoint>(a, "SliderJoint", groundMarker, bodyMarker);

    auto limit = ASMTTranslationLimit::With();
    limit->setName("LowerStop");
    limit->setMarkerI(groundMarker->fullName(""));
    limit->setMarkerJ(bodyMarker->fullName(""));
    limit->settype("=>");
    limit->setlimit("-0.02");
    limit->settol("1e-9");
    a->addLimit(limit);

    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    EXPECT_NEAR(b->xs->back(), -0.02, 1e-7);
    EXPECT_NEAR(b->vxs->back(), 0, 1e-7);
    EXPECT_GE(*std::min_element(b->xs->begin(), b->xs->end()), -0.0200001);
}

TEST(Dynamics, RotationalJointLimitsStopInitialVelocity)
{
    auto a = model(0.1, 1e-8, 0.01);
    auto b = body(a, "Rotor", 1, 0.1);
    b->setOmega3D(0, 0, 10);
    auto groundMarker = marker(a, "Pivot");
    auto bodyMarker = marker(b, "Pivot");
    joint<ASMTRevoluteJoint>(a, "Hinge", groundMarker, bodyMarker);

    for (const auto& [name, comparison, value] : {
             std::tuple{"LowerStop", "=>", "-0.5"},
             std::tuple{"UpperStop", "=<", "0.5"},
         }) {
        auto limit = ASMTRotationLimit::With();
        limit->setName(name);
        limit->setMarkerI(groundMarker->fullName(""));
        limit->setMarkerJ(bodyMarker->fullName(""));
        limit->settype(comparison);
        limit->setlimit(value);
        limit->settol("1e-9");
        a->addLimit(limit);
    }

    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    EXPECT_NEAR(std::abs(b->bryzs->back()), 0.5, 1e-7);
    EXPECT_NEAR(b->omezs->back(), 0, 1e-7);
}

TEST(Dynamics, HardStopsReleaseButDoNotPull)
{
    for (bool rotational : {false, true}) {
        for (double sign : {-1.0, 1.0}) {
            auto a = model(.1, 1e-9, .01);
            a->simulationParameters->sethmax(.0005);
            auto b = body(a, "Body", 1, 1);
            auto i = marker(a, "Guide"), j = marker(b, "Guide");
            std::shared_ptr<ASMTLimit> limit;
            if (rotational) {
                joint<ASMTRevoluteJoint>(a, "Hinge", i, j);
                b->setOmega3D(0, 0, sign*10);
                limit = ASMTRotationLimit::With();
            }
            else {
                joint<ASMTTranslationalJoint>(a, "Slider", i, j);
                b->setVelocity3D(0, 0, sign*10);
                limit = ASMTTranslationLimit::With();
            }
            limit->setName("Stop");
            limit->setMarkerI(i->fullName(""));
            limit->setMarkerJ(j->fullName(""));
            limit->settype(sign > 0 ? "=<" : "=>");
            limit->setlimit(sign > 0 ? "0.35" : "-0.35");
            limit->settol("1e-9");
            a->addLimit(limit);
            auto load = ASMTForceTorque::With();
            load->setName("SeparatingLoad");
            load->setMarkerI(j->fullName(""));
            load->setMarkerJ(i->fullName(""));
            if (rotational) load->setTorque3D(0, 0, -sign*100);
            else load->setForce3D(0, 0, -sign*100);
            a->addForceTorque(load);
            run(a, true);
            ASSERT_NO_FATAL_FAILURE(complete(a, b));
            EXPECT_NEAR(rotational ? b->bryzs->back() : b->zs->back(), sign*.2, .0001);
            EXPECT_NEAR(rotational ? b->omezs->back() : b->vzs->back(),
                        -sign*std::sqrt(30), .001);
        }
    }
}

TEST(Dynamics, CompliantTranslationalLimitReboundsAndReleases)
{
    auto a = model(0.16, 1e-8, 0.002);
    a->simulationParameters->sethmax(0.002);
    auto b = body(a, "Slider", 1, 0.1);
    b->setVelocity3D(-1, 0, 0);
    auto groundMarker = marker(a, "Guide", 0, 0, true);
    auto bodyMarker = marker(b, "Guide", 0, 0, true);
    joint<ASMTTranslationalJoint>(a, "SliderJoint", groundMarker, bodyMarker);

    auto limit = ASMTTranslationLimit::With();
    limit->setName("CompliantLowerStop");
    limit->setMarkerI(groundMarker->fullName(""));
    limit->setMarkerJ(bodyMarker->fullName(""));
    limit->settype("=>");
    limit->setlimit("-0.02");
    limit->settol("1e-9");
    limit->setcompliance("1000", "0");
    a->addLimit(limit);

    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    const auto minimum = *std::min_element(b->xs->begin(), b->xs->end());
    EXPECT_NEAR(minimum, -0.02 - 1.0/std::sqrt(1000.0), 8e-4);
    EXPECT_GT(b->vxs->back(), 0.98);
    EXPECT_GT(b->xs->back(), -0.019);
    EXPECT_GT(*std::max_element(limit->storedEnergies->begin(), limit->storedEnergies->end()), 0.45);
    for (size_t k = 1; k < a->times->size(); ++k) {
        EXPECT_NEAR(0.5*b->vxs->at(k)*b->vxs->at(k) + limit->storedEnergies->at(k), 0.5, 2e-3);
        EXPECT_DOUBLE_EQ(limit->dissipatedPowers->at(k), 0);
    }
}

TEST(Dynamics, CompliantRotationalLimitReboundsAndReleases)
{
    auto a = model(0.5, 1e-8, 0.005);
    a->simulationParameters->sethmax(0.002);
    auto b = body(a, "Rotor", 1, 0.1);
    b->setOmega3D(0, 0, 5);
    auto groundMarker = marker(a, "Pivot");
    auto bodyMarker = marker(b, "Pivot");
    joint<ASMTRevoluteJoint>(a, "Hinge", groundMarker, bodyMarker);

    auto limit = ASMTRotationLimit::With();
    limit->setName("CompliantUpperStop");
    limit->setMarkerI(groundMarker->fullName(""));
    limit->setMarkerJ(bodyMarker->fullName(""));
    limit->settype("=<");
    limit->setlimit("0.2");
    limit->settol("1e-9");
    limit->setcompliance("10", "0");
    a->addLimit(limit);

    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    EXPECT_GT(*std::max_element(b->bryzs->begin(), b->bryzs->end()), 0.65);
    EXPECT_LT(b->omezs->back(), -4.8);
    EXPECT_LT(b->bryzs->back(), 0.19);
}

TEST(Dynamics, CompliantDistanceContactReboundsAndReleases)
{
    auto a = model(1.0, 1e-9, 0.005);
    a->simulationParameters->sethmax(0.001);
    auto first = body(a, "FirstBall", 1, 0.1, -1, 0);
    auto second = body(a, "SecondBall", 1, 0.1, 1, 0);
    first->setVelocity3D(1, 0, 0);
    second->setVelocity3D(-1, 0, 0);

    auto contact = ASMTDistanceLimit::With();
    contact->setName("SphereContact");
    contact->setMarkerI(marker(first, "Centre")->fullName(""));
    contact->setMarkerJ(marker(second, "Centre")->fullName(""));
    contact->settype("=>");
    contact->setlimit("1.0");
    contact->settol("1.0e-9");
    // Nonzero damping exercises the distance constraint's velocity Jacobian,
    // not only the elastic position Jacobian used by the undamped case.
    contact->setcompliance("1000", "5");
    a->addLimit(contact);

    run(a, true);
    double minimumDistance = std::numeric_limits<double>::infinity();
    for (size_t index = 1; index < a->times->size(); ++index) {
        minimumDistance = std::min(
            minimumDistance,
            std::abs(second->xs->at(index) - first->xs->at(index))
        );
    }
    EXPECT_LT(minimumDistance, 1.0);
    EXPECT_LT(first->vxs->back(), 0.0);
    EXPECT_GT(second->vxs->back(), 0.0);
    EXPECT_LT(std::abs(first->vxs->back()), 1.0);
    EXPECT_LT(std::abs(second->vxs->back()), 1.0);
    EXPECT_GT(second->xs->back() - first->xs->back(), 1.0);
    EXPECT_GT(*std::max_element(contact->storedEnergies->begin(), contact->storedEnergies->end()), 0);
    EXPECT_GT(*std::max_element(contact->dissipatedPowers->begin(), contact->dissipatedPowers->end()), 0);
}

TEST(Dynamics, RigidDistanceContactPushesAnUndraggedBody)
{
    auto a = model();
    auto dragged = body(a, "DraggedBall", 1, 0.1, -1, 0);
    auto pushed = body(a, "PushedBall", 1, 0.1, 1, 0);

    auto contact = ASMTDistanceLimit::With();
    contact->setName("SphereContact");
    contact->setMarkerI(marker(dragged, "Centre")->fullName(""));
    contact->setMarkerJ(marker(pushed, "Centre")->fullName(""));
    contact->settype("=>");
    contact->setlimit("2.0");
    contact->settol("1.0e-9");
    a->addLimit(contact);

    a->runPreDrag();
    dragged->setPosition3D(0, 0, 0);
    auto dragParts = std::make_shared<std::vector<Body>>(std::initializer_list<Body>{dragged});
    a->runDragStep(dragParts);

    double draggedX, draggedY, draggedZ;
    double pushedX, pushedY, pushedZ;
    dragged->getPosition3D(draggedX, draggedY, draggedZ);
    pushed->getPosition3D(pushedX, pushedY, pushedZ);
    EXPECT_NEAR(draggedX, 0.0, 2.0e-3);
    EXPECT_NEAR(pushedX - draggedX, 2.0, 1.0e-8);
}

TEST(Dynamics, CompliantLimitDampingReducesReboundWithoutSticking)
{
    auto a = model(0.2, 1e-8, 0.002);
    a->simulationParameters->sethmax(0.002);
    auto b = body(a, "Slider", 1, 0.1);
    b->setVelocity3D(-1, 0, 0);
    auto groundMarker = marker(a, "Guide", 0, 0, true);
    auto bodyMarker = marker(b, "Guide", 0, 0, true);
    joint<ASMTTranslationalJoint>(a, "SliderJoint", groundMarker, bodyMarker);

    auto limit = ASMTTranslationLimit::With();
    limit->setName("DampedLowerStop");
    limit->setMarkerI(groundMarker->fullName(""));
    limit->setMarkerJ(bodyMarker->fullName(""));
    limit->settype("=>");
    limit->setlimit("-0.02");
    limit->settol("1e-9");
    limit->setcompliance("1000", "10");
    a->addLimit(limit);

    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    EXPECT_GT(b->vxs->back(), 0.1);
    EXPECT_LT(b->vxs->back(), 0.95);
    EXPECT_GT(b->xs->back(), -0.019);
}

// Independent RK4 reference for a finite-amplitude *physical* pendulum.
// It does not use the solver's ExactPendulum helper or small-angle approximation.
struct PendulumState { double theta, omega; };
PendulumState reference(double t, double step)
{
    PendulumState q {0.5, 0};
    const double frequencySquared = 2*9.81 / (0.2 + 2*1.0*1.0);
    auto acc = [=](double x) { return -frequencySquared*std::sin(x); };
    const int n = std::max(1, static_cast<int>(std::ceil(t / step)));
    const double h = t/n;
    for (int i = 0; i < n; ++i) {
        const double a1 = acc(q.theta), v1 = q.omega;
        const double a2 = acc(q.theta + h*v1/2), v2 = q.omega + h*a1/2;
        const double a3 = acc(q.theta + h*v2/2), v3 = q.omega + h*a2/2;
        const double a4 = acc(q.theta + h*v3), v4 = q.omega + h*a3;
        q.theta += h*(v1 + 2*v2 + 2*v3 + v4)/6;
        q.omega += h*(a1 + 2*a2 + 2*a3 + a4)/6;
    }
    return q;
}

struct PendulumError { double angle = 0, energy = 0, balance = 0; };
PendulumError pendulum(double tolerance, double output, double maxStep = 0.02)
{
    auto a = model(2.0, tolerance, output);
    a->simulationParameters->sethmax(maxStep);
    auto b = body(a, "Pendulum", 2, 0.2, std::sin(0.5), -std::cos(0.5), 0.5);
    auto j = joint<ASMTRevoluteJoint>(a, "Hinge", marker(a, "Pivot"), marker(b, "Pivot", 0, 1));
    a->constantGravity->setg(0, -9.81, 0);
    run(a, true);
    complete(a, b);
    if (::testing::Test::HasFatalFailure()) throw std::runtime_error("Incomplete pendulum output");
    PendulumError error;
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        const auto q = reference(t, 0.0005);
        const double theta = b->bryzs->at(k), omega = b->omezs->at(k);
        error.angle = std::max(error.angle, std::abs(theta-q.theta));
        const double energy = 0.5*2.2*omega*omega + 2*9.81*(1-std::cos(theta));
        error.energy = std::max(error.energy, std::abs(energy-2*9.81*(1-std::cos(0.5))));
        EXPECT_NEAR(b->xs->at(k), std::sin(theta), 2e-5);
        EXPECT_NEAR(b->ys->at(k), -std::cos(theta), 2e-5);
        // Ground-side hinge reaction is the negative of force on the body.
        error.balance = std::max({error.balance,
            std::abs(reaction(j, k, false, 0) + 2*b->axs->at(k)),
            std::abs(reaction(j, k, false, 1) + 2*(b->ays->at(k)+9.81))});
    }
    return error;
}

TEST(Dynamics, PendulumReferenceConvergence)
{
    for (double t : {0.1, 0.5, 1.0, 2.0}) {
        const auto coarse = reference(t, 0.0005), fine = reference(t, 0.00025);
        EXPECT_NEAR(coarse.theta, fine.theta, 1e-11);
        EXPECT_NEAR(coarse.omega, fine.omega, 1e-11);
    }
}

TEST(Dynamics, PendulumAccuracyEnergyAndToleranceConvergence)
{
    const auto coarse = pendulum(1e-5, 0.02);
    const auto fine = pendulum(1e-9, 0.02);
    const auto resampled = pendulum(1e-9, 0.01);
    RecordProperty("coarse_max_angle_error_rad", number(coarse.angle));
    RecordProperty("fine_max_angle_error_rad", number(fine.angle));
    RecordProperty("fine_max_energy_error_J", number(fine.energy));
    RecordProperty("coarse_max_energy_error_J", number(coarse.energy));
    RecordProperty("resampled_max_angle_error_rad", number(resampled.angle));
    EXPECT_LT(fine.angle, 2e-5);
    EXPECT_LT(fine.energy, 2e-4);
    EXPECT_LE(fine.angle, std::max(2e-7, coarse.angle*0.5));
    EXPECT_LE(fine.energy, std::max(2e-6, coarse.energy*0.5));
    EXPECT_LT(resampled.angle, 2e-5);
}

TEST(Dynamics, PendulumReactionBalanceAndStepConvergence)
{
    const auto normal = pendulum(1e-9, 0.02);
    const auto smallerStep = pendulum(1e-9, 0.02, 0.005);
    RecordProperty("max_force_balance_residual_N", number(normal.balance));
    RecordProperty("small_step_force_balance_residual_N", number(smallerStep.balance));
    // Qualification uses the refined integration step. The coarser run is a
    // convergence probe, not an accepted setting: at 20 ms it exceeds this
    // force-balance bound in FreeCADMbD 8c5cf85. Output spacing remains 20 ms.
    EXPECT_LT(smallerStep.balance, 2e-4);
    EXPECT_LE(smallerStep.balance, std::max(1e-6, normal.balance*0.5));
}

TEST(Dynamics, TorqueDrivenRotor)
{
    auto a = model();
    auto b = body(a, "Rotor", 2, 0.4);
    auto i = marker(a, "Pivot"), j = marker(b, "Pivot");
    joint<ASMTRevoluteJoint>(a, "Hinge", i, j);
    auto load = ASMTForceTorque::With();
    load->setName("Torque");
    load->setMarkerI(i->fullName(""));
    load->setMarkerJ(j->fullName(""));
    load->setTorque3D(0, 0, -0.8);
    a->addForceTorque(load);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        EXPECT_NEAR(b->bryzs->at(k), t*t, 2e-5);
        EXPECT_NEAR(b->omezs->at(k), 2*t, 2e-5);
        EXPECT_NEAR(b->alpzs->at(k), 2, 2e-4);
        EXPECT_NEAR(0.2*std::pow(b->omezs->at(k), 2), 0.8*b->bryzs->at(k), 2e-5);
        EXPECT_NEAR(reaction(load, k, true, 2), -0.8, 1e-12);
        EXPECT_NEAR(load->powers->at(k), 0.8*b->omezs->at(k), 2e-5);
        EXPECT_DOUBLE_EQ(load->storedEnergies->at(k), 0);
        EXPECT_DOUBLE_EQ(load->dissipatedPowers->at(k), 0);
    }
}

TEST(Dynamics, SpringMassDamper)
{
    auto a = model(2, 1e-9);
    auto b = body(a, "Mass", 2, 0.1, 1.2);
    joint<ASMTTranslationalJoint>(a, "Guide", marker(a, "Guide", 0, 0, true),
        marker(b, "Guide", 0, 0, true));
    auto load = ASMTForceTorque::With();
    load->setName("Spring");
    load->setMarkerI(marker(a, "Anchor")->fullName(""));
    load->setMarkerJ(marker(b, "Anchor")->fullName(""));
    load->setSpringDamper(18, 1.2, 1);
    a->addForceTorque(load);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    double previousEnergy = 0.36;
    const double decay = 0.3, wd = std::sqrt(9-decay*decay);
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        const double x = 1+0.2*std::exp(-decay*t)*(std::cos(wd*t)+decay/wd*std::sin(wd*t));
        const double v = -0.2*std::exp(-decay*t)*(wd+decay*decay/wd)*std::sin(wd*t);
        EXPECT_NEAR(b->xs->at(k), x, 2e-5);
        EXPECT_NEAR(b->vxs->at(k), v, 2e-5);
        EXPECT_NEAR(b->axs->at(k), -9*(x-1)-0.6*v, 2e-4);
        EXPECT_NEAR(reaction(load, k, false, 0), 18*(b->xs->at(k)-1)+1.2*b->vxs->at(k), 2e-5);
        EXPECT_NEAR(load->storedEnergies->at(k), 9*std::pow(x-1, 2), 2e-5);
        EXPECT_NEAR(load->dissipatedPowers->at(k), 1.2*v*v, 2e-5);
        EXPECT_NEAR(load->powers->at(k), -(18*(x-1)+1.2*v)*v, 2e-5);
        const double energy = std::pow(b->vxs->at(k), 2)+9*std::pow(b->xs->at(k)-1, 2);
        EXPECT_LE(energy, previousEnergy+2e-6);
        previousEnergy = energy;
    }
}

TEST(Dynamics, TorsionalSpringTracksMultipleTurnsInBothDirections)
{
    for (const double speed : {-20.0, 20.0}) {
        SCOPED_TRACE(speed);
        auto a = model(1, 1e-9, 0.2);
        auto b = body(a, "Rotor", 2, 0.4);
        auto hinge = joint<ASMTRevoluteJoint>(a, "Hinge",
            marker(a, "Pivot"), marker(b, "Pivot"));
        drive(a, hinge, number(speed) + "*time");
        auto load = ASMTForceTorque::With();
        load->setName("Spring");
        load->setMarkerI(marker(a, "SpringGround")->fullName(""));
        load->setMarkerJ(marker(b, "SpringBody")->fullName(""));
        load->setTorsionalSpringDamper(0.1, 0.01, 0);
        a->addForceTorque(load);
        run(a, true);
        ASSERT_NO_FATAL_FAILURE(complete(a, b));
        for (size_t k = 0; k < a->times->size(); ++k) {
            const double t = a->times->at(k);
            if (t <= 1e-12) continue;
            EXPECT_NEAR(reaction(load, k, true, 2), 0.1 * speed * t + 0.01 * speed, 1e-6);
        }
    }
}

TEST(Dynamics, CoupledBushingRejectsActiveAndInvalidMatrices)
{
    std::array<double, 36> k{}, c{};
    k[0] = 2;
    k[35] = 3;
    k[5] = k[30] = 0.5;
    auto load = ASMTForceTorque::With();
    EXPECT_NO_THROW(load->setCoupledBushing(k, c));
    k[5] = k[30] = 3;
    EXPECT_THROW(load->setCoupledBushing(k, c), std::invalid_argument);
    k[5] = 0;
    EXPECT_THROW(load->setCoupledBushing(k, c), std::invalid_argument);
    k[30] = 0;
    c[7] = -1;
    EXPECT_THROW(load->setCoupledBushing(k, c), std::invalid_argument);
    c[7] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(load->setCoupledBushing(k, c), std::invalid_argument);
}

TEST(Dynamics, CoupledBushingTranslationMatchesDampedOscillator)
{
    auto a = model(0.5, 1e-9);
    auto b = body(a, "Mass", 2, 0.1, 1.2);
    joint<ASMTTranslationalJoint>(a, "Guide", marker(a, "Guide", 0, 0, true),
        marker(b, "Guide", 0, 0, true));
    auto load = ASMTForceTorque::With();
    load->setName("CoupledBushing");
    load->setMarkerI(marker(a, "BushingGround", 1)->fullName(""));
    load->setMarkerJ(marker(b, "BushingBody")->fullName(""));
    std::array<double, 36> k{}, c{};
    k[0] = 18;
    k[35] = 2;
    k[5] = k[30] = 1;
    c[0] = 1.2;
    load->setCoupledBushing(k, c);
    a->addForceTorque(load);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    const double decay = 0.3, wd = std::sqrt(9-decay*decay);
    for (size_t i = 0; i < a->times->size(); ++i) {
        const double t = a->times->at(i);
        if (t <= 1e-12) continue;
        const double x = 1 + 0.2*std::exp(-decay*t)
            *(std::cos(wd*t)+decay/wd*std::sin(wd*t));
        EXPECT_NEAR(b->xs->at(i), x, 2e-5);
        // The guide suppresses rotation, but the cross stiffness still creates
        // a torque that must be balanced by the guide reaction.
        EXPECT_NEAR(reaction(load, i, true, 2), b->xs->at(i)-1, 2e-5);
    }
}

TEST(Dynamics, BushingTranslationMatchesDampedOscillator)
{
    auto a = model(2, 1e-9);
    auto b = body(a, "Mass", 2, 0.1, 1.2);
    joint<ASMTTranslationalJoint>(a, "Guide", marker(a, "Guide", 0, 0, true),
        marker(b, "Guide", 0, 0, true));
    auto load = ASMTForceTorque::With();
    load->setName("Bushing");
    load->setMarkerI(marker(a, "BushingGround", 1)->fullName(""));
    load->setMarkerJ(marker(b, "BushingBody")->fullName(""));
    load->setBushing({18, 0, 0}, {1.2, 0, 0}, {0, 0, 0}, {0, 0, 0});
    a->addForceTorque(load);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    double previousEnergy = 0.36;
    const double decay = 0.3, wd = std::sqrt(9-decay*decay);
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        const double x = 1+0.2*std::exp(-decay*t)*(std::cos(wd*t)+decay/wd*std::sin(wd*t));
        const double v = -0.2*std::exp(-decay*t)*(wd+decay*decay/wd)*std::sin(wd*t);
        EXPECT_NEAR(b->xs->at(k), x, 2e-5);
        EXPECT_NEAR(b->vxs->at(k), v, 2e-5);
        EXPECT_NEAR(b->axs->at(k), -9*(x-1)-0.6*v, 2e-4);
        EXPECT_NEAR(reaction(load, k, false, 0), 18*(b->xs->at(k)-1)+1.2*b->vxs->at(k), 2e-5);
        const double energy = std::pow(b->vxs->at(k), 2)+9*std::pow(b->xs->at(k)-1, 2);
        EXPECT_LE(energy, previousEnergy+2e-6);
        previousEnergy = energy;
    }
}

TEST(Dynamics, BushingRotationMatchesDampedOscillator)
{
    auto a = model(2, 1e-9);
    auto b = body(a, "Rotor", 2, 0.4, 0, 0, 0.2);
    joint<ASMTRevoluteJoint>(a, "Hinge", marker(a, "Pivot"), marker(b, "Pivot"));
    auto load = ASMTForceTorque::With();
    load->setName("Bushing");
    load->setMarkerI(marker(a, "BushingGround")->fullName(""));
    load->setMarkerJ(marker(b, "BushingBody")->fullName(""));
    load->setBushing({0, 0, 0}, {0, 0, 0}, {0, 0, 3.6}, {0, 0, 0.24});
    a->addForceTorque(load);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    double previousEnergy = 0.072;
    const double decay = 0.3, wd = std::sqrt(9-decay*decay);
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        const double angle = 0.2*std::exp(-decay*t)
            *(std::cos(wd*t)+decay/wd*std::sin(wd*t));
        const double omega = -0.2*std::exp(-decay*t)
            *(wd+decay*decay/wd)*std::sin(wd*t);
        EXPECT_NEAR(b->bryzs->at(k), angle, 2e-5);
        EXPECT_NEAR(b->omezs->at(k), omega, 2e-5);
        EXPECT_NEAR(b->alpzs->at(k), -9*angle-0.6*omega, 2e-4);
        EXPECT_NEAR(reaction(load, k, true, 2), 3.6*angle+0.24*omega, 2e-5);
        const double energy = 0.2*omega*omega+1.8*angle*angle;
        EXPECT_LE(energy, previousEnergy+2e-6);
        previousEnergy = energy;
    }
}

TEST(Dynamics, RotationalAxisFrictionMatchesCoulombViscousDecay)
{
    auto a = model(0.2, 1e-9, 0.01);
    auto b = body(a, "Rotor", 2, 0.4);
    b->setOmega3D(0, 0, 1);
    auto i = marker(a, "Pivot"), j = marker(b, "Pivot");
    joint<ASMTRevoluteJoint>(a, "Hinge", i, j);
    auto load = ASMTForceTorque::With();
    load->setName("JointFriction");
    load->setMarkerI(i->fullName(""));
    load->setMarkerJ(j->fullName(""));
    load->setAxisFriction(0.8, 0.4, 1e-6, 0.2, true);
    a->addForceTorque(load);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        const double omega = 3*std::exp(-0.5*t)-2;
        const double angle = 6*(1-std::exp(-0.5*t))-2*t;
        EXPECT_NEAR(b->bryzs->at(k), angle, 2e-5);
        EXPECT_NEAR(b->omezs->at(k), omega, 2e-5);
        EXPECT_NEAR(b->alpzs->at(k), -1-0.5*omega, 2e-4);
        EXPECT_NEAR(reaction(load, k, true, 2), 0.4+0.2*omega, 2e-5);
        EXPECT_NEAR(load->powers->at(k), -(0.4+0.2*omega)*omega, 2e-5);
        EXPECT_NEAR(load->dissipatedPowers->at(k), (0.4+0.2*omega)*omega, 2e-5);
        EXPECT_DOUBLE_EQ(load->storedEnergies->at(k), 0);
    }
}

TEST(Dynamics, ReactionBasedJointFrictionUsesTransverseReaction)
{
    {
        auto a = model(0.2, 1e-9, 0.01);
        a->constantGravity->setg(-10, 0, 0);
        auto b = body(a, "Slider", 2, 0.4);
        b->setVelocity3D(0, 0, 1);
        auto i = marker(a, "Rail"), j = marker(b, "Carriage");
        auto guide = joint<ASMTTranslationalJoint>(a, "Guide", i, j);
        auto friction = ASMTForceTorque::With();
        friction->setName("GuideFriction");
        friction->setMarkerI(i->fullName(""));
        friction->setMarkerJ(j->fullName(""));
        friction->setReactionAxisFriction(0.25, 0.25, 1e-6, 0, 1, false, guide);
        a->addForceTorque(friction);
        run(a, true);
        ASSERT_NO_FATAL_FAILURE(complete(a, b));
        for (size_t k = 1; k < a->times->size(); ++k) {
            const double t = a->times->at(k);
            EXPECT_NEAR(b->zs->at(k), t - 1.25*t*t, 3e-5);
            EXPECT_NEAR(b->vzs->at(k), 1 - 2.5*t, 3e-5);
            EXPECT_NEAR(b->azs->at(k), -2.5, 3e-4);
            EXPECT_NEAR(reaction(friction, k, false, 2), 5, 3e-4);
        }
    }
    {
        auto a = model(0.2, 1e-9, 0.01);
        a->constantGravity->setg(-10, 0, 0);
        auto b = body(a, "Rotor", 2, 0.4);
        b->setOmega3D(0, 0, 1);
        auto i = marker(a, "Bearing"), j = marker(b, "Shaft");
        auto bearing = joint<ASMTRevoluteJoint>(a, "Hinge", i, j);
        auto friction = ASMTForceTorque::With();
        friction->setName("BearingFriction");
        friction->setMarkerI(i->fullName(""));
        friction->setMarkerJ(j->fullName(""));
        friction->setReactionAxisFriction(0.2, 0.2, 1e-6, 0, 0.1, true, bearing);
        a->addForceTorque(friction);
        run(a, true);
        ASSERT_NO_FATAL_FAILURE(complete(a, b));
        for (size_t k = 1; k < a->times->size(); ++k) {
            const double t = a->times->at(k);
            EXPECT_NEAR(b->bryzs->at(k), t - 0.5*t*t, 3e-5);
            EXPECT_NEAR(b->omezs->at(k), 1 - t, 3e-5);
            EXPECT_NEAR(b->alpzs->at(k), -1, 3e-4);
            EXPECT_NEAR(reaction(friction, k, true, 2), 0.4, 3e-4);
        }
    }
}

TEST(Dynamics, PrescribedRotorInverseDynamics)
{
    auto a = model();
    auto b = body(a, "Rotor", 2, 0.4);
    auto j = joint<ASMTRevoluteJoint>(a, "Hinge", marker(a, "Pivot"), marker(b, "Pivot"));
    auto d = drive(a, j, "time*time");
    run(a, false);
    ASSERT_NO_FATAL_FAILURE(complete(a, b));
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        EXPECT_NEAR(b->bryzs->at(k), t*t, 2e-6);
        EXPECT_NEAR(b->omezs->at(k), 2*t, 2e-6);
        EXPECT_NEAR(b->alpzs->at(k), 2, 2e-6);
        EXPECT_NEAR(reaction(d, k, true, 2), -0.8, 2e-5);
        EXPECT_NEAR(d->powers->at(k), 1.6*t, 3e-5);
    }
}

TEST(Dynamics, InvalidInputsAreRejected)
{
    auto a = model();
    auto b = body(a, "Body", 2, 0.4);
    for (double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()}) {
        a->simulationParameters->hout = bad;
        EXPECT_THROW(a->runDYNAMIC(), std::invalid_argument);
    }
    a->simulationParameters->hout = 0.02;
    a->simulationParameters->tend = 0;
    EXPECT_THROW(a->runDYNAMIC(), std::invalid_argument);
    a->simulationParameters->tend = 1;
    b->principalMassMarker->mass = 0;
    EXPECT_THROW(a->runDYNAMIC(), std::invalid_argument);
    b->principalMassMarker->mass = 2;
    b->principalMassMarker->setMomentOfInertias(0.4, -0.1, 0.4);
    EXPECT_THROW(a->runDYNAMIC(), std::invalid_argument);
    auto load = ASMTForceTorque::With();
    EXPECT_THROW(load->setSpringDamper(-1, 1, 1), std::invalid_argument);
    EXPECT_THROW(load->setSpringDamper(1, -1, 1), std::invalid_argument);
    EXPECT_THROW(load->setBushing({-1, 0, 0}, {}, {}, {}), std::invalid_argument);
    EXPECT_THROW(load->setBushing({}, {}, {},
        {0, std::numeric_limits<double>::quiet_NaN(), 0}), std::invalid_argument);
    EXPECT_THROW(load->setAxisFriction(0.2, 0.3, 0.1, 0, true), std::invalid_argument);
    EXPECT_THROW(load->setAxisFriction(0.3, 0.2, 0, 0, false), std::invalid_argument);
    EXPECT_THROW(load->setReactionAxisFriction(0.3, 0.2, 0.1, 0, 0, true, nullptr),
        std::invalid_argument);
    EXPECT_THROW(load->setForce3D(0, std::numeric_limits<double>::quiet_NaN(), 0), std::invalid_argument);
    EXPECT_THROW(a->addForceTorque(nullptr), std::invalid_argument);
}

TEST(Dynamics, PrescribedMotionAlongsideFreeDynamics)
{
    auto a = model(0.2, 1e-9);
    auto rotor = body(a, "Rotor", 2, 0.4);
    auto falling = body(a, "Falling", 1, 0.1, 3, 2);
    auto j = joint<ASMTRevoluteJoint>(a, "Hinge", marker(a, "Pivot"), marker(rotor, "Pivot"));
    auto d = drive(a, j, "time*time");
    a->constantGravity->setg(0, -9.81, 0);
    run(a, true);
    ASSERT_NO_FATAL_FAILURE(complete(a, rotor));
    ASSERT_NO_FATAL_FAILURE(complete(a, falling));
    for (size_t k = 2; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        EXPECT_NEAR(rotor->bryzs->at(k), t*t, 2e-6);
        EXPECT_NEAR(rotor->omezs->at(k), 2*t, 2e-5);
        EXPECT_NEAR(reaction(d, k, true, 2), -0.8, 2e-4);
        EXPECT_NEAR(falling->ys->at(k), 2-4.905*t*t, 2e-6);
    }
}

TEST(Dynamics, FinalSampleAndNonzeroStart)
{
    for (double duration : {0.005, 0.055}) {
        auto a = model(0.3+duration, 1e-9);
        a->simulationParameters->tstart = 0.3;
        auto b = body(a, "Body", 2, 0.4, 0, 2);
        a->constantGravity->setg(0, -9.81, 0);
        run(a, true);
        ASSERT_GE(a->times->size(), 3u);
        EXPECT_NEAR(a->times->back(), 0.3+duration, 1e-14);
        EXPECT_NEAR(b->ys->back(), 2-4.905*duration*duration, 2e-7);
        EXPECT_NEAR(b->vys->back(), -9.81*duration, 2e-7);
        EXPECT_NEAR(b->ays->back(), -9.81, 2e-6);
        for (size_t k = 2; k < a->times->size(); ++k) {
            EXPECT_GT(a->times->at(k), a->times->at(k-1));
            EXPECT_LE(a->times->at(k), 0.3+duration);
        }
    }
}

TEST(Dynamics, RerunReplacesResults)
{
    auto a = model(0.1);
    auto b = body(a, "Body", 2, 0.4, 0, 2);
    a->constantGravity->setg(0, -9.81, 0);
    run(a, true);
    const auto firstTimes = *a->times;
    const auto firstY = *b->ys;
    b->setPosition3D(0, 2, 0);
    b->setVelocity3D(0, 0, 0);
    run(a, true);
    EXPECT_EQ(*a->times, firstTimes);
    EXPECT_EQ(*b->ys, firstY);
}

TEST(Dynamics, LoadJacobiansMatchIndependentFiniteDifferences)
{
    auto a = model(0.02);
    auto bi = body(a, "I", 2, 0.4, -0.4, 0.2, 0.3);
    auto bj = body(a, "J", 3, 0.7, 1.1, -0.3, -0.2);
    // A substantial swing exposes an incorrect axis-only torsional tangent.
    const double swing = .8;
    bj->setRotationMatrix(1, 0, 0, 0, std::cos(swing), -std::sin(swing),
                          0, std::sin(swing), std::cos(swing));
    bi->setVelocity3D(0.3, -0.2, 0.4);
    bj->setVelocity3D(-0.1, 0.5, -0.3);
    bi->setOmega3D(0.2, -0.4, 0.7);
    bj->setOmega3D(-0.3, 0.6, 0.2);
    auto mi = marker(bi, "Load", 0.3, 0.2), mj = marker(bj, "Load", -0.2, 0.4);
    run(a, true);
    auto pi = std::static_pointer_cast<Part>(bi->mbdObject);
    auto pj = std::static_pointer_cast<Part>(bj->mbdObject);
    pi->iqX(0); pi->iqE(3); pj->iqX(7); pj->iqE(10);
    auto fi = std::static_pointer_cast<EndFramec>(mi->mbdObject);
    auto fj = std::static_pointer_cast<EndFramec>(mj->mbdObject);
    using V = AppliedForceTorque::Vector;
    const auto wrench = std::make_shared<AppliedForceTorque>(fi, fj, V{1.2, -0.7, 0.9}, V{0.4, 0.8, -0.3});
    const auto spring = std::make_shared<AppliedForceTorque>(fi, fj, 18, 1.2, 1);
    const auto torsion = std::make_shared<AppliedForceTorque>(fi, fj, 18, 1.2, 0.1, true);
    const auto linearFriction = std::make_shared<AppliedForceTorque>(
        fi, fj, AppliedForceTorque::AxisFrictionParameters {0.8, 0.4, 0.2, 0.1, false}
    );
    const auto rotationalFriction = std::make_shared<AppliedForceTorque>(
        fi, fj, AppliedForceTorque::AxisFrictionParameters {0.8, 0.4, 0.2, 0.1, true}
    );
    for (auto load : {wrench, spring, torsion, linearFriction, rotationalFriction}) {
        for (bool velocity : {false, true}) {
            if (!velocity && (load == linearFriction || load == rotationalFriction))
                continue; // Axis direction uses a deliberately lagged position tangent.
            auto jac = std::make_shared<SparseMatrix<double>>(14, 14);
            if (velocity) load->fillpFpydot(jac); else load->fillpFpy(jac);
            for (size_t column = 0; column < 14; ++column) {
                auto p = column < 7 ? pi : pj;
                const size_t local = column % 7;
                FColDsptr q = velocity
                    ? (local < 3 ? p->partFrame->qXdot : p->partFrame->qEdot)
                    : (local < 3 ? p->partFrame->qX : p->partFrame->qE);
                auto& value = q->at(local < 3 ? local : local-3);
                const double original = value, h = 1e-6;
                auto residual = [&]() {
                    p->postDynCorrectorIteration();
                    auto f = std::make_shared<FullColumn<double>>(14, 0.0);
                    load->fillDynError(f);
                    return f;
                };
                value = original+h;
                auto plus = residual();
                value = original-h;
                auto minus = residual();
                value = original;
                p->postDynCorrectorIteration();
                for (size_t row = 0; row < 14; ++row) {
                    const auto found = jac->at(row)->find(column);
                    const double exact = found == jac->at(row)->end() ? 0 : found->second;
                    SCOPED_TRACE("row=" + std::to_string(row) + " col=" + std::to_string(column));
                    EXPECT_NEAR(exact, (plus->at(row)-minus->at(row))/(2*h), 2e-6);
                }
            }
        }
    }
}

TEST(Dynamics, TypedLoadsRoundTrip)
{
    const std::filesystem::path file = "typed-load-roundtrip.asmt";
    for (int kind : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}) {
        auto load = ASMTForceTorque::With();
        load->setName("Load");
        load->setMarkerI("/Benchmark/I");
        load->setMarkerJ("/Benchmark/J");
        if (kind == 1) load->setSpringDamper(18, 1.2, 1);
        else if (kind == 2) load->setTorsionalSpringDamper(18, 1.2, 0.1);
        else if (kind == 3) load->setBushing(
            {1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}
        );
        else if (kind == 4) load->setAxisFriction(4, 3, 0.2, 0.1, false);
        else if (kind == 5) load->setAxisFriction(4, 3, 0.2, 0.1, true);
        else if (kind == 6) {
            std::array<double, 36> k{}, c{};
            k[0] = k[35] = 2;
            k[5] = k[30] = 0.5;
            load->setCoupledBushing(k, c);
        }
        else if (kind == 8 || kind == 9) {
            if (kind == 8) load->setForceFormula3D(1, 0, 0, "2*sin(time)");
            else load->setTorqueFormula3D(0, 1, 0, "3*time");
            load->setFollower(true);
        }
        else {
            load->setForce3D(1.2, -3.4, 5.6);
            load->setTorque3D(-0.8, 0.2, 0.3);
            load->setFollower(kind == 7);
        }
        auto serialize = [&](std::shared_ptr<ASMTForceTorque> item) {
            { std::ofstream stream(file); item->storeOnLevel(stream, 0); }
            std::ifstream stream(file);
            std::vector<std::string> lines;
            std::string line;
            while (std::getline(stream, line)) lines.push_back(line);
            return lines;
        };
        const auto saved = serialize(load);
        auto input = saved;
        input.erase(input.begin());
        auto restored = ASMTForceTorque::With();
        restored->parseASMT(input);
        EXPECT_TRUE(input.empty());
        EXPECT_EQ(serialize(restored), saved);
    }
    std::filesystem::remove(file);
}

TEST(Dynamics, CompliantLimitRoundTrip)
{
    const std::filesystem::path file = "compliant-limit-roundtrip.asmt";
    auto limit = ASMTTranslationLimit::With();
    limit->setName("Stop");
    limit->setMarkerI("/Benchmark/I");
    limit->setMarkerJ("/Benchmark/J");
    limit->settype("=>");
    limit->setlimit("-12.5");
    limit->settol("1e-8");
    limit->setcompliance("2500", "3.5");
    auto serialize = [&](const std::shared_ptr<ASMTTranslationLimit>& item) {
        { std::ofstream stream(file); item->storeOnLevel(stream, 0); }
        std::ifstream stream(file);
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(stream, line)) lines.push_back(line);
        return lines;
    };
    const auto saved = serialize(limit);
    auto input = saved;
    input.erase(input.begin());
    auto restored = ASMTTranslationLimit::With();
    restored->parseASMT(input);
    EXPECT_TRUE(input.empty());
    EXPECT_EQ(serialize(restored), saved);

    auto legacy = saved;
    legacy.resize(legacy.size() - 6);
    legacy.erase(legacy.begin());
    auto restoredLegacy = ASMTTranslationLimit::With();
    restoredLegacy->parseASMT(legacy);
    EXPECT_TRUE(legacy.empty());
    EXPECT_EQ(restoredLegacy->behavior, "Rigid");
    std::filesystem::remove(file);
}

TEST(Dynamics, DistanceLimitRoundTrip)
{
    const std::filesystem::path file = "distance-limit-roundtrip.asmt";
    auto limit = ASMTDistanceLimit::With();
    limit->setName("SphereContact");
    limit->setMarkerI("/Benchmark/I");
    limit->setMarkerJ("/Benchmark/J");
    limit->settype("=>");
    limit->setlimit("12.5");
    limit->settol("1e-8");
    limit->setcompliance("2500", "3.5");
    auto serialize = [&](const std::shared_ptr<ASMTDistanceLimit>& item) {
        { std::ofstream stream(file); item->storeOnLevel(stream, 0); }
        std::ifstream stream(file);
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(stream, line)) lines.push_back(line);
        return lines;
    };
    const auto saved = serialize(limit);
    auto input = saved;
    input.erase(input.begin());
    auto restored = ASMTDistanceLimit::With();
    restored->parseASMT(input);
    EXPECT_TRUE(input.empty());
    EXPECT_EQ(serialize(restored), saved);
    std::filesystem::remove(file);
}

TEST(Dynamics, SmallStepTaylorOperatorAndSingularNodes)
{
    auto bdf = CREATE<StableBackwardDifference>::With();
    bdf->order = 5;
    bdf->time = 0.001;
    bdf->timeNodes = std::make_shared<FullRow<double>>();
    auto past = std::make_shared<std::vector<FColDsptr>>();
    for (int i = 1; i <= 5; ++i) {
        const double t = bdf->time-i*0.0001;
        bdf->timeNodes->push_back(t);
        past->push_back(std::make_shared<FullColumn<double>>(ListD{std::pow(t, 5)}));
    }
    bdf->calcOperatorMatrix();
    auto present = std::make_shared<FullColumn<double>>(ListD{std::pow(bdf->time, 5)});
    EXPECT_NEAR(bdf->derivativepresentpast(5, present, past)->at(0), 120, 1e-7);
    auto copy = bdf->derivativepresentpast(0, present, past);
    EXPECT_NE(copy.get(), present.get());
    EXPECT_EQ(*copy, *present);
    bdf->timeNodes->at(1) = bdf->timeNodes->at(0);
    EXPECT_THROW(bdf->calcOperatorMatrix(), SingularMatrixError);
}

TEST(Dynamics, TwoBodyForceBalanceAndLengthUnits)
{
    for (double scale : {1.0, 1000.0}) {
        auto a = model(0.1, 1e-9);
        auto bi = body(a, "I", 2, 0.4*scale*scale);
        auto bj = body(a, "J", 3, 0.7*scale*scale, 2*scale);
        auto load = ASMTForceTorque::With();
        load->setName("ForcePair");
        load->setMarkerI(marker(bi, "Load")->fullName(""));
        load->setMarkerJ(marker(bj, "Load")->fullName(""));
        load->setForce3D(2*scale, 0, 0);
        a->addForceTorque(load);
        run(a, true);
        ASSERT_NO_FATAL_FAILURE(complete(a, bi));
        ASSERT_NO_FATAL_FAILURE(complete(a, bj));
        for (size_t k = 2; k < a->times->size(); ++k) {
            const double t = a->times->at(k);
            EXPECT_NEAR(bi->xs->at(k)/scale, 0.5*t*t, 2e-7);
            EXPECT_NEAR(bj->xs->at(k)/scale, 2-t*t/3, 2e-7);
            EXPECT_NEAR(bi->axs->at(k)/scale, 1, 2e-6);
            EXPECT_NEAR(bj->axs->at(k)/scale, -2.0/3, 2e-6);
            EXPECT_NEAR((2*bi->vxs->at(k)+3*bj->vxs->at(k))/scale, 0, 2e-7);
        }
    }
}

TEST(Dynamics, CoincidentSpringAttachmentsAreRejected)
{
    auto a = model(0.1);
    auto b = body(a, "Body", 2, 0.4);
    auto load = ASMTForceTorque::With();
    load->setName("Spring");
    load->setMarkerI(marker(a, "Anchor")->fullName(""));
    load->setMarkerJ(marker(b, "Anchor")->fullName(""));
    load->setSpringDamper(18, 1.2, 1);
    a->addForceTorque(load);
    EXPECT_THROW(a->runDYNAMIC(), SimulationStoppingError);
}

TEST(Dynamics, InvalidMotionIsReported)
{
    auto a = model();
    auto b = body(a, "Rotor", 2, 0.4);
    auto j = joint<ASMTRevoluteJoint>(a, "Hinge", marker(a, "Pivot"), marker(b, "Pivot"));
    drive(a, j, "unknown_dynamics_function(time)");
    EXPECT_ANY_THROW(run(a, false));
}

TEST(Dynamics, MetreAndMillimetreInverseDynamics)
{
    // Same eccentric rotor in SI and kg-mm-s. This verifies length-squared
    // inertia/torque scaling through the solver, not the future CAD adapter.
    for (double lengthScale : {1.0, 1000.0}) {
        SCOPED_TRACE("length scale = " + number(lengthScale));
        auto a = model();
        auto b = body(a, "EccentricRotor", 2, 0.4*lengthScale*lengthScale, lengthScale);
        auto j = joint<ASMTRevoluteJoint>(a, "Hinge", marker(a, "Pivot"),
            marker(b, "Pivot", -lengthScale));
        auto d = drive(a, j, "time*time");
        run(a, false);
        ASSERT_NO_FATAL_FAILURE(complete(a, b));
        for (size_t k = 0; k < a->times->size(); ++k) {
            const double t = a->times->at(k);
            if (t <= 1e-12) continue;
            const double theta = t*t, omega = 2*t;
            const double ax = -omega*omega*std::cos(theta)-2*std::sin(theta);
            const double ay = -omega*omega*std::sin(theta)+2*std::cos(theta);
            EXPECT_NEAR(b->xs->at(k)/lengthScale, std::cos(theta), 2e-6);
            EXPECT_NEAR(b->ys->at(k)/lengthScale, std::sin(theta), 2e-6);
            EXPECT_NEAR(reaction(j, k, false, 0)/lengthScale, -2*ax, 2e-4);
            EXPECT_NEAR(reaction(j, k, false, 1)/lengthScale, -2*ay, 2e-4);
            EXPECT_NEAR(reaction(d, k, true, 2)/(lengthScale*lengthScale), -4.8, 2e-4);
        }
    }
}

TEST(Dynamics, SliderCrankKinematicsAndReactions)
{
    const double r = 0.2, l = 0.6, theta0 = 0.3, w = 1.0;
    const double x0 = r*std::cos(theta0)+std::sqrt(l*l-r*r*std::pow(std::sin(theta0), 2));
    const double phi0 = -std::asin(r/l*std::sin(theta0));
    auto a = model();
    auto crank = body(a, "Crank", 1, 0.01, 0, 0, theta0);
    auto rod = body(a, "Rod", 0.8, 0.024,
        (x0+r*std::cos(theta0))/2, r*std::sin(theta0)/2, phi0);
    auto slider = body(a, "Slider", 1.5, 0.01, x0);
    auto hinge = joint<ASMTRevoluteJoint>(a, "CrankHinge", marker(a, "GroundPivot"), marker(crank, "Pivot"));
    joint<ASMTRevoluteJoint>(a, "CrankPin", marker(crank, "Pin", r), marker(rod, "Left", -l/2));
    auto pin = joint<ASMTRevoluteJoint>(a, "SliderPin", marker(slider, "Pin"), marker(rod, "Right", l/2));
    joint<ASMTTranslationalJoint>(a, "Guide", marker(a, "Guide", 0, 0, true), marker(slider, "Guide", 0, 0, true));
    auto d = drive(a, hinge, number(theta0)+"+time");
    run(a, false);
    ASSERT_NO_FATAL_FAILURE(complete(a, crank));
    ASSERT_NO_FATAL_FAILURE(complete(a, rod));
    ASSERT_NO_FATAL_FAILURE(complete(a, slider));
    for (size_t k = 0; k < a->times->size(); ++k) {
        const double t = a->times->at(k);
        if (t <= 1e-12) continue;
        const double s = std::sin(theta0+w*t), c = std::cos(theta0+w*t);
        const double h = std::sqrt(l*l-r*r*s*s);
        const double x = r*c+h, xp = -r*s-r*r*s*c/h;
        const double xpp = -r*c-r*r*(c*c-s*s)/h-r*r*r*r*s*s*c*c/(h*h*h);
        EXPECT_NEAR(slider->xs->at(k), x, 2e-6);
        EXPECT_NEAR(slider->vxs->at(k), w*xp, 2e-6);
        EXPECT_NEAR(slider->axs->at(k), w*w*xpp, 2e-5);
        EXPECT_NEAR(reaction(pin, k, false, 0), 1.5*w*w*xpp, 2e-4);
        // Generalized input torque from d(kinetic energy)/d(theta), constant w.
        const double cxp = (xp-r*s)/2, cyp = r*c/2;
        const double cxpp = (xpp-r*c)/2, cypp = -r*s/2;
        const double phip = -r*c/h, phipp = r*s/h-r*r*r*s*c*c/(h*h*h);
        const double torque = w*w*(0.8*(cxp*cxpp+cyp*cypp)+0.024*phip*phipp+1.5*xp*xpp);
        EXPECT_NEAR(reaction(d, k, true, 2), -torque, 2e-4);
    }
}
} // namespace
