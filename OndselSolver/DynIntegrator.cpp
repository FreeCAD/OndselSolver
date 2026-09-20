/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <string>
#include <cassert>

#include "DynIntegrator.h"
#include "BasicIntegrator.h"
#include "SystemSolver.h"
#include "Solver.h"
#include "Item.h"
#include "SingularMatrixError.h"
#include "SimulationStoppingError.h"
#include "TooSmallStepSizeError.h"
#include "TooManyTriesError.h"
#include "Constraint.h"
#include "Part.h"
#include "DiscontinuityError.h"
#include "BasicDAEIntegrator.h"
#include "System.h"

using namespace MbD;

std::shared_ptr<DynIntegrator> DynIntegrator::With()
{
    auto inst = std::make_shared<DynIntegrator>();
    inst->initialize();
    return inst;
}

void DynIntegrator::assignEquationNumbers()
{
    //"Equation order is p,q,s,u,w,mubar,mu."

    auto parts = system->parts();
    //auto contactEndFrames = system->contactEndFrames();
    //auto uHolders = system->uHolders();
    auto constraints = system->allConstraints();
    ncon = constraints->size();
    auto eqnNo = 0;
    for (auto& part : *parts) {
        part->ipX = eqnNo;
        eqnNo = eqnNo + 3;
        part->ipE = eqnNo;
        eqnNo = eqnNo + 4;
    }
    for (auto& part : *parts) {
        part->iqX(eqnNo);
        eqnNo = eqnNo + 3;
        part->iqE(eqnNo);
        eqnNo = eqnNo + 4;
    }
    //for (auto& endFrm : *contactEndFrames) {
    //    endFrm->is(eqnNo);
    //    eqnNo = eqnNo + endFrm->sSize();
    //}
    //for (auto& uHolder : *uHolders) {
    //    uHolder->iu(eqnNo);
    //    eqnNo += 1;
    //}
    for (auto& con : *constraints) {
        con->iG = eqnNo;
        eqnNo += 1;
    }
    neqn = eqnNo;
}

void DynIntegrator::checkForDiscontinuity()
{
    //"Check for discontinuity in (tpast,t] or [t,tpast) if integrating
    //backward."

    auto t = integrator->t;
    auto tprevious = integrator->tprevious();
    auto epsilon = std::numeric_limits<double>::epsilon();
    double tstartNew;
    if (direction == 0) {
        tstartNew = epsilon;
    }
    else {
        epsilon = std::abs(t) * epsilon;
        tstartNew = ((direction * t) + epsilon) / direction;
    }
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { tstartNew = item->checkForDynDiscontinuityBetweenand(tprevious, tstartNew); });
    if (auto events = system->system->dynamicEvents) {
        const double end = std::min({t, tend, tstartNew});
        const double eventTime = events->locate(tprevious, end, [this](double sample) { interpolateAt(sample); });
        interpolateAt(t);
        if (std::isfinite(eventTime) && eventTime <= end + events->tolerance) {
            checkForOutputThrough(std::nextafter(eventTime, -std::numeric_limits<double>::infinity()));
            events->commit(tprevious, eventTime, [this](double sample) { interpolateAt(sample); });
            system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynStep(); });
            auto types = std::make_shared<std::vector<DiscontinuityType>>(1, EVENT);
            if (tstartNew <= t && std::abs(eventTime-tstartNew) <= events->tolerance) {
                // Do not lose a joint-limit impact coinciding with an event.
                system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) {
                    item->discontinuityAtaddTypeTo(tstartNew, types);
                });
            }
            events->apply(eventTime);
            system->tstartPastsAddFirst(tstart);
            system->tstart = eventTime;
            system->toutFirst = tout <= eventTime + events->tolerance ? std::min(tout + hout, tend) : tout;
            throwDiscontinuityError("Event boundary", types);
        }
        events->commit(tprevious, end, [this](double sample) { interpolateAt(sample); });
        interpolateAt(t);
    }
    if ((direction * tstartNew) > (direction * t)) {
        //"No discontinuity in step"
        return;
    }
    else {
        checkForOutputThrough(tstartNew);
        interpolateAt(tstartNew);
        if (system->system->dynamicEvents)
            system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynStep(); });
        system->tstartPastsAddFirst(tstart);
        system->tstart = tstartNew;
        system->toutFirst = tout;
        auto discontinuityTypes = std::make_shared<std::vector<DiscontinuityType>>();
        system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { item->discontinuityAtaddTypeTo(tstartNew, discontinuityTypes); });
        throwDiscontinuityError("", discontinuityTypes);
    }
}

void DynIntegrator::checkForOutputThrough(double t)
{
    //"Inclusive of t."

    auto tlimit = std::min(t, tend);
    auto thereIsOutput = false;
    if (direction * tend <= (direction * tlimit)) {
        integrator->_continue = false;
    }
    while (direction * tout <= (direction * tlimit)) {
        thereIsOutput = true;
        auto yout = integrator->yDerivat(0, tout);
        auto ydotout = integrator->yDerivat(1, tout);
        auto yddotout = integrator->yDerivat(2, tout);
        system->time(tout);
        system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
            item->setpqsumu(yout);
            item->setpqsumudot(ydotout);
            item->setpqsumuddot(yddotout);
            });
        preDAEOutput();
        system->output();
        if (tout == tend) {
            tout = std::numeric_limits<double>::infinity();
        }
        else {
            const double next = tout + hout;
            const double roundoff = 16*std::numeric_limits<double>::epsilon()
                * std::max({std::abs(tstart), std::abs(tend), std::abs(hout)});
            // Avoid a near-identical regular sample followed by the end sample.
            tout = next >= tend-roundoff ? tend : next;
        }
    }
    if (thereIsOutput) {
        //"Reset system to integrator time."
        system->time(integrator->t);
        auto integ = std::static_pointer_cast<BasicDAEIntegrator>(integrator);
        auto& y = integ->y;
        auto& ydot = integ->ydot;
        system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
            item->setpqsumu(y);
            item->setpqsumudot(ydot);
            });
        postDAEOutput();
    }
}

void DynIntegrator::fillF(FColDsptr col)
{
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { item->fillDynError(col); });
}

void DynIntegrator::fillpFpy(SpMatDsptr mat)
{
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { item->fillpFpy(mat); });
}

void DynIntegrator::fillpFpydot(SpMatDsptr mat)
{
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { item->fillpFpydot(mat); });
}

void DynIntegrator::preRun()
{
    std::string str("MbD: Starting dynamic analysis.");
    system->logString(str);
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { item->preDyn(); });
}

void DynIntegrator::run()
{
    try {
        try {
            try {
                DAEIntegrator::run();
            }
            catch (const SingularMatrixError&) {
                std::stringstream ss;
                ss << "MbD: Solver has encountered a singular matrix." << std::endl;
                ss << "MbD: Check to see if a massless or a very low mass part is under constrained." << std::endl;
                ss << "MbD: Check to see if the system is in a locked position." << std::endl;
                ss << "MbD: Check to see if the error tolerance is too demanding." << std::endl;
                ss << "MbD: Check to see if a curve-curve is about to have multiple contact points." << std::endl;
                auto str = ss.str();
                logString(str);
                throw SimulationStoppingError(str);
            }
        }
        catch (const TooSmallStepSizeError&) {
            std::stringstream ss;
            ss << "MbD: Step size is prevented from going below the user specified minimum." << std::endl;
            ss << "MbD: Check to see if the system is in a locked position." << std::endl;
            ss << "MbD: Check to see if a curve-curve is about to have multiple contact points." << std::endl;
            ss << "MbD: If they are not, lower the permitted minimum step size." << std::endl;
            auto str = ss.str();
            logString(str);
            throw SimulationStoppingError(str);
        }
    }
    catch (const TooManyTriesError&) {
        std::stringstream ss;
        ss << "MbD: Check to see if the error tolerance is too demanding." << std::endl;
        auto str = ss.str();
        logString(str);
        throw SimulationStoppingError(str);
    }
}

void DynIntegrator::fillY(FColDsptr y)
{
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->fillpqsumu(y); });
}

void DynIntegrator::fillYdot(FColDsptr ydot)
{
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->fillpqsumudot(ydot); });
}

void DynIntegrator::incrementTime(double aDouble)
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->storeDynState(); });
    system->time(aDouble);
}

void DynIntegrator::interpolateAt(double t)
{
    auto yout = integrator->yDerivat(0, t);
    auto ydotout = integrator->yDerivat(1, t);
    system->time(t);
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
        item->setpqsumu(yout);
        item->setpqsumudot(ydotout);
    });
    system->partsJointsMotionsLimitsForcesTorquesDo(
        [](std::shared_ptr<Item> item) { item->postDynPredictor(); }
    );
}

size_t DynIntegrator::iterMax()
{
    return system->iterMaxDyn;
}

void DynIntegrator::postDAECorrector()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynCorrector(); });
}

bool DynIntegrator::acceptTrial() const
{
    bool accepted = true;
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) {
        if (!item->acceptDynTrial()) accepted = false;
    });
    return accepted;
}

void DynIntegrator::postDAECorrectorIteration()
{
    throw SimulationStoppingError("To be implemented.");
}

void DynIntegrator::postDAEFirstStep()
{
    if (system->system->dynamicEvents) {
        checkForDiscontinuity();
        system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynFirstStep(); });
        checkForOutputThrough(integrator->t);
        return;
    }
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynFirstStep(); });
    if (integrator->istep >= 0) {
        //"Noise make checking at the start unreliable."
        checkForDiscontinuity();
    }
    checkForOutputThrough(integrator->t);
}

void DynIntegrator::postDAEOutput()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynOutput(); });
}

double DynIntegrator::suggestSmallerOrAcceptFirstStepSize(double hnew)
{
    auto hnew2 = hnew;
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { hnew2 = item->suggestSmallerOrAcceptDynFirstStepSize(hnew2); });
    hnew2 = std::min(hnew2, hmax);
    if (hnew2 < hmin) throw TooSmallStepSizeError("Initial dynamics step is below the permitted minimum");
    return hnew2;
}

double DynIntegrator::suggestSmallerOrAcceptStepSize(double hnew)
{
    auto hnew2 = hnew;
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { hnew2 = item->suggestSmallerOrAcceptDynStepSize(hnew2); });
    if (hnew2 > hmax) {
        hnew2 = hmax;
        system->logString("MbD: Step size is at user specified maximum.");
    }
    if (hnew2 < hmin) {
        std::stringstream ss;
        ss << "MbD: Step size " << hnew2 << " < " << hmin << "  user specified minimum.";
        auto str = ss.str();
        system->logString(str);
        throw TooSmallStepSizeError("");
    }
    return hnew2;
}

void DynIntegrator::updateForDAECorrector()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynCorrectorIteration(); });
}

void DynIntegrator::y(FColDsptr col)
{
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setpqsumu(col); });
}

void DynIntegrator::ydot(FColDsptr col)
{
    system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setpqsumudot(col); });
}

void DynIntegrator::throwDiscontinuityError(const std::string& str, std::shared_ptr<std::vector<DiscontinuityType>> discontinuityTypes)
{
    throw DiscontinuityError(str, discontinuityTypes);
}

void DynIntegrator::useTrialStepStats(std::shared_ptr<SolverStatistics> stats)
{
    system->useDynTrialStepStats(stats);
}

void DynIntegrator::useDAEStepStats(std::shared_ptr<SolverStatistics> stats)
{
    //Do nothing.
}

void DynIntegrator::reportStats()
{
    //Do nothing.
}

void DynIntegrator::postDAEPredictor()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynPredictor(); });
}

void DynIntegrator::postDAEStep()
{
    if (system->system->dynamicEvents) {
        checkForDiscontinuity();
        system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynStep(); });
        checkForOutputThrough(integrator->t);
        return;
    }
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->postDynStep(); });
    if (integrator->istep >= 0) {
        //"Noise make checking at the start unreliable."
        checkForDiscontinuity();
    }
    checkForOutputThrough(integrator->t);
}

void DynIntegrator::postRun()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { item->postDyn(); });
}

void DynIntegrator::preDAECorrector()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->preDynCorrector(); });
}

void DynIntegrator::preDAECorrectorIteration()
{
    throw SimulationStoppingError("To be implemented.");
}

void DynIntegrator::preDAEFirstStep()
{
    throw SimulationStoppingError("To be implemented.");
}

void DynIntegrator::preDAEOutput()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->preDynOutput(); });
}

void DynIntegrator::preDAEPredictor()
{
    system->partsJointsMotionsLimitsForcesTorquesDo([](std::shared_ptr<Item> item) { item->preDynPredictor(); });
}

void DynIntegrator::preDAEStep()
{
    throw SimulationStoppingError("To be implemented.");
}
