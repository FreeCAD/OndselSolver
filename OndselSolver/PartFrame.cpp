/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include<algorithm>
#include <cstdint>

#include "PartFrame.h"
#include "Part.h"
#include "EulerConstraint.h"
#include "AbsConstraint.h"
#include "MarkerFrame.h"
#include "EulerParameters.h"
#include "EulerParametersDot.h"
#include "CREATE.h"
#include "RedundantConstraint.h"
#include "System.h"

using namespace MbD;

PartFrame::PartFrame()
{
}
PartFrame::PartFrame(const std::string& str) : CartesianFrame(str)
{
}
System* PartFrame::root()
{
	return part->root();
}
void PartFrame::initialize()
{
	aGeu = CREATE<EulerConstraint>::With();
	aGeu->owner = this;
	aGabs = std::make_shared<std::vector<std::shared_ptr<Constraint>>>();
	markerFrames = std::make_shared<std::vector<std::shared_ptr<MarkerFrame>>>();
}

void PartFrame::initializeLocally()
{
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->initializeLocally(); });
	aGeu->initializeLocally();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->initializeLocally(); });
}

void PartFrame::initializeGlobally()
{
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->initializeGlobally(); });
	aGeu->initializeGlobally();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->initializeGlobally(); });
}
void PartFrame::setqX(FColDsptr x) {
	qX->copyFrom(x);
}

FColDsptr PartFrame::getqX() {
	return qX;
}

void PartFrame::setqE(FColDsptr x) {
	qE->copyFrom(x);
}

void MbD::PartFrame::setaAap(FMatDsptr mat)
{
	qE = mat->asEulerParameters();
}

FColDsptr PartFrame::getqE() {
	return qE;
}
void PartFrame::setqXdot(FColDsptr x) {
	qXdot = x;
}

FColDsptr PartFrame::getqXdot() {
	return qXdot;
}

void PartFrame::setomeOpO(FColDsptr omeOpO) {
	qEdot = EulerParametersDot<double>::FromqEOpAndOmegaOpO(qE, omeOpO);
}

FColDsptr PartFrame::getomeOpO() {
	return qEdot->omeOpO();
}

void PartFrame::setqXddot(FColDsptr x)
{
	qXddot = x;
}

FColDsptr PartFrame::getqXddot()
{
	return qXddot;
}

void PartFrame::setqEddot(FColDsptr x)
{
	qEddot = x;
}

FColDsptr PartFrame::getqEddot()
{
	return qEddot;
}

FColDsptr PartFrame::omeOpO()
{
	return qEdot->omeOpO();
}

void PartFrame::setPart(Part* x) {
	part = x;
}

Part* PartFrame::getPart() {
	return part;
}

void PartFrame::addMarkerFrame(std::shared_ptr<MarkerFrame> markerFrame)
{
	markerFrame->setPartFrame(this);
	markerFrames->push_back(markerFrame);
}

EndFrmsptr PartFrame::endFrame(const std::string& name)
{
	auto match = std::find_if(markerFrames->begin(), markerFrames->end(), [&](auto& mkr) {return mkr->name == name; });
	return (*match)->endFrames->at(0);
}

void PartFrame::aGabsDo(const std::function<void(std::shared_ptr<Constraint>)>& f)
{
	std::for_each(aGabs->begin(), aGabs->end(), f);
}

void PartFrame::markerFramesDo(const std::function<void(std::shared_ptr<MarkerFrame>)>& f)
{
	std::for_each(markerFrames->begin(), markerFrames->end(), f);
}

void PartFrame::removeRedundantConstraints(std::shared_ptr<std::vector<size_t>> redundantEqnNos)
{
	if (std::find(redundantEqnNos->begin(), redundantEqnNos->end(), aGeu->iG) != redundantEqnNos->end()) {
		auto redunCon = CREATE<RedundantConstraint>::With();
		redunCon->constraint = aGeu;
		aGeu = redunCon;
	}
	for (size_t i = 0; i < aGabs->size(); i++)
	{
		auto& constraint = aGabs->at(i);
		if (std::find(redundantEqnNos->begin(), redundantEqnNos->end(), constraint->iG) != redundantEqnNos->end()) {
			auto redunCon = CREATE<RedundantConstraint>::With();
			redunCon->constraint = constraint;
			aGabs->at(i) = redunCon;
		}
	}
}

void PartFrame::reactivateRedundantConstraints()
{
	if (aGeu->isRedundant()) aGeu = std::dynamic_pointer_cast<RedundantConstraint>(aGeu)->constraint;
	for (size_t i = 0; i < aGabs->size(); i++)
	{
		auto& con = aGabs->at(i);
		if (con->isRedundant()) {
			aGabs->at(i) = std::static_pointer_cast<RedundantConstraint>(con)->constraint;
		}
	}
}

void PartFrame::constraintsReport()
{
	auto redunCons = std::make_shared<std::vector<std::shared_ptr<Constraint>>>();
	aGabsDo([&](std::shared_ptr<Constraint> con) {
		if (con->isRedundant()) {
			redunCons->push_back(con);
		}
		});
	if (aGeu->isRedundant()) redunCons->push_back(aGeu);
	if (redunCons->size() > 0) {
		std::string str = "MbD: " + part->classname() + std::string(" ") + part->name + " has the following constraint(s) removed: ";
		this->logString(str);
		std::for_each(redunCons->begin(), redunCons->end(), [&](auto& con) {
			str = "MbD: " + std::string("    ") + std::string(typeid(*con).name());
			this->logString(str);
			});
	}
}

void PartFrame::prePosIC()
{
	iqX = SIZE_MAX;
	iqE = SIZE_MAX;
	Item::prePosIC();
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrm) { markerFrm->prePosIC(); });
	aGeu->prePosIC();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->prePosIC(); });
}

void PartFrame::prePosKine()
{
	iqX = SIZE_MAX;
	iqE = SIZE_MAX;
	this->calcPostDynCorrectorIteration();
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrm) { markerFrm->prePosKine(); });
	aGeu->prePosKine();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->prePosKine(); });
}

FColDsptr PartFrame::rOpO()
{
	return qX;
}

FMatDsptr PartFrame::aAOp()
{
	return qE->aA;
}

FMatDsptr PartFrame::aC()
{
	return qE->aC;
}

FMatDsptr PartFrame::aCdot()
{
	return qEdot->aCdot;
}

FColDsptr PartFrame::alpOpO()
{
	auto& aB = qE->aB;
	auto& aBdot = qEdot->aBdot;
	return aBdot->timesFullColumn(qEdot)->plusFullColumn(aB->timesFullColumn(qEddot))->times(2.0);
}

FColFMatDsptr PartFrame::pAOppE()
{
	return qE->pApE;
}

void PartFrame::fillEssenConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> essenConstraints)
{
	aGeu->fillEssenConstraints(aGeu, essenConstraints);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillEssenConstraints(con, essenConstraints); });
}

void PartFrame::fillRedundantConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>>)
{
}

void PartFrame::fillConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> allConstraints)
{
	aGeu->fillConstraints(aGeu, allConstraints);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillConstraints(con, allConstraints); });
}

void PartFrame::fillqsu(FColDsptr col)
{
	col->atiputFullColumn(iqX, qX);
	col->atiputFullColumn(iqE, qE);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillqsu(col); });
}

void PartFrame::fillqsuWeights(DiagMatDsptr diagMat)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillqsuWeights(diagMat); });
}

void MbD::PartFrame::fillqsuddotlam(FColDsptr col)
{
	col->atiputFullColumn(iqX, qXddot);
	col->atiputFullColumn(iqE, qEddot);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillqsuddotlam(col); });
	aGeu->fillqsuddotlam(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillqsuddotlam(col); });
}

void PartFrame::fillqsulam(FColDsptr col)
{
	col->atiputFullColumn(iqX, qX);
	col->atiputFullColumn(iqE, qE);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillqsulam(col); });
	aGeu->fillqsulam(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillqsulam(col); });
}

void PartFrame::fillqsudot(FColDsptr col)
{
	col->atiputFullColumn(iqX, qXdot);
	col->atiputFullColumn(iqE, qEdot);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillqsudot(col); });
}

void PartFrame::fillqsudotWeights(DiagMatDsptr diagMat)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillqsudotWeights(diagMat); });
}

void PartFrame::useEquationNumbers()
{
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->useEquationNumbers(); });
	aGeu->useEquationNumbers();
	aGabsDo([](std::shared_ptr<Constraint> con) { con->useEquationNumbers(); });
}

void PartFrame::setqsu(FColDsptr col)
{
	qX->equalFullColumnAt(col, iqX);
	qE->equalFullColumnAt(col, iqE);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->setqsu(col); });
	aGeu->setqsu(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->setqsu(col); });
}

void PartFrame::setqsulam(FColDsptr col)
{
	qX->equalFullColumnAt(col, iqX);
	qE->equalFullColumnAt(col, iqE);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->setqsulam(col); });
	aGeu->setqsulam(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->setqsulam(col); });
}

void PartFrame::setqsudotlam(FColDsptr col)
{
	qXdot->equalFullColumnAt(col, iqX);
	qEdot->equalFullColumnAt(col, iqE);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->setqsudotlam(col); });
	aGeu->setqsudotlam(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->setqsudotlam(col); });
}

void PartFrame::setqsudot(FColDsptr col)
{
	qXdot->equalFullColumnAt(col, iqX);
	qEdot->equalFullColumnAt(col, iqE);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->setqsudot(col); });
}

void PartFrame::postPosICIteration()
{
	Item::postPosICIteration();
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->postPosICIteration(); });
	aGeu->postPosICIteration();
	aGabsDo([](std::shared_ptr<Constraint> con) { con->postPosICIteration(); });
}

void PartFrame::fillPosICError(FColDsptr col)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillPosICError(col); });
	aGeu->fillPosICError(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillPosICError(col); });
}

void PartFrame::fillPosICJacob(SpMatDsptr mat)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillPosICJacob(mat); });
	aGeu->fillPosICJacob(mat);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillPosICJacob(mat); });
}

void PartFrame::postPosIC()
{
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->postPosIC(); });
	aGeu->postPosIC();
	aGabsDo([](std::shared_ptr<Constraint> con) { con->postPosIC(); });
}

void PartFrame::preDyn()
{
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->preDyn(); });
	aGeu->preDyn();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->preDyn(); });
}

void PartFrame::storeDynState()
{
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->storeDynState(); });
	aGeu->storeDynState();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->storeDynState(); });
}

void PartFrame::fillPosKineError(FColDsptr col)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillPosKineError(col); });
	aGeu->fillPosKineError(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillPosKineError(col); });
}

void PartFrame::preVelIC()
{
	Item::preVelIC();
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->preVelIC(); });
	aGeu->preVelIC();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->preVelIC(); });
}

void PartFrame::postVelIC()
{
	qEdot->calcAdotBdotCdot();
	qEdot->calcpAdotpE();
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->postVelIC(); });
	aGeu->postVelIC();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->postVelIC(); });
}

void PartFrame::fillVelICError(FColDsptr col)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillVelICError(col); });
	aGeu->fillVelICError(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillVelICError(col); });
}

void PartFrame::fillVelICJacob(SpMatDsptr mat)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillVelICJacob(mat); });
	aGeu->fillVelICJacob(mat);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillVelICJacob(mat); });
}

void PartFrame::preAccIC()
{
	qXddot = std::make_shared<FullColumn<double>>(3, 0.0);
	qEddot = std::make_shared<FullColumn<double>>(4, 0.0);
	Item::preAccIC();
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->preAccIC(); });
	aGeu->preAccIC();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->preAccIC(); });
}

void PartFrame::fillAccICIterError(FColDsptr col)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillAccICIterError(col); });
	aGeu->fillAccICIterError(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillAccICIterError(col); });
}

void PartFrame::fillAccICIterJacob(SpMatDsptr mat)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillAccICIterJacob(mat); });
	aGeu->fillAccICIterJacob(mat);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillAccICIterJacob(mat); });
}

void PartFrame::setqsuddotlam(FColDsptr col)
{
	qXddot->equalFullColumnAt(col, iqX);
	qEddot->equalFullColumnAt(col, iqE);
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->setqsuddotlam(col); });
	aGeu->setqsuddotlam(col);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->setqsuddotlam(col); });
}

FMatDsptr PartFrame::aBOp()
{
	return qE->aB;
}

void PartFrame::fillPosKineJacob(SpMatDsptr mat)
{
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillPosKineJacob(mat); });
	aGeu->fillPosKineJacob(mat);
	aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillPosKineJacob(mat); });
}

double PartFrame::suggestSmallerOrAcceptDynStepSize(double hnew)
{
	auto hnew2 = hnew;
	auto speed = qXdot->length();
	double htran;
	if (speed < 1.0e-15) {
		htran = 1.0e99;
	}
	else {
		htran = this->root()->translationLimit() / speed;
	}
	if (hnew2 > htran) {
		this->logString("MbD: Time step limited by translation limit per step.");
		hnew2 = htran;
	}
	auto omegaMagnitude = qEdot->omeOpO()->length();
	double hrot;
	if (omegaMagnitude < 1.0e-15) {
		hrot = 1.0e99;
	}
	else {
		hrot = this->root()->rotationLimit() / omegaMagnitude;
	}
	if (hnew2 > hrot) {
		this->logString("MbD: Time step limited by rotation limit per step.");
		hnew2 = hrot;
	}
	markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { hnew2 = markerFrame->suggestSmallerOrAcceptDynStepSize(hnew2); });
	hnew2 = aGeu->suggestSmallerOrAcceptDynStepSize(hnew2);
	aGabsDo([&](std::shared_ptr<Constraint> aGab) { hnew2 = aGab->suggestSmallerOrAcceptDynStepSize(hnew2); });
	return hnew2;
}

void PartFrame::postDynStep()
{
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->postDynStep(); });
	aGeu->postDynStep();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->postDynStep(); });
}

void PartFrame::asFixed()
{
    // Ground at the supplied pose, not at zero coordinates/identity rotation.
    // The Euler normalization constraint supplies the fourth quaternion equation.
    // Leave the largest component unconstrained so that this remains full rank
    // at half-turns as well as near the identity orientation.
    size_t dependent = 0;
    for (size_t i = 1; i < 4; ++i) {
        if (std::abs(qE->at(i)) > std::abs(qE->at(dependent))) {
            dependent = i;
        }
    }
    aGabs->clear();
    for (size_t i = 0; i < 7; ++i) {
        if (i == 3 + dependent) {
            continue;
        }
        auto con = CREATE<AbsConstraint>::With(i);
        con->owner = this;
        con->setConstant(i < 3 ? qX->at(i) : qE->at(i - 3));
        aGabs->push_back(con);
    }
}

void PartFrame::postInput()
{
	qXddot = std::make_shared<FullColumn<double>>(3, 0.0);
	qEddot = std::make_shared<FullColumn<double>>(4, 0.0);
	Item::postInput();
	markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->postInput(); });
	aGeu->postInput();
	aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->postInput(); });
}

void PartFrame::calcPostDynCorrectorIteration()
{
	qE->calcABC();
	qE->calcpApE();
	qEdot->calcAdotBdotCdot();
	qEdot->calcpAdotpE();
}

void PartFrame::fillpqsumu(FColDsptr col)
{
    //"Fill q, s and lam into col."
    col->atiputFullColumn(iqX, qX);
    col->atiputFullColumn(iqE, qE);
    markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillpqsumu(col); });
    aGeu->fillpqsumu(col);
    aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillpqsumu(col); });
}

void PartFrame::fillpqsumudot(FColDsptr col)
{
    col->atiputFullColumn(iqX, qXdot);
    col->atiputFullColumn(iqE, qEdot);
    markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillpqsumudot(col); });
    aGeu->fillpqsumudot(col);
    aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillpqsumudot(col); });
}

void PartFrame::setpqsumu(FColDsptr col)
{
    qX->equalFullColumnAt(col, iqX);
    qE->equalFullColumnAt(col, iqE);
    markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->setpqsumu(col); });
    aGeu->setpqsumu(col);
    aGabsDo([&](std::shared_ptr<Constraint> con) { con->setpqsumu(col); });
}

void PartFrame::setpqsumudot(FColDsptr col)
{
    qXdot->equalFullColumnAt(col, iqX);
    qEdot->equalFullColumnAt(col, iqE);
    markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->setpqsumudot(col); });
    aGeu->setpqsumudot(col);
    aGabsDo([&](std::shared_ptr<Constraint> con) { con->setpqsumudot(col); });
}

void PartFrame::setpqsumuddot(FColDsptr col)
{
    qXddot->equalFullColumnAt(col, iqX);
    qEddot->equalFullColumnAt(col, iqE);
    markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->setpqsumuddot(col); });
    aGeu->setpqsumuddot(col);
    aGabsDo([&](std::shared_ptr<Constraint> con) { con->setpqsumuddot(col); });
}

void PartFrame::postDynPredictor()
{
    CartesianFrame::postDynPredictor();
    markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->postDynPredictor(); });
    aGeu->postDynPredictor();
    aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->postDynPredictor(); });
}

void PartFrame::fillDynError(FColDsptr col)
{
    markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillDynError(col); });
    aGeu->fillDynError(col);
    aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillDynError(col); });
}

void PartFrame::fillpFpy(SpMatDsptr mat)
{
    //markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillpFpy(mat); });
    aGeu->fillpFpy(mat);
    aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillpFpy(mat); });
}

void PartFrame::fillpFpydot(SpMatDsptr mat)
{
    //markerFramesDo([&](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->fillpFpydot(mat); });
    aGeu->fillpFpydot(mat);
    aGabsDo([&](std::shared_ptr<Constraint> con) { con->fillpFpydot(mat); });
}

void PartFrame::postDynCorrectorIteration()
{
    CartesianFrame::postDynCorrectorIteration();
    markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->postDynCorrectorIteration(); });
    aGeu->postDynCorrectorIteration();
    aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->postDynCorrectorIteration(); });
}

void PartFrame::preDynOutput()
{
    CartesianFrame::preDynOutput();
    markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->preDynOutput(); });
    aGeu->preDynOutput();
    aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->preDynOutput(); });
}

void PartFrame::postDynOutput()
{
    CartesianFrame::postDynOutput();
    markerFramesDo([](std::shared_ptr<MarkerFrame> markerFrame) { markerFrame->postDynOutput(); });
    aGeu->postDynOutput();
    aGabsDo([](std::shared_ptr<Constraint> aGab) { aGab->postDynOutput(); });
}
