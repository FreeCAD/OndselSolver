/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
//#include <windows.h>
#include <assert.h>
//#include <debugapi.h>
#include <sstream> 
#include <chrono>

#include "Item.h"
#include "System.h"
#include "Symbolic.h"

using namespace MbD;

Item::Item() {
	auto now = std::chrono::high_resolution_clock::now();
	auto nanoseconds = now.time_since_epoch().count();
	name = std::to_string(nanoseconds);
}

Item::Item(const std::string& str) : name(str)
{
}

System* Item::root()
{
	return owner->root();
}

void MbD::Item::noop()
{
	//No Operations
}

void Item::initialize()
{
}

std::ostream& Item::printOn(std::ostream& s) const
{
	std::string str = typeid(*this).name();
	auto classname = str.substr(11, str.size() - 11);
	s << classname << std::endl;
	return s;
}

void Item::initializeLocally()
{
}

bool MbD::Item::isJointForce()
{
	throw SimulationStoppingError("To be implemented.");
	return false;
}

bool MbD::Item::isJointTorque()
{
	throw SimulationStoppingError("To be implemented.");
	return false;
}

bool MbD::Item::isKinedotIJ()
{
	throw SimulationStoppingError("To be implemented.");
	return false;
}

bool MbD::Item::isKineIJ()
{
	throw SimulationStoppingError("To be implemented.");
	return false;
}

void Item::initializeGlobally()
{
}

void Item::postInput()
{
	//Called once after input
	calcPostDynCorrectorIteration();
}

void Item::calcPostDynCorrectorIteration()
{
}

void MbD::Item::checkForCollisionDiscontinuityBetweenand(double, double)
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::removeRedundantConstraints(std::shared_ptr<std::vector<size_t>>)
{
}

void Item::setpqsumu(FColDsptr)
{
    //Do nothing.
}

void Item::setpqsumuddot(FColDsptr)
{
    //Do nothing.
}

void Item::setpqsumudot(FColDsptr)
{
    //Do nothing.
}

void Item::reactivateRedundantConstraints()
{
}

void MbD::Item::registerName()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::fillPosKineError(FColDsptr)
{
}

void Item::fillPosKineJacob(SpMatDsptr)
{
}

void Item::fillpqsumu(FColDsptr)
{
    //Do nothing.
}

void Item::fillpqsumudot(FColDsptr)
{
    //Do nothing.
}

void Item::fillEssenConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>>)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::fillPerpenConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>>)
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::fillpFpy(SpMatDsptr)
{
    //Do nothing.
}

void Item::fillpFpydot(SpMatDsptr)
{
    //Do nothing.
}

void Item::fillRedundantConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>>)
{
}

void MbD::Item::fillStaticError(FColDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::fillStaticJacob(FMatDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::fillConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>>)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::fillDispConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>>)
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::fillDynError(FColDsptr)
{
    //Do nothing.
}

void Item::fillqsu(FColDsptr)
{
}

void Item::fillqsuWeights(DiagMatDsptr)
{
}

void Item::fillqsulam(FColDsptr)
{
}

void Item::setqsulam(FColDsptr)
{
}

void MbD::Item::simUpdateAll()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::preDyn()
{
	//"Assume positions, velocities and accelerations are valid."
	//"Called once before solving for dynamic solution."
	//"Update all variable dependent instance variables needed for runDYNAMICS even if they 
	//have been calculated in postPosIC, postVelIC and postAccIC."
	//"Calculate p, pdot."
	//"Default is do nothing."
}

void MbD::Item::preDynCorrector()
{
    // Default: no corrector preparation is required.
}

void MbD::Item::preDynCorrectorIteration()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::postDyn()
{
	//"Assume runDYNAMICS ended successfully."
	//"Called once at the end of runDYNAMICS."
	//"Update all instance variables dependent on p,q,s,u,mu,pdot,qdot,sdot,udot,mudot (lam) 
	//regardless of whether they are needed."
	//"This is a subset of update."
	//"Default is do nothing."
}

void MbD::Item::postDynCorrector()
{
    // The corrector iteration already updated dependent quantities.
}

void Item::postDynCorrectorIteration()
{
    //"Called after the predictor stage in the dynamic solution."
    //"Update only instance variables dependent on p,q,s,u,mu,pdot,qdot,sdot,udot,mudot (lam) that are needed for the corrector stage."

    calcPostDynCorrectorIteration();
}

std::string Item::classname()
{
	std::string str = typeid(*this).name();
	auto answer = str.substr(11, str.size() - 11);
	return answer;
}

void Item::preDynFirstStep()
{
	//"Called before the start of the first step in the dynamic solution."
	this->preDynStep();
}

void Item::preDynOutput()
{
    //"Calculate all instance variables just before output."
    calcPostDynCorrectorIteration();
}

void MbD::Item::preDynPredictor()
{
    // Default: no predictor preparation is required.
}

void Item::postDynFirstStep()
{
	this->postDynStep();
}

void Item::postDynOutput()
{
    //"Calculate all instance variables just after output."
    calcPostDynCorrectorIteration();
}

void Item::postDynPredictor()
{
    //"Called after the predictor stage in the dynamic solution."
    //"Update only instance variables dependent on p,q,s,u,mu,pdot,qdot,sdot,udot,mudot (lam)
    //that are needed for the corrector stage."
    //"Needless updating can be expensive here."
    //"This is a subset of update."
    //"Default is do nothing."
    //"updateInSimulation is the interface to the old system."

    calcPostDynCorrectorIteration();
}

void Item::preDynStep()
{
}

void MbD::Item::preICRestart()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::postDynStep()
{
	//"Called after the end of a complete step in the dynamic solution."
	//"Update info before checking for discontinuities."
	//"Default is do nothing."
}

void Item::storeDynState()
{
}

double MbD::Item::suggestSmallerOrAcceptCollisionFirstStepSize(double)
{
	throw SimulationStoppingError("To be implemented.");
	return 0.0;
}

double MbD::Item::suggestSmallerOrAcceptCollisionStepSize(double)
{
	throw SimulationStoppingError("To be implemented.");
	return 0.0;
}

double Item::suggestSmallerOrAcceptDynFirstStepSize(double hnew)
{
	//"Default is return hnew."
	//"Best to do nothing so as not to disrupt the starting algorithm."
	return hnew;
}

double Item::suggestSmallerOrAcceptDynStepSize(double hnew)
{
	//"Default is return hnew."
	return hnew;
}

void Item::preVelIC()
{
	//"Assume positions are valid."
	//"Called once before solving for velocity initial conditions."
	//"Update all variable dependent instance variables needed for velIC even if they have 
	//been calculated in postPosIC."
	//"Variables dependent on t are updated."

	this->calcPostDynCorrectorIteration();
}

void Item::postVelIC()
{
}

void Item::fillqsudot(FColDsptr)
{
}

void MbD::Item::fillqsudotPlam(FColDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::fillqsudotPlamDeriv(FColDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::fillqsudotWeights(DiagMatDsptr)
{
}

void Item::fillVelICError(FColDsptr)
{
}

void Item::fillVelICJacob(SpMatDsptr)
{
}

void Item::setqsudotlam(FColDsptr)
{
}

void MbD::Item::setqsudotPlam(FColDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::setqsudotPlamDeriv(FColDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::preAccIC()
{
	this->calcPostDynCorrectorIteration();
}

void MbD::Item::preCollision()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::preCollisionCorrector()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::preCollisionCorrectorIteration()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::preCollisionDerivativeIC()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::preCollisionPredictor()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::preCollisionStep()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::postAccIC()
{
}

void Item::postAccICIteration()
{
}

void MbD::Item::postCollisionCorrector()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::postCollisionCorrectorIteration()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::postCollisionDerivativeIC()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::postCollisionPredictor()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::postCollisionStep()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::fillqsuddotlam(FColDsptr)
{
}

void Item::fillAccICIterError(FColDsptr)
{
}

void Item::fillAccICIterJacob(SpMatDsptr)
{
}

void MbD::Item::fillCollisionDerivativeICError(FColDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::fillCollisionDerivativeICJacob(SpMatDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::fillCollisionError(FColDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::fillCollisionpFpy(SpMatDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::fillCollisionpFpydot(SpMatDsptr)
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::setqsudot(FColDsptr)
{
}

void Item::setqsuddotlam(FColDsptr)
{
}

std::shared_ptr<StateData> Item::stateData()
{
	throw SimulationStoppingError("To be implemented.");
	return std::shared_ptr<StateData>();
}

void MbD::Item::storeCollisionState()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::discontinuityAtaddTypeTo(double, std::shared_ptr<std::vector<DiscontinuityType>>)
{
}

void MbD::Item::discontinuityAtICAddTo(std::shared_ptr<std::vector<DiscontinuityType>>)
{
	throw SimulationStoppingError("To be implemented.");
}

double Item::checkForDynDiscontinuityBetweenand(double, double t)
{
	//"Check for discontinuity in the last step defined by the interval (tprevious,t]."
	//"Default is assume no discontinuity and return t."

	return t;
}

void Item::constraintsReport()
{
}

void Item::setqsu(FColDsptr)
{
}

void Item::useEquationNumbers()
{
}

void Item::logString(const std::string& str)
{
	this->root()->logString(str);
}

void MbD::Item::logStringwithArgument(const std::string&, const std::string&)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::logStringwithArguments(const std::string&, std::shared_ptr<std::vector<std::string>>)
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::normalImpulse(double)
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::prePosIC()
{
	//"Called once before solving for position initial conditions."
	//"Update all variable dependent instance variables needed for posIC."
	//"This is a subset of update."

	calcPostDynCorrectorIteration();
}

void Item::prePosKine()
{
	this->prePosIC();
}

void MbD::Item::preStatic()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::postPosIC()
{
}

void Item::postPosICIteration()
{
	this->calcPostDynCorrectorIteration();
}

void MbD::Item::postStatic()
{
	throw SimulationStoppingError("To be implemented.");
}

void MbD::Item::postStaticIteration()
{
	throw SimulationStoppingError("To be implemented.");
}

void Item::fillPosICError(FColDsptr)
{
}

void Item::fillPosICJacob(FMatDsptr)
{
}

void Item::fillPosICJacob(SpMatDsptr)
{
}
