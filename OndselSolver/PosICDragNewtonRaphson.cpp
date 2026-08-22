/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include "PosICDragNewtonRaphson.h"
#include "SystemSolver.h"
#include "Part.h"
#include "Constraint.h"

using namespace MbD;

std::shared_ptr<PosICDragNewtonRaphson> MbD::PosICDragNewtonRaphson::With()
{
	auto newtonRaphson = std::make_shared<PosICDragNewtonRaphson>();
	newtonRaphson->initialize();
	return newtonRaphson;
}

void MbD::PosICDragNewtonRaphson::preRun()
{
	std::string str("MbD: Assembling system. ");
	system->logString(str);
	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->prePosIC(); });
}

void MbD::PosICDragNewtonRaphson::askSystemToUpdate()
{
	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->postPosICIteration(); });
}

void MbD::PosICDragNewtonRaphson::postRun()
{
	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->postPosIC(); });
}

void MbD::PosICDragNewtonRaphson::initializeGlobally()
{
	this->assignEquationNumbers();
	system->partsJointsMotionsLimitsForcesTorquesDo([&](std::shared_ptr<Item> item) { item->useEquationNumbers(); });
	this->createVectorsAndMatrices();
	matrixSolver = this->matrixSolverClassNew();
	iterMax = system->iterMaxPosKine;
	dxTol = system->errorTolPosKine;
	for (size_t i = 0; i < qsuWeights->size(); i++)
	{
		qsuWeights->at(i) = 1.0e3;	//minimum weight
	}
	for (auto& part : *dragParts) {
		auto iqX = part->iqX();
		for (size_t i = 0; i < 3; i++)
		{
			qsuWeights->at((size_t)iqX + i) = 1.0e6;	//maximum weight
		}
	}
	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
		item->fillqsu(qsuOld);
		item->fillqsuWeights(qsuWeights);
		item->fillqsulam(x);
		});
}

void MbD::PosICDragNewtonRaphson::fillY()
{
	auto newMinusOld = qsuOld->negated();
	newMinusOld->equalSelfPlusFullColumnAt(x, 0);
	y->zeroSelf();
	y->atiminusFullColumn(0, (qsuWeights->timesFullColumn(newMinusOld)));
	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
		item->fillPosICError(y);
		//std::cout << item->name << *y << std::endl;
		//noop();
		});
	//std::cout << "Final" << *y << std::endl;
}

void MbD::PosICDragNewtonRaphson::fillPyPx()
{
	pypx->zeroSelf();
	pypx->atijminusDiagonalMatrix(0, 0, qsuWeights);
	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) {
		item->fillPosICJacob(pypx);
		//std::cout << *(pypx->at(3)) << std::endl;
		});
	//std::cout << *pypx << std::endl;
}

void MbD::PosICDragNewtonRaphson::passRootToSystem()
{
	system->partsJointsMotionsLimitsDo([&](std::shared_ptr<Item> item) { item->setqsulam(x); });
}

void MbD::PosICDragNewtonRaphson::assignEquationNumbers()
{
	auto parts = system->parts();
	//auto contactEndFrames = system->contactEndFrames();
	//auto uHolders = system->uHolders();
	auto constraints = system->allConstraintsLimits();
	size_t eqnNo = 0;
	for (auto& part : *parts) {
		part->iqX(eqnNo);
		eqnNo = eqnNo + 3;
		part->iqE(eqnNo);
		eqnNo = eqnNo + 4;
	}
	//for (auto& endFrm : *contactEndFrames) {
	//	endFrm->is(eqnNo);
	//	eqnNo = eqnNo + endFrm->sSize();
	//}
	//for (auto& uHolder : *uHolders) {
	//	uHolder->iu(eqnNo);
	//	eqnNo += 1;
	//}
	auto nEqns = eqnNo;	//C++ uses index 0.
	nqsu = nEqns;
	for (auto& con : *constraints) {
		con->iG = eqnNo;
		eqnNo += 1;
	}
	//auto lastEqnNo = eqnNo - 1;
	nEqns = eqnNo;	//C++ uses index 0.
	n = nEqns;
}

void MbD::PosICDragNewtonRaphson::setdragParts(std::shared_ptr<std::vector<std::shared_ptr<Part>>> _dragParts)
{
	dragParts = _dragParts;
}
