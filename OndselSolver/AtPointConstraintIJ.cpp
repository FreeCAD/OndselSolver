/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "AtPointConstraintIJ.h"
#include "DispCompIecJecO.h"
#include "CREATE.h"

using namespace MbD;

AtPointConstraintIJ::AtPointConstraintIJ(EndFrmsptr frmi, EndFrmsptr frmj, size_t axisi) :
	ConstraintIJ(frmi, frmj), axis(axisi)
{
}

void AtPointConstraintIJ::initialize()
{
	ConstraintIJ::initialize();
	initriIeJeO();
}

void AtPointConstraintIJ::initializeLocally()
{
	riIeJeO->initializeLocally();
}

void AtPointConstraintIJ::initializeGlobally()
{
	riIeJeO->initializeGlobally();
}

void AtPointConstraintIJ::initriIeJeO()
{
	riIeJeO = CREATE<DispCompIecJecO>::With(frmI, frmJ, axis);
}

void AtPointConstraintIJ::postInput()
{
	riIeJeO->postInput();
	Constraint::postInput();
}

void AtPointConstraintIJ::calcPostDynCorrectorIteration()
{
	aG = riIeJeO->riIeJeO - aConstant;
}

void AtPointConstraintIJ::prePosIC()
{
	riIeJeO->prePosIC();
	Constraint::prePosIC();
}

ConstraintType AtPointConstraintIJ::type()
{
	return displacement;
}

void AtPointConstraintIJ::postPosICIteration()
{
	riIeJeO->postPosICIteration();
	Item::postPosICIteration();
}

void AtPointConstraintIJ::preVelIC()
{
	riIeJeO->preVelIC();
	Item::preVelIC();
}

void AtPointConstraintIJ::preAccIC()
{
	riIeJeO->preAccIC();
	Constraint::preAccIC();
}

void AtPointConstraintIJ::postDynPredictor()
{
    riIeJeO->postDynPredictor();
    ConstraintIJ::postDynPredictor();
}

void AtPointConstraintIJ::postDynCorrectorIteration()
{
    riIeJeO->postDynCorrectorIteration();
    ConstraintIJ::postDynCorrectorIteration();
}

void AtPointConstraintIJ::preDynOutput()
{
    riIeJeO->preDynOutput();
    ConstraintIJ::preDynOutput();
}

void AtPointConstraintIJ::postDynOutput()
{
    riIeJeO->postDynOutput();
    ConstraintIJ::postDynOutput();
}
