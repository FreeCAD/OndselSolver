/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "LineInPlaneJoint.h"
#include "CREATE.h"
#include "System.h"

using namespace MbD;

MbD::LineInPlaneJoint::LineInPlaneJoint()
{
}

MbD::LineInPlaneJoint::LineInPlaneJoint(const std::string& str) : InPlaneJoint(str)
{
}

void MbD::LineInPlaneJoint::initializeGlobally()
{
	if (constraints->empty())
	{
		this->createInPlaneConstraint();
		addConstraint(CREATE<DirectionCosineConstraintIqcJqc>::ConstraintWith(frmI, frmJ, 2, 2));
		this->root()->hasChanged = true;
	}
	else {
		Joint::initializeGlobally();
	}
}
