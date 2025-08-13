/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "ParallelAxesJoint.h"
#include "CREATE.h"
#include "System.h"

using namespace MbD;

MbD::ParallelAxesJoint::ParallelAxesJoint()
{
}

MbD::ParallelAxesJoint::ParallelAxesJoint(const std::string& str) : Joint(str)
{
}

void MbD::ParallelAxesJoint::initializeGlobally()
{
	if (constraints->empty())
	{
		addConstraint(CREATE<DirectionCosineConstraintIqcJqc>::ConstraintWith(frmI, frmJ, 2, 0));
		addConstraint(CREATE<DirectionCosineConstraintIqcJqc>::ConstraintWith(frmI, frmJ, 2, 1));
		this->root()->hasChanged = true;
	}
	else {
		Joint::initializeGlobally();
	}
}
