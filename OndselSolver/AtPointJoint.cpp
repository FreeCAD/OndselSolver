/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "AtPointJoint.h"
#include "System.h"
#include "CREATE.h"

using namespace MbD;

MbD::AtPointJoint::AtPointJoint()
{
}

MbD::AtPointJoint::AtPointJoint(const std::string& str) : Joint(str)
{
}

std::shared_ptr<AtPointJoint> MbD::AtPointJoint::With()
{
	auto inst = std::make_shared<AtPointJoint>();
	inst->initialize();
	return inst;
}

std::shared_ptr<AtPointJoint> MbD::AtPointJoint::With(const char* str)
{
	auto inst = std::make_shared<AtPointJoint>(str);
	inst->initialize();
	return inst;
}

void MbD::AtPointJoint::createAtPointConstraints()
{
	addConstraint(CREATE<AtPointConstraintIqcJqc>::ConstraintWith(frmI, frmJ, 0));
	addConstraint(CREATE<AtPointConstraintIqcJqc>::ConstraintWith(frmI, frmJ, 1));
	addConstraint(CREATE<AtPointConstraintIqcJqc>::ConstraintWith(frmI, frmJ, 2));
}
