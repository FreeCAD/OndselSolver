/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "MbDMath.h"

void MbD::MbDMath::noop()
{
	//No Operations
}

std::string MbD::MbDMath::xyzFromInt(int i)
{
	switch (i) {
		case 0:
			return "x";	
		case 1:
			return "y";
		case 2:
			return "z";
		default:
			return std::to_string(i);
	}
}

std::string MbD::MbDMath::XYZFromInt(int i)
{
	switch (i) {
	case 0:
		return "X";
	case 1:
		return "Y";
	case 2:
		return "Z";
	default:
		return std::to_string(i);
	}
}
