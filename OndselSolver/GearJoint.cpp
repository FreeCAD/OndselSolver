/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "CREATE.h"
#include "CoaxialGearConstraintIJ.h"
#include "EndFramec.h"
#include "FullColumn.h"
#include "GearConstraintIJ.h"
#include "GearJoint.h"
#include "MarkerFrame.h"
#include "PartFrame.h"
#include "SimulationStoppingError.h"
#include "System.h"

#include <cmath>
#include <limits>

using namespace MbD;

namespace
{

constexpr double coaxialTolerance = 1.0e-9;

FColDsptr frameOriginInO(const EndFrmsptr& frame)
{
	auto markerFrame = frame->getMarkerFrame();
	auto partFrame = markerFrame->getPartFrame();
	// Gear constraints are installed before end-frame caches are refreshed.
	return partFrame->rOpO()->plusFullColumn(partFrame->aAOp()->timesFullColumn(markerFrame->rpmp));
}

FColDsptr frameZAxisInO(const EndFrmsptr& frame)
{
	auto markerFrame = frame->getMarkerFrame();
	auto partFrame = markerFrame->getPartFrame();
	return partFrame->aAOp()->timesFullMatrix(markerFrame->aApm)->column(2);
}

bool framesAreCoaxial(const EndFrmsptr& frameI, const EndFrmsptr& frameJ)
{
	auto axisI = frameZAxisInO(frameI);
	auto axisJ = frameZAxisInO(frameJ);
	const double axisILength = axisI->length();
	const double axisJLength = axisJ->length();
	if (axisILength <= coaxialTolerance || axisJLength <= coaxialTolerance) {
		return false;
	}

	const double cosAngle = axisI->dot(axisJ) / (axisILength * axisJLength);
	if (1.0 - std::abs(cosAngle) > coaxialTolerance) {
		return false;
	}

	auto delta = frameOriginInO(frameJ)->minusFullColumn(frameOriginInO(frameI));
	return delta->cross(axisI)->length() / axisILength <= coaxialTolerance;
}

}  // namespace

MbD::GearJoint::GearJoint()
{
}

MbD::GearJoint::GearJoint(const std::string& str) : Joint(str)
{
}
//
//void MbD::GearJoint::initializeLocally()
//{
//	if (!constraints->empty())
//	{
//		auto constraint = std::static_pointer_cast<GearConstraintIJ>(constraints->back());
//		constraint->initorbitsIJ();
//	}
//	Joint::initializeLocally();
//}

void MbD::GearJoint::setCarrierFrame(EndFrmsptr carrierFrame)
{
	frmK = carrierFrame;
}

void MbD::GearJoint::initializeGlobally()
{
	if (constraints->empty())
	{
		if (framesAreCoaxial(frmI, frmJ)) {
			if (!frmK) {
				throw SimulationStoppingError(
					"Coaxial GearJoint requires a carrier marker."
				);
			}

			auto gearIJK = CoaxialGearConstraintIJ::With(frmI, frmJ, frmK);
			gearIJK->radiusI = radiusI;
			gearIJK->radiusJ = radiusJ;
			gearIJK->setConstant(std::numeric_limits<double>::min());
			addConstraint(gearIJK);
		}
		else {
			auto gearIJ = GearConstraintIJ::With(frmI, frmJ);
			gearIJ->radiusI = radiusI;
			gearIJ->radiusJ = radiusJ;
			gearIJ->setConstant(std::numeric_limits<double>::min());
			addConstraint(gearIJ);
		}
		this->root()->hasChanged = true;
	}
	else {
		Joint::initializeGlobally();
	}
}
