/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include <cstdint>

#include "AngleZConstraintIJ.h"

namespace MbD {
	class AngleZConstraintIqcJc : public AngleZConstraintIJ
	{
		//pGpEI ppGpEIpEI iqEI 
	public:
		AngleZConstraintIqcJc(EndFrmsptr frmi, EndFrmsptr frmj);

		void initthezIeJe() override;
		void addToJointTorqueI(FColDsptr col) override;
		void calc_pGpEI();
		void calc_ppGpEIpEI();
		void calcPostDynCorrectorIteration() override;
		void fillAccICIterError(FColDsptr col) override;
		void fillPosICError(FColDsptr col) override;
		void fillPosICJacob(SpMatDsptr mat) override;
		void fillPosKineJacob(SpMatDsptr mat) override;
		void fillVelICJacob(SpMatDsptr mat) override;
		void fillpFpy(SpMatDsptr mat) override;
		void fillpFpydot(SpMatDsptr mat) override;
		double constraintVelocity() const override;
		void fillGeneralizedForce(FColDsptr col, double multiplier) override;
		void fillGeneralizedForcePositionJacobian(
			SpMatDsptr mat,
			double multiplier,
			double derivative
		) override;
		void fillGeneralizedForceVelocityJacobian(SpMatDsptr mat, double derivative) override;
		void useEquationNumbers() override;

		FRowDsptr pGpEI;
		FMatDsptr ppGpEIpEI;
		size_t iqEI = SIZE_MAX;

	};
}
