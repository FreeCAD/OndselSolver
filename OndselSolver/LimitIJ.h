/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "ConstraintSet.h"

namespace MbD {
	class LimitIJ : public ConstraintSet
	{
		//
	public:
		LimitIJ() = default;
		void fillConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> allConstraints) override;
		void fillDispConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> constraints) override;
		void fillEssenConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> constraints) override;
		void fillPerpenConstraints(std::shared_ptr<std::vector<std::shared_ptr<Constraint>>> constraints) override;
		void fillPosICError(FColDsptr col) override;
		void fillPosICJacob(SpMatDsptr mat) override;
		void fillqsudot(FColDsptr col) override;
		void fillqsulam(FColDsptr col) override;
		void fillqsuddotlam(FColDsptr col) override;
		void fillVelICError(FColDsptr col) override;
		void fillVelICJacob(SpMatDsptr mat) override;
		void fillAccICIterError(FColDsptr col) override;
		void fillAccICIterJacob(SpMatDsptr mat) override;
		void fillpqsumu(FColDsptr col) override;
		void fillpqsumudot(FColDsptr col) override;
		void fillDynError(FColDsptr col) override;
		void fillpFpy(SpMatDsptr mat) override;
		void fillpFpydot(SpMatDsptr mat) override;
		void setqsulam(FColDsptr col) override;
		void setqsudotlam(FColDsptr col) override;
		void setqsuddotlam(FColDsptr col) override;
		void setpqsumu(FColDsptr col) override;
		void setpqsumudot(FColDsptr col) override;
		void useEquationNumbers() override;
		void preDyn() override;
		void preDynStep() override;
		double checkForDynDiscontinuityBetweenand(double tprev, double t) override;
		void discontinuityAtaddTypeTo(
			double t,
			std::shared_ptr<std::vector<DiscontinuityType>> disconTypes
		) override;

		bool satisfied() const;
		void deactivate();
		void activate();
		void setCompliant(double stiffness, double damping);
		double normalReaction() const;
		bool hasTensileReaction() const;
		struct EnergyState {
			double power = 0;
			double storedEnergy = 0;
			double dissipatedPower = 0;
		};
		EnergyState energyState() const;

		double limit = std::numeric_limits<double>::max();
		double tol = std::numeric_limits<double>::max();
		std::string type;
		bool active = false;
		bool compliant = false;
		double stiffness = 0.0;
		double damping = 0.0;
		double previousClearance = std::numeric_limits<double>::max();
		double transitionTime = std::numeric_limits<double>::quiet_NaN();
		double previousReaction = 0;
		bool releasing = false;
	};
}
