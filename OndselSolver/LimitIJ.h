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
		void fillPosICError(FColDsptr col) override;
		void fillPosICJacob(SpMatDsptr mat) override;
		void fillqsulam(FColDsptr col) override;
		void setqsulam(FColDsptr col) override;
		void useEquationNumbers() override;

		bool satisfied() const;
		void deactivate();
		void activate();

		double limit = std::numeric_limits<double>::max();
		double tol = std::numeric_limits<double>::max();
		std::string type;
		bool active = false;
	};
}
