/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#pragma once
#include <string>

namespace MbD {
	class MbDMath
	{
	public:
		void noop();
		static std::string xyzFromInt(int i);
		static std::string XYZFromInt(int i);
	};
}

