/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once
#include <limits>
#include <memory>

namespace MbD {
    class SolverStatistics
    {
    public:
        static std::shared_ptr<SolverStatistics> With();

        size_t iterNo = SIZE_MAX;
        size_t corIterNo = SIZE_MAX;
        double h = std::numeric_limits<double>::min();
        size_t istep = SIZE_MAX;
        size_t order = SIZE_MAX;
        double t = std::numeric_limits<double>::min();
        double truncError = std::numeric_limits<double>::min();
    };
}
