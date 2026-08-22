/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "PosICDragNewtonRaphson.h"

namespace MbD {
    class Part;

    class PosICDragLimitNewtonRaphson : public PosICDragNewtonRaphson
    {
        //Dragging with limits of constrained or under constrained system
        //Assume no redundant constraints
    public:
        static std::shared_ptr<PosICDragLimitNewtonRaphson> With();
        void preRun() override;
        void run() override;
    };
}
