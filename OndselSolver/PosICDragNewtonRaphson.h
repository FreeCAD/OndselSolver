/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "AnyPosICNewtonRaphson.h"

namespace MbD {
    class Part;

    class PosICDragNewtonRaphson : public AnyPosICNewtonRaphson
    {
        //Dragging of constrained or under constrained system
        //Assume no redundant constraints
        //Limits are used in dragging only
    public:
        static std::shared_ptr<PosICDragNewtonRaphson> With();

        void initializeGlobally() override;
        void fillY() override;
        void fillPyPx() override;
        void passRootToSystem() override;
        void assignEquationNumbers() override;

        void preRun() override;
        void askSystemToUpdate() override;
        void postRun() override;

        void setdragParts(std::shared_ptr<std::vector<std::shared_ptr<Part>>> _dragParts);

        std::shared_ptr<std::vector<std::shared_ptr<Part>>> dragParts;
    };
}

