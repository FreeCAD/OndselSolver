/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#pragma once

#include "Joint.h"

namespace MbD {
    class AtPointJoint : public Joint
    {
        //
    public:
        AtPointJoint();
        AtPointJoint(const std::string& str);
        static std::shared_ptr<AtPointJoint> With();
        static std::shared_ptr<AtPointJoint> With(const char* str);

        void createAtPointConstraints();


    };
}

