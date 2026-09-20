// SPDX-License-Identifier: LGPL-2.1-or-later
#include "DistanceLimitIJ.h"

#include "DistanceConstraintIJ.h"
#include "System.h"

using namespace MbD;

std::shared_ptr<DistanceLimitIJ> DistanceLimitIJ::With()
{
    auto result = std::make_shared<DistanceLimitIJ>();
    result->initialize();
    return result;
}

void DistanceLimitIJ::initializeGlobally()
{
    if (constraints->empty()) {
        auto constraint = DistanceConstraintIJ::With(frmI, frmJ);
        constraint->setConstant(limit);
        addConstraint(constraint);
        root()->hasChanged = true;
    }
    else {
        LimitIJ::initializeGlobally();
    }
}
