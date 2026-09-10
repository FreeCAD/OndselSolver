// SPDX-License-Identifier: LGPL-2.1-or-later
#include "ASMTDistanceLimit.h"

#include "DistanceLimitIJ.h"
#include "Units.h"

using namespace MbD;

std::shared_ptr<ASMTDistanceLimit> ASMTDistanceLimit::With()
{
    auto result = std::make_shared<ASMTDistanceLimit>();
    result->initialize();
    return result;
}

std::shared_ptr<ItemIJ> ASMTDistanceLimit::mbdClassNew()
{
    return DistanceLimitIJ::With();
}

void ASMTDistanceLimit::storeOnLevel(std::ofstream& os, size_t level)
{
    storeOnLevelString(os, level, "DistanceLimit");
    ASMTLimit::storeOnLevel(os, level);
}

double ASMTDistanceLimit::coordinateUnit(const Units& units) const
{
    return units.length;
}

double ASMTDistanceLimit::stiffnessUnit(const Units& units) const
{
    return units.length / units.force;
}

double ASMTDistanceLimit::dampingUnit(const Units& units) const
{
    return units.velocity / units.force;
}
