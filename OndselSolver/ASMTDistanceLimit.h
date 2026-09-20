// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include "ASMTLimit.h"

namespace MbD
{
class ASMTDistanceLimit: public ASMTLimit
{
public:
    static std::shared_ptr<ASMTDistanceLimit> With();
    std::shared_ptr<ItemIJ> mbdClassNew() override;
    void storeOnLevel(std::ofstream& os, size_t level) override;
    double coordinateUnit(const Units& units) const override;
    double stiffnessUnit(const Units& units) const override;
    double dampingUnit(const Units& units) const override;
};
}  // namespace MbD
