// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include "LimitIJ.h"

namespace MbD
{
class DistanceLimitIJ: public LimitIJ
{
public:
    static std::shared_ptr<DistanceLimitIJ> With();
    void initializeGlobally() override;
};
}  // namespace MbD
