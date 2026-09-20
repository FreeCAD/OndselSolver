/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "StableBackwardDifference.h"

namespace MbD {


#pragma once
    class StableStartingBDF : public StableBackwardDifference
    {
        //
    public:
        static std::shared_ptr<StableStartingBDF> With();
        void initialize() override;

        void initializeLocally() override;
        double pvdotpv() override;
        void formTaylorMatrix() override;
        void setorder(size_t o) override;
        void setiStep(size_t i) override;
        FColDsptr derivativepresentpastpresentDerivativepastDerivative(size_t n,
            FColDsptr y, std::shared_ptr<std::vector<FColDsptr>> ypast,
            FColDsptr ydot, std::shared_ptr<std::vector<FColDsptr>> ydotpast);
        FColDsptr derivativewith(size_t deriv, std::shared_ptr<std::vector<FColDsptr>> series);


    };
}
