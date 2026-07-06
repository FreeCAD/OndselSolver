/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "Constraint.h"

namespace MbD {
    class AngleZIeqcJeqc;
    class EndFramec;
    using EndFrmsptr = std::shared_ptr<EndFramec>;

    class CoaxialGearConstraintIJ : public Constraint
    {
    public:
        CoaxialGearConstraintIJ(EndFrmsptr frmi, EndFrmsptr frmj, EndFrmsptr frmk);

        static std::shared_ptr<CoaxialGearConstraintIJ> With(
            EndFrmsptr frmi,
            EndFrmsptr frmj,
            EndFrmsptr frmk
        );

        void calcPostDynCorrectorIteration() override;
        void fillAccICIterError(FColDsptr col) override;
        void fillPosICError(FColDsptr col) override;
        void fillPosICJacob(SpMatDsptr mat) override;
        void fillPosKineJacob(SpMatDsptr mat) override;
        void fillVelICJacob(SpMatDsptr mat) override;
        void initialize() override;
        void initializeGlobally() override;
        void initializeLocally() override;
        void postInput() override;
        void postPosICIteration() override;
        void preAccIC() override;
        void prePosIC() override;
        void preVelIC() override;
        double ratio() const;
        void simUpdateAll() override;
        void useEquationNumbers() override;
        std::string constraintSpec() override;

        EndFrmsptr frmI, frmJ, frmK;
        std::shared_ptr<AngleZIeqcJeqc> angleKI, angleKJ;
        double radiusI = 0.0, radiusJ = 0.0;

        FRowDsptr pGpEI, pGpEJ, pGpEK;
        FMatDsptr ppGpEIpEI, ppGpEJpEJ, ppGpEKpEK, ppGpEKpEI, ppGpEKpEJ;
        size_t iqEI = SIZE_MAX, iqEJ = SIZE_MAX, iqEK = SIZE_MAX;
    };
}
