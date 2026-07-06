/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include "CoaxialGearConstraintIJ.h"

#include <cassert>
#include <limits>
#include <utility>

#include "AngleZIeqcJeqc.h"
#include "EndFrameqc.h"
#include "FullColumn.h"
#include "FullMatrix.h"
#include "FullRow.h"
#include "SimulationStoppingError.h"

using namespace MbD;

MbD::CoaxialGearConstraintIJ::CoaxialGearConstraintIJ(
    EndFrmsptr frmi,
    EndFrmsptr frmj,
    EndFrmsptr frmk
) : frmI(std::move(frmi)), frmJ(std::move(frmj)), frmK(std::move(frmk))
{}

std::shared_ptr<CoaxialGearConstraintIJ> MbD::CoaxialGearConstraintIJ::With(
    EndFrmsptr frmi,
    EndFrmsptr frmj,
    EndFrmsptr frmk
)
{
    assert(frmi->isEndFrameqc());
    assert(frmj->isEndFrameqc());
    assert(frmk->isEndFrameqc());

    auto constraint = std::make_shared<CoaxialGearConstraintIJ>(frmi, frmj, frmk);
    constraint->initialize();
    return constraint;
}

void MbD::CoaxialGearConstraintIJ::calcPostDynCorrectorIteration()
{
    angleKI->calcPostDynCorrectorIteration();
    angleKJ->calcPostDynCorrectorIteration();

    const double gearRatio = ratio();
    aG = angleKJ->value() + (gearRatio * angleKI->value()) - aConstant;

    pGpEI = angleKI->pvaluepEJ()->times(gearRatio);
    pGpEJ = angleKJ->pvaluepEJ();
    pGpEK = angleKJ->pvaluepEI()->plusFullRow(angleKI->pvaluepEI()->times(gearRatio));

    ppGpEIpEI = angleKI->ppvaluepEJpEJ()->times(gearRatio);
    ppGpEJpEJ = angleKJ->ppvaluepEJpEJ();
    ppGpEKpEK = angleKJ->ppvaluepEIpEI()->plusFullMatrix(
        angleKI->ppvaluepEIpEI()->times(gearRatio)
    );
    ppGpEKpEI = angleKI->ppvaluepEIpEJ()->times(gearRatio);
    ppGpEKpEJ = angleKJ->ppvaluepEIpEJ();
}

void MbD::CoaxialGearConstraintIJ::fillAccICIterError(FColDsptr col)
{
    col->atiplusFullVectortimes(iqEI, pGpEI, lam);
    col->atiplusFullVectortimes(iqEJ, pGpEJ, lam);
    col->atiplusFullVectortimes(iqEK, pGpEK, lam);

    auto frmIeqc = std::static_pointer_cast<EndFrameqc>(frmI);
    auto frmJeqc = std::static_pointer_cast<EndFrameqc>(frmJ);
    auto frmKeqc = std::static_pointer_cast<EndFrameqc>(frmK);

    auto qEdotI = frmIeqc->qEdot();
    auto qEdotJ = frmJeqc->qEdot();
    auto qEdotK = frmKeqc->qEdot();

    double sum = 0.0;
    sum += pGpEI->timesFullColumn(frmIeqc->qEddot());
    sum += pGpEJ->timesFullColumn(frmJeqc->qEddot());
    sum += pGpEK->timesFullColumn(frmKeqc->qEddot());
    sum += qEdotI->transposeTimesFullColumn(ppGpEIpEI->timesFullColumn(qEdotI));
    sum += qEdotJ->transposeTimesFullColumn(ppGpEJpEJ->timesFullColumn(qEdotJ));
    sum += qEdotK->transposeTimesFullColumn(ppGpEKpEK->timesFullColumn(qEdotK));
    sum += 2.0
        * qEdotK->transposeTimesFullColumn(ppGpEKpEI->timesFullColumn(qEdotI));
    sum += 2.0
        * qEdotK->transposeTimesFullColumn(ppGpEKpEJ->timesFullColumn(qEdotJ));
    col->atiplusNumber(iG, sum);
}

void MbD::CoaxialGearConstraintIJ::fillPosICError(FColDsptr col)
{
    Constraint::fillPosICError(col);
    col->atiplusFullVectortimes(iqEI, pGpEI, lam);
    col->atiplusFullVectortimes(iqEJ, pGpEJ, lam);
    col->atiplusFullVectortimes(iqEK, pGpEK, lam);
}

void MbD::CoaxialGearConstraintIJ::fillPosICJacob(SpMatDsptr mat)
{
    mat->atijplusFullRow(iG, iqEI, pGpEI);
    mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
    mat->atijplusFullRow(iG, iqEJ, pGpEJ);
    mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
    mat->atijplusFullRow(iG, iqEK, pGpEK);
    mat->atijplusFullColumn(iqEK, iG, pGpEK->transpose());

    mat->atijplusFullMatrixtimes(iqEI, iqEI, ppGpEIpEI, lam);
    mat->atijplusFullMatrixtimes(iqEJ, iqEJ, ppGpEJpEJ, lam);
    mat->atijplusFullMatrixtimes(iqEK, iqEK, ppGpEKpEK, lam);

    auto ppGpEKpEIlam = ppGpEKpEI->times(lam);
    mat->atijplusFullMatrix(iqEK, iqEI, ppGpEKpEIlam);
    mat->atijplusTransposeFullMatrix(iqEI, iqEK, ppGpEKpEIlam);

    auto ppGpEKpEJlam = ppGpEKpEJ->times(lam);
    mat->atijplusFullMatrix(iqEK, iqEJ, ppGpEKpEJlam);
    mat->atijplusTransposeFullMatrix(iqEJ, iqEK, ppGpEKpEJlam);
}

void MbD::CoaxialGearConstraintIJ::fillPosKineJacob(SpMatDsptr mat)
{
    mat->atijplusFullRow(iG, iqEI, pGpEI);
    mat->atijplusFullRow(iG, iqEJ, pGpEJ);
    mat->atijplusFullRow(iG, iqEK, pGpEK);
}

void MbD::CoaxialGearConstraintIJ::fillVelICJacob(SpMatDsptr mat)
{
    mat->atijplusFullRow(iG, iqEI, pGpEI);
    mat->atijplusFullColumn(iqEI, iG, pGpEI->transpose());
    mat->atijplusFullRow(iG, iqEJ, pGpEJ);
    mat->atijplusFullColumn(iqEJ, iG, pGpEJ->transpose());
    mat->atijplusFullRow(iG, iqEK, pGpEK);
    mat->atijplusFullColumn(iqEK, iG, pGpEK->transpose());
}

void MbD::CoaxialGearConstraintIJ::initialize()
{
    Constraint::initialize();
    angleKI = AngleZIeqcJeqc::With(frmK, frmI);
    angleKJ = AngleZIeqcJeqc::With(frmK, frmJ);
}

void MbD::CoaxialGearConstraintIJ::initializeGlobally()
{
    angleKI->initializeGlobally();
    angleKJ->initializeGlobally();
}

void MbD::CoaxialGearConstraintIJ::initializeLocally()
{
    angleKI->initializeLocally();
    angleKJ->initializeLocally();
}

void MbD::CoaxialGearConstraintIJ::postInput()
{
    angleKI->postInput();
    angleKJ->postInput();
    if (aConstant == std::numeric_limits<double>::min()) {
        aConstant = angleKJ->value() + (ratio() * angleKI->value());
    }
    Constraint::postInput();
}

void MbD::CoaxialGearConstraintIJ::postPosICIteration()
{
    angleKI->postPosICIteration();
    angleKJ->postPosICIteration();
    Constraint::postPosICIteration();
}

void MbD::CoaxialGearConstraintIJ::preAccIC()
{
    angleKI->preAccIC();
    angleKJ->preAccIC();
    Constraint::preAccIC();
}

void MbD::CoaxialGearConstraintIJ::prePosIC()
{
    angleKI->prePosIC();
    angleKJ->prePosIC();
    Constraint::prePosIC();
}

void MbD::CoaxialGearConstraintIJ::preVelIC()
{
    angleKI->preVelIC();
    angleKJ->preVelIC();
    Constraint::preVelIC();
}

double MbD::CoaxialGearConstraintIJ::ratio() const
{
    if (radiusJ == 0.0) {
        throw SimulationStoppingError("Gear radius is zero.");
    }
    return radiusI / radiusJ;
}

void MbD::CoaxialGearConstraintIJ::simUpdateAll()
{
    angleKI->simUpdateAll();
    angleKJ->simUpdateAll();
    Constraint::simUpdateAll();
}

void MbD::CoaxialGearConstraintIJ::useEquationNumbers()
{
    auto frmIeqc = std::static_pointer_cast<EndFrameqc>(frmI);
    auto frmJeqc = std::static_pointer_cast<EndFrameqc>(frmJ);
    auto frmKeqc = std::static_pointer_cast<EndFrameqc>(frmK);

    iqEI = frmIeqc->iqE();
    iqEJ = frmJeqc->iqE();
    iqEK = frmKeqc->iqE();
}

std::string MbD::CoaxialGearConstraintIJ::constraintSpec()
{
    return "CoaxialGearConstraintIJK";
}
