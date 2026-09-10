/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include "ASMTConstraintSet.h"

#include <array>

#include "ASMTAssembly.h"
#include "ASMTMarker.h"
#include "AppliedForceTorque.h"
#include "FullMatrix.h"
#include "Joint.h"

using namespace MbD;

void MbD::ASMTConstraintSet::initialize()
{
    ASMTItemIJ::initialize();
    powers = std::make_shared<FullRow<double>>();
}

void MbD::ASMTConstraintSet::clearResults()
{
    ASMTItemIJ::clearResults();
    if (powers) {
        powers->clear();
    }
}

void MbD::ASMTConstraintSet::updateFromMbD()
{
    // MbD returns force and torque at connector I. Convert them to the marker-I
    // moment expected by the ASMT result format before saving the reaction.
    auto mbdUnts = mbdUnits();
    auto mbdJoint = std::static_pointer_cast<Joint>(mbdObject);
    auto aFIeO = mbdJoint->aFX()->times(mbdUnts->force);
    auto aTIeO = mbdJoint->aTX()->times(mbdUnts->torque);
    auto rImIeO = mbdJoint->frmI->rmeO()->times(mbdUnts->length);
    auto aFIO = aFIeO;
    auto aTIO = aTIeO->plusFullColumn(rImIeO->cross(aFIeO));
    fxs->push_back(aFIO->at(0));
    fys->push_back(aFIO->at(1));
    fzs->push_back(aFIO->at(2));
    txs->push_back(aTIO->at(0));
    tys->push_back(aTIO->at(1));
    tzs->push_back(aTIO->at(2));

    const auto i = AppliedForceTorque::frameState(mbdJoint->frmI);
    const auto j = AppliedForceTorque::frameState(mbdJoint->frmJ);
    const auto dot = [](const auto& first, const auto& second) {
        return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
    };
    const std::array<double, 3> force {
        aFIeO->at(0) / mbdUnts->force,
        aFIeO->at(1) / mbdUnts->force,
        aFIeO->at(2) / mbdUnts->force,
    };
    const std::array<double, 3> torque {
        aTIeO->at(0) / mbdUnts->torque,
        aTIeO->at(1) / mbdUnts->torque,
        aTIeO->at(2) / mbdUnts->torque,
    };
    const std::array<double, 3> separation {
        j.position[0] - i.position[0],
        j.position[1] - i.position[1],
        j.position[2] - i.position[2],
    };
    const std::array<double, 3> oppositeTorque {
        separation[1] * force[2] - separation[2] * force[1] - torque[0],
        separation[2] * force[0] - separation[0] * force[2] - torque[1],
        separation[0] * force[1] - separation[1] * force[0] - torque[2],
    };
    const double power = dot(force, i.velocity) + dot(torque, i.omega)
        - dot(force, j.velocity) + dot(oppositeTorque, j.omega);
    powers->push_back(power * mbdUnts->torque / mbdUnts->angle / mbdUnts->time);
}

void MbD::ASMTConstraintSet::compareResults(AnalysisType)
{
    if (infxs == nullptr || infxs->empty()) {
        return;
    }
}

void MbD::ASMTConstraintSet::outputResults(AnalysisType)
{}
