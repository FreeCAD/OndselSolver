/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include "Extrapolator.h"
#include "FullColumn.h"

using namespace MbD;

std::shared_ptr<Extrapolator> Extrapolator::With()
{
    auto inst = std::make_shared<Extrapolator>();
    inst->initialize();
    return inst;
}

void Extrapolator::formTaylorMatrix()
{
    //"
    //For method order 3:
    //|    1    (t1 - t)    (t1 - t)^2/2!        (t1 - t)^3/3!    |    |q(t)    |    =    |q(t1)    |
    //|    1    (t2 - t)    (t2 - t)^2/2!        (t2 - t)^3/3!    |    |qd(t)    |        |q(t2)    |
    //|    1    (t3 - t)    (t3 - t)^2/2!        (t3 - t)^3/3!    |    |qdd(t)    |        |q(t3)    |
    //|    1    (t4 - t)    (t4 - t)^2/2!        (t4 - t)^3/3!    |    |qddd(t)|        |q(t4)    |
    //"

    instantiateTaylorMatrix();
    for (size_t i = 0; i < order + 1; i++)
    {
        formTaylorRowwithTimeNodederivative(i, i, 0);
    }
}
