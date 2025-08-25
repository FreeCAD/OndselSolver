/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <algorithm>
#include "FunctionXY.h"
#include "Sum.h"
#include "Constant.h"

using namespace MbD;

MbD::FunctionXY::FunctionXY()
{
}

MbD::FunctionXY::FunctionXY(Symsptr base, Symsptr exp) : x(base), y(exp)
{
}

Symsptr MbD::FunctionXY::copyWith(Symsptr argx, Symsptr argy)
{
	return Symsptr();
}

void MbD::FunctionXY::arguments(Symsptr args)
{
	//args is a Sum with "terms" containing the actual arguments
	auto sum = std::static_pointer_cast<Sum>(args);
	assert(sum->terms->size() == 2);
	x = sum->terms->at(0);
	y = sum->terms->at(1);
}

Symsptr MbD::FunctionXY::expandUntil(Symsptr sptr, std::shared_ptr<std::unordered_set<Symsptr>> set)
{
	auto itr = std::find_if(set->begin(), set->end(), [sptr](Symsptr sym) {return sptr.get() == sym.get(); });
	if (itr != set->end()) return sptr;
	auto newx = x->expandUntil(x, set);
	auto newy = y->expandUntil(y, set);
	auto copy = copyWith(newx, newy);
	return copy;
}

Symsptr MbD::FunctionXY::simplifyUntil(Symsptr sptr, std::shared_ptr<std::unordered_set<Symsptr>> set)
{
	auto itr = std::find_if(set->begin(), set->end(), [sptr](Symsptr sym) {return sptr.get() == sym.get(); });
	if (itr != set->end()) return sptr;
	auto newx = x->simplifyUntil(x, set);
	auto newy = y->simplifyUntil(y, set);
	auto copy = copyWith(newx, newy);
	return copy;
}

void MbD::FunctionXY::createMbD(std::shared_ptr<System> mbdSys, std::shared_ptr<Units> mbdUnits)
{
	x->createMbD(mbdSys, mbdUnits);
	y->createMbD(mbdSys, mbdUnits);
}

bool MbD::FunctionXY::isConstant()
{
	return x->isConstant() && y->isConstant();
}

Symsptr MbD::FunctionXY::differentiateWRT(Symsptr var)
{
	if (this == var.get()) return sptrConstant(1.0);
	auto dfdx = differentiateWRTx();
	auto dfdy = differentiateWRTy();
	auto dxdvar = x->differentiateWRT(var);
	auto dydvar = y->differentiateWRT(var);
	return Symbolic::sum(Symbolic::times(dfdx, dxdvar), Symbolic::times(dfdy, dydvar));
}

