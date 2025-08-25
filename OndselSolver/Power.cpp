/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#include <algorithm>
#include "Power.h"
#include "Constant.h"
#include "Ln.h"

using namespace MbD;

MbD::Power::Power()
{
}

MbD::Power::Power(Symsptr bse, Symsptr ex) : FunctionXY(bse, ex)
{
}

Symsptr MbD::Power::differentiateWRTx()
{
	auto yminus1 = Symbolic::sum(y, sptrConstant(-1.0));
	auto power = Symbolic::raisedTo(x, yminus1);
	auto deriv = Symbolic::times(y, power);
	return deriv->simplified(deriv);
}

Symsptr MbD::Power::differentiateWRTy()
{
	auto lnterm = std::make_shared<Ln>(x);
	auto deriv = Symbolic::times(clonesptr(), lnterm);
	return deriv->simplified();
}

Symsptr MbD::Power::copyWith(Symsptr argx, Symsptr argy)
{
	return std::make_shared<Power>(argx, argy);
}

double MbD::Power::getValue()
{
	return std::pow(x->getValue(), y->getValue());
}

Symsptr MbD::Power::clonesptr()
{
	return std::make_shared<Power>(*this);
}

Symsptr MbD::Power::simplifyUntil(Symsptr sptr, std::shared_ptr<std::unordered_set<Symsptr>> set)
{
	auto itr = std::find_if(set->begin(), set->end(), [sptr](Symsptr sym) {return sptr.get() == sym.get(); });
	if (itr != set->end()) return sptr;
	auto newx = x->simplifyUntil(x, set);
	auto newy = y->simplifyUntil(y, set);
	if (y->isConstant() && y->getValue() == 1) {
		return newx;
	}
	auto copy = copyWith(newx, newy);
	return copy;
}

std::ostream& MbD::Power::printOn(std::ostream& s) const
{
	s << "pow(" << *x << "," << *y << ")";
	return s;
}
