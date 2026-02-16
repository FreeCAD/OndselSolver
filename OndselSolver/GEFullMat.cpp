/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include <cassert>

#include "GEFullMat.h"

using namespace MbD;

void GEFullMat::forwardEliminateWithPivot(size_t)
{
	throw SimulationStoppingError("To be implemented.");
}

void GEFullMat::backSubstituteIntoDU()
{
	answerX = std::make_shared<FullColumn<double>>(n);
	answerX->at(n - 1) = rightHandSideB->at(m - 1) / matrixA->at(m - 1)->at(n - 1);
	for (ssize_t i = (ssize_t)n - 2; i >= 0; i--)	//Use ssize_t because of decrement
	{
		auto rowi = matrixA->at(i);
		double sum = answerX->at(n - 1) * rowi->at(n - 1);
		for (size_t j = i + 1; j < n - 1; j++)
		{
			sum += answerX->at(j) * rowi->at(j);
		}
		answerX->at(i) = (rightHandSideB->at(i) - sum) / rowi->at(i);
	}
}

void GEFullMat::postSolve()
{
	throw SimulationStoppingError("To be implemented.");
}

void GEFullMat::preSolvewithsaveOriginal(FMatDsptr, FColDsptr, bool)
{
	throw SimulationStoppingError("To be implemented.");
}

void GEFullMat::preSolvewithsaveOriginal(SpMatDsptr, FColDsptr, bool)
{
	throw SimulationStoppingError("To be implemented.");
}

double GEFullMat::getmatrixArowimaxMagnitude(size_t i)
{
	return matrixA->at(i)->maxMagnitude();
}

FColDsptr GEFullMat::basicSolvewithsaveOriginal(FMatDsptr fullMat, FColDsptr fullCol, bool saveOriginal)
{
	this->preSolvewithsaveOriginal(fullMat, fullCol, saveOriginal);
	for (size_t p = 0; p < m; p++)
	{
		this->doPivoting(p);
		this->forwardEliminateWithPivot(p);
	}
	this->backSubstituteIntoDU();
	this->postSolve();
	return answerX;
}

FColDsptr GEFullMat::basicSolvewithsaveOriginal(SpMatDsptr, FColDsptr, bool)
{
	throw SimulationStoppingError("To be implemented.");
	return FColDsptr();
}
