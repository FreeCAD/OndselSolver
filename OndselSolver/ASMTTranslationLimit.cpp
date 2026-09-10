#include "ASMTTranslationLimit.h"
#include "SymbolicParser.h"
#include "BasicUserFunction.h"
#include "Constant.h"
#include "TranslationLimitIJ.h"
#include "Units.h"

using namespace MbD;

std::shared_ptr<ASMTTranslationLimit> MbD::ASMTTranslationLimit::With()
{
	auto translationLimit = std::make_shared<ASMTTranslationLimit>();
	translationLimit->initialize();
	return translationLimit;
}

double MbD::ASMTTranslationLimit::coordinateUnit(const Units& units) const { return units.length; }
double MbD::ASMTTranslationLimit::stiffnessUnit(const Units& units) const { return units.length / units.force; }
double MbD::ASMTTranslationLimit::dampingUnit(const Units& units) const { return units.velocity / units.force; }

std::shared_ptr<ItemIJ> MbD::ASMTTranslationLimit::mbdClassNew()
{
	return TranslationLimitIJ::With();
}

void MbD::ASMTTranslationLimit::storeOnLevel(std::ofstream& os, size_t level)
{
	storeOnLevelString(os, level, "TranslationLimit");
	ASMTLimit::storeOnLevel(os, level);
}
