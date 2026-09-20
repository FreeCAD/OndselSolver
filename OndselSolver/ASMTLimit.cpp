#include "ASMTLimit.h"
#include "ASMTAssembly.h"
#include "ASMTJoint.h"
#include "SymbolicParser.h"
#include "BasicUserFunction.h"
#include "Constant.h"
#include "Units.h"
#include "LimitIJ.h"

#include <sstream>

using namespace MbD;

void ASMTLimit::initialize()
{
	ASMTConstraintSet::initialize();
	storedEnergies = std::make_shared<FullRow<double>>();
	dissipatedPowers = std::make_shared<FullRow<double>>();
}

void ASMTLimit::clearResults()
{
	ASMTConstraintSet::clearResults();
	if (storedEnergies) storedEnergies->clear();
	if (dissipatedPowers) dissipatedPowers->clear();
}

void ASMTLimit::updateFromMbD()
{
	const auto state = std::static_pointer_cast<LimitIJ>(mbdObject)->energyState();
	const auto units = mbdUnits();
	const double energyUnit = units->torque / units->angle;
	const double powerUnit = energyUnit / units->time;
	powers->push_back(state.power * powerUnit);
	storedEnergies->push_back(state.storedEnergy * energyUnit);
	dissipatedPowers->push_back(state.dissipatedPower * powerUnit);
}

void MbD::ASMTLimit::initMarkers()
{
	if (motionJoint == "") {
		assert(markerI != "");
		assert(markerJ != "");
	}
	else {
		auto jt = root()->jointAt(motionJoint);
		markerI = jt->markerI;
		markerJ = jt->markerJ;
	}
}

void MbD::ASMTLimit::storeOnLevel(std::ofstream& os, size_t level)
{
	ASMTItemIJ::storeOnLevel(os, level);
	storeOnLevelString(os, level + 1, "MotionJoint");
	storeOnLevelString(os, level + 2, motionJoint);
	storeOnLevelString(os, level + 1, "Limit");
	storeOnLevelString(os, level + 2, limit);
	storeOnLevelString(os, level + 1, "Type");
	storeOnLevelString(os, level + 2, type);
	storeOnLevelString(os, level + 1, "Tol");
	storeOnLevelString(os, level + 2, tol);
	storeOnLevelString(os, level + 1, "Behavior");
	storeOnLevelString(os, level + 2, behavior);
	storeOnLevelString(os, level + 1, "Stiffness");
	storeOnLevelString(os, level + 2, stiffness);
	storeOnLevelString(os, level + 1, "Damping");
	storeOnLevelString(os, level + 2, damping);
}

void MbD::ASMTLimit::readMotionJoint(std::vector<std::string>& lines)
{
	if (readStringOffTop(lines) != "MotionJoint") {
		throw SimulationStoppingError("Expected MotionJoint record.");
	}
	motionJoint = readStringOffTop(lines);
}

void MbD::ASMTLimit::readLimit(std::vector<std::string>& lines)
{
	if (readStringOffTop(lines) != "Limit") {
		throw SimulationStoppingError("Expected Limit record.");
	}
	limit = readStringOffTop(lines);
}

void MbD::ASMTLimit::readType(std::vector<std::string>& lines)
{
	if (readStringOffTop(lines) != "Type") {
		throw SimulationStoppingError("Expected Type record.");
	}
	type = readStringOffTop(lines);
}

void MbD::ASMTLimit::readTol(std::vector<std::string>& lines)
{
	if (readStringOffTop(lines) != "Tol") {
		throw SimulationStoppingError("Expected Tol record.");
	}
	tol = readStringOffTop(lines);
}

void MbD::ASMTLimit::parseASMT(std::vector<std::string>& lines)
{
	ASMTConstraintSet::parseASMT(lines);
	readMotionJoint(lines);
	readLimit(lines);
	readType(lines);
	readTol(lines);
	if (!lines.empty()) {
		std::istringstream next(lines.front());
		std::string record;
		next >> record;
		if (record == "Behavior") readCompliance(lines);
	}
}

void MbD::ASMTLimit::readCompliance(std::vector<std::string>& lines)
{
	if (readStringOffTop(lines) != "Behavior") throw SimulationStoppingError("Expected Behavior record.");
	behavior = readStringOffTop(lines);
	if (readStringOffTop(lines) != "Stiffness") throw SimulationStoppingError("Expected Stiffness record.");
	stiffness = readStringOffTop(lines);
	if (readStringOffTop(lines) != "Damping") throw SimulationStoppingError("Expected Damping record.");
	damping = readStringOffTop(lines);
}

void MbD::ASMTLimit::createMbD(std::shared_ptr<System> mbdSys, std::shared_ptr<Units> mbdUnits)
{
	ASMTConstraintSet::createMbD(mbdSys, mbdUnits);
	auto limitIJ = std::static_pointer_cast<LimitIJ>(mbdObject);
	mbdSys->addLimit(limitIJ);
	//
	auto parser = std::make_shared<SymbolicParser>();
	parser->owner = this;
	std::shared_ptr<BasicUserFunction> userFunc;
	//
	userFunc = std::make_shared<BasicUserFunction>(limit, 1.0);
	parser->parseUserFunction(userFunc);
	auto& geolimit = parser->stack->top();
	geolimit = Symbolic::times(geolimit, sptrConstant(1.0 / coordinateUnit(*mbdUnits)));
	geolimit->createMbD(mbdSys, mbdUnits);
	geolimit = geolimit->simplified(geolimit);
	limitIJ->limit = geolimit->getValue();
	//
	limitIJ->type = type;
	//
	userFunc = std::make_shared<BasicUserFunction>(tol, 1.0);
	parser->parseUserFunction(userFunc);
	auto& geotol = parser->stack->top();
	geotol = Symbolic::times(geotol, sptrConstant(1.0 / coordinateUnit(*mbdUnits)));
	geotol->createMbD(mbdSys, mbdUnits);
	geotol = geotol->simplified(geotol);
	limitIJ->tol = geotol->getValue();
	if (behavior == "Compliant") {
		auto parseConstant = [&](const std::string& expression, double unit) {
			auto function = std::make_shared<BasicUserFunction>(expression, 1.0);
			parser->parseUserFunction(function);
			auto value = parser->stack->top();
			value = Symbolic::times(value, sptrConstant(unit));
			value->createMbD(mbdSys, mbdUnits);
			return value->simplified(value)->getValue();
		};
		limitIJ->setCompliant(
			parseConstant(stiffness, stiffnessUnit(*mbdUnits)),
			parseConstant(damping, dampingUnit(*mbdUnits))
		);
	}
	else if (behavior != "Rigid") {
		throw SimulationStoppingError("Unknown joint-limit behavior.");
	}
}

void MbD::ASMTLimit::setmotionJoint(const std::string& _motionJoint)
{
	motionJoint = _motionJoint;
}

void MbD::ASMTLimit::settype(const std::string& _type)
{
	type = _type;
}

void MbD::ASMTLimit::setlimit(const std::string& _limit)
{
	limit = _limit;
}

void MbD::ASMTLimit::settol(const std::string& _tol)
{
	tol = _tol;
}

void MbD::ASMTLimit::setcompliance(const std::string& newStiffness, const std::string& newDamping)
{
	behavior = "Compliant";
	stiffness = newStiffness;
	damping = newDamping;
}

double MbD::ASMTLimit::coordinateUnit(const Units& units) const { return units.angle; }
double MbD::ASMTLimit::stiffnessUnit(const Units& units) const { return units.angle / units.torque; }
double MbD::ASMTLimit::dampingUnit(const Units& units) const { return units.omega / units.torque; }
