/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/

#pragma once

#include "ASMTConstraintSet.h"
#include "ForceTorqueData.h"

namespace MbD {
	class Units;
    class ASMTLimit : public ASMTConstraintSet
    {
        //
    public:
		void initialize() override;
		void clearResults();
		void updateFromMbD() override;
        virtual void initMarkers();
        void storeOnLevel(std::ofstream& os, size_t level) override;
        void readMotionJoint(std::vector<std::string>& lines);
        void readLimit(std::vector<std::string>& lines);
        void readType(std::vector<std::string>& lines);
        void readTol(std::vector<std::string>& lines);
		void readCompliance(std::vector<std::string>& lines);
        void parseASMT(std::vector<std::string>& lines) override;
        void createMbD(std::shared_ptr<System> mbdSys, std::shared_ptr<Units> mbdUnits) override;
        void setmotionJoint(const std::string& _motionJoint);
        void settype(const std::string& _type);
        void setlimit(const std::string& _limit);
        void settol(const std::string& _tol);
		void setcompliance(const std::string& _stiffness, const std::string& _damping);
		virtual double coordinateUnit(const Units& units) const;
		virtual double stiffnessUnit(const Units& units) const;
		virtual double dampingUnit(const Units& units) const;

		std::string motionJoint, type, limit, tol;
		std::string behavior = "Rigid", stiffness = "0.0", damping = "0.0";
		FRowDsptr storedEnergies, dissipatedPowers;

    };
}

