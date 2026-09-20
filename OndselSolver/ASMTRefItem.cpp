/***************************************************************************
 *   Copyright (c) 2023 Ondsel, Inc.                                       *
 *                                                                         *
 *   This file is part of OndselSolver.                                    *
 *                                                                         *
 *   See LICENSE file for details about copyright.                         *
 ***************************************************************************/
 
#include "ASMTRefItem.h"
#include "CREATE.h"

using namespace MbD;

void MbD::ASMTRefItem::addMarker(std::shared_ptr<ASMTMarker> marker)
{
	markers->push_back(marker);
	marker->owner = this;
}

void MbD::ASMTRefItem::readMarkers(std::vector<std::string>& lines)
{
	assert(lines[0].find("Markers") != std::string::npos);
	lines.erase(lines.begin());
	markers->clear();
	while (!lines.empty() && readString(lines.front()) == "Marker") {
		readMarker(lines);
	}
}

void MbD::ASMTRefItem::readMarker(std::vector<std::string>& lines)
{
	assert(lines[0].find("Marker") != std::string::npos);
	lines.erase(lines.begin());
	auto marker = CREATE<ASMTMarker>::With();
	marker->parseASMT(lines);
	markers->push_back(marker);
	marker->owner = this;
}

void MbD::ASMTRefItem::storeOnLevel(std::ofstream& os, size_t level)
{
	storeOnLevelString(os, level, "RefPoints");
	ASMTSpatialItem::storeOnLevel(os, level+1);
	for (auto& marker : *markers) {
		marker->storeOnLevel(os, level);
	}
}
