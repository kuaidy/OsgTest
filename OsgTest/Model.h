#pragma once
#include <string>
#include <TopoDS_Shape.hxx>
#include <TopLoc_Location.hxx>
#include <TDF_Label.hxx>

struct Model
{
	std::string name;
	TopoDS_Shape shape;
	TopLoc_Location location;
	TDF_Label label;
};