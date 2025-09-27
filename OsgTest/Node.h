#pragma once
#include <string>
#include <TopoDS_Shape.hxx>
#include <TDF_Label.hxx>
#include <vector>

struct Node
{
	std::string name;
	TopoDS_Shape shape;
	TopLoc_Location location;
	TDF_Label label;
	std::vector<Node> children;
};