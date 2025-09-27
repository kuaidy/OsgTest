#pragma once
#include <string>
#include <XCAFDoc_ShapeTool.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <vector>
#include "Node.h"

class StepManager {
public :
	void Read(std::string fileName);
	void TraverseXDEComponent(const TDF_Label& label, 
		const Handle(XCAFDoc_ShapeTool)& shapeTool, 
		const Handle(XCAFDoc_ColorTool)& colorTool, 
		int depth,
		std::vector<Node>& nodes,
		Node* parentNode);
	std::vector<Node> nodes;
};