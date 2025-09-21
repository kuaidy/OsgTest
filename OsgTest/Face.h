#pragma once
#include <TopoDS_Shape.hxx>
// 每个面记录它的3个顶点
struct Face {
	size_t v[3];
	gp_Vec normal;
};