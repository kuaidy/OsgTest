#pragma once
#include <TopoDS_Shape.hxx>
// ÿ�����¼����3������
struct Face {
	size_t v[3];
	gp_Vec normal;
};