#pragma once

#include "tinygltf/tiny_gltf.h"
#include <TopoDS_Shape.hxx>

struct Edge {
    size_t v1, v2;
    Edge(size_t a, size_t b) {
        // ʼ�ձ���С����ǰ����֤�����Ψһ
        if (a < b) { v1 = a; v2 = b; }
        else { v1 = b; v2 = a; }
    }
    bool operator<(const Edge& other) const noexcept {
        // �ȱȽ� v1���ٱȽ� v2
        return std::tie(v1, v2) < std::tie(other.v1, other.v2);
    }
};

// ÿ�����¼����3������
struct Face {
    size_t v[3];
    gp_Vec normal;
};


class GltfManager {
public:
	bool ReadFile(const std::string& fileName);
	gp_Vec GetFaceNormal(const TopoDS_Shape& faceShape);
	void WriteFile();
};