#pragma once

#include "tinygltf/tiny_gltf.h"
#include <TopoDS_Shape.hxx>
#include "Mesh.h"
#include <osg/Matrixd>


class GltfManager {
public:
	bool ReadFile(const std::string& fileName);
	gp_Vec GetFaceNormal(const TopoDS_Shape& faceShape);
	void WriteFile();
	bool GetWorldPositions();
	void ComputeWorldMatrix(const tinygltf::Model& model, int nodeIndex, const osg::Matrixd& parentWorld, std::vector<osg::Matrixd>& worldMatrices);
	osg::Matrixd GetLocalMatrix(const tinygltf::Node& node);
	std::vector<OsgTest::Mesh> meshes;
	std::vector<osg::Matrixd> _worldMatrices;
private:
	tinygltf::Model m_model;
};