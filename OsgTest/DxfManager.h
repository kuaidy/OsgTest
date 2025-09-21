#pragma once
#include <string>
#include <vector>
#include "Mesh.h"
class DxfManager {
public:
	bool Write(const std::string& path,std::vector<OsgTest::Mesh> meshes);
};