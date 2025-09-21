#pragma once
#include <vector>

namespace OsgTest {
	class Mesh {
	public:
		std::vector<float> positions;
		std::vector<int> indices;
		std::vector<float> normals;
		std::vector<float> worldPositions;
		int nodeIndex;
	};
}
