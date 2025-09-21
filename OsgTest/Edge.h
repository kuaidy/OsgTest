#pragma once

#include <tuple>

struct Edge {
	size_t v1, v2;
	Edge(size_t a, size_t b) {
		// 始终保持小的在前，保证无向边唯一
		if (a < b) { v1 = a; v2 = b; }
		else { v1 = b; v2 = a; }
	}
	bool operator<(const Edge& other) const noexcept {
		// 先比较 v1，再比较 v2
		return std::tie(v1, v2) < std::tie(other.v1, other.v2);
	}
};