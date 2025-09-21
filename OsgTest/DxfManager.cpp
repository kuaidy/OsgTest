#include "DxfManager.h"
#include "Edge.h"
#include <map>
#include <unordered_map>
#include "Face.h"
#include <unordered_set>

#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <gp_Pln.hxx>
#include <TopoDS_Face.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <HLRBRep_PolyAlgo.hxx>
#include <HLRBRep_PolyHLRToShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <BRep_Tool.hxx>
#include <HLRBRep_Algo.hxx>
#include <HLRBRep_HLRToShape.hxx>
#include <GeomLProp_SLProps.hxx>

// 顶点 Key
struct PntKey
{
	gp_Pnt p;
	bool operator==(const PntKey& other) const noexcept
	{
		const double eps = 1e-3;
		return std::fabs(p.X() - other.p.X()) < eps &&
			std::fabs(p.Y() - other.p.Y()) < eps &&
			std::fabs(p.Z() - other.p.Z()) < eps;
	}
};


// 自定义哈希
namespace std {
	template<>
	struct hash<PntKey>
	{
		size_t operator()(const PntKey& k) const noexcept
		{
			auto h1 = std::hash<long long>()(llround(k.p.X() * 1e3));
			auto h2 = std::hash<long long>()(llround(k.p.Y() * 1e3));
			auto h3 = std::hash<long long>()(llround(k.p.Z() * 1e3));
			return h1 ^ (h2 << 1) ^ (h3 << 2);
		}
	};
}

struct EdgeKey {
	gp_Pnt p0;
	gp_Pnt p1;
	double eps = 1e-4;

	EdgeKey(const gp_Pnt& a, const gp_Pnt& b) {
		// 保证 p0 < p1，用坐标顺序固定顺序，忽略方向
		if (a.X() < b.X() || (a.X() == b.X() && a.Y() < b.Y()) ||
			(a.X() == b.X() && a.Y() == b.Y() && a.Z() < b.Z())) {
			p0 = a; p1 = b;
		}
		else {
			p0 = b; p1 = a;
		}
	}

	bool operator==(const EdgeKey& other) const noexcept {
		auto eq = [this](const gp_Pnt& a, const gp_Pnt& b) {
			return a.Distance(b) < eps;
			};
		return eq(p0, other.p0) && eq(p1, other.p1);
	}
};

namespace std {
	template <>
	struct hash<EdgeKey> {
		size_t operator()(const EdgeKey& k) const noexcept {
			auto h = [](double x) { return std::hash<long long>()(llround(x * 1e6)); };
			size_t h0 = h(k.p0.X()) ^ (h(k.p0.Y()) << 1) ^ (h(k.p0.Z()) << 2);
			size_t h1 = h(k.p1.X()) ^ (h(k.p1.Y()) << 1) ^ (h(k.p1.Z()) << 2);
			return h0 ^ (h1 << 1);
		}
	};
}

class EdgeSet {
	std::unordered_set<EdgeKey> edges;

public:
	bool Contains(const gp_Pnt& a, const gp_Pnt& b) const {
		return edges.find(EdgeKey(a, b)) != edges.end();
	}

	void Add(const gp_Pnt& a, const gp_Pnt& b) {
		edges.insert(EdgeKey(a, b));
	}

	size_t Size() const { return edges.size(); }
};


bool DxfManager::Write(const std::string& path, std::vector<OsgTest::Mesh> meshes) {
	std::ofstream dxf("test.dxf");
	dxf << "0\nSECTION\n2\nENTITIES\n";

	std::ofstream outlinedxf("outline.dxf");
	outlinedxf << "0\nSECTION\n2\nENTITIES\n";

	for (const auto& mesh : meshes) {

		double tolerance = 1.0e-3;
		// sewing 用于把多个 face 合并成一个 shape
		BRepBuilderAPI_Sewing sewing(tolerance, Standard_True, Standard_True, Standard_True, Standard_True);

		//构建边和面的对应
		std::map<Edge, std::vector<int>> edgeToFaces;
		size_t faceId = 0;
		std::vector<Face> facesOut;
		std::vector<int> indices = mesh.indices;
		std::vector<float> positions = mesh.positions;
		std::unordered_map<PntKey, int> tmpIndices;
		std::vector<PntKey> uniquePnts;

		if (indices.size() == 0) {
			for (size_t i = 0; i + 2 < positions.size(); i += 3) {
				float x = positions[i];
				float y = positions[i + 1];
				float z = positions[i + 2];
				gp_Pnt p(x, y, z);
				PntKey pntKey{ p };
				auto it = tmpIndices.find(pntKey);
				if (it != tmpIndices.end())
				{
					// 已存在，直接复用索引
					indices.push_back(it->second);
				}
				else
				{
					uint32_t newIndex = static_cast<uint32_t>(uniquePnts.size());
					uniquePnts.push_back(pntKey);
					tmpIndices[pntKey] = newIndex;
					indices.push_back(newIndex);
				}
			}

			if (indices.size() % 3 != 0) continue;

			for (size_t i = 0; i + 2 < indices.size(); i += 3) {
				try {
					size_t i0 = indices[i];
					size_t i1 = indices[i + 1];
					size_t i2 = indices[i + 2];

					Face f{ {i0, i1, i2} };
					gp_Pnt p1(uniquePnts[i0].p);
					gp_Pnt p2(uniquePnts[i1].p);
					gp_Pnt p3(uniquePnts[i2].p);
					//计算面的法向量
					gp_Vec v1(p1, p2);
					gp_Vec v2(p1, p3);
					gp_Vec n = v1.Crossed(v2);
					if (n.Magnitude() < 1e-10) continue; // 面退化跳过
					f.normal = n;
					facesOut.push_back(f);
					// 三条边
					Edge e1(i0, i1);
					Edge e2(i1, i2);
					Edge e3(i2, i0);
					edgeToFaces[e1].push_back(faceId);
					edgeToFaces[e2].push_back(faceId);
					edgeToFaces[e3].push_back(faceId);
					++faceId;

					//创建三角形 wire
					BRepBuilderAPI_MakePolygon poly;
					poly.Add(p1);
					poly.Add(p2);
					poly.Add(p3);
					poly.Close();
					TopoDS_Wire wire = poly.Wire();

					TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);
					if (face.IsNull()) continue; // 确保成功
					sewing.Add(face);

					//// 计算三角形法线
					//gp_Vec v11(p2.XYZ() - p1.XYZ());
					//gp_Vec v22(p3.XYZ() - p1.XYZ());
					//gp_Dir normal(v11.Crossed(v22));
					//// 创建平面
					//gp_Pln plane(p1, normal);
					//// 创建 face
					//TopoDS_Face face = BRepBuilderAPI_MakeFace(plane, wire);
					//if (face.Orientation() == TopAbs_REVERSED) {
					//	face.Reverse(); // 确保所有面法线方向一致
					//}
					//sewing.Add(face);
				}
				catch (const Standard_Failure& e)
				{
					size_t i0 = indices[i];
					size_t i1 = indices[i + 1];
					size_t i2 = indices[i + 2];

					gp_Pnt p1(uniquePnts[i0].p);
					gp_Pnt p2(uniquePnts[i1].p);
					gp_Pnt p3(uniquePnts[i2].p);

					std::cout << p1.X() << p1.Y() << p1.Z() << std::endl;
					std::cout << p2.X() << p2.Y() << p2.Z() << std::endl;
					std::cout << p3.X() << p3.Y() << p3.Z() << std::endl;

					std::cerr << "OCCT exception: " << e.GetMessageString() << std::endl;
					continue;
				}
			}

			std::cout << "总共多少个面:" << faceId << std::endl;
			//计算轮廓边
			gp_Vec viewDir(1, 1, 0);
			gp_Pnt origin(0, 0, 0);
			gp_Vec uDir(0, 0, 1);
			auto projectToPlane = [viewDir, origin, uDir](const gp_Pnt& p) {
				gp_Vec normal = viewDir;
				gp_Vec v(origin, p);
				double t = normal.Dot(v) / normal.SquareMagnitude();
				gp_Pnt proj = p.Translated(-t * normal);
				gp_Dir vDir = normal.Crossed(uDir);
				double u = (proj.XYZ() - origin.XYZ()).Dot(uDir.XYZ());
				double vCoord = (proj.XYZ() - origin.XYZ()).Dot(vDir.XYZ());
				return std::pair<double, double>{u, vCoord};
				};

			std::vector<Edge> silhouetteEdges;
			for (const auto& kv : edgeToFaces) {
				const Edge& edge = kv.first;
				const std::vector<int>& faces = kv.second;

				if (faces.size() == 1) {
					// 自由边 → 轮廓边
					silhouetteEdges.push_back(edge);
				}
				else if (faces.size() == 2) {
					const gp_Vec& n1 = facesOut[faces[0]].normal;
					const gp_Vec& n2 = facesOut[faces[1]].normal;
					double dot1 = n1.Dot(viewDir);
					double dot2 = n2.Dot(viewDir);
					if (dot1 * dot2 < 0) {
						silhouetteEdges.push_back(edge); // 剪影边
					}
					else if ((std::abs(dot1) == 1 && std::abs(dot2) == 0) || (std::abs(dot1) == 0 && std::abs(dot2) == 1)) {
						silhouetteEdges.push_back(edge); // 剪影边
					}
				}
			}
			std::cout << "轮廓边" << silhouetteEdges.size() << std::endl;

			for (const Edge& e : silhouetteEdges)
			{
				gp_Pnt p1(uniquePnts[e.v1].p);
				gp_Pnt p2(uniquePnts[e.v2].p);

				auto [x1, y1] = projectToPlane(p1);
				auto [x2, y2] = projectToPlane(p2);

				outlinedxf << "0\nLINE\n8\n0\n";
				outlinedxf << "10\n" << x1 << "\n20\n" << y1 << "\n30\n0.0\n";
				outlinedxf << "11\n" << x2 << "\n21\n" << y2 << "\n31\n0.0\n";
			}

			//合并成一个 shape
			sewing.Perform();
			TopoDS_Shape shape = sewing.SewedShape();
			//// 三角化处理
			//BRepMesh_IncrementalMesh mesh(shape, 1.0e-3);
			//mesh.Perform();
			// 加入 XCAF
				// 1️⃣ 创建 XCAF 文档
			Handle(TDocStd_Document) doc;
			XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
			Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
			if (shapeTool.IsNull()) {
				std::cerr << "ShapeTool not found." << std::endl;
				return 1;
			}
			shapeTool->AddShape(shape);

			TDF_LabelSequence freeShapes;
			shapeTool->GetFreeShapes(freeShapes);
			std::cout << "Total shapes added: " << freeShapes.Length() << std::endl;


			for (Standard_Integer i = 1; i <= freeShapes.Length(); ++i) {
				TDF_Label lab = freeShapes.Value(i);
				TopoDS_Shape shape = shapeTool->GetShape(lab);

				Handle(HLRBRep_Algo) hlr = new HLRBRep_Algo();
				hlr->Add(shape);
				//定义投影视角
				gp_Ax2 view(gp::Origin(), gp::DZ());
				hlr->Projector(HLRAlgo_Projector(view));
				hlr->Update();
				hlr->Hide();
				HLRBRep_HLRToShape hlrToShape(hlr);
				auto visibleEdges = hlrToShape.VCompound(); // 可见边

				//轮廓边
				//TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
				//TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);
				//std::vector<TopoDS_Edge> outlineEdges;
				////求轮廓边
				//EdgeSet silhouetteEdges;

				//for (int i = 1; i <= edgeToFaces.Extent(); ++i) {
				//	TopoDS_Edge edge = TopoDS::Edge(edgeToFaces.FindKey(i));
				//	const TopTools_ListOfShape& faces = edgeToFaces.FindFromIndex(i);

				//	if (faces.Extent() == 1) {
				//		outlineEdges.push_back(edge);
				//		// 边界边
				//		TopoDS_Vertex v1, v2;
				//		TopExp::Vertices(edge, v1, v2);          // 获取边的两个顶点
				//		gp_Pnt p1 = BRep_Tool::Pnt(v1);          // 顶点1坐标
				//		gp_Pnt p2 = BRep_Tool::Pnt(v2);          // 顶点2坐标

				//		silhouetteEdges.Add(p1, p2);             // 加入集合
				//	}
				//	else if (faces.Extent() == 2) {
				//		TopoDS_Face f1 = TopoDS::Face(faces.First());
				//		TopoDS_Face f2 = TopoDS::Face(faces.Last());

				//		gp_Vec n1, n2;
				//		// 获取面法向量（用面中心近似）
				//		Handle(Geom_Surface) surf = BRep_Tool::Surface(f1);
				//		GeomLProp_SLProps props(surf, 0.5, 0.5, 1, 1e-6); // u,v 中心点
				//		n1 = props.Normal();
				//		n2 = props.Normal();

				//		double dp1 = n1.Dot(gp_Vec(0, 0, 1));
				//		double dp2 = n2.Dot(gp_Vec(0, 0, 1));

				//		if (dp1 * dp2 < 0) {
				//			outlineEdges.push_back(edge);
				//			// 两面法向量对视方向符号不同 → 轮廓边
				//			TopoDS_Vertex v1, v2;
				//			TopExp::Vertices(edge, v1, v2);          // 获取边的两个顶点
				//			gp_Pnt p1 = BRep_Tool::Pnt(v1);          // 顶点1坐标
				//			gp_Pnt p2 = BRep_Tool::Pnt(v2);          // 顶点2坐标
				//			silhouetteEdges.Add(p1, p2);             // 加入集合
				//		}
				//	}
				//}

				////输出轮廓边
				//for (auto edge : outlineEdges) {
				//	// 获取顶点
				//	Standard_Real x1, y1, z1, x2, y2, z2;
				//	TopoDS_Vertex v1, v2;
				//	TopExp::Vertices(edge, v1, v2);
				//	gp_Pnt p1 = BRep_Tool::Pnt(v1);
				//	gp_Pnt p2 = BRep_Tool::Pnt(v2);
				//	if (!silhouetteEdges.Contains(p1, p2) && !silhouetteEdges.Contains(p2, p1)) continue;

				//	outlinedxf << "0\nLINE\n8\n0\n";
				//	outlinedxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
				//	outlinedxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
				//}

				//std::cout << silhouetteEdges.size() << std::endl;

		/*		Handle(HLRBRep_Algo) hlrAlgo = new HLRBRep_Algo();
				gp_Ax2 view(gp::Origin(), gp::DY());
				hlrAlgo->Projector(HLRAlgo_Projector(view));
				hlrAlgo->Add(shape);
				hlrAlgo->Update();

				HLRBRep_PolyHLRToShape extractor;
				extractor.Update(hlrAlgo);
				TopoDS_Shape visibleEdges = extractor.VCompound();*/

				// 遍历 TopoDS_Edge 输出线段
				for (TopExp_Explorer exp(visibleEdges, TopAbs_EDGE); exp.More(); exp.Next()) {
					TopoDS_Edge edge = TopoDS::Edge(exp.Current());
					// 获取顶点
					Standard_Real x1, y1, z1, x2, y2, z2;
					TopoDS_Vertex v1, v2;
					TopExp::Vertices(edge, v1, v2);
					gp_Pnt p1 = BRep_Tool::Pnt(v1);
					gp_Pnt p2 = BRep_Tool::Pnt(v2);
					//if (!silhouetteEdges.Contains(p1, p2) && !silhouetteEdges.Contains(p2, p1)) continue;

					dxf << "0\nLINE\n8\n0\n";
					dxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
					dxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
				}
			}

			//std::cout << "总共多少个面:" << faceId << std::endl;
			////计算轮廓边
			//gp_Vec viewDir(0, 1, 0);
			//auto projectToPlane = [](const gp_Pnt& p) { return std::pair<double, double>{p.X(), p.Z()}; };
			//std::vector<Edge> silhouetteEdges;
			//for (const auto& kv : edgeToFaces) {
			//	const Edge& edge = kv.first;
			//	const std::vector<int>& faces = kv.second;

			//	if (faces.size() == 1) {
			//		// 自由边 → 轮廓边
			//		silhouetteEdges.push_back(edge);
			//	}
			//	else if (faces.size() == 2) {
			//		const gp_Vec& n1 = facesOut[faces[0]].normal;
			//		const gp_Vec& n2 = facesOut[faces[1]].normal;
			//		double dot1 = n1.Dot(viewDir);
			//		double dot2 = n2.Dot(viewDir);
			//		if (dot1 * dot2 < 0) {
			//			silhouetteEdges.push_back(edge); // 剪影边
			//		}
			//		else if ((std::abs(dot1) == 1 && std::abs(dot2) == 0) || (std::abs(dot1) == 0 && std::abs(dot2) == 1)) {
			//			silhouetteEdges.push_back(edge); // 剪影边
			//		}
			//	}
			//}

			//for (const Edge& e : silhouetteEdges)
			//{
			//	gp_Pnt p1(positions[3 * e.v1 + 0], positions[3 * e.v1 + 1], positions[3 * e.v1 + 2]);
			//	gp_Pnt p2(positions[3 * e.v2 + 0], positions[3 * e.v2 + 1], positions[3 * e.v2 + 2]);

			//	auto [x1, y1] = projectToPlane(p1);
			//	auto [x2, y2] = projectToPlane(p2);

			//	dxf << "0\nLINE\n8\n0\n";
			//	dxf << "10\n" << x1 << "\n20\n" << y1 << "\n30\n0.0\n";
			//	dxf << "11\n" << x2 << "\n21\n" << y2 << "\n31\n0.0\n";
			//}
		}
		//if (indices.size() % 3 != 0) continue;

		//for (size_t i = 0; i + 2 < indices.size(); i += 3) {
		//	size_t i0 = indices[i];
		//	size_t i1 = indices[i + 1];
		//	size_t i2 = indices[i + 2];

		//	Face f{ {i0, i1, i2} };
		//	gp_Pnt p1(positions[3 * indices[i] + 0], positions[3 * indices[i] + 1], positions[3 * indices[i] + 2]);
		//	gp_Pnt p2(positions[3 * indices[i + 1] + 0], positions[3 * indices[i + 1] + 1], positions[3 * indices[i + 1] + 2]);
		//	gp_Pnt p3(positions[3 * indices[i + 2] + 0], positions[3 * indices[i + 2] + 1], positions[3 * indices[i + 2] + 2]);
		//	//计算面的法向量
		//	//gp_Vec v1(p1, p2);
		//	//gp_Vec v2(p1, p3);
		//	//gp_Vec n = v1.Crossed(v2);
		//	//if (n.Magnitude() > 1e-10)
		//	//	n.Normalize();
		//	//f.normal = n;
		//	//facesOut.push_back(f);

		//	//// 三条边
		//	//Edge e1(i0, i1);
		//	//Edge e2(i1, i2);
		//	//Edge e3(i2, i0);

		//	//edgeToFaces[e1].push_back(faceId);
		//	//edgeToFaces[e2].push_back(faceId);
		//	//edgeToFaces[e3].push_back(faceId);

		//	//++faceId;

		//	//创建三角形 wire
		//	BRepBuilderAPI_MakePolygon poly;
		//	poly.Add(p1);
		//	poly.Add(p2);
		//	poly.Add(p3);
		//	poly.Close();
		//	TopoDS_Wire wire = poly.Wire();

		//	// 计算三角形法线
		//	gp_Vec v1(p2.XYZ() - p1.XYZ());
		//	gp_Vec v2(p3.XYZ() - p1.XYZ());
		//	gp_Dir normal(v1.Crossed(v2));
		//	// 创建平面
		//	gp_Pln plane(p1, normal);
		//	// 创建 face
		//	TopoDS_Face face = BRepBuilderAPI_MakeFace(plane, wire);
		//	//if (face.Orientation() == TopAbs_REVERSED) {
		//	//	face.Reverse(); // 确保所有面法线方向一致
		//	//}
		//	sewing.Add(face);
		//}
		////合并成一个 shape
		//sewing.Perform();
		//TopoDS_Shape shape = sewing.SewedShape();
		////// 三角化处理
		////BRepMesh_IncrementalMesh mesh(shape, 1.0e-3);
		////mesh.Perform();
		//// 加入 XCAF
		//	// 1️⃣ 创建 XCAF 文档
		//Handle(TDocStd_Document) doc;
		//XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
		//Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
		//if (shapeTool.IsNull()) {
		//	std::cerr << "ShapeTool not found." << std::endl;
		//	return 1;
		//}
		//shapeTool->AddShape(shape);

		//Handle(HLRBRep_PolyAlgo) hlrAlgo = new HLRBRep_PolyAlgo();
		//gp_Ax2 view(gp::Origin(), gp::DY());
		//hlrAlgo->Projector(HLRAlgo_Projector(view));
		//hlrAlgo->Load(shape);
		//hlrAlgo->Update();

		//HLRBRep_PolyHLRToShape extractor;
		//extractor.Update(hlrAlgo);
		//TopoDS_Shape visibleEdges = extractor.VCompound();

		//// 遍历 TopoDS_Edge 输出线段
		//for (TopExp_Explorer exp(visibleEdges, TopAbs_EDGE); exp.More(); exp.Next()) {
		//	TopoDS_Edge edge = TopoDS::Edge(exp.Current());

		//	// 获取顶点
		//	Standard_Real x1, y1, z1, x2, y2, z2;
		//	TopoDS_Vertex v1, v2;
		//	TopExp::Vertices(edge, v1, v2);
		//	gp_Pnt p1 = BRep_Tool::Pnt(v1);
		//	gp_Pnt p2 = BRep_Tool::Pnt(v2);

		//	dxf << "0\nLINE\n8\n0\n";
		//	dxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
		//	dxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
		//}

		//std::cout << "总共多少个面:" << faceId << std::endl;
		////计算轮廓边
		//gp_Vec viewDir(0, 1, 0);
		//auto projectToPlane = [](const gp_Pnt& p) { return std::pair<double, double>{p.X(), p.Z()}; };
		//std::vector<Edge> silhouetteEdges;
		//for (const auto& kv : edgeToFaces) {
		//	const Edge& edge = kv.first;
		//	const std::vector<int>& faces = kv.second;

		//	if (faces.size() == 1) {
		//		// 自由边 → 轮廓边
		//		silhouetteEdges.push_back(edge);
		//	}
		//	else if (faces.size() == 2) {
		//		const gp_Vec& n1 = facesOut[faces[0]].normal;
		//		const gp_Vec& n2 = facesOut[faces[1]].normal;
		//		double dot1 = n1.Dot(viewDir);
		//		double dot2 = n2.Dot(viewDir);
		//		if (dot1 * dot2 < 0) {
		//			silhouetteEdges.push_back(edge); // 剪影边
		//		}
		//		else if ((std::abs(dot1) == 1 && std::abs(dot2) == 0) || (std::abs(dot1) == 0 && std::abs(dot2) == 1)) {
		//			silhouetteEdges.push_back(edge); // 剪影边
		//		}
		//	}
		//}

		//for (const Edge& e : silhouetteEdges)
		//{
		//	gp_Pnt p1(positions[3 * e.v1 + 0], positions[3 * e.v1 + 1], positions[3 * e.v1 + 2]);
		//	gp_Pnt p2(positions[3 * e.v2 + 0], positions[3 * e.v2 + 1], positions[3 * e.v2 + 2]);

		//	auto [x1, y1] = projectToPlane(p1);
		//	auto [x2, y2] = projectToPlane(p2);

		//	dxf << "0\nLINE\n8\n0\n";
		//	dxf << "10\n" << x1 << "\n20\n" << y1 << "\n30\n0.0\n";
		//	dxf << "11\n" << x2 << "\n21\n" << y2 << "\n31\n0.0\n";
		//}
	}

	outlinedxf << "0\nENDSEC\n0\nEOF\n";
	outlinedxf.close();

	dxf << "0\nENDSEC\n0\nEOF\n";
	dxf.close();

	return true;
}