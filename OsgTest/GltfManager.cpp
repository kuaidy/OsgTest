#include "GltfManager.h"
#include <gp_Pnt.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <RWGltf_CafReader.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <TDF_LabelSequence.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <HLRBRep_Algo.hxx>
#include <HLRBRep_HLRToShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <Geom_Curve.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Wire.hxx>
#include <TopExp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepGProp_Face.hxx>
#include <GeomLProp_SLProps.hxx>
#include <BRepTools.hxx>
#include <HLRTopoBRep_OutLiner.hxx>
#include <HLRBRep_PolyAlgo.hxx>
#include <HLRBRep_PolyHLRToShape.hxx>
#include <BRep_Builder.hxx>

#include <osg/Matrixd>

//bool GltfManager::ReadFile(const std::string fileName) {
//	//得到文件夹
//	std::filesystem::path path(fileName);
//	std::filesystem::path folderPath = path.parent_path();
//	tinygltf::Model model;
//	tinygltf::TinyGLTF loader;
//	std::string err;
//	std::string warn;
//	bool res = false;
//	if (path.extension() == ".glb") {
//		res = loader.LoadBinaryFromFile(&model, &err, &warn, fileName);
//	}
//	else {
//		res = loader.LoadASCIIFromFile(&model, &err, &warn, fileName);
//	}
//	if (!warn.empty()) {
//		printf("Warn: %s\n", warn.c_str());
//	}
//	if (!err.empty()) {
//		printf("Err: %s\n", err.c_str());
//		return false;
//	}
//	if (!res) {
//		printf("Failed to parse glTF\n");
//		return false;
//	}
//	//获取顶点数据
//	std::vector<float> positions;
//	for (auto &mesh:model.meshes) {
//		for (auto &primitive:mesh.primitives) {
//			auto it = primitive.attributes.find("POSITION");
//			if (it == primitive.attributes.end()) continue;
//			const tinygltf::Accessor& accessor = model.accessors[it->second];
//			const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
//			const tinygltf::Buffer& buffer = model.buffers[view.buffer];
//
//			const float* pos = reinterpret_cast<const float*>(
//				&buffer.data[accessor.byteOffset + view.byteOffset]);
//
//			size_t count = accessor.count; // 顶点个数
//			positions.insert(positions.end(), pos, pos + count * 3);
//		}
//	}
//	//转换成occ对象
//	std::vector<gp_Pnt> vertices;
//	for(TopExp_Exp)
//}
//bool GltfManager::ReadFile(const std::string& fileName) {
//	Handle(TDocStd_Document) doc;
//	Handle(XCAFApp_Application) app= XCAFApp_Application::GetApplication();
//	app->NewDocument("MDTV-XCAF", doc);
//
//	RWGltf_CafReader reader;
//	reader.SetDocument(doc);
//	reader.SetLoadAllScenes(true);     // 如果想加载所有 scenes
//	reader.SetDoublePrecision(false); // 单/双精度顶点（按需）
//	
//	Message_ProgressRange progress;
//	bool ok = reader.Perform(fileName.c_str(), progress);
//	if (!ok) {
//		std::cerr << "read failed" << std::endl;
//		return false;
//	}
//	// 4) 从文档中取出 Shape（XDE/XCAF）
//	Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
//	TDF_LabelSequence seq;
//	shapeTool->GetFreeShapes(seq);
//	for (Standard_Integer i = 1; i <= seq.Length(); ++i) {
//		TDF_Label lab = seq.Value(i);
//		TopoDS_Shape shape = shapeTool->GetShape(lab);
//		//做二维投影
//		Handle(HLRBRep_Algo) hlr = new HLRBRep_Algo();
//		hlr->Add(shape);
//		//定义投影视角
//		gp_Ax2 view(gp::Origin(),gp::DZ());
//		hlr->Projector(HLRAlgo_Projector(view));
//		hlr->Update();
//		hlr->Hide();
//		HLRBRep_HLRToShape hlrToShape(hlr);
//
//		auto visibleEdges = hlrToShape.VCompound(); // 可见边
//		auto hiddenEdges = hlrToShape.HCompound(); // 隐藏边
//
//		//写dxf
//		std::ofstream dxf("test.dxf");
//		dxf << "0\nSECTION\n2\nENTITIES\n";
//
//		for (TopExp_Explorer exp(visibleEdges, TopAbs_EDGE); exp.More();exp.Next()) {
//			TopoDS_Edge edge = TopoDS::Edge(exp.Current());
//			Standard_Real f, l;
//			Handle(Geom_Curve) curve= BRep_Tool::Curve(edge, f, l);
//			if (curve.IsNull()) continue;
//			gp_Pnt p1, p2;
//			curve->D0(f,p1);
//			curve->D0(1,p2);
//			dxf << "0\nLINE\n8\n0\n";
//			dxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
//			dxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
//		}
//
//		dxf << "0\nENDSEC\n0\nEOF\n";
//	}
//}

bool GltfManager::ReadFile(const std::string& fileName) {
	// 2️⃣ 用 tinygltf 解析 GLB
	tinygltf::Model model;
	tinygltf::TinyGLTF loader;
	std::string err, warn;
	bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, fileName);
	if (!ret) {
		ret = loader.LoadBinaryFromFile(&model, &err, &warn, fileName);
		if (!ret) {
			std::cerr << "Failed to load GLB: " << err << std::endl;
			return false;
		}
	}
	m_model = model;

	_worldMatrices.resize(m_model.nodes.size());
	//const tinygltf::Scene& scene = m_model.scenes[m_model.defaultScene];
	//for (int rootNode : scene.nodes) {
	//	osg::Matrixd identity;  // 单位矩阵
	//	identity.makeIdentity();  // 显式设置为单位矩阵
	//	ComputeWorldMatrix(m_model, rootNode, identity, _worldMatrices);
	//}
	for (int i = 0; i < m_model.nodes.size(); ++i) {
		osg::Matrixd identity;  // 单位矩阵
		ComputeWorldMatrix(m_model, i, identity, _worldMatrices);
	}

	GetWorldPositions();

	//for (const auto& mesh : model.meshes) {
	//	OsgTest::Mesh tmpMesh;
	//	for (const auto& prim : mesh.primitives) {
	//		// 获取 POSITION
	//		std::vector<float> positions;
	//		auto accessorIndex = prim.attributes.find("POSITION");
	//		if (accessorIndex != prim.attributes.end()) {
	//			const tinygltf::Accessor& accessor = model.accessors[accessorIndex->second];
	//			const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
	//			const tinygltf::Buffer& buffer = model.buffers[view.buffer];
	//			const unsigned char* dataPtr = buffer.data.data() + view.byteOffset + accessor.byteOffset;
	//			size_t count = accessor.count;
	//			size_t stride = accessor.ByteStride(view);
	//			if (stride == 0) {
	//				stride = tinygltf::GetNumComponentsInType(accessor.type) * tinygltf::GetComponentSizeInBytes(accessor.componentType);
	//			}
	//			for (size_t i = 0; i < count; ++i) {
	//				const float* src = reinterpret_cast<const float*>(dataPtr + i * stride);
	//				for (int j = 0; j < tinygltf::GetNumComponentsInType(accessor.type); ++j) {
	//					positions.push_back(src[j]);
	//				}
	//			}
	//		}
	//		tmpMesh.positions = positions;
	//		//获取索引
	//		std::vector<int> indices;
	//		if (prim.indices >= 0) {
	//			const tinygltf::Accessor& accessor = model.accessors[prim.indices];
	//			const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
	//			const tinygltf::Buffer& buffer = model.buffers[view.buffer];
	//			const unsigned char* dataPtr = buffer.data.data() + view.byteOffset + accessor.byteOffset;
	//			for (size_t k = 0; k < accessor.count; ++k) {
	//				uint32_t index = 0;
	//				switch (accessor.componentType) {
	//				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
	//					index = *(reinterpret_cast<const uint8_t*>(dataPtr + k));
	//					break;
	//				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
	//					index = *(reinterpret_cast<const uint16_t*>(dataPtr + k * sizeof(uint16_t)));
	//					break;
	//				case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
	//					index = *(reinterpret_cast<const uint32_t*>(dataPtr + k * sizeof(uint32_t)));
	//					break;
	//				}
	//				indices.push_back(index);
	//			}
	//		}
	//		tmpMesh.indices = indices;
	//	}
	//	meshes.push_back(tmpMesh);
	//}
}


bool GltfManager::GetWorldPositions() {
	for (int i = 0; i < m_model.nodes.size(); ++i) {
		osg::Matrixd worldMatrixd = _worldMatrices[i];
		if (m_model.nodes[i].mesh >= 0) {
			tinygltf::Mesh mesh = m_model.meshes[m_model.nodes[i].mesh];
			OsgTest::Mesh tmpMesh;
			tmpMesh.nodeIndex = i;
			for (const auto& prim : mesh.primitives) {
				// 获取 POSITION
				std::vector<float> positions;
				auto accessorIndex = prim.attributes.find("POSITION");
				if (accessorIndex != prim.attributes.end()) {
					const tinygltf::Accessor& accessor = m_model.accessors[accessorIndex->second];
					const tinygltf::BufferView& view = m_model.bufferViews[accessor.bufferView];
					const tinygltf::Buffer& buffer = m_model.buffers[view.buffer];
					const unsigned char* dataPtr = buffer.data.data() + view.byteOffset + accessor.byteOffset;
					size_t count = accessor.count;
					size_t stride = accessor.ByteStride(view);
					if (stride == 0) {
						stride = tinygltf::GetNumComponentsInType(accessor.type) * tinygltf::GetComponentSizeInBytes(accessor.componentType);
					}
					for (size_t i = 0; i < count; ++i) {
						const float* src = reinterpret_cast<const float*>(dataPtr + i * stride);
						for (int j = 0; j < tinygltf::GetNumComponentsInType(accessor.type); ++j) {
							positions.push_back(src[j]);
						}
					}
				}
				tmpMesh.positions = positions;
				//得到世界坐标
				std::vector<float> worldPositions;
				for (size_t i = 0; i < positions.size(); i += 3)
				{
					osg::Vec3d localPos(positions[i], positions[i + 1], positions[i + 2]);
					osg::Vec3d worldPos = worldMatrixd.preMult(localPos); // worldMatrix 是节点的世界矩阵
					//osg::Vec3d worldPos = worldMatrixd * localPos; // worldMatrix 是节点的世界矩阵

					//std::swap(worldPos.y(), worldPos.z());
					//worldPos.z() = -worldPos.z(); // 右手保持一致（必要时）

					worldPositions.push_back(worldPos.x());
					worldPositions.push_back(worldPos.y());
					worldPositions.push_back(worldPos.z());
				}
				tmpMesh.worldPositions = worldPositions;

				//获取索引
				std::vector<int> indices;
				if (prim.indices >= 0) {
					const tinygltf::Accessor& accessor = m_model.accessors[prim.indices];
					const tinygltf::BufferView& view = m_model.bufferViews[accessor.bufferView];
					const tinygltf::Buffer& buffer = m_model.buffers[view.buffer];
					const unsigned char* dataPtr = buffer.data.data() + view.byteOffset + accessor.byteOffset;
					for (size_t k = 0; k < accessor.count; ++k) {
						uint32_t index = 0;
						switch (accessor.componentType) {
						case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
							index = *(reinterpret_cast<const uint8_t*>(dataPtr + k));
							break;
						case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
							index = *(reinterpret_cast<const uint16_t*>(dataPtr + k * sizeof(uint16_t)));
							break;
						case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
							index = *(reinterpret_cast<const uint32_t*>(dataPtr + k * sizeof(uint32_t)));
							break;
						}
						indices.push_back(index);
					}
				}
				tmpMesh.indices = indices;
			}
			meshes.push_back(tmpMesh);
		}
	}
	return true;
}

void GltfManager::ComputeWorldMatrix(const tinygltf::Model& model, int nodeIndex, const osg::Matrixd& parentWorld, std::vector<osg::Matrixd>& worldMatrices)
{
	const tinygltf::Node& node = model.nodes[nodeIndex];
	osg::Matrixd local = GetLocalMatrix(node);
	osg::Matrixd world = local * parentWorld;
	worldMatrices[nodeIndex] = world;
	for (int child : node.children)
		ComputeWorldMatrix(model, child, world, worldMatrices);
}

osg::Matrixd GltfManager::GetLocalMatrix(const tinygltf::Node& node)
{
	osg::Matrixd m;
	if (!node.matrix.empty()) {
		m.set(node.matrix[0], node.matrix[1], node.matrix[2], node.matrix[3],
			node.matrix[4], node.matrix[5], node.matrix[6], node.matrix[7],
			node.matrix[8], node.matrix[9], node.matrix[10], node.matrix[11],
			node.matrix[12], node.matrix[13], node.matrix[14], node.matrix[15]);
	}
	else {
		// scale
		double sx = node.scale.size() > 0 ? node.scale[0] : 1.0;
		double sy = node.scale.size() > 1 ? node.scale[1] : 1.0;
		double sz = node.scale.size() > 2 ? node.scale[2] : 1.0;
		m.preMultScale(osg::Vec3d(sx, sy, sz));
		// rotation (quaternion w,x,y,z)
		if (node.rotation.size() == 4) {
			osg::Quat q(node.rotation[0], node.rotation[1],
				node.rotation[2], node.rotation[3]);
			m.preMultRotate(q);
		}
		// translation
		if (node.translation.size() == 3) {
			m.setTrans(osg::Vec3d(node.translation[0],
				node.translation[1],
				node.translation[2]));
		}
	}
	return m;
}



//bool GltfManager::ReadFile(const std::string& fileName) {
//	// 1️⃣ 创建 XCAF 文档
//	Handle(TDocStd_Document) doc;
//	XCAFApp_Application::GetApplication()->NewDocument("MDTV-XCAF", doc);
//	Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
//	if (shapeTool.IsNull()) {
//		std::cerr << "ShapeTool not found." << std::endl;
//		return 1;
//	}
//
//	// 2️⃣ 用 tinygltf 解析 GLB
//	tinygltf::Model model;
//	tinygltf::TinyGLTF loader;
//	std::string err, warn;
//	bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, fileName);
//	if (!ret) {
//		ret = loader.LoadBinaryFromFile(&model,&err,&warn,fileName);
//		if (!ret) {
//			std::cerr << "Failed to load GLB: " << err << std::endl;
//			return false;
//		}
//	}
//
//	//顶点
//	//std::vector<float> positions;
//	//索引
//	//std::vector<float> indices;
//	// 写入 DXF
//	std::ofstream dxf("outline.dxf");
//	dxf << "0\nSECTION\n2\nENTITIES\n";
//
//	std::ofstream polydxf("poly.dxf");
//	polydxf << "0\nSECTION\n2\nENTITIES\n";
//
//	for (const auto& mesh : model.meshes) {
//		double tolerance = 1.0e-4;
//		// sewing 用于把多个 face 合并成一个 shape
//		BRepBuilderAPI_Sewing sewing(tolerance, Standard_True, Standard_True, Standard_True, Standard_True);
//		OsgTest::Mesh tmpMesh;
//		for (const auto& prim : mesh.primitives) {
//			// 获取 POSITION
//			std::vector<float> positions;
//			auto accessorIndex = prim.attributes.find("POSITION");
//			if (accessorIndex != prim.attributes.end()) {
//				const tinygltf::Accessor& accessor = model.accessors[accessorIndex->second];
//				const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
//				const tinygltf::Buffer& buffer = model.buffers[view.buffer];
//				const unsigned char* dataPtr = buffer.data.data() + view.byteOffset + accessor.byteOffset;
//				size_t count = accessor.count;
//				size_t stride = accessor.ByteStride(view);
//				if (stride == 0) {
//					stride = tinygltf::GetNumComponentsInType(accessor.type) * tinygltf::GetComponentSizeInBytes(accessor.componentType);
//				}
//				for (size_t i = 0; i < count; ++i) {
//					const float* src = reinterpret_cast<const float*>(dataPtr + i * stride);
//					for (int j = 0; j < tinygltf::GetNumComponentsInType(accessor.type); ++j) {
//						positions.push_back(src[j]);
//					}
//				}
//			}
//			tmpMesh.positions = positions;
//			//获取索引
//			std::vector<int> indices;
//			if (prim.indices >= 0) {
//				const tinygltf::Accessor& accessor = model.accessors[prim.indices];
//				const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
//				const tinygltf::Buffer& buffer = model.buffers[view.buffer];
//				const unsigned char* dataPtr = buffer.data.data() + view.byteOffset + accessor.byteOffset;
//				for (size_t k = 0; k < accessor.count; ++k) {
//					uint32_t index = 0;
//					switch (accessor.componentType) {
//					case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
//						index = *(reinterpret_cast<const uint8_t*>(dataPtr + k));
//						break;
//					case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
//						index = *(reinterpret_cast<const uint16_t*>(dataPtr + k * sizeof(uint16_t)));
//						break;
//					case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
//						index = *(reinterpret_cast<const uint32_t*>(dataPtr + k * sizeof(uint32_t)));
//						break;
//					}
//					indices.push_back(index);
//				}
//			}
//			tmpMesh.indices = indices;
//
//			//构建边和面的对应
//			std::map<Edge, std::vector<int>> edgeToFaces;
//			size_t faceId = 0;
//			std::vector<Face> facesOut;
//			for (size_t i = 0; i + 2 < indices.size(); i += 3) {
//				size_t i0 = indices[i];
//				size_t i1 = indices[i + 1];
//				size_t i2 = indices[i + 2];
//
//				Face f{ {i0, i1, i2} };
//				gp_Pnt p1(positions[3 * indices[i] + 0], positions[3 * indices[i] + 1], positions[3 * indices[i] + 2]);
//				gp_Pnt p2(positions[3 * indices[i + 1] + 0], positions[3 * indices[i + 1] + 1], positions[3 * indices[i + 1] + 2]);
//				gp_Pnt p3(positions[3 * indices[i + 2] + 0], positions[3 * indices[i + 2] + 1], positions[3 * indices[i + 2] + 2]);
//				//计算面的法向量
//				gp_Vec v1(p1, p2);
//				gp_Vec v2(p1, p3);
//				gp_Vec n = v1.Crossed(v2);
//				if (n.Magnitude() > 1e-10)
//					n.Normalize();
//				f.normal = n;
//				facesOut.push_back(f);
//
//				// 三条边
//				Edge e1(i0, i1);
//				Edge e2(i1, i2);
//				Edge e3(i2, i0);
//
//				edgeToFaces[e1].push_back(faceId);
//				edgeToFaces[e2].push_back(faceId);
//				edgeToFaces[e3].push_back(faceId);
//
//				++faceId;
//			}
//			std::cout << "总共多少个面:" << faceId << std::endl;
//			//计算轮廓边
//			gp_Vec viewDir(0, 0, 1);
//			auto projectToPlane = [](const gp_Pnt& p) { return std::pair<double, double>{p.X(), p.Y()}; };
//			std::vector<Edge> silhouetteEdges;
//			for (const auto& kv : edgeToFaces) {
//				const Edge& edge = kv.first;
//				const std::vector<int>& faces = kv.second;
//
//				if (faces.size() == 1) {
//					// 自由边 → 轮廓边
//					silhouetteEdges.push_back(edge);
//				}
//				else if (faces.size() == 2) {
//					const gp_Vec& n1 = facesOut[faces[0]].normal;
//					const gp_Vec& n2 = facesOut[faces[1]].normal;
//					double dot1 = n1.Dot(viewDir);
//					double dot2 = n2.Dot(viewDir);
//					if (dot1 * dot2 < 0) {
//						silhouetteEdges.push_back(edge); // 剪影边
//					}
//					else if ((std::abs(dot1) == 1 && std::abs(dot2) == 0) || (std::abs(dot1) == 0 && std::abs(dot2) == 1)) {
//						silhouetteEdges.push_back(edge); // 剪影边
//					}
//				}
//			}
//
//			auto projectToXY = [](const gp_Pnt& p) -> std::pair<double, double>
//				{
//					return { p.X(), p.Z() };  // 直接丢掉 Z
//				};
//
//			for (const Edge& e : silhouetteEdges)
//			{
//				gp_Pnt p1(positions[3 * e.v1 + 0], positions[3 * e.v1 + 1], positions[3 * e.v1 + 2]);
//				gp_Pnt p2(positions[3 * e.v2 + 0], positions[3 * e.v2 + 1], positions[3 * e.v2 + 2]);
//
//				auto [x1, y1] = projectToPlane(p1);
//				auto [x2, y2] = projectToPlane(p2);
//
//				dxf << "0\nLINE\n8\n0\n";
//				dxf << "10\n" << x1 << "\n20\n" << y1 << "\n30\n0.0\n";
//				dxf << "11\n" << x2 << "\n21\n" << y2 << "\n31\n0.0\n";
//			}
//
//			// 按三角形生成 face
//			if (size(indices) > 0) {
//				for (size_t i = 0; i < size(indices); i += 3) {
//					gp_Pnt p1(positions[3 * indices[i] + 0],
//						positions[3 * indices[i] + 1],
//						positions[3 * indices[i] + 2]);
//					gp_Pnt p2(positions[3 * indices[i + 1] + 0],
//						positions[3 * indices[i + 1] + 1],
//						positions[3 * indices[i + 1] + 2]);
//					gp_Pnt p3(positions[3 * indices[i + 2] + 0],
//						positions[3 * indices[i + 2] + 1],
//						positions[3 * indices[i + 2] + 2]);
//
//					// 创建三角形 wire
//					BRepBuilderAPI_MakePolygon poly;
//					poly.Add(p1);
//					poly.Add(p2);
//					poly.Add(p3);
//					poly.Close();
//					TopoDS_Wire wire = poly.Wire();
//
//					// 计算三角形法线
//					gp_Vec v1(p2.XYZ() - p1.XYZ());
//					gp_Vec v2(p3.XYZ() - p1.XYZ());
//					gp_Dir normal(v1.Crossed(v2));
//					// 创建平面
//					gp_Pln plane(p1, normal);
//					// 创建 face
//					TopoDS_Face face = BRepBuilderAPI_MakeFace(plane, wire);
//					//if (face.Orientation() == TopAbs_REVERSED) {
//					//	face.Reverse(); // 确保所有面法线方向一致
//					//}
//					sewing.Add(face);
//				}
//			}
//
//			//直接使用网格三角形
//			TColgp_Array1OfPnt nodes(1, positions.size() / 3);
//			for (size_t i = 0; i < positions.size() / 3; ++i)
//				nodes.SetValue(i + 1, gp_Pnt(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]));
//
//			Poly_Array1OfTriangle triangles(1, indices.size() / 3);
//			for (size_t i = 0; i < indices.size() / 3; ++i)
//				triangles.SetValue(i + 1, Poly_Triangle(
//					indices[3 * i] + 1,
//					indices[3 * i + 1] + 1,
//					indices[3 * i + 2] + 1));
//
//			//Handle(Poly_Triangulation) triangulation = new Poly_Triangulation(nodes, triangles, Standard_False);
//			//// ==== 3. 创建一个占位 Face，并关联 Triangulation ====
//			//gp_Pln plane(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
//			//BRepBuilderAPI_MakeFace faceMaker(plane);
//			//TopoDS_Face face = faceMaker.Face();
//			//BRep_Builder builder;
//			//builder.UpdateFace(face, triangulation);
//
//			// 创建面并附加三角网格
//	/*		BRep_Builder builder;
//			TopoDS_Face face;
//			builder.MakeFace(face);
//			builder.UpdateFace(face, triangulation, Standard_True);*/
//
//			// ==== 4. HLR PolyAlgo 计算隐藏线 ====
//			//Handle(HLRBRep_PolyAlgo) hlrAlgo = new HLRBRep_PolyAlgo();
//			//gp_Ax2 view(gp::Origin(), gp::DZ());
//			//hlrAlgo->Projector(HLRAlgo_Projector(view));
//			//hlrAlgo->Load(face);
//			//hlrAlgo->Update();
//
//			//HLRBRep_PolyHLRToShape extractor;
//			//extractor.Update(hlrAlgo);
//			//TopoDS_Shape visibleEdges = extractor.VCompound();
//
//			// ==== 5. 输出 DXF ====
//
//			// 遍历 TopoDS_Edge 输出线段
//			//for (TopExp_Explorer exp(visibleEdges, TopAbs_EDGE); exp.More(); exp.Next()) {
//			//	TopoDS_Edge edge = TopoDS::Edge(exp.Current());
//
//			//	// 获取顶点
//			//	Standard_Real x1, y1, z1, x2, y2, z2;
//			//	TopoDS_Vertex v1, v2;
//			//	TopExp::Vertices(edge, v1, v2);
//			//	gp_Pnt p1 = BRep_Tool::Pnt(v1);
//			//	gp_Pnt p2 = BRep_Tool::Pnt(v2);
//
//			//	polydxf << "0\nLINE\n8\n0\n";
//			//	polydxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
//			//	polydxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
//			//}
//
//		}
//		meshes.push_back(tmpMesh);
//		// 合并成一个 shape
//		sewing.Perform();
//		TopoDS_Shape shape = sewing.SewedShape();
//		//// 三角化处理
//		//BRepMesh_IncrementalMesh mesh(shape, 1.0e-3);
//		//mesh.Perform();
//		// 加入 XCAF
//		shapeTool->AddShape(shape);
//	}
//
//	dxf << "0\nENDSEC\n0\nEOF\n";
//	dxf.close();
//
//	polydxf << "0\nENDSEC\n0\nEOF\n";
//	polydxf.close();
//
//	// 4️⃣ 输出统计
//	TDF_LabelSequence freeShapes;
//	shapeTool->GetFreeShapes(freeShapes);
//	std::cout << "Total shapes added: " << freeShapes.Length() << std::endl;
//
//
//	gp_Vec viewDir(1, 0, 0);
//	double eps = 1.0e-4;
//	for (Standard_Integer i = 1; i <= freeShapes.Length(); ++i) {
//		TDF_Label lab = freeShapes.Value(i);
//		TopoDS_Shape shape = shapeTool->GetShape(lab);
//
//		TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
//		TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);
//
//		////做二维投影
//		//Handle(HLRBRep_Algo) hlr = new HLRBRep_Algo();
//		//hlr->Add(shape);
//		////定义投影视角
//		//gp_Ax2 view(gp::Origin(), gp::DZ());
//		//hlr->Projector(HLRAlgo_Projector(view));
//		//hlr->Update();
//		//hlr->Hide();
//		//HLRBRep_HLRToShape hlrToShape(hlr);
//
//		//auto visibleEdges = hlrToShape.VCompound(); // 可见边
//		//auto hiddenEdges = hlrToShape.HCompound(); // 隐藏边
//		//auto rg1lineEdges = hlrToShape.Rg1LineVCompound();
//		//auto outlineEdges = hlrToShape.OutLineVCompound();//轮廓线
//
//
//		Handle(HLRBRep_PolyAlgo) polyAlgo = new HLRBRep_PolyAlgo();
//		gp_Ax2 view(gp::Origin(), gp::DZ());
//		polyAlgo->Projector(HLRAlgo_Projector(view));
//		polyAlgo->Load(shape);   // 加载模型
//		polyAlgo->Update();      // 计算隐藏线
//
//		// 提取可见边
//		HLRBRep_PolyHLRToShape extractor;
//		extractor.Update(polyAlgo);
//		TopoDS_Shape visibleEdges = extractor.VCompound(); // 可见线集合
//
//		//写dxf
//		std::ofstream visibledxf("visible.dxf");
//		visibledxf << "0\nSECTION\n2\nENTITIES\n";
//		for (TopExp_Explorer exp(visibleEdges, TopAbs_EDGE); exp.More(); exp.Next()) {
//			TopoDS_Edge edge = TopoDS::Edge(exp.Current());
//			gp_Pnt p1;
//			gp_Pnt p2;
//
//			TopoDS_Vertex v1, v2;
//			TopExp::Vertices(edge, v1, v2); // 获取 edge 两端顶点
//			p1 = BRep_Tool::Pnt(v1);
//			p2 = BRep_Tool::Pnt(v2);
//
//			//Standard_Real f, l;
//			//Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, f, l);
//			//if (!curve.IsNull()) {
//			//	p1 = curve->Value(f);
//			//	p2 = curve->Value(l);
//			//}
//			//else {
//			//	TopoDS_Vertex v1, v2;
//			//	TopExp::Vertices(edge, v1, v2); // 获取 edge 两端顶点
//			//	p1 = BRep_Tool::Pnt(v1);
//			//	p2 = BRep_Tool::Pnt(v2);
//			//}
//
//			visibledxf << "0\nLINE\n8\n0\n";
//			visibledxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
//			visibledxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
//		}
//		visibledxf << "0\nENDSEC\n0\nEOF\n";
//
//
//
//		//std::vector<TopoDS_Edge> contourEdges;
//		//for (int i = 1; i <= edgeToFaces.Extent(); i++) {
//		//	const TopoDS_Edge& edge = TopoDS::Edge(edgeToFaces.FindKey(i));
//		//	const TopTools_ListOfShape& faces = edgeToFaces.FindFromIndex(i);
//
//		//	if (faces.Extent() == 2) {
//		//		// 计算两个面的法向量
//		//		gp_Vec normal1 = GetFaceNormal(faces.First());
//		//		gp_Vec normal2 = GetFaceNormal(faces.Last());
//
//		//		// 检查法向量与视图方向的点积符号是否不同
//		//		double dot1 = normal1.Dot(viewDir);
//		//		double dot2 = normal2.Dot(viewDir);
//
//		//		if ((dot1 * dot2 < 0)) {
//		//			contourEdges.push_back(edge);
//		//		}
//		//		else if ((std::abs(dot1) == 1 && std::abs(dot2) == 0) || (std::abs(dot1) == 0 && std::abs(dot2) == 1)) {
//		//			contourEdges.push_back(edge);
//		//		}
//		//	}
//		//	else if (faces.Extent() == 1) {
//		//		//自由边
//		//		contourEdges.push_back(edge);
//		//	}
//		//}
//		//////写dxf
//		//std::ofstream dxf("hlr.dxf");
//		//dxf << "0\nSECTION\n2\nENTITIES\n";
//		//for (auto contourEdge : contourEdges) {
//		//	gp_Pnt p1;
//		//	gp_Pnt p2;
//
//		//	TopoDS_Vertex v1, v2;
//		//	TopExp::Vertices(contourEdge, v1, v2); // 获取 edge 两端顶点
//		//	p1 = BRep_Tool::Pnt(v1);
//		//	p2 = BRep_Tool::Pnt(v2);
//
//		//	dxf << "0\nLINE\n8\n0\n";
//		//	dxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
//		//	dxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
//		//}
//		//dxf << "0\nENDSEC\n0\nEOF\n";
//
//		////获取轮廓边
//
//		//检测是否为轮廓边
//		//TopTools_IndexedDataMapOfShapeListOfShape edgeToFaces;
//		//TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edgeToFaces);
//		//std::vector<TopoDS_Edge> contourEdges;
//		//for (int i = 1; i <= edgeToFaces.Extent(); i++) {
//		//	const TopoDS_Edge& edge = TopoDS::Edge(edgeToFaces.FindKey(i));
//		//	const TopTools_ListOfShape& faces = edgeToFaces.FindFromIndex(i);
//
//		//	if (faces.Extent() == 2) {
//		//		// 计算两个面的法向量
//		//		gp_Vec normal1 = GetFaceNormal(faces.First());
//		//		gp_Vec normal2 = GetFaceNormal(faces.Last());
//
//		//		// 检查法向量与视图方向的点积符号是否不同
//		//		double dot1 = normal1.Dot(viewDir);
//		//		double dot2 = normal2.Dot(viewDir);
//
//		//		if ((dot1 * dot2 <= 0)) {
//		//			contourEdges.push_back(edge);
//		//		}
//		//		else if ((std::abs(dot1) == 1 && std::abs(dot2) == 0) || (std::abs(dot1) == 0 && std::abs(dot2) == 1)) {
//		//			contourEdges.push_back(edge);
//		//		}
//		//	}
//		//	else if (faces.Extent() == 1) {
//		//		//自由边
//		//		contourEdges.push_back(edge);
//		//	}
//		//}
//		////写dxf
//		//std::ofstream dxf("test.dxf");
//		//dxf << "0\nSECTION\n2\nENTITIES\n";
//		//for (auto contourEdge : contourEdges) {
//		//	gp_Pnt p1;
//		//	gp_Pnt p2;
//
//		//	TopoDS_Vertex v1, v2;
//		//	TopExp::Vertices(contourEdge, v1, v2); // 获取 edge 两端顶点
//		//	p1 = BRep_Tool::Pnt(v1);
//		//	p2 = BRep_Tool::Pnt(v2);
//
//		//	dxf << "0\nLINE\n8\n0\n";
//		//	dxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
//		//	dxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
//		//}
//		//dxf << "0\nENDSEC\n0\nEOF\n";
//
//		////获取轮廓边
//		//std::vector<TopoDS_Edge> allEdges;
//		//for (int i = 1; i <= edgeToFaces.Extent(); i++) {
//		//	const TopoDS_Edge& edge = TopoDS::Edge(edgeToFaces.FindKey(i));
//		//	const TopTools_ListOfShape& faces = edgeToFaces.FindFromIndex(i);
//		//	allEdges.push_back(edge);
//		//}
//		////写dxf
//		//std::ofstream alldxf("all.dxf");
//		//alldxf << "0\nSECTION\n2\nENTITIES\n";
//		//for (auto contourEdge : allEdges) {
//		//	gp_Pnt p1;
//		//	gp_Pnt p2;
//
//		//	TopoDS_Vertex v1, v2;
//		//	TopExp::Vertices(contourEdge, v1, v2); // 获取 edge 两端顶点
//		//	p1 = BRep_Tool::Pnt(v1);
//		//	p2 = BRep_Tool::Pnt(v2);
//
//		//	alldxf << "0\nLINE\n8\n0\n";
//		//	alldxf << "10\n" << p1.X() << "\n20\n" << p1.Y() << "\n30\n" << p1.Z() << "\n";
//		//	alldxf << "11\n" << p2.X() << "\n21\n" << p2.Y() << "\n31\n" << p2.Z() << "\n";
//		//}
//		//alldxf << "0\nENDSEC\n0\nEOF\n";
//	}
//}








gp_Vec GltfManager::GetFaceNormal(const TopoDS_Shape& faceShape) {
	TopoDS_Face face = TopoDS::Face(faceShape);

	// 方法1：使用面的中心点计算法向量

	// 1. 获取面的几何表面
	Handle(Geom_Surface) surface = BRep_Tool::Surface(face);

	// 2. 获取参数范围
	Standard_Real uMin, uMax, vMin, vMax;
	BRepTools::UVBounds(face, uMin, uMax, vMin, vMax);

	// 3. 计算中心参数
	Standard_Real uCenter = (uMin + uMax) / 2.0;
	Standard_Real vCenter = (vMin + vMax) / 2.0;

	// 4. 使用几何属性计算法向量
	GeomLProp_SLProps props(surface, uCenter, vCenter, 1, 1e-6);

	if (props.IsNormalDefined()) {
		gp_Vec normal = props.Normal();
		// 5. 处理面的方向
		if (face.Orientation() == TopAbs_REVERSED) {
			normal.Reverse();
		}
		// 6. 归一化
		normal.Normalize();
		return normal;
	}
	// 备用方案：返回默认法向量
	return gp_Vec(0, 0, 1);
}