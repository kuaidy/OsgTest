// OsgTest.cpp
#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osgGA/OrbitManipulator>
#include "CameraQuatRotateHandler.h"
#include "MyTrackballManipulator.h"
#include <osg/ShapeDrawable>
#include "GltfManager.h"
//#include <osgDB/Registry>
#include <osg/MatrixTransform>

#include "DxfManager.h"


osg::ref_ptr<osg::Node> CreateMeshNode(const std::vector<float>& verts, const std::vector<int>& idx)
{
	osg::ref_ptr<osg::Vec3Array> osgVerts = new osg::Vec3Array;
	osgVerts->reserve(verts.size() / 3);
	for (size_t i = 0; i < verts.size(); i += 3)
	{
		osgVerts->push_back(osg::Vec3(verts[i], verts[i + 1], verts[i + 2]));
	}

	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
	geom->setVertexArray(osgVerts.get());

	if (!idx.empty()) {
		osg::ref_ptr<osg::DrawElementsUInt> osgIndices = new osg::DrawElementsUInt(GL_TRIANGLES);
		osgIndices->insert(osgIndices->end(), idx.begin(), idx.end());
		geom->addPrimitiveSet(osgIndices.get());
	}
	else {
		geom->addPrimitiveSet(new osg::DrawArrays(GL_TRIANGLES, 0, osgVerts->size()));
	}

	// 可选：法线、UV、颜色
	// geom->setNormalArray(normals, osg::Array::BIND_PER_VERTEX);
	// geom->setTexCoordArray(0, texcoords);

	osg::ref_ptr<osg::Geode> geode = new osg::Geode;
	geode->addDrawable(geom.get());
	return geode.get();
}


int main()
{
	// 创建 Viewer
	osgViewer::Viewer viewer;
	// 设置窗口矩形: x, y, width, height
	osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
	traits->x = 100;       // 窗口左上角X位置
	traits->y = 100;       // 窗口左上角Y位置
	traits->width = 800;   // 窗口宽度
	traits->height = 600;  // 窗口高度
	traits->windowDecoration = true; // 显示窗口边框和标题栏
	traits->doubleBuffer = true;

	osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());

	osg::ref_ptr<osg::Camera> camera = viewer.getCamera();
	camera->setGraphicsContext(gc.get());
	camera->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));
	GLenum buffer = traits->doubleBuffer ? GL_BACK : GL_FRONT;
	camera->setDrawBuffer(buffer);
	camera->setReadBuffer(buffer);

	viewer.setCameraManipulator(new MyTrackballManipulator());


	osg::ref_ptr<osg::Group> root = new osg::Group;
	//// 读取 cow.osg
	//osg::ref_ptr<osg::Node> node = osgDB::readNodeFile("C:/Users/kdyonly/Desktop/cow.osg");
	//if (!node)
	//{
	//	std::cout << "模型加载失败，请检查 OSG 插件和文件路径！" << std::endl;
	//	return 1;
	//}
	//root->addChild(node.get());

	//读取glb
	GltfManager gltfManager;
	gltfManager.ReadFile("C:/Users/kdyonly/Desktop/model.glb");
	for (auto mesh : gltfManager.meshes) {
		//osg::ref_ptr<osg::MatrixTransform> transform = new osg::MatrixTransform();
		//transform->setMatrix(gltfManager._worldMatrices[mesh.nodeIndex]); // 世界矩阵
		//transform->addChild(CreateMeshNode(mesh.positions, mesh.indices));
		root->addChild(CreateMeshNode(mesh.worldPositions, mesh.indices));
	}
	DxfManager dxfManager;
	dxfManager.Write("", gltfManager.meshes);

	// 创建旋转中心标记（小球）
	//osg::ref_ptr<osg::Sphere> sphere = new osg::Sphere(osg::Vec3(10, 0, 0),0.1); // 半径 0.1
	//osg::ref_ptr<osg::ShapeDrawable> sphereDrawable = new osg::ShapeDrawable(sphere);
	//osg::ref_ptr<osg::Geode> sphereGeode = new osg::Geode;
	//sphereGeode->addDrawable(sphereDrawable.get());
	//root->addChild(sphereGeode);

	// 设置场景数据
	viewer.setSceneData(root);
	// 启动渲染循环
	return viewer.run();
}


