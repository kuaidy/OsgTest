// OsgTest.cpp
#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osgGA/OrbitManipulator>
#include "CameraQuatRotateHandler.h"
#include "MyTrackballManipulator.h"
#include <osg/ShapeDrawable>
#include "GltfManager.h"
//#include <osgDB/Registry>

int main()
{
	GltfManager gltfManager;
	gltfManager.ReadFile("C:/Users/kdyonly/Desktop/chahu.glb");
	return 0;
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
	// 读取 cow.osg
	osg::ref_ptr<osg::Node> node = osgDB::readNodeFile("C:/Users/kdyonly/Desktop/cow.osg");
	if (!node)
	{
		std::cout << "模型加载失败，请检查 OSG 插件和文件路径！" << std::endl;
		return 1;
	}
	root->addChild(node.get());

	// 创建旋转中心标记（小球）
	osg::ref_ptr<osg::Sphere> sphere = new osg::Sphere(osg::Vec3(10, 0, 0),0.1); // 半径 0.1
	osg::ref_ptr<osg::ShapeDrawable> sphereDrawable = new osg::ShapeDrawable(sphere);
	osg::ref_ptr<osg::Geode> sphereGeode = new osg::Geode;
	sphereGeode->addDrawable(sphereDrawable.get());
	root->addChild(sphereGeode);

	// 设置场景数据
	viewer.setSceneData(root);
	// 启动渲染循环
	return viewer.run();
}