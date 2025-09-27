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
#include <osg/PositionAttitudeTransform>

#include "DxfManager.h"

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>
#include <libavformat/avformat.h>
}

#include <ixwebsocket/IXWebSocketServer.h>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS_Face.hxx>
#include "StepManager.h"

osg::ref_ptr<osg::Node> createScene() {
	// 读取 cow.osg
	osg::ref_ptr<osg::Node> node = osgDB::readNodeFile("C:/Users/kdyonly/Desktop/cow.osg");
	return node;
}

// 存储所有活跃的 WebSocket 客户端
std::vector<std::weak_ptr<ix::WebSocket>> g_clients;

// 保护客户端列表的互斥锁
std::mutex g_clientsMutex;

std::string packetData;

int mainWebsocket() {
	const int W = 640, H = 480, FPS = 25;

	ix::initNetSystem();
	ix::WebSocketServer server(8008);
	server.setOnClientMessageCallback([](std::shared_ptr<ix::ConnectionState> connectionState, ix::WebSocket& webSocket, const ix::WebSocketMessagePtr& msg) {
		if (msg->type == ix::WebSocketMessageType::Open)
		{
			std::cout << "New connection" << std::endl;
			std::cout << "id: " << connectionState->getId() << std::endl;
			std::cout << "Uri: " << msg->openInfo.uri << std::endl;
			std::cout << "Headers:" << std::endl;
			for (auto it : msg->openInfo.headers)
			{
				std::cout << "\t" << it.first << ": " << it.second << std::endl;
			}
		}
		else if (msg->type == ix::WebSocketMessageType::Message)
		{
			std::cout << "Received: " << msg->str << std::endl;
			webSocket.send(packetData);
		}
		});

	auto res = server.listen();
	std::cout << "Listening on port " << server.getPort() << std::endl;
	// Per message deflate connection is enabled by default. It can be disabled
	// which might be helpful when running on low power devices such as a Rasbery Pi
	server.disablePerMessageDeflate();
	// Run the server in the background. Server can be stoped by calling server.stop()
	server.start();
	// Block until server.stop() is called.
	server.wait();

	// OSG Viewer
	osgViewer::Viewer viewer;
	osg::ref_ptr<osg::Group> root = new osg::Group;
	viewer.setSceneData(root);

	// 设置 FBO 渲染
	osg::ref_ptr<osg::Camera> camera = new osg::Camera;
	camera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);
	camera->setRenderOrder(osg::Camera::PRE_RENDER);
	camera->setViewport(0, 0, 1280, 720);
	root->addChild(camera);

	root->addChild(createScene());

	viewer.setUpViewInWindow(0, 0, 1280, 720);

	// 初始化 FFmpeg
	const AVCodec* codec = avcodec_find_encoder(AV_CODEC_ID_H264);
	AVCodecContext* c = avcodec_alloc_context3(codec);
	c->width = 1280;
	c->height = 720;
	c->pix_fmt = AV_PIX_FMT_YUV420P;
	c->time_base = { 1,30 };
	avcodec_open2(c, codec, nullptr);

	AVFrame* frame = av_frame_alloc();
	frame->format = c->pix_fmt;
	frame->width = c->width;
	frame->height = c->height;
	av_frame_get_buffer(frame, 0);

	SwsContext* sws_ctx = sws_getContext(
		c->width, c->height, AV_PIX_FMT_RGB24,
		c->width, c->height, AV_PIX_FMT_YUV420P,
		SWS_BICUBIC, nullptr, nullptr, nullptr
	);

	while (!viewer.done()) {
		viewer.frame();

		// 读取屏幕像素 RGB
		std::vector<unsigned char> rgbBuffer(1280 * 720 * 3);
		glReadPixels(0, 0, 1280, 720, GL_RGB, GL_UNSIGNED_BYTE, rgbBuffer.data());

		// 转 YUV
		uint8_t* inData[1] = { rgbBuffer.data() };
		int inLineSize[1] = { 3 * 1280 };
		sws_scale(sws_ctx, inData, inLineSize, 0, 720, frame->data, frame->linesize);

		// 编码 H.264
		AVPacket pkt;
		av_init_packet(&pkt);
		pkt.data = nullptr;
		pkt.size = 0;
		avcodec_send_frame(c, frame);
		if (avcodec_receive_packet(c, &pkt) == 0) {
			packetData = std::string(reinterpret_cast<const char*>(pkt.data), pkt.size);
			av_packet_unref(&pkt);
		}
	}
	return 0;
}



osg::ref_ptr<osg::Node> ShapeToOsg(const TopoDS_Shape& shape)
{
	// 1. 先做三角化
	BRepMesh_IncrementalMesh mesh(shape, 0.1); // 0.1 为容差，可调

	osg::ref_ptr<osg::Geode> geode = new osg::Geode();
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array();
	osg::ref_ptr<osg::DrawElementsUInt> indices = new osg::DrawElementsUInt(GL_TRIANGLES);

	// 2. 遍历面
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
		TopoDS_Face face = TopoDS::Face(exp.Current());
		TopLoc_Location loc;
		Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
		if (tri.IsNull()) continue;

		gp_Trsf trsf = loc.Transformation();

		Standard_Integer nb = tri->NbNodes();
		const Poly_Array1OfTriangle& tris = tri->Triangles();
		Standard_Integer vStart = vertices->size();

		for (Standard_Integer i = 1; i <= nb; ++i) {
			gp_Pnt p = tri->Node(i).Transformed(trsf);
			vertices->push_back(osg::Vec3(p.X(), p.Y(), p.Z()));
		}

		// 三角索引
		for (Standard_Integer i = tris.Lower(); i <= tris.Upper(); i++) {
			Standard_Integer n1, n2, n3;
			tris(i).Get(n1, n2, n3);
			indices->push_back(vStart + n1 - 1);
			indices->push_back(vStart + n2 - 1);
			indices->push_back(vStart + n3 - 1);
		}
	}

	geom->setVertexArray(vertices);
	geom->addPrimitiveSet(indices);
	geode->addDrawable(geom);
	return geode.release();
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


	StepManager stepManager;
	stepManager.Read("C:/Users/kdyonly/Desktop/assm1.STEP");
	for (int i = 0; i < stepManager.nodes.size();++i) {
		root->addChild(ShapeToOsg(stepManager.nodes[i].shape));
	}
	//// 读取 cow.osg
	//osg::ref_ptr<osg::Node> node = osgDB::readNodeFile("C:/Users/kdyonly/Desktop/cow.osg");
	//if (!node)
	//{
	//	std::cout << "模型加载失败，请检查 OSG 插件和文件路径！" << std::endl;
	//	return 1;
	//}
	//root->addChild(node.get());

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


