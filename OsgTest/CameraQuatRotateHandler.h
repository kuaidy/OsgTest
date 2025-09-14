#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osgGA/GUIEventHandler>
#include <osg/Vec3>
#include <osg/Quat>

class CameraQuatRotateHandler : public osgGA::GUIEventHandler {
public:
    CameraQuatRotateHandler(osgViewer::Viewer* v);
	bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter&) override;
private:
    osgViewer::Viewer* viewer;
    osg::Vec3d center;
    double sensitivity = 0.5;

    // 上一次鼠标位置
    float lastX = 0.0f;
    float lastY = 0.0f;
    bool dragging = false;

    // 相机初始状态
    osg::Vec3d initialEye;
    osg::Vec3d initialUp;
    osg::Vec3d initialCenter;
};