#include "CameraQuatRotateHandler.h"

CameraQuatRotateHandler::CameraQuatRotateHandler(osgViewer::Viewer* v):viewer(v){
    osg::Matrix viewMat = viewer->getCamera()->getViewMatrix();
    viewMat.getLookAt(initialEye, initialCenter,initialUp);
}

bool CameraQuatRotateHandler::handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa) {
    switch (ea.getEventType()) {
    case osgGA::GUIEventAdapter::PUSH:
        lastX = ea.getX();
        lastY = ea.getY();
        dragging = true;
        return true;

    case osgGA::GUIEventAdapter::RELEASE:
        dragging = false;
        return true;

    case osgGA::GUIEventAdapter::DRAG:
        if (!dragging) return false;
        {
            float deltaX = ea.getX() - lastX;
            float deltaY = ea.getY() - lastY;
            lastX = ea.getX();
            lastY = ea.getY();

            //std::cout << deltaX << std::endl;
            //std::cout << lastY << std::endl;
            
            // 计算相机到中心向量
            osg::Vec3d dir = initialEye - center;

            // 水平旋转（绕 up）
            osg::Quat quatY(deltaX * sensitivity, initialUp);
            dir = quatY * dir;

            // 垂直旋转（绕右向量）
            osg::Vec3d right = dir ^ initialUp;
            right.normalize();
            osg::Quat quatX(deltaY * sensitivity, right);
            dir = quatX * dir;

            osg::Vec3d newEye = center + dir;
            viewer->setCameraManipulator(nullptr);
            viewer->getCamera()->setViewMatrixAsLookAt(newEye, center, initialUp);
            
            return true;
        }

    default:
        return false;
    }
}