#include "MyTrackballManipulator.h"

MyTrackballManipulator::MyTrackballManipulator() {
	m_rotateCenter = osg::Vec3d(10, 0, 0);
}

osg::Matrixd MyTrackballManipulator::getMatrix() const
{
	return osg::Matrixd::translate(-m_rotateCenter)*
		osg::Matrixd::rotate(_rotation) *
		osg::Matrixd::translate(0., 0., _distance) *
		osg::Matrixd::translate(m_rotateCenter + _center);
}

//osg::Matrixd MyTrackballManipulator::getInverseMatrix() const {
//	return osg::Matrixd::translate(-m_rotateCenter) *
//		osg::Matrixd::rotate(_rotation.inverse()) *
//		osg::Matrixd::translate(m_rotateCenter - _center) *
//		osg::Matrixd::translate(0.0, 0.0, -_distance);
//}