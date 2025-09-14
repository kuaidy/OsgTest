#include <osgGA/TrackballManipulator>
class MyTrackballManipulator :public osgGA::TrackballManipulator {
public:
	MyTrackballManipulator();
	osg::Matrixd getMatrix() const override;
	//osg::Matrixd getInverseMatrix() const override;
private:
	osg::Matrix m_rotationMat;
	osg::Vec3d m_rotateCenter;
};