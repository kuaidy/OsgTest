# CMakeList.txt: OsgTest 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
# 查找 OpenSceneGraph
include_directories("C:/MyProject/OsgTest/ThirdParty/osg-3.6.3/include")
link_directories("C:/MyProject/OsgTest/ThirdParty/osg-3.6.3/lib")

# tinyfltf
include_directories("C:/MyProject/OsgTest/ThirdParty/tinygltf-2.9.5/include")
link_directories("C:/MyProject/OsgTest/ThirdParty/tinygltf-2.9.5/libd")

# occt
include_directories("C:/MyProject/OsgTest/ThirdParty/occt/install/inc")
link_directories("C:/MyProject/OsgTest/ThirdParty/occt/install/win64/vc14/libd")

#ffmpeg
include_directories("C:/MyProject/OsgTest/ThirdParty/ffmpeg-n8.0-latest-win64-gpl-shared-8.0/include")
link_directories("C:/MyProject/OsgTest/ThirdParty/ffmpeg-n8.0-latest-win64-gpl-shared-8.0/lib")

include_directories("C:/MyProject/vcpkg/installed/x64-windows/include")
link_directories("C:/MyProject/vcpkg/installed/x64-windows/debug/lib")

# 将源代码添加到此项目的可执行文件。
add_executable (OsgTest "OsgTest.cpp" "OsgTest.h" "CameraQuatRotateHandler.cpp" "CameraQuatRotateHandler.h" "MyTrackballManipulator.cpp" "GltfManager.h" "GltfManager.cpp" "Mesh.h" "Mesh.cpp" "DxfManager.h" "DxfManager.cpp" "Edge.h" "Face.h" "StepManager.h" "StepManager.cpp" "Model.h" "Model.cpp" "Node.h" "Node.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET OsgTest PROPERTY CXX_STANDARD 20)
endif()

# occ
set(OCCT_LIBS
    TKernel 
    TKMath 
    TKBRep 
    TKLCAF 
    TKXCAF 
    TKRWMesh 
    TKHLR 
    TKXDE 
    TKG3d 
    TKTopAlgo 
    TKGeomBase 
    TKMesh
    TKSTEPBase
    TKCAF
    TKSTEP
    TKXSBase
    TKPrim
    TKBO
    TKService
    TKV3d
    TKOpenGl
    TKXDESTEP
    TKSTEPAttr
    TKSTEP209
)


target_link_libraries(
    OsgTest
    debug 
    osgdbd 
    osgGAd 
    osgUtild 
    osgViewerd 
    osgd 
    OpenThreadsd 
    tinygltf
    ${OCCT_LIBS}
    ixwebsocket 
    Ws2_32
    mbedcrypto
    mbedtls
    mbedx509
    crypt32
    zlibd
    bcrypt
    opengl32
    avutil
    avcodec
    swscale
)