#include "StepManager.h"
#include <XCAFApp_Application.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_Tool.hxx>
#include <TDataStd_Name.hxx>
#include <vector>
#include "Node.h"

void StepManager::Read(std::string fileName) {
	// 1. 创建 XDE 应用和文档
	Handle(XCAFApp_Application) app = XCAFApp_Application::GetApplication();
	Handle(TDocStd_Document) doc;
	app->NewDocument("XmlXCAF", doc);

	// 2. 使用 STEPCAFControl_Reader（支持 XDE）
	STEPCAFControl_Reader cafReader;
	//cafReader.SetColorMode(true);
	//cafReader.SetNameMode(true);
	//cafReader.SetLayerMode(true);

	// 3. 读取 STEP 文件
	IFSelect_ReturnStatus status = cafReader.ReadFile(fileName.c_str());
	if (status == IFSelect_RetDone) {
		cafReader.Transfer(doc);
	}
	else {
		std::cerr << "Read STEP failed." << std::endl;
		return;
	}

	// 4. 获取 XDE 工具
	Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
	Handle(XCAFDoc_ColorTool) colorTool = XCAFDoc_DocumentTool::ColorTool(doc->Main());

	// 5. 获取所有根组件（Top-Level Assemblies/Parts）
	TDF_LabelSequence labels;
	shapeTool->GetFreeShapes(labels); // 获取无父级的形状标签

	std::cout << "Found " << labels.Length() << " top-level components:" << std::endl;

	// 6. 递归遍历每个根组件
	for (int i = 1; i <= labels.Length(); i++) {
		TDF_Label compLabel = labels.Value(i);
		TraverseXDEComponent(compLabel, shapeTool, colorTool, 0, nodes, nullptr);
	}

}

void StepManager::TraverseXDEComponent(
	const TDF_Label& label,
	const Handle(XCAFDoc_ShapeTool)& shapeTool,
	const Handle(XCAFDoc_ColorTool)& colorTool,
	int depth,
	std::vector<Node>& nodes,
	Node* parentNode
)
{
	Node node;
	std::string indent(depth * 2, ' ');
	Handle(TDataStd_Name) nameAttr;
	std::string name = "<Unnamed>";
	// 查询标签是否包含 TDataStd_Name 属性
	if (label.FindAttribute(TDataStd_Name::GetID(), nameAttr)) {
		TCollection_ExtendedString extStr = nameAttr->Get();
		// 转换为 C 字符串
		name = Standard_CString(extStr.ToExtString());
	}
	node.name = name;
	// 获取形状
	TopoDS_Shape shape;
	if (shapeTool->GetShape(label, shape) && !shape.IsNull()) {
		std::cout << " ✅";
	}
	else {
		std::cout << " ❌";
	}
	node.shape = shape;
	// 获取位置（变换矩阵）
	TopLoc_Location loc;
	loc = shapeTool->GetLocation(label);
	if (!loc.IsIdentity()) {
		std::cout << indent << "Location: Not Identity" << std::endl;
	}
	node.location = loc;
	if (parentNode) {
		parentNode->children.push_back(node);
	}
	nodes.push_back(node);

	// 递归处理子组件
	TDF_LabelSequence children;
	shapeTool->GetComponents(label, children); // 获取子组件
	for (int i = 1; i <= children.Length(); i++) {
		TraverseXDEComponent(children.Value(i), shapeTool, colorTool, depth + 1, nodes, &node);
	}
}