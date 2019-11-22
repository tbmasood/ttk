#include <ttkTensorMonotonousSubdivision.h>

#define MODULE_ERROR_S MODULE_TMS "Error: "
#ifndef TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)    \
  if(COND) {                            \
    cerr << MODULE_ERROR_S MSG << endl; \
    return RET;                         \
  }
#else // TTK_ENABLE_KAMIKAZE
#define TTK_ABORT_KK(COND, MSG, RET)
#endif // TTK_ENABLE_KAMIKAZE

vtkStandardNewMacro (ttkTensorMonotonousSubdivision);

int ttkTensorMonotonousSubdivision::FillInputPortInformation(int port,
		vtkInformation *info) {
	switch (port) {
	case 0:
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
		break;
	}
	return 1;
}

int ttkTensorMonotonousSubdivision::FillOutputPortInformation(int port,
		vtkInformation *info) {
	switch (port) {
	case 0:
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
		break;
	}
	return 1;
}

vtkSmartPointer<vtkDataArray> ttkTensorMonotonousSubdivision::AllocateField(
		vtkDataArray *const inputField, int ntuples) const {

	vtkSmartPointer < vtkDataArray > res { };

	// allocate the memory for the output scalar field
	switch (inputField->GetDataType()) {
	case VTK_CHAR:
		res = vtkSmartPointer < vtkCharArray > ::New();
		break;
	case VTK_DOUBLE:
		res = vtkSmartPointer < vtkDoubleArray > ::New();
		break;
	case VTK_FLOAT:
		res = vtkSmartPointer < vtkFloatArray > ::New();
		break;
	case VTK_INT:
		res = vtkSmartPointer < vtkIntArray > ::New();
		break;
	case VTK_ID_TYPE:
		res = vtkSmartPointer < vtkIdTypeArray > ::New();
		break;
	case VTK_LONG:
		res = vtkSmartPointer < vtkLongArray > ::New();
		break;
	default:
		std::stringstream msg;
		msg << MODULE_TMS
		"Unsupported data array type" << endl;
		dMsg(std::cout, msg.str(), fatalMsg);
		break;
	}
	res->SetNumberOfComponents(inputField->GetNumberOfComponents());
	res->SetNumberOfTuples(ntuples);
	res->SetName(inputField->GetName());
	return res;
}

int ttkTensorMonotonousSubdivision::InterpolateFields(
		vtkDataSet *const input,
		vtkUnstructuredGrid *const output) const {

	const size_t npointdata = input->GetPointData()->GetNumberOfArrays();
	const size_t ncelldata = input->GetCellData()->GetNumberOfArrays();

	const auto outPointsNumber = baseWorker_.getNumberOfVertices();

	for (size_t i = 0; i < npointdata; ++i) {
		auto inputField = input->GetPointData()->GetArray(i);
		if (inputField == nullptr) {
			return -2;
		}

#define DISPATCH_INTERPOLATE_DIS(CASE, TYPE)                \
  case CASE:                                                \
    baseWorker_.interpolateDiscreteField<TYPE>(             \
      static_cast<TYPE *>(inputField->GetVoidPointer(0)),   \
      static_cast<TYPE *>(outputField->GetVoidPointer(0)),  \
		inputField->GetNumberOfComponents());               \
    break
#define DISPATCH_INTERPOLATE_CONT(CASE, TYPE)               \
  case CASE:                                                \
    baseWorker_.interpolateContinuousField<TYPE>(           \
      static_cast<TYPE *>(inputField->GetVoidPointer(0)),   \
      static_cast<TYPE *>(outputField->GetVoidPointer(0)),  \
		inputField->GetNumberOfComponents());               \
    break

		auto outputField = AllocateField(inputField, outPointsNumber);
		if (outputField == nullptr) {
			return -3;
		}

		switch (inputField->GetDataType()) {
		DISPATCH_INTERPOLATE_DIS(VTK_CHAR, char)
;			DISPATCH_INTERPOLATE_DIS(VTK_INT, int);
			DISPATCH_INTERPOLATE_DIS(VTK_LONG, long);
			DISPATCH_INTERPOLATE_DIS(VTK_ID_TYPE, vtkIdType);
			DISPATCH_INTERPOLATE_CONT(VTK_FLOAT, float);
			DISPATCH_INTERPOLATE_CONT(VTK_DOUBLE, double);
		}
		output->GetPointData()->AddArray(outputField);
	}

	const auto outCellsNumber = baseWorker_.getNumberOfTriangles();

	for (size_t i = 0; i < ncelldata; ++i) {
		auto inputField = input->GetCellData()->GetArray(i);
		if (inputField == nullptr) {
			return -2;
		}

		auto outputField = AllocateField(inputField, outCellsNumber);
		if (outputField == nullptr) {
			return -3;
		}

		switch (inputField->GetDataType()) {
		vtkTemplateMacro(
				baseWorker_.interpolateCellDataField < VTK_TT
						> (static_cast<VTK_TT*>(inputField->GetVoidPointer(0)), static_cast<VTK_TT*>(outputField->GetVoidPointer(
								0)), inputField->GetNumberOfComponents()));
		}
		output->GetCellData()->AddArray(outputField);
	}

	return 0;
}

int ttkTensorMonotonousSubdivision::doIt(std::vector<vtkDataSet*> &inputs,
		std::vector<vtkDataSet*> &outputs) {

	ttk::Memory m;

	vtkDataSet *input = inputs[0];
	auto output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);

	auto triangulation = ttkTriangulation::getTriangulation(input);
	ttk::Triangulation triangulationSubdivision;

	if (triangulation == nullptr) {
		return -1;
	}

	triangulation->setWrapper(this);

	baseWorker_.setupTriangulation(triangulation);
	baseWorker_.setWrapper(this);
	baseWorker_.setOutputTriangulation(&triangulationSubdivision);
	baseWorker_.setTensorDataArray(
			input->GetPointData()->GetTensors()->GetVoidPointer(0));
	baseWorker_.setSubdivisionField(SubdivisionField);
	baseWorker_.setGenerateAnisotropyField(GenerateAnisotropyField);
	baseWorker_.setGenerateDeterminantField(GenerateDeterminantField);
	baseWorker_.setGenerateTraceField(GenerateTraceField);
	baseWorker_.setGenerateEigenValuesField(GenerateEigenValuesField);

	switch (input->GetPointData()->GetTensors()->GetDataType()) {

	case VTK_CHAR:
		baseWorker_.execute<char>();
		break;

	case VTK_DOUBLE:
		baseWorker_.execute<double>();
		break;

	case VTK_FLOAT:
		baseWorker_.execute<float>();
		break;

	case VTK_INT:
		baseWorker_.execute<int>();
		break;

	case VTK_SHORT:
		baseWorker_.execute<short>();
		break;

	case VTK_LONG:
		baseWorker_.execute<long>();
		break;

	default:
		std::stringstream msg;
		msg << MODULE_TMS
		"Unsupported data type :(" << endl;
		dMsg(cerr, msg.str(), fatalMsg);
		return -3;
	}

	// first iteration: interpolate input scalar fields
	int ret = InterpolateFields(input, output);
	TTK_ABORT_KK(ret < 0, "Error interpolating input data array(s)", -1);

	// generated 3D coordinates
	auto points = vtkSmartPointer < vtkPoints > ::New();
	for (size_t i = 0; i < points_.size() / 3; i++) {
		points->InsertNextPoint(&points_[3 * i]);
	}
	output->SetPoints(points);

	// generated triangles
	const size_t dataPerCell = 4;
	auto cells = vtkSmartPointer < vtkCellArray > ::New();
	for (size_t i = 0; i < cells_.size() / dataPerCell; i++) {
		cells->InsertNextCell(3, &cells_[dataPerCell * i + 1]);
	}
	output->SetCells(VTK_TRIANGLE, cells);

	// cell dimension
	auto cellDim = vtkSmartPointer < ttkSimplexIdTypeArray > ::New();
	cellDim->SetName("CellDimension");
	cellDim->SetVoidArray(pointDim_.data(), pointDim_.size(), 1);
	output->GetPointData()->AddArray(cellDim);

	// anisotropy
	if (GenerateAnisotropyField) {
		auto anisotropy = vtkSmartPointer < vtkFloatArray > ::New();
		anisotropy->SetName("Anisotropy");
		anisotropy->SetVoidArray(anisotropy_.data(), anisotropy_.size(), 1);
		output->GetPointData()->AddArray(anisotropy);
	}

	// determinant
	if (GenerateDeterminantField) {
		auto determinant = vtkSmartPointer < vtkFloatArray > ::New();
		determinant->SetName("Determinant");
		determinant->SetVoidArray(determinant_.data(), determinant_.size(), 1);
		output->GetPointData()->AddArray(determinant);
	}

	// trace
	if (GenerateTraceField) {
		auto trace = vtkSmartPointer < vtkFloatArray > ::New();
		trace->SetName("Trace");
		trace->SetVoidArray(trace_.data(), trace_.size(), 1);
		output->GetPointData()->AddArray(trace);
	}

	// lambda1 and lambda2
	if (GenerateEigenValuesField) {
		auto lambda1 = vtkSmartPointer < vtkFloatArray > ::New();
		lambda1->SetName("EigenValue_High");
		lambda1->SetVoidArray(lambda1_.data(), lambda1_.size(), 1);
		output->GetPointData()->AddArray(lambda1);

		auto lambda2 = vtkSmartPointer < vtkFloatArray > ::New();
		lambda2->SetName("EigenValue_Low");
		lambda2->SetVoidArray(lambda2_.data(), lambda2_.size(), 1);
		output->GetPointData()->AddArray(lambda2);
	}

	{
		std::stringstream msg;
		msg << MODULE_TMS
		"Memory usage: " << m.getElapsedUsage() << " MB." << endl;
		dMsg(std::cout, msg.str(), memoryMsg);
	}

	return 0;
}
