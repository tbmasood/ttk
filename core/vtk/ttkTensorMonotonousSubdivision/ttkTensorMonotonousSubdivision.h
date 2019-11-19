/// \ingroup vtk
/// \class ttkTensorMonotonousSubdivision
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the tensorMonotonousSubdivision processing
/// package.
///
/// VTK wrapping code for the @TensorMonotonousSubdivision package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::TensorMonotonousSubdivision
#pragma once

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// TTK code includes
#include <TensorMonotonousSubdivision.h>
#include <ttkWrapper.h>

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkTensorMonotonousSubdivision
#else
class ttkTensorMonotonousSubdivision
#endif
: public vtkDataSetAlgorithm,
public ttk::Wrapper {

public:
	static ttkTensorMonotonousSubdivision *New();
	vtkTypeMacro(ttkTensorMonotonousSubdivision, vtkDataSetAlgorithm);

	// default ttk setters
	vtkSetMacro(debugLevel_, int);

	vtkGetMacro(SubdivisionLevel, unsigned int);
	vtkSetMacro(SubdivisionLevel, unsigned int);

	void SetThreadNumber(int threadNumber) {
		ThreadNumber = threadNumber;
		SetThreads();
	}
	void SetUseAllCores(bool onOff) {
		UseAllCores = onOff;
		SetThreads();
	}

	int FillInputPortInformation(int /*port*/, vtkInformation *info) override {
		info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
		return 1;
	}

protected:
	ttkTensorMonotonousSubdivision() {
		SetNumberOfInputPorts(1);
		SetNumberOfOutputPorts(1);
	}

	TTK_SETUP();

	/**
	 * @brief Allocate an output array of same type that input array
	 */
	vtkSmartPointer<vtkDataArray>
	AllocateScalarField(vtkDataArray *const inputScalarField,
			int ntuples) const;

	int InterpolateScalarFields(vtkUnstructuredGrid *const input,
			vtkUnstructuredGrid *const output) const;

private:
	// number of subdivisions
	unsigned int SubdivisionLevel {1};

	// output 3D coordinates of generated points: old points first, then edge
	// middles, then triangle barycenters
	std::vector<float> points_ {};
	// output triangles
	std::vector<ttk::LongSimplexId> cells_ {};
	// generated point cell id
	std::vector<ttk::SimplexId> pointId_ {};
	// generated points dimension: 0 vertex of parent triangulation, 1 edge
	// middle, 2 triangle barycenter
	std::vector<ttk::SimplexId> pointDim_ {};
	std::vector<float> anisotropy_ {};
	std::vector<float> determinant_ {};
	std::vector<float> trace_ {};
	std::vector<float> lambda1_ {};
	std::vector<float> lambda2_ {};

	// base worker
	ttk::TensorMonotonousSubdivision baseWorker_ {points_, cells_, pointId_, pointDim_, anisotropy_, determinant_, trace_, lambda1_, lambda2_};
};
