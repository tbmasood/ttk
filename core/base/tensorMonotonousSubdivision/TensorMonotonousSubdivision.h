/// \ingroup base
/// \class ttk::TensorMonotonousSubdivision
/// \author Talha Bin Masood (talha.bin.masood@liu.se)
/// \date November 2019
///
/// \brief Subdivise a triangulation according to triangle barycenter
///
/// %TensorMonotonousSubdivision generates a new, finer triangulation from
/// an input triangulation. Every triangle is divided in six new
/// triangles using the 3 edges middles and the triangle barycenter.
///
/// Scalar data on vertices (point data) with continuous values
/// (float/double) can be interpolated on the new
/// triangulation. Scalar data on input triangles can be replicated on
/// the new triangles.
///
/// \sa ttk::Triangulation
/// \sa ttkTensorMonotonousSubdivision.cpp %for a usage example.

#pragma once

#include <type_traits>

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>

#define MODULE_TMS "[TensorMonotonousSubdivision] "

namespace ttk {

class TensorMonotonousSubdivision: public Debug {

public:
	TensorMonotonousSubdivision(std::vector<float> &points,
			std::vector<LongSimplexId> &cells, std::vector<SimplexId> &pointDim,
			std::vector<float> &anisotropy, std::vector<float> &determinant,
			std::vector<float> &trace, std::vector<float> &lambda1,
			std::vector<float> &lambda2) :
			points_ { points }, cells_ { cells }, pointDim_ { pointDim }, anisotropy_ {
					anisotropy }, determinant_ { determinant }, trace_ { trace }, lambda1_ {
					lambda1 }, lambda2_ { lambda2 } {
	}

	inline void setOutputTriangulation(Triangulation *const triangulation) {
		outputTriangl_ = triangulation;
	}
	inline void setupTriangulation(Triangulation *const triangulation) {
		inputTriangl_ = triangulation;
		if (inputTriangl_ != nullptr) {
			inputTriangl_->preprocessVertexNeighbors();
			inputTriangl_->preprocessEdges();
			inputTriangl_->preprocessTriangles();
			inputTriangl_->preprocessTriangleEdges();
		}
		numInputVertices_ = inputTriangl_->getNumberOfVertices();
		numInputEdges_ = inputTriangl_->getNumberOfEdges();
		numInputTriangles_ = inputTriangl_->getNumberOfTriangles();
	}

	/** @brief Return the number of vertices in the output triangulation
	 */
	inline SimplexId getNumberOfVertices() const {
		return numOutputVertices_;
	}

	/** @brief Return the number of triangles in the output triangulation
	 */
	inline LongSimplexId getNumberOfTriangles() const {
		return numOutputTriangles_;
	}

	template<typename dataType>
	int execute() {

		Timer t;

		SimplexId vertexNumber = inputTriangl_->getNumberOfVertices();
		subdivideTriangulation<dataType>();
		buildOutputTriangulation();

		{
			std::stringstream msg;
			msg << MODULE_TMS "Data-set (" << vertexNumber
					<< " points) processed in " << t.getElapsedTime() << " s. ("
					<< threadNumber_ << " thread(s))." << std::endl;
			dMsg(std::cout, msg.str(), timeMsg);
		}

		return 0;
	}

	void setTensorDataArray(void *tensorData) {
		tensorData_ = tensorData;
	}

	int setSubdivisionField(const unsigned int subdivisionField) {
		subdivisionField_ = subdivisionField;
		return 0;
	}

	int setGenerateAnisotropyField(const bool generateAnisotropyField) {
		generateAnisotropyField_ = generateAnisotropyField;
		return 0;
	}

	int setGenerateDeterminantField(const bool generateDeterminantField) {
		generateDeterminantField_ = generateDeterminantField;
		return 0;
	}

	int setGenerateTraceField(const bool generateTraceField) {
		generateTraceField_ = generateTraceField;
		return 0;
	}

	int setGenerateEigenValuesField(const bool generateEigenValuesField) {
		generateEigenValuesField_ = generateEigenValuesField;
		return 0;
	}

	/**
	 * @brief Interpolate floating-point point data on subdivised triangulation
	 *
	 * Copy values on parent vertices, interpolate on edges and barycenters
	 *
	 * @param[in] data Pointer to input data on parent triangulation
	 * @param[out] output Allocated buffer to be filled
	 *
	 * @return 0 in case of success
	 */
	template<typename T>
	int interpolateContinuousField(const T *data, T *output,
			int numComponents) const {
		static_assert(
				std::is_floating_point<T>::value, "Floating point type required.");
		if (inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
			return 1;
		}
		const auto nOutVerts = this->getNumberOfVertices();
		if (nOutVerts < 0 || nOutVerts < numInputVertices_) {
			return 1;
		}

		// copy data on parent vertices
		std::copy(data, data + numInputVertices_ * numComponents, output);

		// interpolate on edges
		for (SimplexId i = 0; i < numDividedEdges_; ++i) {
			SimplexId a { }, b { };
			inputTriangl_->getEdgeVertex(pointEdgeMap_[i], 0, a);
			inputTriangl_->getEdgeVertex(pointEdgeMap_[i], 1, b);
			float t = edgeBarycenters_[i];
			for (int j = 0; j < numComponents; ++j) {
				output[numComponents * (numInputVertices_ + i) + j] = (1 - t)
						* data[numComponents * a + j]
						+ t * data[numComponents * b + j];
			}
		}

		// interpolate on triangle barycenters
		for (SimplexId i = 0; i < numDividedTriangles_; ++i) {
			SimplexId a { }, b { }, c { };
			inputTriangl_->getTriangleVertex(pointTriangleMap_[i], 0, a);
			inputTriangl_->getTriangleVertex(pointTriangleMap_[i], 1, b);
			inputTriangl_->getTriangleVertex(pointTriangleMap_[i], 2, c);
			float alpha = triangleBarycenters_[2 * i + 0];
			float beta = triangleBarycenters_[2 * i + 1];
			float gamma = 1 - alpha - beta;
			for (int j = 0; j < numComponents; ++j) {
				output[numComponents
						* (numInputVertices_ + numDividedEdges_ + i) + j] =
						alpha * data[numComponents * a + j]
								+ beta * data[numComponents * b + j]
								+ gamma * data[numComponents * c + j];
			}
		}
		return 0;
	}

	/**
	 * @brief Interpolate integer point data on subdivised triangulation
	 *
	 * Copy values on parent vertices, put 0 elsewhere.
	 *
	 * @param[in] data Pointer to input data on parent triangulation
	 * @param[out] output Allocated buffer to be filled
	 *
	 * @return 0 in case of success
	 */
	template<typename T>
	int interpolateDiscreteField(const T *data, T *output,
			int numComponents) const {
		static_assert(std::is_integral<T>::value, "Integral type required.");
		if (inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
			return 1;
		}
		const auto nOutVerts = this->getNumberOfVertices();
		if (nOutVerts < 0 || nOutVerts < numInputVertices_) {
			return 1;
		}
		std::fill(output, output + nOutVerts * numComponents, T { 0 });
		std::copy(data, data + numInputVertices_ * numComponents, output);
		return 0;
	}

	/**
	 * @brief Interpolate cell data on subdivised triangulation
	 *
	 * Copy parent triangle value on the six children triangles.
	 *
	 * @param[in] data Pointer to input data on parent triangulation
	 * @param[out] output Allocated buffer to be filled
	 *
	 * @return 0 in case of success
	 */
	template<typename T>
	int interpolateCellDataField(const T *data, T *output,
			int numComponents) const {
		if (inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
			return 1;
		}
		const auto nOutTriangles = this->getNumberOfTriangles();
		if (nOutTriangles < 0) {
			return 1;
		}
		for (LongSimplexId i = 0; i < nOutTriangles; ++i) {
			for (int j = 0; j < numComponents; ++j) {
				output[numComponents * i + j] = data[numComponents
						* originalCellMap_[i] + j];
			}
		}
		return 0;
	}

private:
	template<typename dataType>
	int subdivideTriangulation() {

		// not implemented for dimension >= 3
		if (inputTriangl_->getDimensionality() >= 3) {
			std::stringstream msg;
			msg << MODULE_TMS
			"Not yet implemented for dimension 3 and above" << std::endl;
			dMsg(std::cout, msg.str(), infoMsg);
			return 1;
		}

		dataType *tensorData = (dataType*) tensorData_;
		std::vector<bool> edgeHasCriticalPoint(numInputEdges_, false);
		std::vector<float> edgeBarycenters(numInputEdges_, -1.0);
		std::vector<float> edgeCriticalValues(numInputEdges_, -1.0);
		std::vector<bool> triangleHasCriticalPoint(numInputTriangles_, false);
		std::vector<float> triangleBarycenters(2 * numInputTriangles_, 0.0);

		numDividedEdges_ = 0;
		for (SimplexId i = 0; i < numInputEdges_; ++i) {
			SimplexId v1, v2;
			inputTriangl_->getEdgeVertex(i, 0, v1);
			inputTriangl_->getEdgeVertex(i, 1, v2);

			dataType v1_E = tensorData[9 * v1];
			dataType v1_F = tensorData[9 * v1 + 1];
			dataType v1_G = tensorData[9 * v1 + 4];

			dataType v2_E = tensorData[9 * v2];
			dataType v2_F = tensorData[9 * v2 + 1];
			dataType v2_G = tensorData[9 * v2 + 4];

			float A = (v2_E - v1_E - v2_G + v1_G) * (v2_E - v1_E - v2_G + v1_G)
					+ 4 * (v2_F - v1_F) * (v2_F - v1_F);
			float B = 2 * (v2_E - v1_E - v2_G + v1_G) * (v1_E - v1_G)
					+ 8 * v1_F * (v2_F - v1_F);
			float C = (v1_E - v1_G) * (v1_E - v1_G) + 4 * v1_F * v1_F;
			if (A != 0) {
				float t = -B / (2 * A);
				edgeBarycenters[i] = t;
				if (t > 0 && t < 1) {
					edgeHasCriticalPoint[i] = true;
					numDividedEdges_++;
					edgeCriticalValues[i] = A * t * t + B * t + C;
				}
			} else {
				edgeBarycenters[i] = 0;
			}
		}

		numDividedTriangles_ = 0;
		for (SimplexId i = 0; i < numInputTriangles_; ++i) {
			SimplexId v1, v2, v3;
			inputTriangl_->getTriangleVertex(i, 0, v1);
			inputTriangl_->getTriangleVertex(i, 1, v2);
			inputTriangl_->getTriangleVertex(i, 2, v3);

			dataType v1_E = tensorData[9 * v1];
			dataType v1_F = tensorData[9 * v1 + 1];
			dataType v1_G = tensorData[9 * v1 + 4];

			dataType v2_E = tensorData[9 * v2];
			dataType v2_F = tensorData[9 * v2 + 1];
			dataType v2_G = tensorData[9 * v2 + 4];

			dataType v3_E = tensorData[9 * v3];
			dataType v3_F = tensorData[9 * v3 + 1];
			dataType v3_G = tensorData[9 * v3 + 4];

			float Ex = v1_E - v3_E;
			float Fx = v1_F - v3_F;
			float Gx = v1_G - v3_G;

			float Ey = v2_E - v3_E;
			float Fy = v2_F - v3_F;
			float Gy = v2_G - v3_G;

			float Ec = v3_E;
			float Fc = v3_F;
			float Gc = v3_G;

			float A = (Ex - Gx) * (Ex - Gx) + 4 * Fx * Fx;
			float B = 2 * (Ex - Gx) * (Ey - Gy) + 8 * Fx * Fy;
			float C = (Ey - Gy) * (Ey - Gy) + 4 * Fy * Fy;
			float D = 2 * (Ex - Gx) * (Ec - Gc) + 8 * Fx * Fc;
			float E = 2 * (Ey - Gy) * (Ec - Gc) + 8 * Fy * Fc;

			float H = 4 * A * C - B * B;

			bool hasCrit = false;
			if (H != 0) {
				float alpha = (-2 * C * D + B * E) / H;
				float beta = (-2 * A * E + B * D) / H;
				float gamma = 1 - alpha - beta;
				triangleBarycenters[2 * i + 0] = alpha;
				triangleBarycenters[2 * i + 1] = beta;
				hasCrit = alpha > 0 && beta > 0 && gamma > 0;
				triangleHasCriticalPoint[i] = hasCrit;
				if (hasCrit) {
					numDividedTriangles_++;
				}
			} else {
				triangleBarycenters[2 * i + 0] = 1;
				triangleBarycenters[2 * i + 1] = 0;
			}
		}

		numOutputVertices_ = numInputVertices_ + numDividedEdges_
				+ numDividedTriangles_;
		points_.clear();
		points_.resize(3 * numOutputVertices_);
		pointDim_.clear();
		pointDim_.resize(numOutputVertices_);
		if (generateAnisotropyField_) {
			anisotropy_.clear();
			anisotropy_.resize(numOutputVertices_);
		}
		if (generateDeterminantField_) {
			determinant_.clear();
			determinant_.resize(numOutputVertices_);
		}
		if (generateTraceField_) {
			trace_.clear();
			trace_.resize(numOutputVertices_);
		}
		if (generateEigenValuesField_) {
			lambda1_.clear();
			lambda1_.resize(numOutputVertices_);
			lambda2_.clear();
			lambda2_.resize(numOutputVertices_);
		}

		for (SimplexId i = 0; i < numInputVertices_; ++i) {
			pointDim_[i] = 0;

			dataType E = tensorData[9 * i];
			dataType F = tensorData[9 * i + 1];
			dataType G = tensorData[9 * i + 4];
			float anisotropy = 0;
			if (generateAnisotropyField_ || generateEigenValuesField_) {
				anisotropy = (E - G) * (E - G) + 4 * F * F;
				anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
			}
			if (generateAnisotropyField_) {
				anisotropy_[i] = anisotropy;
			}
			if (generateDeterminantField_) {
				determinant_[i] = E * G - (F * F);
			}
			if (generateTraceField_) {
				trace_[i] = E + G;
			}
			if (generateEigenValuesField_) {
				lambda1_[i] = (E + G + anisotropy) / 2;
				lambda2_[i] = (E + G - anisotropy) / 2;
			}

			float v_x, v_y, v_z;
			inputTriangl_->getVertexPoint(i, v_x, v_y, v_z);
			points_[3 * i + 0] = v_x;
			points_[3 * i + 1] = v_y;
			points_[3 * i + 2] = v_z;
		}

		edgeBarycenters_.clear();
		edgeBarycenters_.resize(numDividedEdges_);
		pointEdgeMap_.clear();
		pointEdgeMap_.resize(numDividedEdges_);
		SimplexId index = numInputVertices_;
		std::vector < SimplexId > edgeIndexMap(numInputEdges_);
		for (SimplexId i = 0; i < numInputEdges_; ++i) {
			if (edgeHasCriticalPoint[i]) {
				SimplexId v1, v2;
				inputTriangl_->getEdgeVertex(i, 0, v1);
				inputTriangl_->getEdgeVertex(i, 1, v2);
				float v1_x, v1_y, v1_z;
				inputTriangl_->getVertexPoint(v1, v1_x, v1_y, v1_z);
				float v2_x, v2_y, v2_z;
				inputTriangl_->getVertexPoint(v2, v2_x, v2_y, v2_z);

				dataType v1_E = tensorData[9 * v1];
				dataType v1_F = tensorData[9 * v1 + 1];
				dataType v1_G = tensorData[9 * v1 + 4];

				dataType v2_E = tensorData[9 * v2];
				dataType v2_F = tensorData[9 * v2 + 1];
				dataType v2_G = tensorData[9 * v2 + 4];

				float t = edgeBarycenters[i];

				dataType E = t * v2_E + (1 - t) * v1_E;
				dataType F = t * v2_F + (1 - t) * v1_F;
				dataType G = t * v2_G + (1 - t) * v1_G;

				float anisotropy = 0;
				if (generateAnisotropyField_ || generateEigenValuesField_) {
					anisotropy = (E - G) * (E - G) + 4 * F * F;
					anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
				}
				if (generateAnisotropyField_) {
					anisotropy_[index] = anisotropy;
				}
				if (generateDeterminantField_) {
					determinant_[index] = E * G - (F * F);
				}
				if (generateTraceField_) {
					trace_[index] = E + G;
				}
				if (generateEigenValuesField_) {
					lambda1_[index] = (E + G + anisotropy) / 2;
					lambda2_[index] = (E + G - anisotropy) / 2;
				}

				pointDim_[index] = 1;
				points_[3 * index + 0] = t * v2_x + (1 - t) * v1_x;
				points_[3 * index + 1] = t * v2_y + (1 - t) * v1_y;
				points_[3 * index + 2] = t * v2_z + (1 - t) * v1_z;
				edgeBarycenters_[index - numInputVertices_] = t;
				pointEdgeMap_[index - numInputVertices_] = i;

				edgeIndexMap[i] = index;
				++index;
			} else {
				edgeIndexMap[i] = -1;
			}
		}

		triangleBarycenters_.clear();
		triangleBarycenters_.resize(2 * numDividedTriangles_);
		pointTriangleMap_.clear();
		pointTriangleMap_.resize(numDividedTriangles_);
		std::vector < SimplexId > triangleIndexMap(numInputTriangles_);
		for (SimplexId i = 0; i < numInputTriangles_; ++i) {
			if (triangleHasCriticalPoint[i]) {
				SimplexId v1, v2, v3;
				inputTriangl_->getTriangleVertex(i, 0, v1);
				inputTriangl_->getTriangleVertex(i, 1, v2);
				inputTriangl_->getTriangleVertex(i, 2, v3);
				float v1_x, v1_y, v1_z;
				inputTriangl_->getVertexPoint(v1, v1_x, v1_y, v1_z);
				float v2_x, v2_y, v2_z;
				inputTriangl_->getVertexPoint(v2, v2_x, v2_y, v2_z);
				float v3_x, v3_y, v3_z;
				inputTriangl_->getVertexPoint(v3, v3_x, v3_y, v3_z);

				dataType v1_E = tensorData[9 * v1];
				dataType v1_F = tensorData[9 * v1 + 1];
				dataType v1_G = tensorData[9 * v1 + 4];

				dataType v2_E = tensorData[9 * v2];
				dataType v2_F = tensorData[9 * v2 + 1];
				dataType v2_G = tensorData[9 * v2 + 4];

				dataType v3_E = tensorData[9 * v3];
				dataType v3_F = tensorData[9 * v3 + 1];
				dataType v3_G = tensorData[9 * v3 + 4];

				float alpha = triangleBarycenters[2 * i + 0];
				float beta = triangleBarycenters[2 * i + 1];
				float gamma = 1 - alpha - beta;

				dataType E = alpha * v1_E + beta * v2_E + gamma * v3_E;
				dataType F = alpha * v1_F + beta * v2_F + gamma * v3_F;
				dataType G = alpha * v1_G + beta * v2_G + gamma * v3_G;

				float anisotropy = 0;
				if (generateAnisotropyField_ || generateEigenValuesField_) {
					anisotropy = (E - G) * (E - G) + 4 * F * F;
					anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
				}
				if (generateAnisotropyField_) {
					anisotropy_[index] = anisotropy;
				}
				if (generateDeterminantField_) {
					determinant_[index] = E * G - (F * F);
				}
				if (generateTraceField_) {
					trace_[index] = E + G;
				}
				if (generateEigenValuesField_) {
					lambda1_[index] = (E + G + anisotropy) / 2;
					lambda2_[index] = (E + G - anisotropy) / 2;
				}
				
				SimplexId temp = 2
						* (index - numInputVertices_ - numDividedEdges_);
				triangleBarycenters_[temp] = alpha;
				triangleBarycenters_[temp + 1] = beta;
				pointTriangleMap_[index - numInputVertices_ - numDividedEdges_] =
						i;

				points_[3 * index + 0] = alpha * v1_x + beta * v2_x
						+ gamma * v3_x;
				points_[3 * index + 1] = alpha * v1_y + beta * v2_y
						+ gamma * v3_y;
				points_[3 * index + 2] = alpha * v1_z + beta * v2_z
						+ gamma * v3_z;

				pointDim_[index] = 2;
				triangleIndexMap[i] = index;
				++index;
			} else {
				triangleIndexMap[i] = -1;
			}
		}

		cells_.clear();
		cells_.reserve(numInputTriangles_ * 4 * 6);
		originalCellMap_.clear();
		originalCellMap_.reserve(numInputTriangles_ * 6);
		SimplexId cellIndex = 0;
		for (LongSimplexId i = 0; i < numInputTriangles_; ++i) {
			SimplexId e1, e2, e3;
			inputTriangl_->getTriangleEdge(i, 0, e1);
			inputTriangl_->getTriangleEdge(i, 1, e2);
			inputTriangl_->getTriangleEdge(i, 2, e3);

			SimplexId e1v1, e1v2;
			inputTriangl_->getEdgeVertex(e1, 0, e1v1);
			inputTriangl_->getEdgeVertex(e1, 1, e1v2);
			SimplexId e2v1, e2v2;
			inputTriangl_->getEdgeVertex(e2, 0, e2v1);
			inputTriangl_->getEdgeVertex(e2, 1, e2v2);
			SimplexId e3v1, e3v2;
			inputTriangl_->getEdgeVertex(e3, 0, e3v1);
			inputTriangl_->getEdgeVertex(e3, 1, e3v2);

			SimplexId v1 = e1v1;
			SimplexId v2 = e1v2;
			SimplexId v3 = (e2v1 == v1 || e2v1 == v2) ? e2v2 : e2v1;
			if (e2v1 == v3) {
				if (e2v2 == v1) {
					SimplexId temp = e2;
					e2 = e3;
					e3 = temp;
				}
			} else {
				if (e2v1 == v1) {
					SimplexId temp = e2;
					e2 = e3;
					e3 = temp;
				}
			}

			SimplexId v12, v23, v13, v123;
			if (triangleHasCriticalPoint[i]) {
				v123 = triangleIndexMap[i];
				if (edgeHasCriticalPoint[e1]) {
					v12 = edgeIndexMap[e1];
					addOutputTriangle(cellIndex, v123, v1, v12, i);
					addOutputTriangle(cellIndex, v123, v12, v2, i);
				} else {
					addOutputTriangle(cellIndex, v123, v1, v2, i);
				}
				if (edgeHasCriticalPoint[e2]) {
					v23 = edgeIndexMap[e2];
					addOutputTriangle(cellIndex, v123, v2, v23, i);
					addOutputTriangle(cellIndex, v123, v23, v3, i);
				} else {
					addOutputTriangle(cellIndex, v123, v2, v3, i);
				}
				if (edgeHasCriticalPoint[e3]) {
					v13 = edgeIndexMap[e3];
					addOutputTriangle(cellIndex, v123, v3, v13, i);
					addOutputTriangle(cellIndex, v123, v13, v1, i);
				} else {
					addOutputTriangle(cellIndex, v123, v3, v1, i);
				}
			} else {
				int numCrits = 0;
				if (edgeHasCriticalPoint[e1]) {
					v12 = edgeIndexMap[e1];
					numCrits++;
				}
				if (edgeHasCriticalPoint[e2]) {
					v23 = edgeIndexMap[e2];
					numCrits++;
				}
				if (edgeHasCriticalPoint[e3]) {
					v13 = edgeIndexMap[e3];
					numCrits++;
				}
				switch (numCrits) {
				case 0:
					addOutputTriangle(cellIndex, v1, v2, v3, i);
					break;
				case 1:
					if (edgeHasCriticalPoint[e1]) {
						addOutputTriangle(cellIndex, v12, v3, v1, i);
						addOutputTriangle(cellIndex, v12, v3, v2, i);
					} else if (edgeHasCriticalPoint[e2]) {
						addOutputTriangle(cellIndex, v23, v1, v2, i);
						addOutputTriangle(cellIndex, v23, v1, v3, i);
					} else if (edgeHasCriticalPoint[e3]) {
						addOutputTriangle(cellIndex, v13, v2, v3, i);
						addOutputTriangle(cellIndex, v13, v2, v1, i);
					}
					break;
				case 2:
					if (!edgeHasCriticalPoint[e1]) {
						addOutputTriangle(cellIndex, v23, v13, v3, i);
						if (edgeCriticalValues[e2] < edgeCriticalValues[e3]) {
							addOutputTriangle(cellIndex, v23, v13, v1, i);
							addOutputTriangle(cellIndex, v23, v2, v1, i);
						} else {
							addOutputTriangle(cellIndex, v13, v23, v2, i);
							addOutputTriangle(cellIndex, v13, v1, v2, i);
						}
					} else if (!edgeHasCriticalPoint[e2]) {
						addOutputTriangle(cellIndex, v12, v13, v1, i);
						if (edgeCriticalValues[e1] < edgeCriticalValues[e3]) {
							addOutputTriangle(cellIndex, v12, v13, v3, i);
							addOutputTriangle(cellIndex, v12, v2, v3, i);
						} else {
							addOutputTriangle(cellIndex, v13, v12, v2, i);
							addOutputTriangle(cellIndex, v13, v3, v2, i);
						}
					} else if (!edgeHasCriticalPoint[e3]) {
						addOutputTriangle(cellIndex, v12, v23, v2, i);
						if (edgeCriticalValues[e1] < edgeCriticalValues[e2]) {
							addOutputTriangle(cellIndex, v12, v23, v3, i);
							addOutputTriangle(cellIndex, v12, v1, v3, i);
						} else {
							addOutputTriangle(cellIndex, v23, v12, v1, i);
							addOutputTriangle(cellIndex, v23, v3, v1, i);
						}
					}
					break;
				case 3:
					float val12 = edgeCriticalValues[e1];
					float val23 = edgeCriticalValues[e2];
					float val13 = edgeCriticalValues[e3];
					if (val12 < val23 && val12 < val13) {
						addOutputTriangle(cellIndex, v12, v1, v13, i);
						addOutputTriangle(cellIndex, v12, v13, v3, i);
						addOutputTriangle(cellIndex, v12, v3, v23, i);
						addOutputTriangle(cellIndex, v12, v23, v2, i);
					} else if (val23 < val12 && val23 < val13) {
						addOutputTriangle(cellIndex, v23, v2, v12, i);
						addOutputTriangle(cellIndex, v23, v12, v1, i);
						addOutputTriangle(cellIndex, v23, v1, v13, i);
						addOutputTriangle(cellIndex, v23, v13, v3, i);
					} else {
						addOutputTriangle(cellIndex, v13, v3, v23, i);
						addOutputTriangle(cellIndex, v13, v23, v2, i);
						addOutputTriangle(cellIndex, v13, v2, v12, i);
						addOutputTriangle(cellIndex, v13, v12, v1, i);
					}
					break;
				}
			}
		}
		cells_.resize(cellIndex);
		numOutputTriangles_ = cellIndex / 4;
		originalCellMap_.resize(numOutputTriangles_);
		return 0;
	}

	int buildOutputTriangulation();
	void addOutputTriangle(int &offset, const SimplexId &v1,
			const SimplexId &v2, const SimplexId &v3,
			const LongSimplexId &inputCellIndex);
	float intersectSegments(const float &x1, const float &y1, const float &x2,
			const float &y2, const float &x3, const float &y3, const float &x4,
			const float &y4);

	unsigned int subdivisionField_ { 0 };
	bool generateAnisotropyField_ { true };
	bool generateDeterminantField_ { false };
	bool generateTraceField_ { false };
	bool generateEigenValuesField_ { false };

	SimplexId numInputVertices_ { };
	SimplexId numInputEdges_ { };
	SimplexId numInputTriangles_ { };

	SimplexId numDividedEdges_ { };
	SimplexId numDividedTriangles_ { };

	SimplexId numOutputVertices_ { };
	SimplexId numOutputTriangles_ { };

	Triangulation *inputTriangl_ { };
	void *tensorData_ { };

	std::vector<float> &points_;
	std::vector<LongSimplexId> &cells_;
	std::vector<SimplexId> &pointDim_;
	std::vector<float> &anisotropy_;
	std::vector<float> &determinant_;
	std::vector<float> &trace_;
	std::vector<float> &lambda1_;
	std::vector<float> &lambda2_;

	std::vector<float> edgeBarycenters_ { };
	std::vector<float> triangleBarycenters_ { };
	std::vector<SimplexId> pointEdgeMap_ { };
	std::vector<LongSimplexId> pointTriangleMap_ { };
	std::vector<LongSimplexId> originalCellMap_ { };

	// output triangulation built on output points & output cells
	Triangulation *outputTriangl_ { };
};
} // namespace ttk
