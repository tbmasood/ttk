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

#define MODULE_S "[TensorMonotonousSubdivision] "

namespace ttk {

class TensorMonotonousSubdivision: public Debug {

public:
	TensorMonotonousSubdivision(std::vector<float> &points,
			std::vector<LongSimplexId> &cells, std::vector<SimplexId> &pointId,
			std::vector<SimplexId> &pointDim, std::vector<float> &anisotropy,
			std::vector<float> &determinant, std::vector<float> &trace,
			std::vector<float> &lambda1, std::vector<float> &lambda2) :
			points_ { points }, cells_ { cells }, pointId_ { pointId }, pointDim_ {
					pointDim }, anisotropy_ { anisotropy }, determinant_ {
					determinant }, trace_ { trace }, lambda1_ { lambda1 }, lambda2_ {
					lambda2 } {
	}

	inline void setOutputTriangulation(Triangulation *const triangulation) {
		outputTriangl_ = triangulation;
	}
	inline void setInputPoints(const void *const addr) {
		inputPoints_ = static_cast<const float*>(addr);
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
	inline SimplexId getNumberOfTriangles() const {
		return numOutputTriangles_;
	}

	template<typename dataType>
	int execute() {

		Timer t;

		SimplexId vertexNumber = inputTriangl_->getNumberOfVertices();
		subdiviseTriangulation<dataType>();
		buildOutputTriangulation();

		{
			std::stringstream msg;
			msg << MODULE_S "Data-set (" << vertexNumber
					<< " points) processed in " << t.getElapsedTime() << " s. ("
					<< threadNumber_ << " thread(s))." << std::endl;
			dMsg(std::cout, msg.str(), timeMsg);
		}

		return 0;
	}

	void setTensorDataArray(void *tensorData) {
		tensorData_ = tensorData;
	}

private:
	template<typename dataType>
	int subdiviseTriangulation() {

		// not implemented for dimension >= 3
		if (inputTriangl_->getDimensionality() >= 3) {
			std::stringstream msg;
			msg << MODULE_S
			"Not yet implemented for dimension 3 and above" << std::endl;
			dMsg(std::cout, msg.str(), infoMsg);
			return 1;
		}

		dataType *tensorData = (dataType*) tensorData_;
		std::vector<bool> edgeHasMinima(numInputEdges_, false);
		std::vector<float> edgeBarycenters(numInputEdges_, -1.0);
		std::vector<float> edgeMinValues(numInputEdges_, -1.0);
		std::vector<bool> triangleHasMinima(numInputTriangles_, false);
		std::vector<float> triangleBarycenters(2 * numInputTriangles_, 0.0);

		numDividedEdges_ = 0;
		for (SimplexId i = 0; i < numInputEdges_; ++i) {
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

			float A = (v2_E - v1_E - v2_G + v1_G) * (v2_E - v1_E - v2_G + v1_G)
					+ 4 * (v2_F - v1_F) * (v2_F - v1_F);
			float B = 2 * (v2_E - v1_E - v2_G + v1_G) * (v1_E - v1_G)
					+ 8 * v1_F * (v2_F - v1_F);
			float C = (v1_E - v1_G) * (v1_E - v1_G) + 4 * v1_F * v1_F;
			if (A != 0) {
				float t = -B / (2 * A);
				edgeBarycenters[i] = t;
				if (t > 0 && t < 1) {
					edgeHasMinima[i] = true;
					numDividedEdges_++;
					edgeMinValues[i] = A * t * t + B * t + C;
				}
			}
		}
		std::cout << MODULE_S "Processed Edges, found " << numDividedEdges_
				<< " with minima." << std::endl;

		numDividedTriangles_ = 0;
		for (SimplexId i = 0; i < numInputTriangles_; ++i) {
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
				triangleHasMinima[i] = hasCrit;
				if (hasCrit) {
					numDividedTriangles_++;
				}
			}
		}
		std::cout << MODULE_S "Processed Triangles, found "
				<< numDividedTriangles_ << " with minima." << std::endl;

		points_.clear();
		pointId_.clear();
		pointDim_.clear();
		anisotropy_.clear();
		determinant_.clear();
		trace_.clear();
		lambda1_.clear();
		lambda2_.clear();
		numOutputVertices_ = numInputVertices_ + numDividedEdges_
				+ numDividedTriangles_;
		points_.resize(3 * numOutputVertices_);
		pointId_.resize(numOutputVertices_);
		pointDim_.resize(numOutputVertices_);
		anisotropy_.resize(numOutputVertices_);
		determinant_.resize(numOutputVertices_);
		trace_.resize(numOutputVertices_);
		lambda1_.resize(numOutputVertices_);
		lambda2_.resize(numOutputVertices_);

		std::copy(inputPoints_, inputPoints_ + numInputVertices_ * 3,
				points_.begin());

		for (SimplexId i = 0; i < numInputVertices_; ++i) {
			pointId_[i] = i;
			pointDim_[i] = 0;

			dataType E = tensorData[9 * i];
			dataType F = tensorData[9 * i + 1];
			dataType G = tensorData[9 * i + 4];
			float anisotropy = (E - G) * (E - G) + 4 * F * F;
			anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
			float trace = E + G;
			float determinant = E * G - (F * F);
			float discriminant = (trace * trace - 4 * determinant) / 2;
			discriminant = (discriminant < 0) ? 0 : sqrt(discriminant);
			anisotropy_[i] = anisotropy;
			trace_[i] = trace;
			determinant_[i] = determinant;
			lambda1_[i] = (trace + anisotropy) / 2;
			lambda2_[i] = (trace - anisotropy) / 2;
		}
		std::cout << MODULE_S "Added input vertices to output." << std::endl;

		edgeBarycenters_.clear();
		edgeBarycenters_.resize(numDividedEdges_);
		SimplexId index = numInputVertices_;
		std::vector < SimplexId > edgeIndexMap(numInputEdges_);
		for (SimplexId i = 0; i < numInputEdges_; ++i) {
			if (edgeHasMinima[i]) {
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

				float anisotropy = (E - G) * (E - G) + 4 * F * F;
				anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
				float trace = E + G;
				float determinant = E * G - (F * F);
				float discriminant = (trace * trace - 4 * determinant) / 2;
				discriminant = (discriminant < 0) ? 0 : sqrt(discriminant);
				anisotropy_[index] = anisotropy;
				trace_[index] = trace;
				determinant_[index] = determinant;
				lambda1_[index] = (trace + anisotropy) / 2;
				lambda2_[index] = (trace - anisotropy) / 2;

				pointId_[index] = index;
				pointDim_[index] = 1;

				points_[3 * index + 0] = t * v2_x + (1 - t) * v1_x;
				points_[3 * index + 1] = t * v2_y + (1 - t) * v1_y;
				points_[3 * index + 2] = t * v2_z + (1 - t) * v1_z;
				edgeBarycenters_[index - numInputVertices_] = t;

				edgeIndexMap[i] = index;
				++index;
			} else {
				edgeIndexMap[i] = -1;
			}
		}
		std::cout << MODULE_S "Added divided edges to output." << std::endl;

		triangleBarycenters_.clear();
		triangleBarycenters_.resize(2 * numDividedTriangles_);
		std::vector < SimplexId > triangleIndexMap(numInputTriangles_);
		for (SimplexId i = 0; i < numInputTriangles_; ++i) {
			if (triangleHasMinima[i]) {
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

				float anisotropy = (E - G) * (E - G) + 4 * F * F;
				anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
				float trace = E + G;
				float determinant = E * G - (F * F);
				anisotropy_[index] = anisotropy;
				trace_[index] = trace;
				determinant_[index] = determinant;
				lambda1_[index] = (trace + anisotropy) / 2;
				lambda2_[index] = (trace - anisotropy) / 2;
				SimplexId temp = 2
						* (index - numInputVertices_ - numDividedEdges_);
				triangleBarycenters_[temp] = alpha;
				triangleBarycenters_[temp + 1] = beta;

				points_[3 * index + 0] = alpha * v1_x + beta * v2_x
						+ gamma * v3_x;
				points_[3 * index + 1] = alpha * v1_y + beta * v2_y
						+ gamma * v3_y;
				points_[3 * index + 2] = alpha * v1_z + beta * v2_z
						+ gamma * v3_z;

				pointId_[index] = index;
				pointDim_[index] = 2;
				triangleIndexMap[i] = index;
				++index;
			} else {
				triangleIndexMap[i] = -1;
			}
		}
		std::cout << MODULE_S "Added divided triangles to output." << std::endl;

		cells_.clear();
		cells_.reserve(numInputTriangles_ * 4 * 6);
		SimplexId cellIndex = 0;
		for (SimplexId i = 0; i < numInputTriangles_; ++i) {
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
			if (triangleHasMinima[i]) {
				v123 = triangleIndexMap[i];
				if (edgeHasMinima[e1]) {
					v12 = edgeIndexMap[e1];
					addOutputTriangle(cellIndex, v123, v1, v12);
					addOutputTriangle(cellIndex, v123, v12, v2);
				} else {
					addOutputTriangle(cellIndex, v123, v1, v2);
				}
				if (edgeHasMinima[e2]) {
					v23 = edgeIndexMap[e2];
					addOutputTriangle(cellIndex, v123, v2, v23);
					addOutputTriangle(cellIndex, v123, v23, v3);
				} else {
					addOutputTriangle(cellIndex, v123, v2, v3);
				}
				if (edgeHasMinima[e3]) {
					v13 = edgeIndexMap[e3];
					addOutputTriangle(cellIndex, v123, v3, v13);
					addOutputTriangle(cellIndex, v123, v13, v1);
				} else {
					addOutputTriangle(cellIndex, v123, v3, v1);
				}
			} else {
				int numCrits = 0;
				if (edgeHasMinima[e1]) {
					v12 = edgeIndexMap[e1];
					numCrits++;
				}
				if (edgeHasMinima[e2]) {
					v23 = edgeIndexMap[e2];
					numCrits++;
				}
				if (edgeHasMinima[e3]) {
					v13 = edgeIndexMap[e3];
					numCrits++;
				}
				switch (numCrits) {
				case 0:
					addOutputTriangle(cellIndex, v1, v2, v3);
					break;
				case 1:
					if (edgeHasMinima[e1]) {
						addOutputTriangle(cellIndex, v12, v3, v1);
						addOutputTriangle(cellIndex, v12, v3, v2);
					} else if (edgeHasMinima[e2]) {
						addOutputTriangle(cellIndex, v23, v1, v2);
						addOutputTriangle(cellIndex, v23, v1, v3);
					} else if (edgeHasMinima[e3]) {
						addOutputTriangle(cellIndex, v13, v2, v3);
						addOutputTriangle(cellIndex, v13, v2, v1);
					}
					break;
				case 2:
					if (!edgeHasMinima[e1]) {
						addOutputTriangle(cellIndex, v23, v13, v3);
						if (edgeMinValues[e2] < edgeMinValues[e3]) {
							addOutputTriangle(cellIndex, v23, v13, v1);
							addOutputTriangle(cellIndex, v23, v2, v1);
						} else {
							addOutputTriangle(cellIndex, v13, v23, v2);
							addOutputTriangle(cellIndex, v13, v1, v2);
						}
					} else if (!edgeHasMinima[e2]) {
						addOutputTriangle(cellIndex, v12, v13, v1);
						if (edgeMinValues[e1] < edgeMinValues[e3]) {
							addOutputTriangle(cellIndex, v12, v13, v3);
							addOutputTriangle(cellIndex, v12, v2, v3);
						} else {
							addOutputTriangle(cellIndex, v13, v12, v2);
							addOutputTriangle(cellIndex, v13, v3, v2);
						}
					} else if (!edgeHasMinima[e3]) {
						addOutputTriangle(cellIndex, v12, v23, v2);
						if (edgeMinValues[e1] < edgeMinValues[e2]) {
							addOutputTriangle(cellIndex, v12, v23, v3);
							addOutputTriangle(cellIndex, v12, v1, v3);
						} else {
							addOutputTriangle(cellIndex, v23, v12, v1);
							addOutputTriangle(cellIndex, v23, v3, v1);
						}
					}
					break;
				case 3:
					float val12 = edgeMinValues[e1];
					float val23 = edgeMinValues[e2];
					float val13 = edgeMinValues[e3];
					if (val12 < val23 && val12 < val13) {
						addOutputTriangle(cellIndex, v12, v1, v13);
						addOutputTriangle(cellIndex, v12, v13, v3);
						addOutputTriangle(cellIndex, v12, v3, v23);
						addOutputTriangle(cellIndex, v12, v23, v2);
					} else if (val23 < val12 && val23 < val13) {
						addOutputTriangle(cellIndex, v23, v2, v12);
						addOutputTriangle(cellIndex, v23, v12, v1);
						addOutputTriangle(cellIndex, v23, v1, v13);
						addOutputTriangle(cellIndex, v23, v13, v3);
					} else {
						addOutputTriangle(cellIndex, v13, v3, v23);
						addOutputTriangle(cellIndex, v13, v23, v2);
						addOutputTriangle(cellIndex, v13, v2, v12);
						addOutputTriangle(cellIndex, v13, v12, v1);
					}
					break;
				}
			}
		}
		cells_.resize(cellIndex);
		numOutputTriangles_ = cellIndex / 4;
		std::cout << MODULE_S "Added all the new traingles to the output."
				<< std::endl;
		return 0;
	}

	int buildOutputTriangulation();
	void addOutputTriangle(int &offset, const SimplexId &v1,
			const SimplexId &v2, const SimplexId &v3);
	float intersectSegments(const float &x1, const float &y1, const float &x2,
			const float &y2, const float &x3, const float &y3, const float &x4,
			const float &y4);

	SimplexId numInputVertices_ { };
	SimplexId numInputEdges_ { };
	SimplexId numInputTriangles_ { };

	SimplexId numDividedEdges_ { };
	SimplexId numDividedTriangles_ { };

	SimplexId numOutputVertices_ { };
	SimplexId numOutputTriangles_ { };

	Triangulation *inputTriangl_ { };
	const float *inputPoints_ { };
	void *tensorData_ { };

	std::vector<float> &points_;
	std::vector<LongSimplexId> &cells_;
	std::vector<SimplexId> &pointId_;
	std::vector<SimplexId> &pointDim_;
	std::vector<float> &anisotropy_;
	std::vector<float> &determinant_;
	std::vector<float> &trace_;
	std::vector<float> &lambda1_;
	std::vector<float> &lambda2_;

	std::vector<float> edgeBarycenters_ { };
	std::vector<float> triangleBarycenters_ { };

	// output triangulation built on output points & output cells
	Triangulation *outputTriangl_ { };
};
} // namespace ttk
