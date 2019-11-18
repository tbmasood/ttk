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
			std::vector<SimplexId> &pointDim) :
			points_ { points }, cells_ { cells }, pointId_ { pointId }, pointDim_ {
					pointDim } {
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

		SimplexId maxPoints = numInputVertices_ + numInputEdges_
				+ numInputTriangles_;
		points_.clear();
		points_.reserve(maxPoints * 3);

		cells_.clear();
		cells_.reserve(numInputTriangles_ * 4 * 6);

		pointId_.clear();
		pointId_.reserve(maxPoints);

		pointDim_.clear();
		pointDim_.reserve(maxPoints);

		//std::copy(inputPoints_, inputPoints_ + numInputVertices_ * 3,
		//		points_.begin());

		for (SimplexId i = 0; i < numInputVertices_; ++i) {
			points_.push_back(inputPoints_[3 * i + 0]);
			points_.push_back(inputPoints_[3 * i + 1]);
			points_.push_back(inputPoints_[3 * i + 2]);
			pointId_.push_back(i);
			pointDim_.push_back(0);
		}

		dataType *tensorData = (dataType*) tensorData_;
		std::vector<bool> edgeHasMinima(numInputEdges_, false);
		std::vector<float> edgeCenters(3 * numInputEdges_, 0.0);
		std::vector<float> edgeMinValues(numInputEdges_, -1.0);
		std::vector<bool> triangleHasMinima(numInputTriangles_, false);
		std::vector<float> triangleH(numInputTriangles_, -1.0);
		std::vector<float> triangleCenters(3 * numInputTriangles_, 0.0);

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

			dataType v1_E = tensorData[9 * v1];
			dataType v1_F = tensorData[9 * v1 + 1];
			dataType v1_G = tensorData[9 * v1 + 4];
			float v1_x, v1_y, v1_z;
			inputTriangl_->getVertexPoint(v1, v1_x, v1_y, v1_z);

			dataType v2_E = tensorData[9 * v2];
			dataType v2_F = tensorData[9 * v2 + 1];
			dataType v2_G = tensorData[9 * v2 + 4];
			float v2_x, v2_y, v2_z;
			inputTriangl_->getVertexPoint(v2, v2_x, v2_y, v2_z);

			dataType v3_E = tensorData[9 * v3];
			dataType v3_F = tensorData[9 * v3 + 1];
			dataType v3_G = tensorData[9 * v3 + 4];
			float v3_x, v3_y, v3_z;
			inputTriangl_->getVertexPoint(v3, v3_x, v3_y, v3_z);

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
			float F = (Ec - Gc) * (Ec - Gc) + 4 * Fc * Fc;

			float H = 4 * A * C - B * B;
			triangleH[i] = H;

			bool hasCrit = false;
			float critX, critY, critZ;
			if (H != 0) {
				float alpha = (-2 * C * D + B * E) / H;
				float beta = (-2 * A * E + B * D) / H;
				float gamma = 1 - alpha - beta;
				hasCrit = alpha > 0 && beta > 0 && gamma > 0;
				critX = alpha * v1_x + beta * v2_x + gamma * v3_x;
				critY = alpha * v1_y + beta * v2_y + gamma * v3_y;
				critZ = alpha * v1_z + beta * v2_z + gamma * v3_z;
				triangleHasMinima[i] = hasCrit;
				if (hasCrit) {
					triangleCenters[3 * i + 0] = critX;
					triangleCenters[3 * i + 1] = critY;
					triangleCenters[3 * i + 2] = critZ;
				}

				float t = (2 * A - B + D - E) / (2 * (A - B + C));
				if (t > 0 && t < 1) {
					edgeHasMinima[e1] = true;
					alpha = 1 - t;
					beta = t;
					gamma = 0;
					edgeCenters[3 * e1 + 0] = alpha * v1_x + beta * v2_x;
					edgeCenters[3 * e1 + 1] = alpha * v1_y + beta * v2_y;
					edgeCenters[3 * e1 + 2] = alpha * v1_z + beta * v2_z;
					edgeMinValues[e1] = A * alpha * alpha + B * alpha * beta
							+ C * beta * beta + D * alpha + E * beta + F;
				}

				t = -E / (2 * C);
				if (t > 0 && t < 1) {
					edgeHasMinima[e2] = true;
					alpha = 0;
					beta = t;
					gamma = 1 - t;
					edgeCenters[3 * e2 + 0] = beta * v2_x + gamma * v3_x;
					edgeCenters[3 * e2 + 1] = beta * v2_y + gamma * v3_y;
					edgeCenters[3 * e2 + 2] = beta * v2_z + gamma * v3_z;
					edgeMinValues[e2] = C * beta * beta + E * beta + F;
				}

				t = -D / (2 * A);
				if (t > 0 && t < 1) {
					edgeHasMinima[e3] = true;
					alpha = t;
					beta = 0;
					gamma = 1 - t;
					edgeCenters[3 * e3 + 0] = alpha * v1_x + gamma * v3_x;
					edgeCenters[3 * e3 + 1] = alpha * v1_y + gamma * v3_y;
					edgeCenters[3 * e3 + 2] = alpha * v1_z + gamma * v3_z;
					edgeMinValues[e3] = A * alpha * alpha + D * alpha + F;
				}
			} else {
				float symmX1 = -D / (2 * A);
				float symmY1 = 0;
				float symmX2 = 0;
				float symmY2 = -D / B;
				float symmIntersect12 = intersectSegments(1, 0, 0, 1, symmX1,
						symmY1, symmX2, symmY2);
				float symmIntersect23 = intersectSegments(0, 1, 0, 0, symmX1,
						symmY1, symmX2, symmY2);
				float symmIntersect13 = intersectSegments(1, 0, 0, 0, symmX1,
						symmY1, symmX2, symmY2);
				bool intersects12 = symmIntersect12 >= 0
						&& symmIntersect12 <= 1;
				bool intersects23 = symmIntersect23 >= 0
						&& symmIntersect23 <= 1;
				bool intersects13 = symmIntersect13 >= 0
						&& symmIntersect13 <= 1;
				if (intersects12) {
					edgeHasMinima[e1] = true;
					float alpha = 1 - symmIntersect12;
					float beta = symmIntersect12;
					edgeCenters[3 * e1 + 0] = alpha * v1_x + beta * v2_x;
					edgeCenters[3 * e1 + 1] = alpha * v1_y + beta * v2_y;
					edgeCenters[3 * e1 + 2] = alpha * v1_z + beta * v2_z;
					edgeMinValues[e1] = A * alpha * alpha + B * alpha * beta
							+ C * beta * beta + D * alpha + E * beta + F;
				}
				if (intersects13) {
					edgeHasMinima[e3] = true;
					float alpha = 1 - symmIntersect13;
					float gamma = symmIntersect13;
					edgeCenters[3 * e3 + 0] = alpha * v1_x + gamma * v3_x;
					edgeCenters[3 * e3 + 1] = alpha * v1_y + gamma * v3_y;
					edgeCenters[3 * e3 + 2] = alpha * v1_z + gamma * v3_z;
					edgeMinValues[e3] = A * alpha * alpha + D * alpha + F;
				}
				if (intersects23) {
					edgeHasMinima[e2] = true;
					float beta = 1 - symmIntersect23;
					float gamma = symmIntersect23;
					edgeCenters[3 * e2 + 0] = beta * v2_x + gamma * v3_x;
					edgeCenters[3 * e2 + 1] = beta * v2_y + gamma * v3_y;
					edgeCenters[3 * e2 + 2] = beta * v2_z + gamma * v3_z;
					edgeMinValues[e2] = C * beta * beta + E * beta + F;
				}
			}
		}

		SimplexId index = numInputVertices_;
		std::vector < SimplexId > edgeIndexMap(numInputEdges_);
		for (SimplexId i = 0; i < numInputEdges_; ++i) {
			if (edgeHasMinima[i]) {
				points_.push_back(edgeCenters[3 * i + 0]);
				points_.push_back(edgeCenters[3 * i + 1]);
				points_.push_back(edgeCenters[3 * i + 2]);
				pointId_.push_back(index);
				pointDim_.push_back(1);
				edgeIndexMap[i] = index;
				++index;
			} else {
				edgeIndexMap[i] = -1;
			}
		}

		std::vector < SimplexId > triangleIndexMap(numInputTriangles_);
		for (SimplexId i = 0; i < numInputTriangles_; ++i) {
			if (triangleHasMinima[i]) {
				points_.push_back(triangleCenters[3 * i + 0]);
				points_.push_back(triangleCenters[3 * i + 1]);
				points_.push_back(triangleCenters[3 * i + 2]);
				pointId_.push_back(index);
				pointDim_.push_back(2);
				triangleIndexMap[i] = index;
				++index;
			} else {
				triangleIndexMap[i] = -1;
			}
		}
		points_.resize(3 * index);
		pointId_.resize(index);
		pointDim_.resize(index);
		numOutputVertices_ = index;

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
			if (triangleH[i] == 0) {
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
						addOutputTriangle(cellIndex, v3, v13, v23);
						addOutputTriangle(cellIndex, v13, v23, v1);
						addOutputTriangle(cellIndex, v23, v1, v2);
					} else if (!edgeHasMinima[e2]) {
						addOutputTriangle(cellIndex, v1, v12, v13);
						addOutputTriangle(cellIndex, v12, v13, v2);
						addOutputTriangle(cellIndex, v13, v2, v3);
					} else if (!edgeHasMinima[e3]) {
						addOutputTriangle(cellIndex, v2, v12, v23);
						addOutputTriangle(cellIndex, v12, v23, v1);
						addOutputTriangle(cellIndex, v23, v1, v3);
					}
					break;
				case 3:
					addOutputTriangle(cellIndex, v12, v23, v13);
					addOutputTriangle(cellIndex, v12, v13, v1);
					addOutputTriangle(cellIndex, v12, v23, v1);
					addOutputTriangle(cellIndex, v23, v13, v3);
					break;
				}
			} else if (triangleHasMinima[i]) {
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
		return 0;
	}

	int buildOutputTriangulation();
	void addOutputTriangle(int &offset, const SimplexId &v1,
			const SimplexId &v2, const SimplexId &v3);
	float intersectSegments(const float &x1, const float &y1, const float &x2,
			const float &y2, const float &x3, const float &y3, const float &x4,
			const float &y4);

	// input triangulation properties
	SimplexId numInputVertices_ { };
	SimplexId numInputEdges_ { };
	SimplexId numInputTriangles_ { };

	SimplexId numOutputVertices_ { };
	SimplexId numOutputTriangles_ { };

	// input triangulation
	Triangulation *inputTriangl_ { };
	// array of input points coordinates
	const float *inputPoints_ { };
	// list of input point data
	void *tensorData_ { };

	// output 3D coordinates of generated points: old points first, then edge
	// middles, then triangle barycenters
	std::vector<float> &points_;
	// output triangles
	std::vector<LongSimplexId> &cells_;
	// generated point cell id
	std::vector<SimplexId> &pointId_;
	// generated points dimension: 0 vertex of parent triangulation, 1 edge
	// middle, 2 triangle barycenter
	std::vector<SimplexId> &pointDim_;

	std::vector<float> anisotropy_;
	std::vector<float> determinant_;
	std::vector<float> trace_;
	std::vector<float> lambda1_;
	std::vector<float> lambda2_;

	// output triangulation built on output points & output cells
	Triangulation *outputTriangl_ { };
};
} // namespace ttk
