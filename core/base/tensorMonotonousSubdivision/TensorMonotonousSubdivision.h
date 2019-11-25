/// \ingroup base
/// \class ttk::TensorMonotonousSubdivision
/// \author Talha Bin Masood (talha.bin.masood@liu.se)
/// \date November 2019
///
/// \brief Subdivide a triangulation according to tensor invariant.
///
/// Under linear interpolation of tensor components, tensor invariants
/// like anisotropy and determinant are quadratic. To compute the contour
/// tree correctly for these invariants, it is thus important to subdivide
/// the input mesh into set of triangles where these invariants have
/// monotonous behaviour.
///
/// We assume 2D symmteric tensors are provided at the mesh vertices:
/// ( e f )
/// ( f g )
///
/// This class allows subdivison of 2D tensor meshes into monotonous
/// triangle meshes based on two invariants:
/// 1. Anisotropy: defined as the difference between Eigen values,
///    or sqrt((e-g)^2 + 4*f^2)
/// 2. Determinant: the determinant of the tensor (e*g - f^2)
///
/// Data on vertices (point data) with continuous values
/// (float/double) can be interpolated on the new
/// triangulation. Scalar data on input triangles can be replicated on
/// the new triangles.
///
/// \b Related \b publications \n
/// "Topology-guided Tessellation of Quadratic Elements" \n
/// Scott Dillard, Vijay Natarajan, Gunther Weber, Valerio Pascucci and Bernd
/// Hamann\n J. Computational Geometry and Applications, 19(2), 2009 \n
///
/// \sa ttk::Triangulation

#pragma once

// base code includes
#include <Triangulation.h>
#include <Wrapper.h>

#define MODULE_TMS "[TensorMonotonousSubdivision] "

namespace ttk {

  class TensorMonotonousSubdivision : public Debug {

  public:
    TensorMonotonousSubdivision(std::vector<float> &points,
                                std::vector<LongSimplexId> &cells,
                                std::vector<SimplexId> &pointDim,
                                std::vector<float> &anisotropy,
                                std::vector<float> &determinant,
                                std::vector<float> &trace,
                                std::vector<float> &lambda1,
                                std::vector<float> &lambda2)
      : points_{points}, cells_{cells}, pointDim_{pointDim},
        anisotropy_{anisotropy}, determinant_{determinant}, trace_{trace},
        lambda1_{lambda1}, lambda2_{lambda2} {
    }

    inline void setOutputTriangulation(Triangulation *const triangulation) {
      outputTriangl_ = triangulation;
    }
    inline void setupTriangulation(Triangulation *const triangulation) {
      inputTriangl_ = triangulation;
      if(inputTriangl_ != nullptr) {
        inputTriangl_->preprocessVertexNeighbors();
        inputTriangl_->preprocessEdges();
        inputTriangl_->preprocessTriangles();
        inputTriangl_->preprocessTriangleEdges();
      }
      numInputVertices_ = inputTriangl_->getNumberOfVertices();
      numInputEdges_ = inputTriangl_->getNumberOfEdges();
      numInputTriangles_ = inputTriangl_->getNumberOfTriangles();
    }

    /**
     * @brief Return the number of vertices in the output triangulation
     */
    inline SimplexId getNumberOfVertices() const {
      return numOutputVertices_;
    }

    /**
     * Return the number of triangles in the output triangulation
     */
    inline LongSimplexId getNumberOfTriangles() const {
      return numOutputTriangles_;
    }

    template <typename dataType>
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

    /**
     * @brief Set the tensor field for the input points. Although each tensor is
     * given as a 3x3 tensor in this array, but during computation we assume 2x2
     * symmteric tensor is given for each point. So:
     *
     * Input tensor ( a b c )  is treated as ( a b 0 )
     *              ( d e f )                ( b e 0 )
     *              ( g h i )                ( 0 0 0 )
     */
    inline void setTensorDataArray(void *tensorData) {
      tensorData_ = tensorData;
    }

    /**
     * @brief Sets the criteria for mesh subdivision.
     *
     * @param[subdivisionField] If set to 0, then subdivision is done based on
     * anisotropy. If set to 1, then subdivision is done based on determinant.
     */
    inline void setSubdivisionField(const unsigned int subdivisionField) {
      subdivisionField_ = subdivisionField;
    }

    /**
     * @brief A flag to set if tensor anisotropy should be computed for each
     * vertex in the output mesh.
     */
    inline void setGenerateAnisotropyField(const bool generateAnisotropyField) {
      generateAnisotropyField_ = generateAnisotropyField;
    }

    /**
     * @brief A flag to set if tensor determinant should be computed for each
     * vertex in the output mesh.
     */
    inline void
      setGenerateDeterminantField(const bool generateDeterminantField) {
      generateDeterminantField_ = generateDeterminantField;
    }

    /**
     * @brief A flag to set if tensor trace should be computed for each
     * vertex in the output mesh.
     */
    inline void setGenerateTraceField(const bool generateTraceField) {
      generateTraceField_ = generateTraceField;
    }

    /**
     * @brief A flag to set if tensor eigen values should be computed for each
     * vertex in the output mesh.
     */
    inline void
      setGenerateEigenValuesField(const bool generateEigenValuesField) {
      generateEigenValuesField_ = generateEigenValuesField;
    }

    /**
     * @brief Interpolate floating-point point data on subdivised triangulation.
     *
     * Copy values on parent vertices, interpolate on edges and triangle centers
     *
     * @param[in] data Pointer to input data on parent triangulation
     * @param[out] output Allocated buffer to be filled
     *
     * @return 0 in case of success
     */
    template <typename T>
    int interpolateContinuousField(const T *data,
                                   T *output,
                                   int numComponents) const {
      static_assert(
        std::is_floating_point<T>::value, "Floating point type required.");
      if(inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
        return 1;
      }
      const auto nOutVerts = this->getNumberOfVertices();
      if(nOutVerts < 0 || nOutVerts < numInputVertices_) {
        return 1;
      }

      // copy data on parent vertices
      std::copy(data, data + numInputVertices_ * numComponents, output);

      // interpolate on edges
      for(SimplexId i = 0; i < numDividedEdges_; ++i) {
        SimplexId a{}, b{};
        inputTriangl_->getEdgeVertex(pointEdgeMap_[i], 0, a);
        inputTriangl_->getEdgeVertex(pointEdgeMap_[i], 1, b);
        float t = edgeBarycenters_[i];
        for(int j = 0; j < numComponents; ++j) {
          output[numComponents * (numInputVertices_ + i) + j]
            = (1 - t) * data[numComponents * a + j]
              + t * data[numComponents * b + j];
        }
      }

      // interpolate on triangle barycenters
      for(SimplexId i = 0; i < numDividedTriangles_; ++i) {
        SimplexId a{}, b{}, c{};
        inputTriangl_->getTriangleVertex(pointTriangleMap_[i], 0, a);
        inputTriangl_->getTriangleVertex(pointTriangleMap_[i], 1, b);
        inputTriangl_->getTriangleVertex(pointTriangleMap_[i], 2, c);
        float alpha = triangleBarycenters_[2 * i + 0];
        float beta = triangleBarycenters_[2 * i + 1];
        float gamma = 1 - alpha - beta;
        for(int j = 0; j < numComponents; ++j) {
          output[numComponents * (numInputVertices_ + numDividedEdges_ + i) + j]
            = alpha * data[numComponents * a + j]
              + beta * data[numComponents * b + j]
              + gamma * data[numComponents * c + j];
        }
      }
      return 0;
    }

    /**
     * @brief Interpolate integer point data on subdivided triangulation.
     *
     * Copy values on parent vertices, put 0 elsewhere.
     *
     * @param[in] data Pointer to input data on parent triangulation
     * @param[out] output Allocated buffer to be filled
     *
     * @return 0 in case of success
     */
    template <typename T>
    int interpolateDiscreteField(const T *data,
                                 T *output,
                                 int numComponents) const {
      static_assert(std::is_integral<T>::value, "Integral type required.");
      if(inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
        return 1;
      }
      const auto nOutVerts = this->getNumberOfVertices();
      if(nOutVerts < 0 || nOutVerts < numInputVertices_) {
        return 1;
      }
      std::fill(output, output + nOutVerts * numComponents, T{0});
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
    template <typename T>
    int interpolateCellDataField(const T *data,
                                 T *output,
                                 int numComponents) const {
      if(inputTriangl_ == nullptr || outputTriangl_ == nullptr) {
        return 1;
      }
      const auto nOutTriangles = this->getNumberOfTriangles();
      if(nOutTriangles < 0) {
        return 1;
      }
      for(LongSimplexId i = 0; i < nOutTriangles; ++i) {
        for(int j = 0; j < numComponents; ++j) {
          output[numComponents * i + j]
            = data[numComponents * originalCellMap_[i] + j];
        }
      }
      return 0;
    }

  private:
    /**
     * @brief The key method which generates the subdivision of a tensor mesh
     * into a mesh with monotonous invariant behavior.
     *
     * @return 0 in case of success
     */
    template <typename dataType>
    int subdivideTriangulation() {

      // not implemented for dimension >= 3
      if(inputTriangl_->getDimensionality() >= 3) {
        std::stringstream msg;
        msg << MODULE_TMS "Not implemented for dimension 3 and above."
            << std::endl;
        dMsg(std::cout, msg.str(), infoMsg);
        return 1;
      }

      if(tensorData_ == nullptr) {
        std::stringstream msg;
        msg << MODULE_TMS "Tensor data not available for this mesh."
            << std::endl;
        dMsg(std::cout, msg.str(), infoMsg);
        return 1;
      }

      dataType *tensorData = (dataType *)tensorData_;
      std::vector<bool> edgeHasCriticalPoint(numInputEdges_, false);
      std::vector<float> edgeBarycenters(numInputEdges_, -1.0);
      std::vector<float> edgeCriticalValues(numInputEdges_, -1.0);
      std::vector<bool> triangleHasCriticalPoint(numInputTriangles_, false);
      std::vector<float> triangleBarycenters(2 * numInputTriangles_, 0.0);

      for(SimplexId i = 0; i < numInputEdges_; ++i) {
        SimplexId v1, v2;
        inputTriangl_->getEdgeVertex(i, 0, v1);
        inputTriangl_->getEdgeVertex(i, 1, v2);

        dataType v1_E = tensorData[9 * v1];
        dataType v1_F = tensorData[9 * v1 + 1];
        dataType v1_G = tensorData[9 * v1 + 4];

        dataType v2_E = tensorData[9 * v2];
        dataType v2_F = tensorData[9 * v2 + 1];
        dataType v2_G = tensorData[9 * v2 + 4];

        float A = 0, B = 0, C = 0;
        if(subdivisionField_ == 0) {
          // anisotropy
          A = (v2_E - v1_E - v2_G + v1_G) * (v2_E - v1_E - v2_G + v1_G)
              + 4 * (v2_F - v1_F) * (v2_F - v1_F);
          B = 2 * (v2_E - v1_E - v2_G + v1_G) * (v1_E - v1_G)
              + 8 * v1_F * (v2_F - v1_F);
          C = (v1_E - v1_G) * (v1_E - v1_G) + 4 * v1_F * v1_F;
        } else {
          // determinant
          A = (v2_E - v1_E) * (v2_G - v1_G) - (v2_F - v1_F) * (v2_F - v1_F);
          B = v1_E * (v2_G - v1_G) + v1_G * (v2_E - v1_E)
              - 2 * v1_F * (v2_F - v1_F);
          C = v1_E * v1_G - v1_F * v1_F;
        }
        if(A != 0) {
          float t = -B / (2 * A);
          edgeBarycenters[i] = t;
          if(t > 0 && t < 1) {
            edgeHasCriticalPoint[i] = true;
            edgeCriticalValues[i] = A * t * t + B * t + C;
          }
        } else {
          edgeBarycenters[i] = 0.5;
        }
      }

      for(SimplexId i = 0; i < numInputTriangles_; ++i) {
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

        float A = 0, B = 0, C = 0, D = 0, E = 0;
        float Ex = v1_E - v3_E;
        float Fx = v1_F - v3_F;
        float Gx = v1_G - v3_G;

        float Ey = v2_E - v3_E;
        float Fy = v2_F - v3_F;
        float Gy = v2_G - v3_G;

        float Ec = v3_E;
        float Fc = v3_F;
        float Gc = v3_G;

        if(subdivisionField_ == 0) {
          // anisotropy
          A = (Ex - Gx) * (Ex - Gx) + 4 * Fx * Fx;
          B = 2 * (Ex - Gx) * (Ey - Gy) + 8 * Fx * Fy;
          C = (Ey - Gy) * (Ey - Gy) + 4 * Fy * Fy;
          D = 2 * (Ex - Gx) * (Ec - Gc) + 8 * Fx * Fc;
          E = 2 * (Ey - Gy) * (Ec - Gc) + 8 * Fy * Fc;
        } else {
          // determinant
          A = Ex * Gx - Fx * Fx;
          B = Ex * Gy + Ey * Gx - 2 * Fx * Fy;
          C = Ey * Gy - Fy * Fy;
          D = Ex * Gc + Ec * Gx - 2 * Fx * Fc;
          E = Ey * Gc + Ec * Gy - 2 * Fy * Fc;
        }

        float H = 4 * A * C - B * B;

        bool hasCrit = false;
        if(H != 0) {
          float alpha = (-2 * C * D + B * E) / H;
          float beta = (-2 * A * E + B * D) / H;
          float gamma = 1 - alpha - beta;
          triangleBarycenters[2 * i + 0] = alpha;
          triangleBarycenters[2 * i + 1] = beta;
          hasCrit = alpha > 0 && beta > 0 && gamma > 0;
        } else {
          triangleBarycenters[2 * i + 0] = 1 / 3.0;
          triangleBarycenters[2 * i + 1] = 1 / 3.0;
        }
        triangleHasCriticalPoint[i] = hasCrit;
      }

      std::vector<SimplexId> edgeIndexMap(numInputEdges_);
      SimplexId index = numInputVertices_;
      numDividedEdges_ = 0;
      for(SimplexId i = 0; i < numInputEdges_; ++i) {
        if(edgeHasCriticalPoint[i]) {
          ++numDividedEdges_;
          edgeIndexMap[i] = index;
          ++index;
        } else {
          edgeIndexMap[i] = -1;
        }
      }
      std::vector<SimplexId> triangleIndexMap(numInputTriangles_);
      numDividedTriangles_ = 0;
      for(SimplexId i = 0; i < numInputTriangles_; ++i) {
        if(triangleHasCriticalPoint[i]) {
          ++numDividedTriangles_;
          triangleIndexMap[i] = index;
          ++index;
        } else {
          triangleIndexMap[i] = -1;
        }
      }

      numOutputVertices_
        = numInputVertices_ + numDividedEdges_ + numDividedTriangles_;
      points_.clear();
      points_.resize(3 * numOutputVertices_);
      pointDim_.clear();
      pointDim_.resize(numOutputVertices_);
      anisotropy_.clear();
      determinant_.clear();
      trace_.clear();
      lambda1_.clear();
      lambda2_.clear();
      if(generateAnisotropyField_) {
        anisotropy_.resize(numOutputVertices_);
      }
      if(generateDeterminantField_) {
        determinant_.resize(numOutputVertices_);
      }
      if(generateTraceField_) {
        trace_.resize(numOutputVertices_);
      }
      if(generateEigenValuesField_) {
        lambda1_.resize(numOutputVertices_);
        lambda2_.resize(numOutputVertices_);
      }

      for(SimplexId i = 0; i < numInputVertices_; ++i) {
        pointDim_[i] = 0;
        dataType E = tensorData[9 * i];
        dataType F = tensorData[9 * i + 1];
        dataType G = tensorData[9 * i + 4];
        float anisotropy = 0;
        if(generateAnisotropyField_ || generateEigenValuesField_) {
          anisotropy = (E - G) * (E - G) + 4 * F * F;
          anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
        }
        if(generateAnisotropyField_) {
          anisotropy_[i] = anisotropy;
        }
        if(generateDeterminantField_) {
          determinant_[i] = E * G - (F * F);
        }
        if(generateTraceField_) {
          trace_[i] = E + G;
        }
        if(generateEigenValuesField_) {
          lambda1_[i] = (E + G + anisotropy) / 2;
          lambda2_[i] = (E + G - anisotropy) / 2;
        }

        float v_x, v_y, v_z;
        inputTriangl_->getVertexPoint(i, v_x, v_y, v_z);
        pointDim_[i] = 0;
        points_[3 * i + 0] = v_x;
        points_[3 * i + 1] = v_y;
        points_[3 * i + 2] = v_z;
      }

      edgeBarycenters_.clear();
      edgeBarycenters_.resize(numDividedEdges_);
      pointEdgeMap_.clear();
      pointEdgeMap_.resize(numDividedEdges_);

      for(SimplexId i = 0; i < numInputEdges_; ++i) {
        if(edgeHasCriticalPoint[i]) {
          SimplexId edgeIndex = edgeIndexMap[i];
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
          if(generateAnisotropyField_ || generateEigenValuesField_) {
            anisotropy = (E - G) * (E - G) + 4 * F * F;
            anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
          }
          if(generateAnisotropyField_) {
            anisotropy_[edgeIndex] = anisotropy;
          }
          if(generateDeterminantField_) {
            determinant_[edgeIndex] = E * G - (F * F);
          }
          if(generateTraceField_) {
            trace_[edgeIndex] = E + G;
          }
          if(generateEigenValuesField_) {
            lambda1_[edgeIndex] = (E + G + anisotropy) / 2;
            lambda2_[edgeIndex] = (E + G - anisotropy) / 2;
          }

          pointDim_[edgeIndex] = 1;
          points_[3 * edgeIndex + 0] = t * v2_x + (1 - t) * v1_x;
          points_[3 * edgeIndex + 1] = t * v2_y + (1 - t) * v1_y;
          points_[3 * edgeIndex + 2] = t * v2_z + (1 - t) * v1_z;
          edgeBarycenters_[edgeIndex - numInputVertices_] = t;
          pointEdgeMap_[edgeIndex - numInputVertices_] = i;
        }
      }

      triangleBarycenters_.clear();
      triangleBarycenters_.resize(2 * numDividedTriangles_);
      pointTriangleMap_.clear();
      pointTriangleMap_.resize(numDividedTriangles_);

      for(SimplexId i = 0; i < numInputTriangles_; ++i) {
        if(triangleHasCriticalPoint[i]) {
          SimplexId triangleIndex = triangleIndexMap[i];
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
          if(generateAnisotropyField_ || generateEigenValuesField_) {
            anisotropy = (E - G) * (E - G) + 4 * F * F;
            anisotropy = anisotropy < 0 ? 0 : sqrt(anisotropy);
          }
          if(generateAnisotropyField_) {
            anisotropy_[triangleIndex] = anisotropy;
          }
          if(generateDeterminantField_) {
            determinant_[triangleIndex] = E * G - (F * F);
          }
          if(generateTraceField_) {
            trace_[triangleIndex] = E + G;
          }
          if(generateEigenValuesField_) {
            lambda1_[triangleIndex] = (E + G + anisotropy) / 2;
            lambda2_[triangleIndex] = (E + G - anisotropy) / 2;
          }

          SimplexId temp = triangleIndex - numInputVertices_ - numDividedEdges_;
          triangleBarycenters_[2 * temp + 0] = alpha;
          triangleBarycenters_[2 * temp + 1] = beta;
          pointTriangleMap_[temp] = i;

          pointDim_[triangleIndex] = 2;
          points_[3 * triangleIndex + 0]
            = alpha * v1_x + beta * v2_x + gamma * v3_x;
          points_[3 * triangleIndex + 1]
            = alpha * v1_y + beta * v2_y + gamma * v3_y;
          points_[3 * triangleIndex + 2]
            = alpha * v1_z + beta * v2_z + gamma * v3_z;
        }
      }

      cells_.clear();
      originalCellMap_.clear();
      // Worst case, all the input triangles have a critical point and they need
      // to be divided into 6 triangles each.
      cells_.reserve(numInputTriangles_ * 4 * 6);
      originalCellMap_.reserve(numInputTriangles_ * 6);
      SimplexId cellIndex = 0;
      for(LongSimplexId i = 0; i < numInputTriangles_; ++i) {
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
        if(e2v1 == v3) {
          if(e2v2 == v1) {
            SimplexId temp = e2;
            e2 = e3;
            e3 = temp;
          }
        } else {
          if(e2v1 == v1) {
            SimplexId temp = e2;
            e2 = e3;
            e3 = temp;
          }
        }
        // After above steps, e1 is between v1 and v2,
        // e2 is between v2 and v3, and e3 is between v1 and v3

        SimplexId v12, v23, v13, v123;
        if(triangleHasCriticalPoint[i]) {
          v123 = triangleIndexMap[i];
          if(edgeHasCriticalPoint[e1]) {
            v12 = edgeIndexMap[e1];
            addOutputTriangle(cellIndex, v123, v1, v12, i);
            addOutputTriangle(cellIndex, v123, v12, v2, i);
          } else {
            addOutputTriangle(cellIndex, v123, v1, v2, i);
          }
          if(edgeHasCriticalPoint[e2]) {
            v23 = edgeIndexMap[e2];
            addOutputTriangle(cellIndex, v123, v2, v23, i);
            addOutputTriangle(cellIndex, v123, v23, v3, i);
          } else {
            addOutputTriangle(cellIndex, v123, v2, v3, i);
          }
          if(edgeHasCriticalPoint[e3]) {
            v13 = edgeIndexMap[e3];
            addOutputTriangle(cellIndex, v123, v3, v13, i);
            addOutputTriangle(cellIndex, v123, v13, v1, i);
          } else {
            addOutputTriangle(cellIndex, v123, v3, v1, i);
          }
        } else {
          int numCrits = 0;
          if(edgeHasCriticalPoint[e1]) {
            v12 = edgeIndexMap[e1];
            numCrits++;
          }
          if(edgeHasCriticalPoint[e2]) {
            v23 = edgeIndexMap[e2];
            numCrits++;
          }
          if(edgeHasCriticalPoint[e3]) {
            v13 = edgeIndexMap[e3];
            numCrits++;
          }
          switch(numCrits) {
            case 0:
              addOutputTriangle(cellIndex, v1, v2, v3, i);
              break;
            case 1:
              if(edgeHasCriticalPoint[e1]) {
                addOutputTriangle(cellIndex, v12, v3, v1, i);
                addOutputTriangle(cellIndex, v12, v3, v2, i);
              } else if(edgeHasCriticalPoint[e2]) {
                addOutputTriangle(cellIndex, v23, v1, v2, i);
                addOutputTriangle(cellIndex, v23, v1, v3, i);
              } else if(edgeHasCriticalPoint[e3]) {
                addOutputTriangle(cellIndex, v13, v2, v3, i);
                addOutputTriangle(cellIndex, v13, v2, v1, i);
              }
              break;
            case 2:
              // This case is correctly handled for anisotropy subdivision
              // TO DO: I am not sure if determinant subdivision is also correct
              if(!edgeHasCriticalPoint[e1]) {
                addOutputTriangle(cellIndex, v23, v13, v3, i);
                if(edgeCriticalValues[e2] < edgeCriticalValues[e3]) {
                  addOutputTriangle(cellIndex, v23, v13, v1, i);
                  addOutputTriangle(cellIndex, v23, v2, v1, i);
                } else {
                  addOutputTriangle(cellIndex, v13, v23, v2, i);
                  addOutputTriangle(cellIndex, v13, v1, v2, i);
                }
              } else if(!edgeHasCriticalPoint[e2]) {
                addOutputTriangle(cellIndex, v12, v13, v1, i);
                if(edgeCriticalValues[e1] < edgeCriticalValues[e3]) {
                  addOutputTriangle(cellIndex, v12, v13, v3, i);
                  addOutputTriangle(cellIndex, v12, v2, v3, i);
                } else {
                  addOutputTriangle(cellIndex, v13, v12, v2, i);
                  addOutputTriangle(cellIndex, v13, v3, v2, i);
                }
              } else if(!edgeHasCriticalPoint[e3]) {
                addOutputTriangle(cellIndex, v12, v23, v2, i);
                if(edgeCriticalValues[e1] < edgeCriticalValues[e2]) {
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
              if(val12 < val23 && val12 < val13) {
                addOutputTriangle(cellIndex, v12, v1, v13, i);
                addOutputTriangle(cellIndex, v12, v13, v3, i);
                addOutputTriangle(cellIndex, v12, v3, v23, i);
                addOutputTriangle(cellIndex, v12, v23, v2, i);
              } else if(val23 < val12 && val23 < val13) {
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
      // The size of cells_ is 4 times the number of generated output triangles.
      cells_.resize(cellIndex);
      numOutputTriangles_ = cellIndex / 4;
      originalCellMap_.resize(numOutputTriangles_);
      return 0;
    }

    /*
     * @brief Builds the output triangulation after the points_ and cells_ lists
     * have been already generated by the subdivideTriangulation()
     */
    int buildOutputTriangulation();

    /*
     * @brief Adds an output triangle in the cells_ list.
     */
    void addOutputTriangle(int &offset,
                           const SimplexId &v1,
                           const SimplexId &v2,
                           const SimplexId &v3,
                           const LongSimplexId &inputCellIndex);

    unsigned int subdivisionField_{0};
    bool generateAnisotropyField_{true};
    bool generateDeterminantField_{false};
    bool generateTraceField_{false};
    bool generateEigenValuesField_{false};

    /*
     * @brief Number of vertices in the input triangulation.
     */
    SimplexId numInputVertices_{};

    /*
     * @brief Number of edges in the input triangulation.
     */
    SimplexId numInputEdges_{};

    /*
     * @brief Number of triangles in the input triangulation.
     */
    SimplexId numInputTriangles_{};

    /*
     * @brief Number of edges which have been divided because of presence of
     * an edge critical point.
     */
    SimplexId numDividedEdges_{};

    /*
     * @brief Number of triangles which have been divided
     * because of presence of a triangle critical point.
     */
    SimplexId numDividedTriangles_{};

    /*
     * @brief Number of vertices in the output triangulation, it should be equal
     * to (numInputVertices_ + numDividedEdges_ + numDividedTriangles_)
     */
    SimplexId numOutputVertices_{};

    /*
     * @brief Number of triangles in the output.
     */
    SimplexId numOutputTriangles_{};

    /*
     * @brief Pointer to input triangulation.
     */
    Triangulation *inputTriangl_{};

    /*
     * @brief Pointer to input tensor data array.
     */
    void *tensorData_{};

    /*
     * @brief Points in the output triangulation, saved as
     * (x1, y1, z1, x2, y2, z2, ...)
     */
    std::vector<float> &points_;

    /*
     * @brief Triangles in the output triangulation, saved as
     * (3, t1v1, t1v2, t1v3, 3, t2v1, t2v2, t2v3, ...)
     */
    std::vector<LongSimplexId> &cells_;

    /*
     * @brief The list containing the type of the output point:
     *  0 - input vertex
     *  1 - edge critical point
     *  2 - tringle critical point
     */
    std::vector<SimplexId> &pointDim_;

    /*
     * @brief Tensor anisotropy value computed for each point of the output.
     */
    std::vector<float> &anisotropy_;

    /*
     * @brief Tensor determinant computed for each point of the output.
     */
    std::vector<float> &determinant_;

    /*
     * @brief Tensor trace computed for each point of the output.
     */
    std::vector<float> &trace_;

    /*
     * @brief The larger of the Eigen values of the tensor computed for each
     * point of the output.
     */
    std::vector<float> &lambda1_;

    /*
     * @brief The smaller of the Eigen values of the tensor computed for each
     * point of the output.
     */
    std::vector<float> &lambda2_;

    /*
     * @brief A value between 0 and 1 is saved for edges which have an edge
     * critical point, this value is used for interpolation of fields.
     */
    std::vector<float> edgeBarycenters_{};

    /*
     * @brief The barycentric coordinate of the critical point within the
     * triangle is saved which is used for interpolation.
     */
    std::vector<float> triangleBarycenters_{};

    /*
     * @brief For the edge critical points in the output, a mapping from point
     * index to edge index in the input triangulation.
     */
    std::vector<SimplexId> pointEdgeMap_{};

    /*
     * @brief For the triangle critical points in the output, a mapping from
     * point index to the triangle index in the input triangulation.
     */
    std::vector<LongSimplexId> pointTriangleMap_{};

    /*
     * @brief Map from traingle index in the output to the triangle index in the
     * input.
     */
    std::vector<LongSimplexId> originalCellMap_{};

    /*
     * @brief Output triangulation built on output points and output cells.
     */
    Triangulation *outputTriangl_{};
  };
} // namespace ttk
