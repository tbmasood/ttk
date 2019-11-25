#include "TensorMonotonousSubdivision.h"

void ttk::TensorMonotonousSubdivision::addOutputTriangle(
  int &offset,
  const SimplexId &v1,
  const SimplexId &v2,
  const SimplexId &v3,
  const LongSimplexId &inputCellIndex) {
  cells_.push_back(3);
  cells_.push_back(v1);
  cells_.push_back(v2);
  cells_.push_back(v3);
  originalCellMap_.push_back(inputCellIndex);
  offset += 4;
}

int ttk::TensorMonotonousSubdivision::buildOutputTriangulation() {
  // ensure subdivision is already performed
  if(points_.empty() || cells_.empty()) {
    return 1;
  }

  // ensure output triangulation allocated by caller
  if(outputTriangl_ == nullptr) {
    return 2;
  }

  outputTriangl_->setInputPoints(points_.size() / 3, points_.data());
  outputTriangl_->setInputCells(cells_.size() / 4, cells_.data());

  return 0;
}
