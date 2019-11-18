#include "TensorMonotonousSubdivision.h"
#include "Geometry.h"
#include <array>

void ttk::TensorMonotonousSubdivision::addOutputTriangle(int &offset,
		const SimplexId &v1, const SimplexId &v2, const SimplexId &v3) {
	cells_.push_back(3);
	cells_.push_back(v1);
	cells_.push_back(v2);
	cells_.push_back(v3);
	offset += 4;
}

float ttk::TensorMonotonousSubdivision::intersectSegments(const float &x1,
		const float &y1, const float &x2, const float &y2, const float &x3,
		const float &y3, const float &x4, const float &y4) {
	return (((x1 - x3) * (y3 - y4)) - ((y1 - y3) * (x3 - x4)))
			/ (((x1 - x2) * (y3 - y4)) - ((y1 - y2) * (x3 - x4)));
}

int ttk::TensorMonotonousSubdivision::buildOutputTriangulation() {
// ensure subdivision is already performed
	if (points_.empty() || cells_.empty()) {
		return 1;
	}

// ensure output triangulation allocated by caller
	if (outputTriangl_ == nullptr) {
		return 2;
	}

	outputTriangl_->setInputPoints(points_.size() / 3, points_.data());
	outputTriangl_->setInputCells(cells_.size() / 4, cells_.data());

	return 0;
}
