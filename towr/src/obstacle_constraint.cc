/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <towr/constraints/obstacle_constraint.h>


namespace towr {


ObstacleConstraint::ObstacleConstraint (const HeightMap::Ptr& terrain)
    :ConstraintSet(kSpecifyLater, "obstacle-")
{
  terrain_ = terrain;
}

void
ObstacleConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(ee_motion_id_);

  // skip first node, b/c already constrained by initial stance
  for (int id=1; id<ee_motion_->GetNodes().size(); ++id)
    node_ids_.push_back(id);

  int constraint_count = node_ids_.size();
  SetRows(constraint_count);
}

Eigen::VectorXd
ObstacleConstraint::GetValues () const
{
  VectorXd g(GetRows());

  auto nodes = ee_motion_->GetNodes();
  int row = 0;
  for (int id : node_ids_) {
    Vector3d p = nodes.at(id).p();
    auto obstacle_list = terrain_->GetObstacles();
    g(row++) = CalculateDistanceToNearestObstacle(p.x(), p.y());
  }

  return g;
}

double ObstacleConstraint::CalculateDistanceToNearestObstacle(double x, double y) const
{
  auto segments = terrain_->GetObstacles();
  double min_distance = std::numeric_limits<double>::max();
  for (const auto& segment : segments) {
    double distance = DistanceToLineSegment({x, y}, segment);
    min_distance = std::min(min_distance, distance);
  }
  return min_distance;
}

double ObstacleConstraint::DistanceToLineSegment(const Eigen::Vector2d& point, const std::pair<Eigen::Vector2d, Eigen::Vector2d>& segment) const
{
  Eigen::Vector2d p = point;
  Eigen::Vector2d v = segment.first;
  Eigen::Vector2d w = segment.second;

  // Vector from v to w
  Eigen::Vector2d vw = w - v;

  // Vector from v to p
  Eigen::Vector2d vp = p - v;

  // Calculate squared distances
  double vw_length_squared = vw.squaredNorm();
  if (vw_length_squared == 0.0) {
    // Segment is a point
    return (p - v).norm();
  }

  // Calculate the projection of point onto the line segment
  double projection = vp.dot(vw) / vw_length_squared;
  if (projection < 0.0) {
      // Before segment
      return (p - v).norm();
  } else if (projection > 1.0) {
      // After segment
      return (p - w).norm();
  }

  // Projection is within the segment, calculate perpendicular distance
  Eigen::Vector2d projection_point = v + projection * vw;
  return (p - projection_point).norm();
}

std::tuple<double, double, double> ObstacleConstraint::DistanceAndDerivativesToLineSegment(const Eigen::Vector2d& point, const std::pair<Eigen::Vector2d, Eigen::Vector2d>& segment) const {
    Eigen::Vector2d p = point;
    Eigen::Vector2d a = segment.first;
    Eigen::Vector2d b = segment.second;
    Eigen::Vector2d ab = b - a;
    Eigen::Vector2d ap = p - a;
    Eigen::Vector2d bp = p - b;

    double ab_length_squared = ab.squaredNorm();
    if (ab_length_squared == 0.0) {
        // Segment is a point
        double p_x = (p - a).normalized().x();
        double p_y = (p - a).normalized().y();
        return {ap.norm(), p_x, p_y};
    }

    double t = ap.dot(ab) / ab_length_squared;
    Eigen::Vector2d projection_point;
    if (t < 0.0) {
        projection_point = a;  // Closest to point a
    } else if (t > 1.0) {
        projection_point = b;  // Closest to point b
    } else {
        projection_point = a + t * ab;  // On the segment
    }

    Eigen::Vector2d d = p - projection_point;
    double distance = d.norm();
    Eigen::Vector2d d_normalized = (distance > 0.0) ? d.normalized() : Eigen::Vector2d(0.0, 0.0);

    // Calculate derivatives
    double derivative_x = -d_normalized.x();
    double derivative_y = -d_normalized.y();

    return {distance, derivative_x, derivative_y};
}


ObstacleConstraint::VecBound
ObstacleConstraint::GetBounds () const
{
  VecBound bounds(GetRows());

  int row = 0;
  for (int id : node_ids_) {
    // if (ee_motion_->IsConstantNode(id))
    //   bounds.at(row) = ifopt::BoundZero;
    // else
    bounds.at(row) = ifopt::Bounds(min_distance_to_obstacle_, ifopt::inf);
    row++;
  }

  return bounds;
}

void ObstacleConstraint::FillJacobianBlock (std::string var_set, Jacobian& jac) const {
  if (var_set == ee_motion_->GetName()) {
    auto nodes = ee_motion_->GetNodes();
    int row = 0;
    for (int id : node_ids_) {
      Vector3d p = nodes.at(id).p();
      auto segments = terrain_->GetObstacles();

      // Initialize derivatives to zero
      double derivative_x = 0.0;
      double derivative_y = 0.0;

      for (const auto& segment : segments) {
        double distance, temp_derivative_x, temp_derivative_y;
        std::tie(distance, temp_derivative_x, temp_derivative_y) = 
            DistanceAndDerivativesToLineSegment({p.x(), p.y()}, segment);

        if (distance < min_distance_to_obstacle_) {
          derivative_x = temp_derivative_x;
          derivative_y = temp_derivative_y;
          break; // Break if this is the nearest obstacle
        }
      }

      // Fill the Jacobian entries for this node
      int idx_x = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kPos, X));
      int idx_y = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(id, kPos, Y));
      jac.coeffRef(row, idx_x) = derivative_x;
      jac.coeffRef(row, idx_y) = derivative_y;
      row++;
    }
  }
}

} /* namespace towr */
