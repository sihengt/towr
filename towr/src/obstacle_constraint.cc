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


ObstacleConstraint::ObstacleConstraint (const HeightMap::Ptr& terrain, double T, double dt, 
                                        const SplineHolder& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "obstacle")
{
  // We need the position of the base to calculate our constraint
  base_linear_ = spline_holder.base_linear_;
  
  // Terrain holds obstacle information.
  terrain_ = terrain;
  // TODO: Not all terrains have obstacles. For those without, we could essentially not use this constraint.
  obstacle_list = terrain_->GetObstacles();

  // Pretty much only one constraint (distance of base to obstacle) for each dt.
  SetRows(GetNumberOfNodes());
}

void 
ObstacleConstraint::UpdateConstraintAtInstance(double t, int k, VectorXd& g) const
{
  double x = base_linear_->GetPoint(t).p().x();
  double y = base_linear_->GetPoint(t).p().x();
  g.Row(k) = CalculateDistanceToNearestObstacle(x, y);
}

void
ObstacleConstraint::UpdateBoundsAtInstance(double t, int k, VecBound& bounds) const
{
  bounds.at(k) = ifopt::Bounds(min_distance_to_obstacle_, ifopt::inf);
}

void 
ObstacleConstraint::UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const
{
  if (var_set == id::base_lin_nodes) {
    // Get the base position at time t
    Eigen::Vector2d base_pos(base_linear_->GetPoint(t).p().x(), base_linear_->GetPoint(t).p().y());

    // Find the closest obstacle and get the derivatives
    double min_distance = std::numeric_limits<double>::max();
    double derivative_x = 0.0;
    double derivative_y = 0.0;

    for (const auto& segment : obstacle_list) {
      auto [distance, dx, dy] = DistanceAndDerivativesToLineSegment(base_pos, segment);
      if (distance < min_distance) {
        min_distance = distance;
        derivative_x = dx;
        derivative_y = dy;
      }
    }

    // Fill the Jacobian matrix
    // Assuming jac is already of correct size and initialized to zero
    int idx_x = k * 3; // Assuming 3 variables (x, y, z) per node
    int idx_y = k * 3 + 1;

    jac.coeffRef(0, idx_x) = derivative_x;
    jac.coeffRef(0, idx_y) = derivative_y;
  }
}

double ObstacleConstraint::CalculateDistanceToNearestObstacle(double x, double y) const
{
  double min_distance = std::numeric_limits<double>::max();
  for (const auto& segment : obstacle_list) {
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



// TODO: OBSOLETE!
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
