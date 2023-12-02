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
#include <iostream>

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
  double y = base_linear_->GetPoint(t).p().y();
  g[k] = CalculateDistanceToNearestObstacle(x, y);
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

    if (obstacle_list.size() == 0)
      return;     

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
  double min_distance = 100.0;

  for (const auto& segment : obstacle_list) {
    double distance = DistanceToLineSegment({x, y}, segment);
    min_distance = std::min(min_distance, distance);
  }
  return min_distance;
}

Eigen::Vector2d ObstacleConstraint::calculateEdgeDirection(const Eigen::Vector2d& point, const std::pair<Eigen::Vector2d, Eigen::Vector2d>& segment) const {
    Eigen::Vector2d edgeVector = segment.second - segment.first; // Vector along the edge of the obstacle
    Eigen::Vector2d normalVector(-edgeVector.y(), edgeVector.x()); // Normal vector to the edge

    // Check the side of the point relative to the obstacle edge
    if ((point - segment.first).dot(normalVector) > 0) {
        // Point is on the left side of the obstacle (assuming standard coordinate system)
        // Rotate the edge vector 90 degrees counter-clockwise to guide along the edge
        return Eigen::Vector2d(-normalVector.x(), -normalVector.y());
    } else {
        // Point is on the right side of the obstacle
        // Rotate the edge vector 90 degrees clockwise
        return normalVector;
    }
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

    double ab_length_squared = ab.squaredNorm();
    double epsilon = 1e-6;  // small value to avoid division by zero

    // Check if line segment is actually a point (tiny obstacle)
    if (ab_length_squared < epsilon) {
        // Segment is effectively a point
        double norm_ap = ap.norm();
        if (norm_ap < 0.1) {  // bufferZoneWidth is a predefined constant for the buffer zone width
            Eigen::Vector2d edgeDirection = calculateEdgeDirection(point, segment);
            double derivative_x = edgeDirection.x();
            double derivative_y = edgeDirection.y();
            return {norm_ap, derivative_x, derivative_y};
        } else {
          if (norm_ap < epsilon) {
              // Point is on the segment point, derivatives are zero
              return {0.0, 0.0, 0.0};
          } else {
              double p_x = ap.x() / norm_ap;
              double p_y = ap.y() / norm_ap;
              return {norm_ap, p_x, p_y};
          }
        }
    }

    double t = ap.dot(ab) / (ab_length_squared + epsilon);
    Eigen::Vector2d projection_point = (t < 0.0) ? a : (t > 1.0) ? b : a + t * ab;

    Eigen::Vector2d d = p - projection_point;
    double distance = d.norm();

    if (distance < epsilon) {
        // Point is very close or on the segment, derivatives are zero
        return {0.0, 0.0, 0.0};
    }

    Eigen::Vector2d d_normalized = d / distance;
    double derivative_x = -d_normalized.x();
    double derivative_y = -d_normalized.y();

    return {distance, derivative_x, derivative_y};
}

} /* namespace towr */
