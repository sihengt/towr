#ifndef TOWR_CONSTRAINTS_OBSTACLE_CONSTRAINT_H_
#define TOWR_CONSTRAINTS_OBSTACLE_CONSTRAINT_H_

#include <ifopt/constraint_set.h>
#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/terrain/height_map.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>

#include <towr/variables/spline.h>
#include <towr/variables/spline_holder.h>

#include "time_discretization_constraint.h"


namespace towr {

/**
 * @brief 
 *
 * @ingroup Constraints
 */
class ObstacleConstraint : public TimeDiscretizationConstraint {
public:
  using Vector3d = Eigen::Vector3d;

  /**
   * @brief Constructs a obstacle constraint.
   * @param terrain  The terrain height value and slope for each position x,y.
   * @param ee_motion_id The name of the endeffector variable set.
   */
  ObstacleConstraint (const HeightMap::Ptr& terrain, double T, double dt, const SplineHolder& spline_holder);
  virtual ~ObstacleConstraint () = default;
  
  /**
   * @brief Sets the constraint value a specific time t, corresponding to node k.
   * @param t  The time along the trajectory to set the constraint.
   * @param k  The index of the time t, so t=k*dt
   * @param[in/out] g  The complete vector of constraint values, for which the
   *                   corresponding row must be filled.
   */
  void UpdateConstraintAtInstance(double t, int k, VectorXd& g) const override;
  
    /**
   * @brief Sets upper/lower bound a specific time t, corresponding to node k.
   * @param t  The time along the trajectory to set the bounds.
   * @param k  The index of the time t, so t=k*dt
   * @param[in/out] b The complete vector of bounds, for which the corresponding
   *                  row must be set.
   */
  void UpdateBoundsAtInstance(double t, int k, VecBound& b) const override;
  
  /**
   * @brief Sets Jacobian rows at a specific time t, corresponding to node k.
   * @param t  The time along the trajcetory to set the bounds.
   * @param k  The index of the time t, so t=k*dt
   * @param var_set The name of the ifopt variables currently being queried for.
   * @param[in/out] jac  The complete Jacobian, for which the corresponding
   *                     row and columns must be set.
   */
  void UpdateJacobianAtInstance(double t, int k, std::string var_set, Jacobian& jac) const override;
  
private:
  NodesVariablesPhaseBased::Ptr ee_motion_; ///< the position of the endeffector.
  HeightMap::Ptr terrain_;    ///< the height map of the current terrain.
  NodeSpline::Ptr base_linear_;

  std::string ee_motion_id_;  ///< the name of the endeffector variable set.
  std::vector<int> node_ids_; ///< the indices of the nodes constrained.
  std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> obstacle_list;

  double CalculateDistanceToNearestObstacle(double x, double y) const;
  double DistanceToLineSegment(const Eigen::Vector2d& point, const std::pair<Eigen::Vector2d, Eigen::Vector2d>& segment) const;
  Eigen::Vector2d calculateEdgeDirection(const Eigen::Vector2d& point, const std::pair<Eigen::Vector2d, Eigen::Vector2d>& segment) const;
  std::tuple<double, double, double> DistanceAndDerivativesToLineSegment(const Eigen::Vector2d& point, const std::pair<Eigen::Vector2d, Eigen::Vector2d>& segment) const;
  double min_distance_to_obstacle_ = 0.1; // [m]
};

} /* namespace towr */

#endif /* TOWR_CONSTRAINTS_OBSTACLE_CONSTRAINT_H_ */
