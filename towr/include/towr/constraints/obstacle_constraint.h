#ifndef TOWR_CONSTRAINTS_OBSTACLE_CONSTRAINT_H_
#define TOWR_CONSTRAINTS_OBSTACLE_CONSTRAINT_H_

#include <ifopt/constraint_set.h>

#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/terrain/height_map.h>

namespace towr {

/**
 * @brief 
 *
 * @ingroup Constraints
 */
class ObstacleConstraint : public ifopt::ConstraintSet {
public:
  using Vector3d = Eigen::Vector3d;

  /**
   * @brief Constructs a obstacle constraint.
   * @param terrain  The terrain height value and slope for each position x,y.
   * @param ee_motion_id The name of the endeffector variable set.
   */
  ObstacleConstraint (const HeightMap::Ptr& terrain, std::string ee_motion_id);
  virtual ~ObstacleConstraint () = default;

  void InitVariableDependedQuantities(const VariablesPtr& x) override;

  VectorXd GetValues() const override;
  VecBound GetBounds() const override;
  void FillJacobianBlock (std::string var_set, Jacobian&) const override;

private:
  NodesVariablesPhaseBased::Ptr ee_motion_; ///< the position of the endeffector.
  HeightMap::Ptr terrain_;    ///< the height map of the current terrain.

  std::string ee_motion_id_;  ///< the name of the endeffector variable set.
  std::vector<int> node_ids_; ///< the indices of the nodes constrained.

  double CalculateDistanceToNearestObstacle(double x, double y);
  double DistanceToLineSegment(const Eigen::Vector2d& point, const std::pair<Eigen::Vector2d, Eigen::Vector2d>& segment);
  std::tuple<double, double, double> DistanceAndDerivativesToLineSegment(const Eigen::Vector2d& point, const std::pair<Eigen::Vector2d, Eigen::Vector2d>& segment);
};

} /* namespace towr */

#endif /* TOWR_CONSTRAINTS_OBSTACLE_CONSTRAINT_H_ */
