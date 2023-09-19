
#ifndef POS_HPP
#define POS_HPP

#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>

class POS
{

public:
  POS();
  POS(std::vector<Eigen::Vector3d> points,
      std::pair<Eigen::Vector4d, Eigen::Vector3d> pose);
  POS(std::vector<std::vector<float>> points,
      std::pair<Eigen::Vector4d, Eigen::Vector3d> pose);

  void Transform(Eigen::Vector4d r_im,
                 Eigen::Vector3d t,
                 double s);

  // Accessors
  // Return latest by default
  // User can request a specific set: 0, 1, 2, ...
  std::vector<Eigen::Vector3d> GetPoints(int idx = -1);
  std::pair<Eigen::Vector4d, Eigen::Vector3d> GetPose(int idx = -1);
  int LenHistory();

  void AddData(std::vector<Eigen::Vector3d> points,
               std::pair<Eigen::Vector4d, Eigen::Vector3d> pose);

  Eigen::Vector4d ComputeQuaternionRotation(const Eigen::Vector4d& v1, 
                                            const Eigen::Vector4d& v2);

  Eigen::Vector4d EulerToQuaternion(const Eigen::Vector3d& euler);

  Eigen::Vector3d QuaternionToEuler(const Eigen::Quaterniond& quat);

  Eigen::Vector4d QuaternionFromAngleAxis(const Eigen::Vector3d& axis,
                                          const double ang);

  double DegToRad(double deg);

  double RadToDeg(double rad);

private:
  Eigen::Vector4d ConcatenateQuaternions(const Eigen::Vector4d &qvec1,
                                         const Eigen::Vector4d &qvec2);

  Eigen::Vector4d RotateQuaternion(const Eigen::Vector4d &base,
                                        const Eigen::Vector4d &rot);

  Eigen::Vector4d InvertQuaternion(const Eigen::Vector4d &qvec);

  Eigen::Vector4d NormalizeQuaternion(const Eigen::Vector4d &qvec);

  Eigen::Vector3d QuaternionRotatePoint(const Eigen::Vector4d &qvec,
                                        const Eigen::Vector3d &point);

  Eigen::Vector3d QuaternionRotatePointAngle(const Eigen::Vector3d &axis,
                                             const Eigen::Vector3d &p,
                                             const double angle);

  Eigen::Vector4d VectorToQuaternion(const Eigen::Vector3d& euler_vector);



  std::vector<std::vector<Eigen::Vector3d>> points_;
  std::vector<std::pair<Eigen::Vector4d, Eigen::Vector3d>> pose_;
};



#endif