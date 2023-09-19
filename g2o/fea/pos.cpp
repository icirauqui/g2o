#include "pos.hpp"

POS::POS() {}

POS::POS(std::vector<Eigen::Vector3d> points,
         std::pair<Eigen::Vector4d, Eigen::Vector3d> pose){
  points_.push_back(points);
  pose_.push_back(pose);
}

POS::POS(std::vector<std::vector<float>> points,
         std::pair<Eigen::Vector4d, Eigen::Vector3d> pose){
  std::vector<Eigen::Vector3d> eigen_points;
  for (unsigned int i = 0; i < points.size(); i++) {
    eigen_points.push_back(Eigen::Vector3d(points[i][0], points[i][1], points[i][2]));
  }
  points_.push_back(eigen_points);
  pose_.push_back(pose);
}

void POS::Transform(Eigen::Vector4d r_im,
                    Eigen::Vector3d t,
                    double s) {
  // Get latest
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose = GetPose();
  std::vector<Eigen::Vector3d> points = GetPoints();

  // Translate
  pose.second += t;
  for (unsigned int i = 0; i < points.size(); i++) {
    points[i] += t;
  }

  // Rotate
  pose.first = ConcatenateQuaternions(pose.first, r_im);
  //Eigen::Vector4d r_pt = InvertQuaternion(r_im);
  for (unsigned int i = 0; i < points.size(); i++) {
    points[i] -= pose.second;
    points[i] = QuaternionRotatePoint(r_im, points[i]);
    points[i] += pose.second;
  }

  // Scale
  pose.second *= s;
  for (unsigned int i = 0; i < points.size(); i++) {
    points[i] -= pose.second;
    points[i] *= s;
    points[i] += pose.second;
  }

  // Store new
  pose_.push_back(pose);
  points_.push_back(points);
}

std::vector<Eigen::Vector3d> POS::GetPoints(int idx)
{
  if (idx == -1)
  {
    idx = points_.size() - 1;
  }
  return points_[idx];
}

std::pair<Eigen::Vector4d, Eigen::Vector3d> POS::GetPose(int idx)
{
  if (idx == -1)
  {
    idx = pose_.size() - 1;
  }
  return pose_[idx];
}

int POS::LenHistory()
{
  return points_.size();
}

void POS::AddData(std::vector<Eigen::Vector3d> points,
                  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose)
{
  points_.push_back(points);
  pose_.push_back(pose);
}

Eigen::Vector4d POS::ComputeQuaternionRotation(const Eigen::Vector4d& v1, 
                                               const Eigen::Vector4d& v2) {

    Eigen::Quaterniond q1(v1(0), v1(1), v1(2), v1(3));
    Eigen::Quaterniond q2(v2(0), v2(1), v2(2), v2(3));
    
    // Compute the quaternion difference between q1 and q2
    Eigen::Quaterniond rot = q2 * q1.inverse();

    // Normalize the resulting quaternion to ensure it's a unit quaternion
    rot.normalize();

    return Eigen::Vector4d(rot.w(), rot.x(), rot.y(), rot.z());
}


Eigen::Vector4d POS::EulerToQuaternion(const Eigen::Vector3d& euler) {
    // Convert Euler angles (roll, pitch, yaw) to quaternion
    Eigen::Quaterniond quat = Eigen::AngleAxisd(euler(2), Eigen::Vector3d::UnitZ()) *
                              Eigen::AngleAxisd(euler(1), Eigen::Vector3d::UnitY()) *
                              Eigen::AngleAxisd(euler(0), Eigen::Vector3d::UnitX());
    return Eigen::Vector4d(quat.w(), quat.x(), quat.y(), quat.z());
}

Eigen::Vector3d POS::QuaternionToEuler(const Eigen::Quaterniond& quat) {
    // Convert quaternion to Euler angles (roll, pitch, yaw)
    Eigen::Vector3d euler = quat.toRotationMatrix().eulerAngles(2, 1, 0);
    return euler;
}

Eigen::Vector4d POS::QuaternionFromAngleAxis(const Eigen::Vector3d& axis,
                                             const double ang) {
  return Eigen::Vector4d(cos(ang/2), axis(0)*sin(ang/2), axis(1)*sin(ang/2), axis(2)*sin(ang/2));
}

double POS::DegToRad(double deg) {
    return deg * M_PI / 180.0;
}

double POS::RadToDeg(double rad) {
    return rad * 180.0 / M_PI;
}

Eigen::Vector4d POS::ConcatenateQuaternions(const Eigen::Vector4d &qvec1,
                                            const Eigen::Vector4d &qvec2)
{
  const Eigen::Vector4d normalized_qvec1 = NormalizeQuaternion(qvec1);
  const Eigen::Vector4d normalized_qvec2 = NormalizeQuaternion(qvec2);
  const Eigen::Quaterniond quat1(normalized_qvec1(0), normalized_qvec1(1),
                                 normalized_qvec1(2), normalized_qvec1(3));
  const Eigen::Quaterniond quat2(normalized_qvec2(0), normalized_qvec2(1),
                                 normalized_qvec2(2), normalized_qvec2(3));
  const Eigen::Quaterniond cat_quat = quat2 * quat1;
  return NormalizeQuaternion(
      Eigen::Vector4d(cat_quat.w(), cat_quat.x(), cat_quat.y(), cat_quat.z()));
}

Eigen::Vector4d POS::RotateQuaternion(const Eigen::Vector4d &base,
                                      const Eigen::Vector4d &rot)
{
  Eigen::Quaterniond base_quat(base(0), base(1), base(2), base(3));
  Eigen::Quaterniond rot_quat(rot(0), rot(1), rot(2), rot(3));
  Eigen::Vector4d rot_1 = InvertQuaternion(rot);
  Eigen::Quaterniond rot_quat_1(rot_1(0), rot_1(1), rot_1(2), rot_1(3));
  Eigen::Quaterniond rotated = rot_quat * base_quat * rot_quat_1;
  return Eigen::Vector4d(rotated.w(), rotated.x(), rotated.y(), rotated.z());
}

Eigen::Vector4d POS::InvertQuaternion(const Eigen::Vector4d &qvec)
{
  return Eigen::Vector4d(qvec(0), -qvec(1), -qvec(2), -qvec(3));
}

Eigen::Vector4d POS::NormalizeQuaternion(const Eigen::Vector4d &qvec) {
  const double norm = qvec.norm();
  if (norm == 0) {
    return Eigen::Vector4d(1.0, qvec(1), qvec(2), qvec(3));
  }
  else {
    return qvec / norm;
  }
}

//Eigen::Vector3d POS::QuaternionRotatePoint(const Eigen::Vector4d &qvec,
//                                           const Eigen::Vector3d &point)
//{
//  const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
//  const Eigen::Quaterniond quat(normalized_qvec(0), normalized_qvec(1),
//                                normalized_qvec(2), normalized_qvec(3));
//  return quat * point;
//}

Eigen::Vector3d POS::QuaternionRotatePoint(const Eigen::Vector4d &q,
                                           const Eigen::Vector3d &p)
{
  const Eigen::Vector4d qn = NormalizeQuaternion(q);
  const Eigen::Quaterniond quat(qn(0), qn(1), qn(2), qn(3));

  Eigen::Quaternion p_quat = Eigen::Quaternion(0.0, p(0), p(1), p(2));

  Eigen::Quaternion p_rot = quat * p_quat * quat.inverse();

  Eigen::Vector3d p_rot_vec(p_rot.x(), p_rot.y(), p_rot.z());

  return p_rot_vec;
}


Eigen::Vector3d POS::QuaternionRotatePointAngle(const Eigen::Vector3d &axis,
                                             const Eigen::Vector3d &p,
                                             const double angle) {
  double ang = angle*M_PI/180;
  Eigen::Vector4d q = Eigen::Vector4d(cos(ang/2), axis(0)*sin(ang/2), axis(1)*sin(ang/2), axis(2)*sin(ang/2));
  
  Eigen::Vector4d q_norm = NormalizeQuaternion(q);
  Eigen::Quaterniond q1(q_norm(0), q_norm(1), q_norm(2), q_norm(3));
  
  Eigen::Quaternion p_quat = Eigen::Quaternion(0.0, p(0), p(1), p(2));
  Eigen::Quaternion p1 = q1 * p_quat * q1.inverse();
  Eigen::Vector3d p1_vec(p1.x(), p1.y(), p1.z());

  return p1_vec;
}

Eigen::Vector4d POS::VectorToQuaternion(const Eigen::Vector3d& euler_vector) {
    // Convert the Euler vector to a quaternion
    Eigen::Vector4d quaternion;
    quaternion.head<3>() = euler_vector.normalized() * std::sin(euler_vector.norm() / 2.0);
    quaternion(3) = std::cos(euler_vector.norm() / 2.0);
    return quaternion;
}

