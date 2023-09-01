#include "fem.hpp"

FEM::FEM(std::string element): element_(element) {}

void FEM::AddPoint(Eigen::Vector3d point) {
  points_.push_back(point);
  points_alive_.push_back(true);
}


bool FEM::InitCloud() {
  pc_.width = points_.size();
  pc_.height = 1;
  pc_.points.resize(pc_.width * pc_.height);

  for (size_t i = 0; i < pc_.points.size(); ++i) {
    pc_.points[i].x = points_[i](0);
    pc_.points[i].y = points_[i](1);
    pc_.points[i].z = points_[i](2);
  }

  return pc_.width > 0;
}

bool FEM::MovingLeastSquares() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc_));

  // Create a KD-Tree
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);

  // Output has the PointNormal type in order to store the normals calculated by MLS
  pcl::PointCloud<pcl::PointNormal> mls_points;

  // Init object (second point type is for the normals, even if unused)
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
 
  mls.setComputeNormals (true);

  // Set parameters
  mls.setInputCloud (cloud);
  mls.setPolynomialOrder (2);
  mls.setSearchMethod (tree);
  mls.setSearchRadius (2);

  // Reconstruct
  mls.process (mls_points);

  if (mls_points.size() == 0) {
    std::cout << "MLS failed" << std::endl;
    return false;
  }

  // Get corresponding indexes: for each output point, returns the index of the
  // input one.
  mls_indices_.clear();

  pcl::PointIndicesPtr pIdx1 = mls.getCorrespondingIndices();
  for (unsigned int i = 0; i < pIdx1->indices.size(); i++)
    mls_indices_.push_back(pIdx1->indices[i]);

  int idxit = 0;
  for (unsigned int i = 0; i < points_.size(); i++) {
    int currentpos = i;
    if (currentpos == mls_indices_[idxit]) {
      idxit++;
    } else {
      points_alive_[i] = false;
    }
  }

  // Replace pc_ with mls_points
  pc_.width = mls_points.size();
  pc_.height = 1;
  pc_.points.resize(pc_.width * pc_.height);

  for (size_t i = 0; i < pc_.points.size(); ++i) {
    pc_.points[i].x = mls_points[i].x;
    pc_.points[i].y = mls_points[i].y;
    pc_.points[i].z = mls_points[i].z;
  }

  return true;
}


bool FEM::Triangulate() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc_));

  // Normal estimation
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud (cloud);
  n.setInputCloud(cloud);
  n.setSearchMethod(tree);
  n.setKSearch(20);
  n.compute(*normals);

  // Store normals in normals_ as Eigen::Vector3d
  normals_.clear();
  for (unsigned int i = 0; i < normals->size(); i++) {
    Eigen::Vector3d normal;
    normal << normals->points[i].normal_x,
              normals->points[i].normal_y,
              normals->points[i].normal_z;
    normals_.push_back(normal);
  }

  // Concatenate the XYZ and normal fields*
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*cloud, *normals, *cloud_with_normals);

    // Create search tree*
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal>);
  tree2->setInputCloud (cloud_with_normals);

  // Initialize objects
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;

  // Set the maximum distance between connected points (maximum edge length)
  gp3.setSearchRadius (1.5);

  // Set typical values for the parameters
  gp3.setMu (2.5);
  gp3.setMaximumNearestNeighbors (100);
  gp3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
  gp3.setMinimumAngle(M_PI/18); // 10 degrees
  gp3.setMaximumAngle(2*M_PI/3); // 120 degrees
  gp3.setNormalConsistency(false);

  // Get result
  gp3.setInputCloud (cloud_with_normals);
  gp3.setSearchMethod (tree2);
  gp3.reconstruct (mesh_);

  // Additional vertex information
  // To which part (cloud) does each vertex belong
  std::vector<int> parts = gp3.getPartIDs();
  // Whether the vertex status is [-1,0,1,2,3] = [NONE,FREE,FRINGE,BOUNDARY,COMPLETED]
  std::vector<int> states = gp3.getPointStates();

  // Get largest part
  std::unordered_map<int, int> map;
  int max_val = -1;
  int max_key = -1;
  for (unsigned int i=0; i<parts.size(); i++) {
    if (map.find(parts[i]) == map.end()) {
      map[parts[i]] = 1;
      if (max_val < 1) {
        max_val = 1;
        max_key = parts[i];
      }
    }
    else {
      map[parts[i]]++;
      if (map[parts[i]] > max_val) {
        max_val = map[parts[i]];
        max_key = parts[i];
      }
    }
  }

  // Get triangles of largest part only
  triangles_.clear();
  for (unsigned int i=0; i<mesh_.polygons.size(); i++) {
    unsigned int nver0 = mesh_.polygons[i].vertices[0];
    unsigned int nver1 = mesh_.polygons[i].vertices[1];
    unsigned int nver2 = mesh_.polygons[i].vertices[2];

    if (parts[nver0] == max_key && parts[nver1] == max_key && parts[nver2] == max_key) {
      std::vector<int> triangle;
      triangle.push_back(nver0);
      triangle.push_back(nver1);
      triangle.push_back(nver2);
      triangles_.push_back(triangle);
    }
  }
  
  return triangles_.size() > 0 ? true : false;
}



bool FEM::Compute(bool moving_least_squares) {
  bool ok = InitCloud();

  if (ok && moving_least_squares) {
    ok = MovingLeastSquares();
  }

  if (ok)
    ok = Triangulate();

  return ok;
}

void FEM::ComputeExtrusion() {
  std::vector<double> distances;
  for (std::vector<int> triangle : triangles_) {
    pcl::PointXYZ p0 = pc_.points[triangle[0]];
    pcl::PointXYZ p1 = pc_.points[triangle[1]];
    pcl::PointXYZ p2 = pc_.points[triangle[2]];

    double d01 = pcl::euclideanDistance(p0, p1);
    double d12 = pcl::euclideanDistance(p1, p2);
    double d20 = pcl::euclideanDistance(p2, p0);

    distances.push_back(d01);
    distances.push_back(d12);
    distances.push_back(d20);
  }

  // Get median
  std::sort(distances.begin(), distances.end());
  element_height_ = distances[distances.size()/2];
  
  // Compute second layer at a distance of 1/2 element height and with the direction of the normal vector
  points2_.clear();
  pc2_.width = points_.size();
  pc2_.height = 1;
  pc2_.points.resize(pc2_.width * pc2_.height);

  // Average normal
  Eigen::Vector3d normal;
  for (Eigen::Vector3d n : normals_) {
    normal -= n;
  }
  normal /= normals_.size();
  
  for (unsigned int i=0; i<points_.size(); i++) {
    Eigen::Vector3d point2 = points_[i] - element_height_/2 * normal;
    points2_.push_back(point2);

    pc2_.points[i].x = point2(0);
    pc2_.points[i].y = point2(1);
    pc2_.points[i].z = point2(2);
  }

  for (std::vector<int> triangle : triangles_) {
    std::vector<int> element;
    element.push_back(triangle[0] + points_.size());
    element.push_back(triangle[1] + points_.size());
    element.push_back(triangle[2] + points_.size());
    element.push_back(triangle[0]);
    element.push_back(triangle[1]);
    element.push_back(triangle[2]);
    elements_.push_back(element);
  }
}

std::vector<Eigen::Vector3d> FEM::GetExtrusionDelta() {
  std::vector<Eigen::Vector3d> delta;
  for (unsigned int i=0; i<points_.size(); i++) {
    delta.push_back(points2_[i] - points_[i]);
  }
  return delta;
}

std::vector<Eigen::Vector3d> FEM::GetExtrusion() {
  return points2_;
}

void FEM::SetExtrusion(std::vector<Eigen::Vector3d> extrusion_delta, double element_height) {
  for (unsigned int i=0; i<points_.size(); i++) {
    Eigen::Vector3d pt2 = points_[i] + extrusion_delta[i];
    points2_.push_back(pt2);
  }
  element_height_ = element_height;
}

double FEM::GetElementHeight() {
  return element_height_;
}


void FEM::ViewMesh(bool extrusion, 
                   std::vector<Eigen::Vector3d> cloud2,
                   std::vector<Eigen::Vector3d> cloud2extrusion,
                   std::pair<Eigen::Vector4d, Eigen::Vector3d> pose1,
                   std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2) {
                    
  pcl::PointCloud<pcl::PointXYZ> pc;
  pc.width = cloud2.size();
  pc.height = 1;
  pc.points.resize(pc_.width * pc_.height);
  for (unsigned int i=0; i<cloud2.size(); i++) {
    pc.points[i].x = cloud2[i](0);
    pc.points[i].y = cloud2[i](1);
    pc.points[i].z = cloud2[i](2);
  }

  ViewMesh(extrusion, pc, cloud2extrusion, pose1, pose2);
}

void FEM::ViewMesh(bool extrusion, 
                   pcl::PointCloud<pcl::PointXYZ> cloud2,
                   std::vector<Eigen::Vector3d> cloud2extrusion,
                   std::pair<Eigen::Vector4d, Eigen::Vector3d> pose1,
                   std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ> (pc_));

  pcl::visualization::PCLVisualizer viewer;
  viewer.setBackgroundColor(1, 1, 1);
  viewer.addPointCloud(cloud);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.0);
  
  for (unsigned int i=0; i<triangles_.size(); i++) {
    pcl::PointXYZ p0 = cloud->points[triangles_[i][0]];
    pcl::PointXYZ p1 = cloud->points[triangles_[i][1]];
    pcl::PointXYZ p2 = cloud->points[triangles_[i][2]];
    std::string name = "triangle_1_" + std::to_string(i);

    // Add line
    viewer.addLine<pcl::PointXYZ>(p0, p1, 0.0, 0.7, 0.0, name + "a");
    viewer.addLine<pcl::PointXYZ>(p1, p2, 0.0, 0.7, 0.0, name + "b");
    viewer.addLine<pcl::PointXYZ>(p2, p0, 0.0, 0.7, 0.0, name + "c");
  }

  if (extrusion) {
    for (unsigned int i=0; i<triangles_.size(); i++) {
      Eigen::Vector3d p0e = points2_[triangles_[i][0]];
      Eigen::Vector3d p1e = points2_[triangles_[i][1]];
      Eigen::Vector3d p2e = points2_[triangles_[i][2]];

      pcl::PointXYZ p0(p0e(0), p0e(1), p0e(2));
      pcl::PointXYZ p1(p1e(0), p1e(1), p1e(2));
      pcl::PointXYZ p2(p2e(0), p2e(1), p2e(2));
      std::string name = "triangle_1_" + std::to_string(i) + "_extrusion";

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name + "_a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.9, 0.7, name + "_b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.9, 0.7, name + "_c");
    }

    for (unsigned int i=0; i<points2_.size(); i++) {
      pcl::PointXYZ p0 = cloud->points[i];
      pcl::PointXYZ p1(points2_[i](0), points2_[i](1), points2_[i](2));
      std::string name = "line_1_" + std::to_string(i);

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name);
    }
  }

  // If pose1 not zero, then add as a thick point
  if (pose1.second.norm() > 0) {
    pcl::PointXYZ p1(pose1.second(0), pose1.second(1), pose1.second(2));
    viewer.addSphere(p1, 0.1, 0.0, 0.7, 0.0, "pose1");

      // Add a line in the direction of the quaternion in pose2.first
      std::pair<Eigen::Vector3d, Eigen::Vector3d> line_pts = QuaternionLine(pose1.first, pose1.second);
      pcl::PointXYZ p1_a(line_pts.first(0), line_pts.first(1), line_pts.first(2));
      pcl::PointXYZ p1_b(line_pts.second(0), line_pts.second(1), line_pts.second(2));
      viewer.addLine<pcl::PointXYZ>(p1_a, p1_b, 0.0, 0.7, 0.0, "pose1_line");
  }

  if (cloud2.points.size() > 0) {
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud2_color(cloud2.makeShared(), 0, 0, 255);
    viewer.addPointCloud(cloud2.makeShared(), cloud2_color, "cloud2");
    viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud2");

    for (unsigned int i=0; i<triangles_.size(); i++) {
      pcl::PointXYZ p0 = cloud2.points[triangles_[i][0]];
      pcl::PointXYZ p1 = cloud2.points[triangles_[i][1]];
      pcl::PointXYZ p2 = cloud2.points[triangles_[i][2]];
      std::string name = "triangle_2_" + std::to_string(i);

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.0, 0.0, name + "a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.0, 0.0, name + "b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.0, 0.0, name + "c");
    }

    if (extrusion) {
      for (unsigned int i=0; i<triangles_.size(); i++) {
        Eigen::Vector3d p0e = cloud2extrusion[triangles_[i][0]];
        Eigen::Vector3d p1e = cloud2extrusion[triangles_[i][1]];
        Eigen::Vector3d p2e = cloud2extrusion[triangles_[i][2]];

        pcl::PointXYZ p0(p0e(0), p0e(1), p0e(2));
        pcl::PointXYZ p1(p1e(0), p1e(1), p1e(2));
        pcl::PointXYZ p2(p2e(0), p2e(1), p2e(2));
        std::string name = "triangle_2_" + std::to_string(i) + "_extrusion";

        // Add line
        viewer.addLine<pcl::PointXYZ>(p0, p1, 0.9, 0.7, 0.7, name + "_a");
        viewer.addLine<pcl::PointXYZ>(p1, p2, 0.9, 0.7, 0.7, name + "_b");
        viewer.addLine<pcl::PointXYZ>(p2, p0, 0.9, 0.7, 0.7, name + "_c");
      }

      for (unsigned int i=0; i<points2_.size(); i++) {
        pcl::PointXYZ p0 = cloud2.points[i];
        pcl::PointXYZ p1(cloud2extrusion[i](0), cloud2extrusion[i](1), cloud2extrusion[i](2));
        std::string name = "line_2_" + std::to_string(i);

        // Add line
        viewer.addLine<pcl::PointXYZ>(p0, p1, 0.9, 0.7, 0.7, name);
      }

    }

    // If pose2 not zero, then add as a thick point
    if (pose2.second.norm() > 0) {
      pcl::PointXYZ p2(pose2.second(0), pose2.second(1), pose2.second(2));
      viewer.addSphere(p2, 0.1, 0.7, 0.0, 0.0, "pose2");

      // Add a line in the direction of the quaternion in pose2.first
      std::pair<Eigen::Vector3d, Eigen::Vector3d> line_pts = QuaternionLine(pose2.first, pose2.second);
      pcl::PointXYZ p2_a(line_pts.first(0), line_pts.first(1), line_pts.first(2));
      pcl::PointXYZ p2_b(line_pts.second(0), line_pts.second(1), line_pts.second(2));
      viewer.addLine<pcl::PointXYZ>(p2_a, p2_b, 0.7, 0.0, 0.0, "pose2_line");

    }
  }

  viewer.addCoordinateSystem(0.5);


  Eigen::Vector3d p(0.0,0.0,-1.0);
  double ang = 30*M_PI/180;
  Eigen::Vector3d axis(1,0,0);
  Eigen::Vector4d q = Eigen::Vector4d(cos(ang/2), axis(0)*sin(ang/2), axis(1)*sin(ang/2), axis(2)*sin(ang/2));
  
  double norm = q.norm();
  Eigen::Vector4d q_norm;
  if (norm == 0) {
    q_norm = Eigen::Vector4d(1.0, q(1), q(2), q(3));
  }
  else {
    q_norm = q / norm;
  }
  Eigen::Quaterniond q1(q_norm(0), q_norm(1), q_norm(2), q_norm(3));
  
  Eigen::Quaternion p_quat = Eigen::Quaternion(0.0, p(0), p(1), p(2));
  Eigen::Quaternion p1 = q1 * p_quat * q1.inverse();
  //std::cout << "Original " << p_quat << std::endl;
  //std::cout << "Rotated  " << p1 << " : " << p1.w() << " " << p1.x() << " " << p1.y() << " " << p1.z() << std::endl;
  Eigen::Vector3d p1_vec(p1.x(), p1.y(), p1.z());

  //std::cout << "Original " << p << std::endl;
  //std::cout << "Rotated  " << p1_vec << std::endl;
  pcl::PointXYZ p_o(0.0, 0.0, 0.0);
  pcl::PointXYZ p_a(p(0), p(1), p(2));
  pcl::PointXYZ p_b(p1_vec(0), p1_vec(1), p1_vec(2));
  // Add both points to viewer
  viewer.addSphere(p_a, 0.1, 0.0, 0.0, 0.0, "p_a");
  viewer.addSphere(p_b, 0.1, 0.0, 1.0, 0.0, "p_b");
  viewer.addLine<pcl::PointXYZ>(p_o, p_a, 0.0, 0.0, 0.0, "p_line_oa");
  viewer.addLine<pcl::PointXYZ>(p_o, p_b, 0.0, 0.0, 0.0, "p_line_ob");


  viewer.spin();
  viewer.close();
}


std::vector<std::vector<float>> FEM::GetNodes() {
  std::vector<std::vector<float>> points;
  for (Eigen::Vector3d pt: points_) {
    std::vector<float> point;
    point.push_back(pt[0]);
    point.push_back(pt[1]);
    point.push_back(pt[2]);
    points.push_back(point);
  }
  for (Eigen::Vector3d pt: points2_) {
    std::vector<float> point;
    point.push_back(pt[0]);
    point.push_back(pt[1]);
    point.push_back(pt[2]);
    points.push_back(point);
  }
  return points;
}

std::vector<Eigen::Vector3d> FEM::GetEigenNodes() {
  return points_;
}

std::vector<std::vector<int>> FEM::GetTriangles() {
  return triangles_;
}

void FEM::SetTriangles(std::vector<std::vector<int>> triangles) {
  triangles_ = triangles;
}

std::vector<std::vector<int>> FEM::GetElements() {
  return elements_;
}

void FEM::SetElements(std::vector<std::vector<int>> elements) {
  elements_ = elements;
}

pcl::PointCloud<pcl::PointXYZ> FEM::GetCloud() {
  return pc_;
}

std::pair<Eigen::Vector3d, Eigen::Vector3d> FEM::QuaternionLine2(
    Eigen::Vector4d qvec, Eigen::Vector3d point, double radius) {

  // Step 1: Normalize the quaternion
  Eigen::Quaterniond normalizedQuaternion(qvec.normalized());

  // Step 2: Convert the quaternion to a rotation matrix
  Eigen::Matrix3d rotationMatrix = normalizedQuaternion.toRotationMatrix();

  // Step 3: Define the line endpoint
  Eigen::Vector3d direction(radius, 0.0, 0.0); // Direction of the line
  Eigen::Vector3d endpoint = point + rotationMatrix * direction;

  return std::make_pair(point, endpoint);
}

std::pair<Eigen::Vector3d, Eigen::Vector3d> FEM::QuaternionLine(
  Eigen::Vector4d qvec, Eigen::Vector3d point, double radius) {

  Eigen::Quaterniond q(qvec(0), qvec(1), qvec(2), qvec(3));
  Eigen::Matrix3d rotation_matrix = q.normalized().toRotationMatrix();

  Eigen::Vector3d reference_vector(1.0, 0.0, 0.0);

  Eigen::Vector3d direction = rotation_matrix * reference_vector;
  direction.normalize();

  Eigen::Vector3d endpoint = point + direction;

  return std::make_pair(point, endpoint);
}