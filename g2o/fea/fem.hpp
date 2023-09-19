#ifndef FE_HPP
#define FE_HPP


#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include <eigen3/Eigen/Dense>


// PCL Libraries
//#include <pcl/io/pcd_io.h>
//#include <pcl/io/vtk_io.h>
//#include <pcl/point_types.h>

//#include <pcl/common/common.h>
//#include <pcl/features/normal_3d_omp.h>
#include <pcl/surface/mls.h>
//#include <pcl/surface/impl/mls.hpp>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/visualization/cloud_viewer.h>

#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <pcl/surface/gp3.h>

#include <pcl/visualization/vtk.h>

#include <pcl/console/parse.h>
#include <pcl/io/vtk_lib_io.h>



class FEM {

public:

  FEM(std::string element);

  void AddPoint(Eigen::Vector3d point);

  bool InitCloud();

  bool MovingLeastSquares();
  
  bool Triangulate();

  bool Compute(bool moving_least_squares = true);

  void ComputeExtrusion();
  std::vector<Eigen::Vector3d> GetExtrusionDelta();
  std::vector<Eigen::Vector3d> GetExtrusion();
  void SetExtrusion(std::vector<Eigen::Vector3d> extrusion_delta, double element_height);
  double GetElementHeight();

  void ViewMesh(bool extrusion = false,
                std::vector<Eigen::Vector3d> cloud2 = std::vector<Eigen::Vector3d>(),
                std::vector<Eigen::Vector3d> cloud2extrusion = std::vector<Eigen::Vector3d>(),
                std::pair<Eigen::Vector4d, Eigen::Vector3d> pose1 = std::make_pair(Eigen::Vector4d(0.0, 0.0, 0.0, 0.0), Eigen::Vector3d(0.0, 0.0, 0.0)),
                std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2 = std::make_pair(Eigen::Vector4d(0.0, 0.0, 0.0, 0.0), Eigen::Vector3d(0.0, 0.0, 0.0)));

  void ViewMesh(bool extrusion = false,
                pcl::PointCloud<pcl::PointXYZ> cloud2 = pcl::PointCloud<pcl::PointXYZ>(),
                std::vector<Eigen::Vector3d> cloud2extrusion = std::vector<Eigen::Vector3d>(),
                std::pair<Eigen::Vector4d, Eigen::Vector3d> pose1 = std::make_pair(Eigen::Vector4d(0.0, 0.0, 0.0, 0.0), Eigen::Vector3d(0.0, 0.0, 0.0)),
                std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2 = std::make_pair(Eigen::Vector4d(0.0, 0.0, 0.0, 0.0), Eigen::Vector3d(0.0, 0.0, 0.0)));

  std::vector<std::vector<float>> GetNodes();
  std::vector<Eigen::Vector3d> GetEigenNodes();

  std::vector<std::vector<int>> GetTriangles();
  void SetTriangles(std::vector<std::vector<int>> triangles);

  std::vector<std::vector<int>> GetElements();
  void SetElements(std::vector<std::vector<int>> elements);

  pcl::PointCloud<pcl::PointXYZ> GetCloud();


private:

  std::pair<Eigen::Vector3d, Eigen::Vector3d> QuaternionLine2(
    Eigen::Vector4d qvec, Eigen::Vector3d point, double radius = 1.0);
  std::pair<Eigen::Vector3d, Eigen::Vector3d> QuaternionLine(
    Eigen::Vector4d qvec, Eigen::Vector3d point, double radius = 1.0);

  std::string element_;

  std::vector<Eigen::Vector3d> points_, points2_;
  std::vector<bool> points_alive_;

  pcl::PointCloud<pcl::PointXYZ> pc_, pc2_;

  std::vector<std::vector<int>> triangles_, elements_;

  std::vector<int> mls_indices_;

  std::vector<Eigen::Vector3d> normals_;

  pcl::PolygonMesh mesh_;

  double element_height_ = 0.0;

  // Interface 
  float mls_search_radius_ = 1.0;
  int mls_polynomial_order_ = 3;
  float mesh_mu_ = 2.5;
  float mesh_search_radius_ = 10.0;
  int mesh_max_neighbours_ = 25;
  int mesh_surf_angle_ = 150;
  int mesh_min_angle_ = 5;
  int mesh_max_angle_ = 85;


};
















#endif