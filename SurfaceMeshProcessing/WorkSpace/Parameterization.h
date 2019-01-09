/*本参数化方法参考论文<Parametrization and smooth approximation of surface triangulations>*/
#pragma once
#ifndef PARAMETERIZATION
#define PARAMETERIZATION

#include "MeshViewer/MeshDefinition.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <Eigen/SparseLU>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <fstream>
using namespace Eigen;

class ShapePreserveParameter
{
public:
	ShapePreserveParameter(Mesh input_mesh_);
	~ShapePreserveParameter(){};
private:
	bool getFirstBoundaryVertex();
	void getBoundaryVertexArray();
	void resetBoundaryVertex();
	std::vector<Mesh::Point> local_parametrization(Mesh::VertexHandle mid_vertex);//只对非边界点处理
	bool judge_intersect(Mesh::Point Line1_start_point, Mesh::Point Line1_end_point, Mesh::Point Line2_start_point, Mesh::Point Line2_end_point);//各用两个点表示一条直线
	float area_of_triangle(Mesh::Point p1, Mesh::Point p2, Mesh::Point p3);
	void cal_barycentric_coordinate(Mesh::Point p1, Mesh::Point p2, Mesh::Point p3, Mesh::Point mid_point, Mesh::Point &coordinate);//用一个三维点coordinate表示重心坐标
	void get_lamda_of_vertex(Mesh::VertexHandle mid_vertex, std::vector<float> &lamda_array);//只对非边界点处理,计算这个点的1-邻域点对应的系数
	void get_lamda_of_vertex_average(Mesh::VertexHandle mid_vertex, std::vector<float> &lamda_array);//采用均匀参数化
	void find_interior_point(std::vector<Mesh::VertexHandle> &interior_vertex_handle,std::vector<int> &new_vertex_idex);//将内点进行重新编号
	void construct_equation(SparseMatrix<float> &L, MatrixXf &b_x, MatrixXf &b_y, std::vector<Mesh::VertexHandle> &interior_vertex_handle);//组装矩阵以及解方程组
	void solve_and_reset_mesh();
public:
	std::vector<Mesh::VertexHandle> getBoundaryVertex();
	Mesh get_output_mesh();
	void parametrization();
private:
	Mesh input_mesh;
	Mesh output_mesh;
	Mesh::VertexHandle first_boundary_vertex;
	std::vector<Mesh::VertexHandle> boundary_vertex_array;
};

#endif // !PARAMETERIZATION