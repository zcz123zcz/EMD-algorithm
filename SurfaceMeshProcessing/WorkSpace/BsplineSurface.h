#pragma once
#ifndef BSPLINE_SURFACE
#define BSPLINE_SURFACE

#include "Parameterization.h"
#include "BasisFuns.h"
#include <memory>
#include <fstream>
#include <ctime>
class BsplineSurface
{
public:
	BsplineSurface(Mesh input_mesh_,int alpha_,int k_u_,int num_u_,int k_v_,int num_v_);
	~BsplineSurface();

public:
	void get_control_point();

private:
	void get_parameter();
	std::vector<float> get_knot(int k, int num);//k是基的阶数, num是划分0-1参数的份数
	void get_all_knot();
	void construct_local_smooth_matrix(int lamda, std::vector<int> &row, std::vector<int> &col, std::vector<float> &value);//构建光滑矩阵的子矩阵
	void construct_smoth_matrix();//构建光滑系数矩阵
private:
	void construct_cofficient_matrix();//构建拟合系数矩阵
	void remesh();//根据节点和参数重新计算网格
public:
	Mesh get_output_mesh();
private:
	Mesh input_mesh;
	float alpha;//光滑系数
	int k_u;//u方向的基阶数
	int num_u;//u方向的控制节点个数为num_u+1
	int k_v;//v方向基阶数
	int num_v;
private:
	Mesh parameter;
	std::vector<float> U;//u方向节点
	std::vector<float> V;//v方向节点
private:
	SparseMatrix<float> BaseMatrix;
	SparseMatrix<float> K_s;//smooth矩阵
	SparseMatrix<float> K_d;//拟合矩阵
	MatrixXf b;
	MatrixXf X;//节点
private:
	Mesh output_mesh;

};

#endif BSPLINE_SURFACE