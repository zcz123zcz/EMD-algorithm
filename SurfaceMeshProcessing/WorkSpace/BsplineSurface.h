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
	std::vector<float> get_knot(int k, int num);//k�ǻ��Ľ���, num�ǻ���0-1�����ķ���
	void get_all_knot();
	void construct_local_smooth_matrix(int lamda, std::vector<int> &row, std::vector<int> &col, std::vector<float> &value);//�����⻬������Ӿ���
	void construct_smoth_matrix();//�����⻬ϵ������
private:
	void construct_cofficient_matrix();//�������ϵ������
	void remesh();//���ݽڵ�Ͳ������¼�������
public:
	Mesh get_output_mesh();
private:
	Mesh input_mesh;
	float alpha;//�⻬ϵ��
	int k_u;//u����Ļ�����
	int num_u;//u����Ŀ��ƽڵ����Ϊnum_u+1
	int k_v;//v���������
	int num_v;
private:
	Mesh parameter;
	std::vector<float> U;//u����ڵ�
	std::vector<float> V;//v����ڵ�
private:
	SparseMatrix<float> BaseMatrix;
	SparseMatrix<float> K_s;//smooth����
	SparseMatrix<float> K_d;//��Ͼ���
	MatrixXf b;
	MatrixXf X;//�ڵ�
private:
	Mesh output_mesh;

};

#endif BSPLINE_SURFACE