/*给定u,v方向的节点，计算基函数*/
#pragma once
#ifndef BASIS_FUNS
#define BASIS_FUNS
#include <vector>
#include <iostream>

class BasisFuns
{
public:
	BasisFuns(std::vector<float> U_, std::vector<float> V_, int k_u_, int k_v_);
	int findSpan(std::vector<float> knot,float t);
	float basis_U(int i, float t);
	float basis_V(int i, float t);
private:
	float basis(std::vector<float> knot, int i, float t, int p);
private:
	std::vector<float> U;//U方向的节点
	std::vector<float> V;//V方向的节点
	int k_u;//U方向基函数的阶
	int k_v;//V方向基函数的阶
};


#endif BASIS_FUNS