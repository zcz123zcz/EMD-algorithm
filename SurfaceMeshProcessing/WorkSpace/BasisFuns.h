/*����u,v����Ľڵ㣬���������*/
#pragma once
#ifndef BASIS_FUNS
#define BASIS_FUNS
#include <vector>
#include <iostream>

class BasisFuns
{
public:
	BasisFuns(std::vector<float> U_, std::vector<float> V_, int k_u_, int k_v_);
	int findSpan_U(float t);
	int findSpan_V(float t);
public:
	float basis_u_v(int r, int c, float u, float v);
private:
	float basis(std::vector<float> knot, int i, float t, int p);
	float basis_U(int i, float t);
	float basis_V(int i, float t);
	void remove_repeat_knot();
	int findSpan(std::vector<float> knot, float t, int k);
private:
	std::vector<float> U;//U����Ľڵ�
	std::vector<float> V;//V����Ľڵ�
	std::vector<float> U_t;//ȥ���ظ��ڵ�
	std::vector<float> V_t;//ȥ���ظ��ڵ�
	int k_u;//U����������Ľ�
	int k_v;//V����������Ľ�
};


#endif BASIS_FUNS