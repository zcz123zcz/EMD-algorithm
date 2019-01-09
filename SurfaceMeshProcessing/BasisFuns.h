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
	int findSpan(std::vector<float> knot,float t);
	float basis_U(int i, float t);
	float basis_V(int i, float t);
private:
	float basis(std::vector<float> knot, int i, float t, int p);
private:
	std::vector<float> U;//U����Ľڵ�
	std::vector<float> V;//V����Ľڵ�
	int k_u;//U����������Ľ�
	int k_v;//V����������Ľ�
};


#endif BASIS_FUNS