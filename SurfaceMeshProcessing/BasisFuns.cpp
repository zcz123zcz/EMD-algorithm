#include "BasisFuns.h"

BasisFuns::BasisFuns(std::vector<float> U_, std::vector<float> V_, int k_u_, int k_v_) :U(U_), V(V_), k_u(k_u_), k_v(k_v_)
{

}

float BasisFuns::basis(std::vector<float> knot, int i, float t, int p)
{
	float delta = 1e-10;
	float y = 0;
	if (p == 0)
	{
		if (t < knot[i + 1] && t >= knot[i]) 
		{
			y = 1;
		}
		return y;
	}
	else
	{
		y = (t - knot[i]) / (knot[i + p] - knot[i] + delta)*basis(knot, i, t, p - 1) + (knot[i + p + 1] - t) / (knot[i + p + 1] - knot[i + 1] + delta)*basis(knot, i + 1, t, p - 1);
		return y;
	}
}

int BasisFuns::findSpan(std::vector<float> knot, float t)
{
	int n = knot.size();
	int low = 0;
	int high = n - 1;
	int mid = (low + high) / 2;

	if (t > knot[high - 1])
	{
		return high - 1;
	}
	while (t < knot[mid] || t >= knot[mid + 1])
	{
		if (t < knot[mid])
		{
			high = mid;
		}
		else
		{
			low = mid;
		}
		mid = (low + high) / 2;
	}
	return mid;
}

float BasisFuns::basis_U(int i, float t)
{
	return basis(U, i, t, k_u);
}

float BasisFuns::basis_V(int i, float t)
{
	return basis(V, i, t, k_v);
}