
#include "BsplineSurface.h"

BsplineSurface::BsplineSurface(Mesh input_mesh_, int alpha_, int k_u_, int num_u_, int k_v_, int num_v_) 
	:input_mesh(input_mesh_), alpha(alpha_), k_u(k_u_), k_v(k_v_), num_u(num_u_), num_v(num_v_)
{

}

BsplineSurface::~BsplineSurface()
{

}

void BsplineSurface::get_parameter()
{
	std::tr1::shared_ptr<ShapePreserveParameter> sInv(new ShapePreserveParameter(input_mesh));
	parameter = (sInv.get())->get_output_mesh();
}

std::vector<float> BsplineSurface::get_knot(int k, int num)
{
	std::vector<float> knot;
	for (int i = 0; i <= k;i++)
	{
		knot.push_back(0.0f);
	}
	for (int i = 1; i <= num - k;i++)
	{
		knot.push_back(float(i) / float(num - k + 1));
	}
// 	int temp_num = (num - k) / 4;
// 	float temp_gap = (float(1)/float(3)) / float(temp_num);
// 	float cur_gap = 0;
// 	for (int i = 0; i < temp_num;i++)
// 	{
// 		knot.push_back(float(i + 1)*temp_gap);
// 		cur_gap = float(i + 1)*temp_gap;
// 	}
// 	temp_gap = float(1.0 / 3.0) / float(2*temp_num);
// 	float new_gap = 0;
// 	for (int i = 0; i < 2*temp_num; i++)
// 	{
// 		knot.push_back(float(i + 1)*temp_gap+cur_gap);
// 		new_gap = float(i + 1)*temp_gap + cur_gap;
// 	}
// 
// 	int left_num = (num - k) - 3 * temp_num;//剩余的节点数
// 	float end = float(num - k) / float(num - k + 1);
// 	float left_gap = (end - new_gap)/float(left_num);
// 	for (int i = 0; i < left_num;i++)
// 	{
// 		knot.push_back(new_gap + float(i + 1)*left_gap);
// 	}

	for (int i = 0; i <= k;i++)
	{
		knot.push_back(1.0f);
	}
	return knot;
}

void BsplineSurface::get_all_knot()
{
	U = get_knot(k_u, num_u);
	V = get_knot(k_v, num_v);
}

void BsplineSurface::construct_local_smooth_matrix(int lamda, std::vector<int> &row, std::vector<int> &col, std::vector<float> &value)
{
	row.push_back(lamda*(num_u + 1)), col.push_back(lamda*(num_u + 1)), value.push_back(1.0);
	row.push_back(lamda*(num_u + 1)), col.push_back(lamda*(num_u + 1)+1), value.push_back(-1.0);
	for (int i = lamda*(num_u + 1)+1; i < lamda*(num_u + 1) + num_u;i++)
	{
		row.push_back(i), col.push_back(i-1), value.push_back(-1.0);
		row.push_back(i), col.push_back(i), value.push_back(2.0);
		row.push_back(i), col.push_back(i + 1), value.push_back(-1.0);
	}
	row.push_back((lamda + 1)*(num_u + 1) - 1), col.push_back((lamda + 1)*(num_u + 1) - 2), value.push_back(-1.0);
	row.push_back((lamda + 1)*(num_u + 1) - 1), col.push_back((lamda + 1)*(num_u + 1) - 1), value.push_back(1.0);
}

void BsplineSurface::construct_smoth_matrix()
{
	std::vector<int> row;
	std::vector<int> col;
	std::vector<float> value;
	for (int lamda = 0; lamda <= num_v;lamda++)
	{
		construct_local_smooth_matrix(lamda, row, col, value);
	}
	std::vector<Triplet<float>> triplets;
	for (int i = 0; i < row.size();i++)
	{
		//std::cout << row[i] << "   " << col[i] << "    " << value[i] << std::endl;
		triplets.emplace_back(row[i], col[i], value[i]);
	}
	K_s = SparseMatrix<float>((num_u + 1)*(num_v + 1), (num_u + 1)*(num_v + 1));
	K_s.setFromTriplets(triplets.begin(), triplets.end());
}
void BsplineSurface::construct_cofficient_matrix()
{
	get_parameter();//先得到参数
	get_all_knot();//计算节点
	std::tr1::shared_ptr<BasisFuns> sInv(new BasisFuns(U,V,k_u,k_v));
	int M = (num_u + 1)*(num_v + 1);
	int N = input_mesh.n_vertices();
	BaseMatrix=SparseMatrix<float>(N, M);
	MatrixXf P(N, 3);
	std::vector<int> row;std::vector<int> col;std::vector<float> value;
	clock_t startTime, endTime;
	startTime = clock();//计时开始
	for (Mesh::VertexIter v_it = input_mesh.vertices_begin(); v_it != input_mesh.vertices_end();v_it++)
	{
		P((*v_it).idx(), 0) = input_mesh.point(*v_it).data()[0], P((*v_it).idx(), 1) = input_mesh.point(*v_it).data()[1], P((*v_it).idx(), 2) = input_mesh.point(*v_it).data()[2];
		Mesh::Point point = parameter.point(*v_it);
		int id_u = (sInv.get())->findSpan_U(point.data()[0]);
		int id_v = (sInv.get())->findSpan_V(point.data()[1]);
		/*确定基函数不为0的区间*/
		int left_v = id_v - k_v >= 0 ? id_v - k_v : 0;
		int right_v = id_v + k_v <= num_v ? id_v + k_v : num_v;
		int low_u = id_u - k_u >= 0 ? id_u - k_u : 0;
		int up_u = id_u + k_u <= num_u ? id_u + k_u : num_u;

		for (int i = low_u; i <= up_u;i++)
		{
			for (int j = left_v; j <= right_v;j++)
			{
				float temp_value = (sInv.get())->basis_u_v(i, j, point.data()[0], point.data()[1]);
				int temp_basis_id = i*(num_v+1) + j;
				row.push_back((*v_it).idx()), col.push_back(temp_basis_id), value.push_back(temp_value);
			}
		}
	}
	endTime = clock();
	std::cout << "矩阵装配的时间为: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
	std::vector<Triplet<float>> triplets;
	for (int i = 0; i < row.size();i++)
	{
		triplets.emplace_back(row[i], col[i], value[i]);
	}
	BaseMatrix.setFromTriplets(triplets.begin(), triplets.end());
	K_d = BaseMatrix.transpose()*BaseMatrix;
	b = BaseMatrix.transpose()*P;
}

void BsplineSurface::get_control_point()
{
	construct_smoth_matrix();
	construct_cofficient_matrix();
 	SparseMatrix<float> K_com = K_d+alpha*K_s;
	std::cout << K_com.rows() << "   " << K_com.cols() << std::endl;
	//SparseLU<SparseMatrix<float>> solver;
	//SparseQR<SparseMatrix<float>,AMDOrdering<int>> solver;
	SimplicialLDLT<SparseMatrix<float>> solver;
	solver.compute(K_com);
	if (solver.info() != Success)
	{
		std::cout << "分解失败" << std::endl;
	}
	std::cout << "矩阵分解完成" << std::endl;
	 X = solver.solve(b);
	std::cout << "求解成功！" << std::endl;
	//std::cout << X<< std::endl;

	std::ofstream outer;
	outer.open("control_point.txt");
	for (int i = 0; i < X.rows();i++)
	{
		outer << X(i, 0) << " " << X(i, 1) << "  " << X(i, 2) << std::endl;
	}
	outer.close();
}

void BsplineSurface::remesh()
{
	get_control_point();
	output_mesh = input_mesh;
	std::tr1::shared_ptr<BasisFuns> sInv(new BasisFuns(U, V, 3, 3));
	clock_t startTime, endTime;
	startTime = clock();//计时开始
	MatrixXf new_Point = BaseMatrix*X;
	for (Mesh::VertexIter v_it = input_mesh.vertices_begin(); v_it != input_mesh.vertices_end(); v_it++)
	{
// 		Mesh::Point point = parameter.point(*v_it);
// 		int id_u = (sInv.get())->findSpan_U(point.data()[0]);
// 		int id_v = (sInv.get())->findSpan_V(point.data()[1]);
// 		/*确定基函数不为0的区间*/
// 		int left_v = id_v - k_v >= 0 ? id_v - k_v : 0;
// 		int right_v = id_v + k_v <= num_v ? id_v + k_v : num_v;
// 		int low_u = id_u - k_u >= 0 ? id_u - k_u : 0;
// 		int up_u = id_u + k_u <= num_u ? id_u + k_u : num_u;
// 		Mesh::Point new_point(0.0, 0.0, 0.0);
// 		for (int i = low_u; i <= up_u; i++)
// 		{
// 			for (int j = left_v; j <= right_v; j++)
// 			{
// 				float temp_value = (sInv.get())->basis_u_v(i, j, point.data()[0], point.data()[1]);
// 				int temp_basis_id = i*(num_v + 1) + j;
// 				Mesh::Point temp_point(X(temp_basis_id, 0), X(temp_basis_id, 1), X(temp_basis_id, 2));
// 				new_point += temp_value*temp_point;
// 			}
// 		}
		Mesh::Point Point(new_Point((*v_it).idx(), 0), new_Point((*v_it).idx(), 1), new_Point((*v_it).idx(), 2));
		output_mesh.set_point(*v_it, Point);
	}
	endTime = clock();
	std::cout << "重新计算网格的时间为: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
}

Mesh BsplineSurface::get_output_mesh()
{
	remesh();
	return output_mesh;
}