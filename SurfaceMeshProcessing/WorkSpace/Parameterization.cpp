
#include "Parameterization.h"

ShapePreserveParameter::ShapePreserveParameter(Mesh input_mesh_) :input_mesh(input_mesh_), first_boundary_vertex(-1)
{
	std::cout << "开始" << std::endl;
	solve_and_reset_mesh();
}

bool ShapePreserveParameter::getFirstBoundaryVertex()
{
	typedef Mesh::Point Point;
	Mesh::VertexIter vIt = input_mesh.vertices_begin();
	Mesh::VertexIter vEnd = input_mesh.vertices_end();
	for (; vIt != vEnd;vIt++)
	{
		if (input_mesh.is_boundary(vIt))
		{
			first_boundary_vertex = vIt;
			return true;
		}
	}
	return false;
}

void ShapePreserveParameter::getBoundaryVertexArray()
{
	bool is_open_mesh = getFirstBoundaryVertex();
	if (is_open_mesh)
	{
		boundary_vertex_array.clear();
		uint size = input_mesh.n_vertices();
		std::vector<int> status(size, -1);
		Mesh::VertexHandle current_boundary_vertex = first_boundary_vertex;
		boundary_vertex_array.push_back(current_boundary_vertex);
		status[current_boundary_vertex.idx()] = 1;
		Mesh::HalfedgeHandle current_halfedge;
		while (true)
		{

			for (Mesh::VertexOHalfedgeIter voh_it = input_mesh.voh_iter(current_boundary_vertex); voh_it; ++voh_it)
			{
				if (input_mesh.is_boundary(*voh_it))
				{
					current_halfedge = *voh_it;
					break;
				}
			}
			Mesh::VertexHandle next_boundary_vertex = input_mesh.to_vertex_handle(current_halfedge);
			if (status[next_boundary_vertex.idx()] == 1)
			{
				break;
			}
			else
			{
				current_boundary_vertex = next_boundary_vertex;
				boundary_vertex_array.push_back(current_boundary_vertex);
				status[current_boundary_vertex.idx()] = 1;
			}

		}

	}
}

std::vector<Mesh::VertexHandle> ShapePreserveParameter::getBoundaryVertex()
{
	getBoundaryVertexArray();
	return boundary_vertex_array;
}

void ShapePreserveParameter::resetBoundaryVertex()
{
	output_mesh = input_mesh;
	getBoundaryVertexArray();
	int num_of_boundary_vertex = boundary_vertex_array.size();//边界点的总数量
	int increment = num_of_boundary_vertex / 4;
	float increment_length = 1.0f / float(increment);
	/*第一条边*/
	for (int i = 0; i < increment;i++)
	{
		float new_x = i*increment_length;
		Mesh::Point new_point(new_x, 0.0, 0.0);
		output_mesh.set_point(boundary_vertex_array[i], new_point);
	}

	/*第二条边*/
	for (int i = increment; i < 2 * increment; i++)
	{
		float new_y = (i-increment)*increment_length;
		Mesh::Point new_point(1.0, new_y, 0.0);
		output_mesh.set_point(boundary_vertex_array[i], new_point);
	}

	/*第三条边*/
	for (int i = 2 * increment; i < 3 * increment; i++)
	{
		float new_x = 1.0 - (i - 2 * increment)*increment_length;
		Mesh::Point new_point(new_x, 1.0, 0.0);
		output_mesh.set_point(boundary_vertex_array[i], new_point);
	}

	/*第四条边*/
	int new_increment = num_of_boundary_vertex - 3 * increment;
	increment_length = 1.0f / float(new_increment);
	for (int i = 3 * increment; i < num_of_boundary_vertex;i++)
	{
		float new_y = 1.0 - (i - 3* increment)*increment_length;
		Mesh::Point new_point(0.0, new_y, 0.0);
		output_mesh.set_point(boundary_vertex_array[i], new_point);
	}
}

Mesh ShapePreserveParameter::get_output_mesh()
{
	//resetBoundaryVertex();
	return output_mesh;
}

std::vector<Mesh::Point> ShapePreserveParameter::local_parametrization(Mesh::VertexHandle mid_vertex)
{
	float epsilon = 1e-9, pi = 3.1415926;
	std::vector<float> angle_array;//每个点对应的内角
	angle_array.push_back(0.0);//特殊处理第一个角
	std::vector<float> length_array;//每条边的长度
	float theta = 0.0f;//三维内角和
	std::vector<Mesh::Point> plane_vertex_array;
	Mesh::Point current_vertex(0.0, 0.0, 0.0);//中心点
	plane_vertex_array.push_back(current_vertex);

	for (Mesh::VertexOHalfedgeIter vh_it = input_mesh.voh_begin(mid_vertex); vh_it != input_mesh.voh_end(mid_vertex);vh_it++)
	{
		float current_length = input_mesh.calc_edge_length(*vh_it);
		length_array.push_back(current_length);
		Mesh::HalfedgeHandle other_halfedge = input_mesh.next_halfedge_handle(input_mesh.next_halfedge_handle(*vh_it));//与当前边成夹角的边
		Mesh::Point p1 = input_mesh.point(mid_vertex);
		Mesh::Point p2 = input_mesh.point(input_mesh.to_vertex_handle(*vh_it));
		Mesh::Point p3 = input_mesh.point(input_mesh.from_vertex_handle(other_halfedge));

		float temp_angle = acos(((p2 - p1) | (p3 - p1)) / ((p2 - p1).norm()*(p3 - p1).norm()));
		angle_array.push_back(temp_angle), theta += temp_angle;
	}

	float temp_theta = 0.0;
	for (int i = 0; i < angle_array.size() - 1;i++)
	{
		angle_array[i] = 2 * pi*angle_array[i]/ (theta + epsilon);
		temp_theta += angle_array[i];
		Mesh::Point temp_point(cos(temp_theta), sin(temp_theta), 0.0);
		plane_vertex_array.push_back(temp_point);
	}
	return plane_vertex_array;
}

void ShapePreserveParameter::parametrization()
{
	for (Mesh::VertexIter v_it = input_mesh.vertices_begin(); v_it != input_mesh.vertices_end();v_it++)
	{
		if (!input_mesh.is_boundary(*v_it))
		{
			std::vector<Mesh::Point> plane_point_array = local_parametrization(*v_it);
			for (int i = 0; i < plane_point_array.size();i++)
			{
				std::cout << plane_point_array[i] << std::endl;
			}
			std::cout << std::endl;
		}
	}
}

/*判断两条平面直线是否相交，且交点位于第二条直线两个顶点之间，利用直线的参数形式来计算*/
/*即判断交点的参数t2=((x2-x1)(y3-y1)-(x3-x1)(y2-y1))/((x4-x3)(y2-y1)-(x2-x1)(y4-y3))是否在0-1之间*/
bool ShapePreserveParameter::judge_intersect(Mesh::Point Line1_start_point, Mesh::Point Line1_end_point, Mesh::Point Line2_start_point, Mesh::Point Line2_end_point)
{
	float epsilon = 1e-8;
	Mesh::Point vec1 = Line1_end_point - Line1_start_point;
	Mesh::Point vec2 = Line2_end_point - Line2_start_point;
	if ((vec1%vec2).norm()<epsilon)//无交点
	{
		//std::cout << "无交点" << std::endl;
		return false;
	}
	if ((Line1_start_point - Line2_start_point).norm()<=1e-8 || (Line1_start_point - Line2_end_point).norm()<=1e-8)//第二种特殊情况，即有点重合
	{
		//std::cout << "顶点重合" << std::endl;
		return false;
	}
	float x1 = Line1_start_point.data()[0], y1 = Line1_start_point.data()[1];
	float x2 = Line1_end_point.data()[0], y2 = Line1_end_point.data()[1];
	float x3 = Line2_start_point.data()[0], y3 = Line2_start_point.data()[1];
	float x4 = Line2_end_point.data()[0], y4 = Line2_end_point.data()[1];
	float t = ((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)) / ((x4 - x3)*(y2 - y1) - (x2 - x1)*(y4 - y3) + epsilon);
	//std::cout << t << std::endl;
	if (t>=0&&t<=1)
	{
		return true;
	}
	return false;
}

float ShapePreserveParameter::area_of_triangle(Mesh::Point p1, Mesh::Point p2, Mesh::Point p3)
{
	Mesh::Point vec1 = p2 - p1, vec2 = p3 - p1;
	float area = 0.5*(vec2%vec1).norm();
	return area;
}

void ShapePreserveParameter::cal_barycentric_coordinate(Mesh::Point p1, Mesh::Point p2, Mesh::Point p3, Mesh::Point mid_point, Mesh::Point &coordinate)
{
	float total_area = area_of_triangle(p1, p2, p3)+(1e-9);
	float area1 = area_of_triangle(mid_point, p2, p3);
	float area2 = area_of_triangle(p1, mid_point, p3);
	float area3 = area_of_triangle(p1, p2, mid_point);
	coordinate = Mesh::Point(area1/total_area, area2/total_area, area3/total_area);
}

void ShapePreserveParameter::get_lamda_of_vertex(Mesh::VertexHandle mid_vertex,std::vector<float> &lamda_array)
{
	lamda_array.clear();
	std::vector<Mesh::Point> plane_vertex_array = local_parametrization(mid_vertex);
	plane_vertex_array.push_back(plane_vertex_array[1]);//把第一个邻域点再此加入，防止数组越界
	Mesh::Point mid_point = input_mesh.point(mid_vertex);
	int numOfNeighbor = plane_vertex_array.size() - 2;
	MatrixXf Miu(numOfNeighbor, numOfNeighbor);
	Miu.setZero();
	for (int i = 1; i <=numOfNeighbor;i++)//对每一个领域点，装配局部系数矩阵Miu的一行
	{
		int k = -1;
		for (int j = 1; j <= numOfNeighbor;j++)//找到该邻域点对应的边，详情见论文
		{
			bool intersect = judge_intersect(plane_vertex_array[i], plane_vertex_array[0], plane_vertex_array[j], plane_vertex_array[j + 1]);
			if (intersect)
			{
				k = j;
				break;
			}
		}
		Mesh::Point coordinate(0.0, 0.0, 0.0);
		//std::cout << numOfNeighbor << "   " << k << "   " << plane_vertex_array.size()<<std::endl;
		cal_barycentric_coordinate(plane_vertex_array[i], plane_vertex_array[k], plane_vertex_array[k + 1], plane_vertex_array[0],coordinate);
		//std::cout << "死于这里"<< std::endl;
		Miu(i - 1, i - 1) = coordinate.data()[0], Miu(i - 1, k - 1) = coordinate.data()[1];
		if (k==numOfNeighbor)
		{
			Miu(i - 1, 0) = coordinate.data()[2];
		}
		else
		{
			Miu(i - 1, k) = coordinate.data()[2];
		}
	}
	//std::cout << Miu << std::endl;
	/*再把每一列相加取平均*/
	for (int i = 0; i < numOfNeighbor;i++)
	{
		float temp_lamda = 0.0;
		for (int j = 0; j < numOfNeighbor;j++)
		{
			temp_lamda += Miu(j, i);
		}
		temp_lamda /= float(numOfNeighbor);
		lamda_array.push_back(temp_lamda);
	}

	float sum = 0;
	for (int it = 0; it < lamda_array.size() - 1; it++)
	{
		sum += lamda_array[it];
	}
	lamda_array.back() = 1.0 - sum;

	for (int i = 0; i < lamda_array.size();i++)
	{
		if (lamda_array[i]<=0.0||lamda_array[i]>=1.0)
		{
			std::cout << "系数越界！"<<std::endl;
		}
	}
}

void ShapePreserveParameter::get_lamda_of_vertex_average(Mesh::VertexHandle mid_vertex, std::vector<float> &lamda_array)
{
	lamda_array.clear();
	int num = input_mesh.valence(mid_vertex);

	for (Mesh::VertexVertexIter it = input_mesh.vv_begin(mid_vertex); it != input_mesh.vv_end(mid_vertex); it++)
	{
		lamda_array.push_back(1.0 / float(num));
	}
	for (int i = 0; i < lamda_array.size(); i++)
	{
		if (lamda_array[i] <= 0.0 || lamda_array[i] >= 1.0)
		{
			std::cerr << "系数越界！" << std::endl;
		}
	}
}

void ShapePreserveParameter::find_interior_point(std::vector<Mesh::VertexHandle> &interior_vertex_handle,std::vector<int> &new_vertex_idex)
{
	int num = input_mesh.n_vertices();
	new_vertex_idex = std::vector<int>(num, -1);
	int temp_id = 0;
	interior_vertex_handle.clear();
	for (Mesh::VertexIter v_it = input_mesh.vertices_begin(); v_it != input_mesh.vertices_end();v_it++)
	{
		if (!input_mesh.is_boundary(*v_it))
		{
			new_vertex_idex[(*v_it).idx()] = temp_id;
			temp_id++;
			interior_vertex_handle.push_back(*v_it);
		}
	}
}

void ShapePreserveParameter::construct_equation(SparseMatrix<float> &L, MatrixXf &b_x, MatrixXf &b_y, std::vector<Mesh::VertexHandle> &interior_vertex_handle)
{
	std::vector<int> new_vertex_idex;
	std::cout << "找边界点" << std::endl;
	find_interior_point(interior_vertex_handle,new_vertex_idex);
	/*先固定边界点*/
	std::cout << "装配矩阵" << std::endl;
	resetBoundaryVertex();
	/*装配矩阵*/
	int size = interior_vertex_handle.size();
	L = SparseMatrix<float>(size, size);
	b_x=MatrixXf(size, 1); b_x.setZero();
	b_y = MatrixXf(size, 1); b_y.setZero();
	std::vector<int> row;
	std::vector<int> col;
	std::vector<float> value;
	std::vector<float> lamda_array;
	for (int i = 0; i < size;i++)
	{
		get_lamda_of_vertex(interior_vertex_handle[i], lamda_array);
		//get_lamda_of_vertex_average(interior_vertex_handle[i], lamda_array);
		row.push_back(i);
		col.push_back(i);
		value.push_back(1.0);
		/*对于每个点的邻点进行判断，如果是边界点就放入右端项，否则放入系数矩阵*/
		int id = 0;
		for (Mesh::VertexOHalfedgeIter voh_it = output_mesh.voh_begin(interior_vertex_handle[i]); voh_it != output_mesh.voh_end(interior_vertex_handle[i]);voh_it++)
		{
 			Mesh::VertexHandle temp_vertex = output_mesh.to_vertex_handle(*voh_it);
			if (output_mesh.is_boundary(temp_vertex))//放入右端项
			{
				b_x(i,0) += lamda_array[id] * output_mesh.point(temp_vertex).data()[0]; 
				b_y(i,0)+= lamda_array[id] * output_mesh.point(temp_vertex).data()[1];
			}
			else
			{
				row.push_back(i);
				col.push_back(new_vertex_idex[temp_vertex.idx()]);
				if (new_vertex_idex[temp_vertex.idx()] == i)
				{
					std::cout << "出错了！" << std::endl;
				}
				value.push_back(-lamda_array[id]);
			}

			id++;
		}
	}
	std::cout << "装配完成" << std::endl;
	std::vector<Triplet<float>> triplets;
	for (int i = 0; i < row.size();i++)
	{
		triplets.emplace_back(row[i], col[i], value[i]);
	}
	L.setFromTriplets(triplets.begin(),triplets.end());
}

void ShapePreserveParameter::solve_and_reset_mesh()
{
	SparseMatrix<float> L;
	MatrixXf b_x, b_y;
	std::vector<Mesh::VertexHandle> interior_vertex_handle;
	construct_equation(L, b_x, b_y, interior_vertex_handle);
	SparseLU<SparseMatrix<float>> solver;
	//SparseQR<SparseMatrix<float>,AMDOrdering<int>> solver;
	//SimplicialLDLT<SparseMatrix<float>> solver;
	solver.compute(L);
	if (solver.info() != Success)
	{
		std::cout << "分解失败" << std::endl;
	}
	std::cout << "参数化矩阵分解完成" << std::endl;
	MatrixXf x = solver.solve(b_x);
	MatrixXf y = solver.solve(b_y);
	std::cout << "参数化求解完成" << std::endl;
	for (int i = 0; i < interior_vertex_handle.size();i++)
	{
		Mesh::Point new_point(x(i, 0), y(i, 0), 0.0);
		output_mesh.set_point(interior_vertex_handle[i], new_point);
	}
}
