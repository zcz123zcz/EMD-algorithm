
#include "Parameterization.h"

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
	resetBoundaryVertex();
	return output_mesh;
}