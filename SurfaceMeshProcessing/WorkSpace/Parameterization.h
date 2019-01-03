
#pragma once
#ifndef PARAMETERIZATION
#define PARAMETERIZATION

#include "MeshViewer/MeshDefinition.h"
#include <vector>

class ShapePreserveParameter
{
public:
	ShapePreserveParameter(Mesh input_mesh_) :input_mesh(input_mesh_), first_boundary_vertex(-1){};
	~ShapePreserveParameter(){};
private:
	bool getFirstBoundaryVertex();
	void getBoundaryVertexArray();
	void resetBoundaryVertex();
public:
	std::vector<Mesh::VertexHandle> getBoundaryVertex();
	Mesh get_output_mesh();
private:
	Mesh input_mesh;
	Mesh output_mesh;
	Mesh::VertexHandle first_boundary_vertex;
	std::vector<Mesh::VertexHandle> boundary_vertex_array;
};

#endif // !PARAMETERIZATION