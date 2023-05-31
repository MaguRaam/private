/*
 * triangulation.hh
 *      Author: sunder
 */

#ifndef TRIANGULATION_HH_
#define TRIANGULATION_HH_

#include "headers.hh"

struct GeometryInfo {

	static const int dim = 2;
	static const int vertices_per_cell = 3;
    static const int vertices_per_face = 3;
    static const int faces_per_cell = 3;

};

/* Structures defining various mesh entities */

struct Vertex {
	double x[2];       // Coordinates of the vertex
	int nc;            // Number of cells sharing the vertex
	int *cell;         // List of cells sharing the vertex
	int bcond;         // Boundary condition at the vertex
	~Vertex() {
		delete[] cell;
	}
};

struct Face {
	int vrtx[2];       // List of vertices that make the face
	int cell[2];       // Two cells sharing the face: Left cell (0), Right cell (1)
	int lf[2];         // Local face index of the cell sharing the face
	double nv[2];      // Unit normal, always from left to right
	double area;       // Area of the face
	int bcond;         // Boundary condition assigned to the face
	bool at_boundary;  // Whether the face is at boundary or not
};

struct Cell {
	int vrtx[3];       // List of vertices that make the cell
	int nghbr[3];      // List of face neighbours
	int face[3];       // List of faces that makeup the cell
	double xc[2];      // Centroid of the cell
	double vol;        // Volume of the cell
};

/* Triangulation structure for holding mesh related data */

struct Triangulation {

	int nVrtx;         // No. of vertices in the mesh
	int nCell;         // No. of cells in the mesh
	int nFace;         // No. of faces in the mesh

	Vertex* vrtx;      // Data related to the vertices in the mesh
	Face*   face;      // Data related to the faces in the mesh
	Cell*   cell;      // Data related to the cells in the mesh

    double heffn;      // Effective spacing based on # of nodes
    double heffc;      // Effective spacing based on # of cells
    double heffv;      // Average of sqrt(volume).
	double heffv_min;  // Minimum sqrt(volume).
    double heffv_max;  // Maximum sqrt(volume).
    double heffi_min;  // Minimum incircle diameter in the mesh.
    double heffi_max;  // Maximum incircle diameter in the mesh.

    /* Constructor and destructor */

    Triangulation ();  // Constructor
    ~Triangulation();  // Destructor


    double incircle_diameter(const int&) const;
    void plot_triangulation(const double*) const;


};

// Additional geometry related functions

double triangle_area(const double&, const double&,
                    const double&, const double&,
                    const double&, const double&);


double area_of_face(const double&, const double&,
                    const double&, const double&);

bool point_lies_in_sector(double,double, double, double, double, double, double, double);


#endif /* TRIANGULATION_HH_ */
