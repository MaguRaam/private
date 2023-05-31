/*
 * triangulation.cc
 *      Author: sunder
 */

#include "triangulation.hh"

//----------------------------------------------------------------------------
// Constructor for the triangulation class. Takes the name of the grid file
// as input. Allocates appropriate memory and forms the grid connectivity
//----------------------------------------------------------------------------

Triangulation::Triangulation() {

	std::cout << "------------------------------------------------------------------------" << std::endl;
	std::cout << "Building triangulation" << std::endl;

	int i,k,vk,v0,v1,vk0,vk1,vk2,ck,counter;
	double x0,x1,y0,y1,xf,yf,dx,dy,temp;

	//------------------------------------------------------
	// 1) Read the vertex file
	//------------------------------------------------------

	nVrtx = 0;

	std::cout << "    >> Reading vertex file...";

	std::ifstream vrtx_file("mesh.node");

	if ( !(vrtx_file.is_open()) ) {
		std::cerr << "Error. Unable to open vertex file" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	std::string line;   // Read a line in the text file
	std::string buffer; // Read a buffer in the text file

	std::getline(vrtx_file, line);
	std::istringstream instr(line);

	// Get the number of vertices in mesh

	instr >> nVrtx;

	// Allocate memory

	vrtx = new Vertex[nVrtx];

	// Read coordinates of the vertex

	for (i = 0; i < nVrtx; ++i) {
		std::getline(vrtx_file, line);
		instr.clear();
		instr.str(line);
		instr >> buffer
			  >> vrtx[i].x[0]
		      >> vrtx[i].x[1]
		      >> vrtx[i].bcond;

		// Set the number of cells surrounding a vertex to zero

		vrtx[i].nc = 0;
	}

	vrtx_file.close();

	std::cout << " nVrtx = " << nVrtx << ". Done.";

	//------------------------------------------------------
	// 2) Read cell file
	//------------------------------------------------------

	nCell = 0;

	std::cout << "\n    >> Reading cell file...";

	std::ifstream cell_file("mesh.ele");

	if ( !(cell_file.is_open()) ) {
		std::cerr << "Error. Unable to open cell file" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// Get the number of cells

	std::getline(cell_file, line);
	instr.clear(); instr.str(line);
	instr >> nCell;

	// Allocate memory

	cell = new Cell[nCell];

	// Read connectivity of the cells

	for (i = 0; i < nCell; ++i) {
		std::getline(cell_file, line);
		instr.clear();
		instr.str(line);
		instr >> buffer
			  >> cell[i].vrtx[0]
		      >> cell[i].vrtx[1]
		      >> cell[i].vrtx[2];
	}

	cell_file.close();

	std::cout << " nCell = " << nCell << ". Done.";

	//------------------------------------------------------
	// 3) Read face file
	//------------------------------------------------------

	nFace = 0;

	std::cout << "\n    >> Reading face file...";

	std::ifstream face_file("mesh.edge");

	if ( !(face_file.is_open()) ) {
		std::cerr << "Error. Unable to open face file" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// Get the number of faces

	std::getline(face_file, line);
	instr.clear(); instr.str(line);
	instr >> nFace;

	// Allocate memory

	face = new Face[nFace];

	// Read vertices of the face and boundary conditions

	for (i = 0; i < nFace; ++i) {
		std::getline(face_file, line);
		instr.clear();
		instr.str(line);
		instr >> buffer
			  >> face[i].vrtx[0]
		      >> face[i].vrtx[1]
		      >> face[i].bcond;
	}

	face_file.close();

	std::cout << " nFace = " << nFace << ". Done.";

	//------------------------------------------------------
	// 4) Find the cells surrounding a vertex
	//------------------------------------------------------

	std::cout << "\n    >> Building vertex to cell connectivity...";

	// Find the number of cells around a vertex

	for (i = 0; i < nCell; ++i) {
		for (k = 0; k < GeometryInfo::vertices_per_cell; ++k) {
			vk = cell[i].vrtx[k];
			vrtx[vk].nc++;
		}
	}

	// Allocate the vertex-to-cell connectivity and reinitialize the number of cells surrounding a vertex

	for (i = 0; i < nVrtx; ++i) {
		vrtx[i].cell = new int[vrtx[i].nc];
		vrtx[i].nc = 0;
	}

	// Fill the cell arrays

	for (i = 0; i < nCell; ++i) {
		for (k = 0; k < GeometryInfo::vertices_per_cell; ++k) {
			vk = cell[i].vrtx[k];
			vrtx[vk].cell[vrtx[vk].nc] = i;
			vrtx[vk].nc++;
		}
	}

	std::cout << " Done.";

	//------------------------------------------------------
	// 5) Find the cells adjacent to a face
	//------------------------------------------------------

	std::cout << "\n    >> Building face to cell connectivity...";

	for (i = 0; i < nFace; ++i) {

		// Start by assuming the face is orphan -> no cells around it.

		face[i].cell[0] = -1;
		face[i].cell[1] = -1;

		v0 = face[i].vrtx[0]; v1 = face[i].vrtx[1];

		counter = 0;

		for (k = 0; k < vrtx[v0].nc; ++k) {

			ck = vrtx[v0].cell[k];

            vk0 = cell[ck].vrtx[0];
            vk1 = cell[ck].vrtx[1];
            vk2 = cell[ck].vrtx[2];

            // Check if any two vertices match (v0,v1)

            if ((v0==vk0 && v1==vk1) ||
                (v0==vk1 && v1==vk0) ||
                (v0==vk1 && v1==vk2) ||
                (v0==vk2 && v1==vk1) ||
                (v0==vk2 && v1==vk0) ||
                (v0==vk0 && v1==vk2) ) {

                face[i].cell[counter] = ck;
                counter = counter+1;
            }

		}

		// Also assign boundary tag to the face. (If the right face is missing =>)

		if (face[i].cell[1] == -1)
			face[i].at_boundary = true;
		else
			face[i].at_boundary = false;
	}

	std::cout << " Done.";

	//------------------------------------------------------
	// 6) Find the faces which make up a cell
	//------------------------------------------------------

	std::cout << "\n    >> Building cell to face connectivity...";

	for (i = 0; i < nFace; ++i) {

		// Get the vertices of the face

        v0 = face[i].vrtx[0];
        v1 = face[i].vrtx[1];

        // First check the left cell

        ck = face[i].cell[0];

        vk0 = cell[ck].vrtx[0];
        vk1 = cell[ck].vrtx[1];
        vk2 = cell[ck].vrtx[2];

        if ((v0 == vk0 && v1 == vk1) ||
            (v0 == vk1 && v1 == vk0) ) {
            cell[ck].face[0] = i;
        }

        if ((v0 == vk1 && v1 == vk2) ||
            (v0 == vk2 && v1 == vk1) ) {
            cell[ck].face[1] = i;
        }

        if ((v0 == vk2 && v1 == vk0) ||
            (v0 == vk0 && v1 == vk2) ) {
            cell[ck].face[2] = i;
        }

        // If the right cell exists, do the same for the right cell

        if (face[i].cell[1] != -1) {

            ck = face[i].cell[1];

            vk0 = cell[ck].vrtx[0];
            vk1 = cell[ck].vrtx[1];
            vk2 = cell[ck].vrtx[2];

            if ((v0 == vk0 && v1 == vk1) ||
                (v0 == vk1 && v1 == vk0) ) {
                cell[ck].face[0] = i;
            }

            if ((v0 == vk1 && v1 == vk2) ||
                (v0 == vk2 && v1 == vk1) ) {
                cell[ck].face[1] = i;
            }

            if ((v0 == vk2 && v1 == vk0) ||
                (v0 == vk0 && v1 == vk2) ) {
                cell[ck].face[2] = i;
            }

        }
	}

	std::cout << " Done.";

	//------------------------------------------------------
	// 7) Find the direct (von-Neumann) neighbours of a cell
	//------------------------------------------------------

	std::cout << "\n    >> Finding direct neighbours of a cell...";

    for (i = 0; i < nCell; ++i) {

        for (counter = 0; counter < GeometryInfo::faces_per_cell; ++counter) {

            k = cell[i].face[counter];

            if (face[k].cell[0] == i)
                cell[i].nghbr[counter] = face[k].cell[1];
            else
                cell[i].nghbr[counter] = face[k].cell[0];
        }
    }

    std::cout << " Done.";

	//------------------------------------------------------
	// 8) For a given face f, find the local index of the
    //    face in its left and right cells
	//------------------------------------------------------

    std::cout << "\n    >> Finding local face index of face in straddling cells...";

    for (i = 0; i < nFace; ++i) {

    	// Left cell

    	ck = face[i].cell[0];

    	for (k = 0; k < GeometryInfo::faces_per_cell; ++k) {
    		if (cell[ck].face[k] == i) {
    			face[i].lf[0] = k;
    			break;
    		}
    	}

    	// Right cell

    	if ( !(face[i].at_boundary) ) {

        	ck = face[i].cell[1];

        	for (k = 0; k < GeometryInfo::faces_per_cell; ++k) {
        		if (cell[ck].face[k] == i) {
        			face[i].lf[1] = k;
        			break;
        		}
        	}
    	}
    }

    std::cout << " Done.";

	//------------------------------------------------------
	// 9) Calculate some metrics related to the grid
	//------------------------------------------------------

    std::cout << "\n    >> Calculating grid metrics...";

    // a) Cell volumes and centroids

    for (i = 0; i < nCell; ++i) {

        vk0 = cell[i].vrtx[0];
        vk1 = cell[i].vrtx[1];
        vk2 = cell[i].vrtx[2];

        cell[i].vol = triangle_area(vrtx[vk0].x[0],vrtx[vk0].x[1],
        		                    vrtx[vk1].x[0],vrtx[vk1].x[1],
									vrtx[vk2].x[0],vrtx[vk2].x[1]);

        cell[i].xc[0] = (1./3.)*(vrtx[vk0].x[0] + vrtx[vk1].x[0] + vrtx[vk2].x[0]);
        cell[i].xc[1] = (1./3.)*(vrtx[vk0].x[1] + vrtx[vk1].x[1] + vrtx[vk2].x[1]);
    }

    // b) Face areas and normals

    for (i = 0; i < nFace; ++i) {

        vk0 = face[i].vrtx[0];
        vk1 = face[i].vrtx[1];

        x0 = vrtx[vk0].x[0]; y0 = vrtx[vk0].x[1];
        x1 = vrtx[vk1].x[0]; y1 = vrtx[vk1].x[1];

        face[i].area = std::hypot( (x0-x1), (y0-y1) );

        face[i].nv[0] = -(y0-y1)/face[i].area;
        face[i].nv[1] =  (x0-x1)/face[i].area;

        // Make sure that normal is always oriented from left cell to right cell

        xf = 0.5*(x0+x1); yf = 0.5*(y0+y1);

        dx = xf - cell[face[i].cell[0]].xc[0];
        dy = yf - cell[face[i].cell[0]].xc[1];

        // If the normal is not in proper direction, negate the normal

        if ( (face[i].nv[0]*dx + face[i].nv[1]*dy) < 0.0) {
            face[i].nv[0] = -face[i].nv[0];
            face[i].nv[1] = -face[i].nv[1];
        }
    }

    std::cout << " Done.";

	//------------------------------------------------------
	// 10) Calculate effective mesh spacings
	//------------------------------------------------------

    std::cout << "\n    >> Calculating effective mesh spacings...";

    heffn     = std::sqrt( 1.0/static_cast<double>(nVrtx) ); // Effective spacing based on # of nodes
    heffc     = std::sqrt( 1.0/static_cast<double>(nCell) ); // Effective spacing based on # of cells
    heffv     = std::sqrt( cell[0].vol );
    heffv_min = std::sqrt( cell[0].vol );
    heffv_max = std::sqrt( cell[0].vol );
    heffi_min = incircle_diameter(0);
    heffi_max = incircle_diameter(0);

    for (i = 1; i < nCell; ++i) {
    	temp = std::sqrt( cell[i].vol );
        heffv     = heffv + temp;
        heffv_min = std::min( heffv_min, temp );
        heffv_max = std::max( heffv_max, temp );
        temp = incircle_diameter(i);
        heffi_min = std::min( heffi_min, temp );
        heffi_max = std::max( heffi_max, temp );
    }

    heffv = heffv/static_cast<double>(nCell);

    std::cout << " heffi_min = " << heffi_min << " .Done.";

	std::cout << "\nDone." << std::endl;
	std::cout << "------------------------------------------------------------------------" << std::endl;
}

//----------------------------------------------------------------------------
// Destructor, free all the memory
//----------------------------------------------------------------------------

Triangulation::~Triangulation() {

	delete[] vrtx;
	delete[] cell;

}

//----------------------------------------------------------------------------
// Calculate in-circle diameter of ith cell
//----------------------------------------------------------------------------

double Triangulation::incircle_diameter(const int& iCell) const {

	double x1 = vrtx[cell[iCell].vrtx[0]].x[0]; double y1 = vrtx[cell[iCell].vrtx[0]].x[1];
	double x2 = vrtx[cell[iCell].vrtx[1]].x[0]; double y2 = vrtx[cell[iCell].vrtx[1]].x[1];
	double x3 = vrtx[cell[iCell].vrtx[2]].x[0]; double y3 = vrtx[cell[iCell].vrtx[2]].x[1];

	// use heron's formula

	double a = std::hypot(x2-x1, y2-y1);
	double b = std::hypot(x3-x2, y3-y2);
	double c = std::hypot(x1-x3, y1-y3);
	double s = 0.5*(a + b + c);

	double r = std::sqrt(s*(s-a)*(s-b)*(s-c))/s;

	return 2.0*r;
}


//----------------------------------------------------------------------------
// Plot triangulation in VTK format.
//----------------------------------------------------------------------------

void Triangulation::plot_triangulation(const double* U) const {
    std::ofstream vtk;
    const std::string filename = "grid.vtk";
    vtk.open (filename);
    vtk.flags( std::ios::dec | std::ios::scientific );
    vtk.precision(6);

    if ( !( vtk.is_open() ) ) {
        std::cerr << "Error. Unable to open file: " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    vtk << "# vtk DataFile Version 3.0" << "\n";
    vtk << "Unstructured Mesh" << "\n";
    vtk << "ASCII" << "\n";
    vtk << "\nDATASET UNSTRUCTURED_GRID" << "\n";

    // Point data

    vtk << "POINTS " << nVrtx << " double" << "\n";

    for (int v = 0; v < nVrtx ; ++v)
        vtk << vrtx[v].x[0] << " " << vrtx[v].x[1] << " " << 0.0 << "\n";

    // Cell connectivity

    vtk << "CELLS " << nCell << " " << 4*nCell  << "\n";

    for (int c = 0; c < nCell; ++c) {

        vtk << 3 << " ";

        for (int local_v = 0; local_v < GeometryInfo::vertices_per_cell; ++local_v )
            vtk << cell[c].vrtx[local_v] << " ";

        vtk << "\n";
    }

    // Cell type - 5 for triangles

    vtk << "CELL_TYPES " << nCell << "\n";

    for (int c = 0; c < nCell; ++c)
        vtk << 5 << "\n";

    vtk << "CELL_DATA " << nCell << "\n";
    vtk << "SCALARS U double 1" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";

    for (int c = 0; c < nCell; ++c)
    	vtk << U[c] << std::endl;

    vtk.close();
}

//----------------------------------------------------------------------------
// Given the three coordinates of a triangle, find its area
//----------------------------------------------------------------------------

double triangle_area(const double& x1, const double& y1,
						const double& x2, const double& y2,
						const double& x3, const double& y3) {

	// use shoe-lace formula

	return 0.5*std::abs(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x3*y2 - x2*y1);
}

//----------------------------------------------------------------------------
// Given a sector with vertex (x0,y0) and arms (x1,y1) and (x2,y2) check if
// a point (x,y) lies inside the sector (angle is always taken to be less
// than 180 degrees)
//----------------------------------------------------------------------------

bool point_lies_in_sector(double x0,double y0, double x1, double y1, double x2, double y2, double x, double y) {

	double startAngle = std::atan2((y1-y0),(x1-x0));
	double endAngle   = std::atan2((y2-y0),(x2-x0));
	double angle      = std::atan2((y-y0),(x-x0));


	if (startAngle < 0.0) {
		startAngle = 2.0*M_PI + startAngle;
	}

	if (endAngle < 0.0) {
		endAngle = 2.0*M_PI + endAngle;
	}

	if (angle < 0.0) {
		angle = 2.0*M_PI + angle;
	}

	// Make them in ascending order


	if (startAngle > endAngle) {
		double temp = startAngle;
		startAngle = endAngle;
		endAngle = temp;
	}

	if (endAngle - startAngle < M_PI) {
		if (angle >= startAngle && angle <= endAngle)
			return true;
		else
			return false;
	}

	else {
		if (angle >= startAngle && angle >= endAngle)
			return true;
		else if (angle <= startAngle && angle <= endAngle)
			return true;
		else
			return false;

	}
}



