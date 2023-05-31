/*
 * hype.cc
 *      Author: sunder
 */

#include "hype.hh"

//----------------------------------------------------------------------------
// Constructor for the HyPE class
//----------------------------------------------------------------------------

HyPE::HyPE() {

	std::cout << "                 A simple 2D unstructured WENO code" << std::endl;
	std::cout << "                      Written by Dasika Sunder " << std::endl;
	std::cout << "\n    >> Allocating memory and initializing solution...";

	int i, q, f, iCell, iVar;
	double X0,X1,X2,Y0,Y1,Y2,x,y;
	Array1D<double> q0(nVar);

	const int NGP_Vol = 6; // STRANG5, order 6, degree of precision 4 (for initializng the solution)

	double xiGP[]  = {0.816847572980459,0.091576213509771,0.091576213509771,0.108103018168070,0.445948490915965,0.445948490915965};
	double etaGP[] = {0.091576213509771,0.816847572980459,0.091576213509771,0.445948490915965,0.108103018168070,0.445948490915965};
	double wGP[]   = {0.109951743655322,0.109951743655322,0.109951743655322,0.223381589678011,0.223381589678011,0.223381589678011};

	N = 2;
	nDOF = 6;

	// Allocate memory

	stencil = new Stencil[tria.nCell];
	uh.reinit(tria.nCell, nVar, nDOF);
	wh.reinit(tria.nVrtx, nVar);
	uh0.reinit(tria.nCell, nVar);
	duh.reinit(tria.nCell,nVar);
	lim.reinit(tria.nCell);

	// Set time related data

	CFL = 0.4;
	time = 0.0;
	tend = 0.2;
	timestep = 0;
	write_interval = 500;
	dt = 0.0;
	rk_stage = 1;

	// Set quadrature related data

	NGP_face = 2;

	xiGP_Face.reinit(GeometryInfo::faces_per_cell, NGP_face);
	etaGP_Face.reinit(GeometryInfo::faces_per_cell, NGP_face);
	wGP_Face.reinit(GeometryInfo::faces_per_cell, NGP_face);
	phiGP_Face.reinit(GeometryInfo::faces_per_cell,NGP_face, nDOF);

	Array1D<double> xGP_1D(NGP_face), wGP_1D(NGP_face);
	QGauss(0.0, 1.0, xGP_1D, wGP_1D, NGP_face); // Gauss quadrature points and weights in [0,1]

	// Face 0 (0,0)--------(1,0)

	f = 0;

	for (q= 0; q < NGP_face; ++q) {
		xiGP_Face(f,q)  = xGP_1D(q);
		etaGP_Face(f,q) = 0.0;
		wGP_Face(f,q)   = wGP_1D(q);
	}

	// Face 1 (1,0) ------ (0,1)

	f = 1;

	for (q= 0; q < NGP_face; ++q) {
		straight_face_to_inclined_face(xGP_1D(q), xiGP_Face(f,q), etaGP_Face(f,q));
		wGP_Face(f,q)   = wGP_1D(q);
	}

	// Face 2 (1,0) ------ (0,0)

	f = 2;

	for (q= 0; q < NGP_face; ++q) {
		xiGP_Face(f,q)  = 0.0;
		etaGP_Face(f,q) = xGP_1D(NGP_face - q -1);
		wGP_Face(f,q)   = wGP_1D(q);
	}

	// Evaluate basis functions at the quadrature points

	for (f = 0; f < GeometryInfo::faces_per_cell; ++f) {
		for (q = 0; q < NGP_face; ++q) {
			for (i = 0; i < nDOF; ++i) {
				phiGP_Face(f,q,i) = basis(xiGP_Face(f,q), etaGP_Face(f,q), i);
			}
		}
	}

	// Initialize nodal values and their basis functions (for positivity preserving)

	N_Node = 6;

	xiNode.reinit(N_Node); etaNode.reinit(N_Node); phiNode.reinit(N_Node,nDOF);

	xiNode[0] = 0.0; etaNode[0] = 0.0;
	xiNode[1] = 0.5; etaNode[1] = 0.0;
	xiNode[2] = 1.0; etaNode[2] = 0.0;
	xiNode[3] = 0.0; etaNode[3] = 0.5;
	xiNode[4] = 0.5; etaNode[4] = 0.5;
	xiNode[5] = 0.0; etaNode[5] = 1.0;

	for (q = 0; q < N_Node; ++q) {
		for (i = 0; i < nDOF; ++i) {
			phiNode(q,i) = basis(xiNode[q], etaNode[q], i);
		}
	}

	// Initialize solution with the initial condition

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

    	// Vertices of the cell

    	X0 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[0]; Y0 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[1];
    	X1 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[0]; Y1 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[1];
    	X2 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[0]; Y2 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[1];

		for (iVar = 0; iVar < nVar; ++iVar)
			uh(iCell,iVar,0) = 0.0;

		for (q = 0; q < NGP_Vol; ++ q) {

    		reference_to_physical(X0, Y0, X1, Y1, X2, Y2, xiGP[q], etaGP[q], x, y);

    		icond(x,y,q0);

    		for (iVar = 0; iVar < nVar; ++iVar)
    			uh(iCell,iVar,0) += wGP[q]*q0[iVar];
    	}
	}

	std::cout << " Done.";
}

//----------------------------------------------------------------------------
// Destructor: Free all the memory in the class
//----------------------------------------------------------------------------

HyPE::~HyPE() {


	delete[] stencil;
}

//----------------------------------------------------------------------------
// Select stencils using sectorial search algorithm
//----------------------------------------------------------------------------

void HyPE::select_stencils() {

	std::cout << "\n    >> Selecting stencils...";

	int i, j, k, c, iCell, iVrtx, count;
	int* int_store = new int[200];
	std::set<int> lstencil;
	bool already_counted;
	double x0,y0,x1,y1,x2,y2,x,y;

	//------------------------------------------------------
	// 1) Select stencils
	//------------------------------------------------------

	// 1a) Find primary neighbours

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

		for (k =0 ; k < GeometryInfo::vertices_per_cell; ++k) {

			iVrtx = tria.cell[iCell].vrtx[k];

			for (j = 0; j < tria.vrtx[iVrtx].nc; ++j) {

				if (tria.vrtx[iVrtx].cell[j] != iCell) {
					lstencil.insert(tria.vrtx[iVrtx].cell[j]);
				}

			}
		}

        stencil[iCell].nPrimary = lstencil.size();
        stencil[iCell].Primary = new int[stencil[iCell].nPrimary];

        count = 0;

        for (std::set<int>::iterator it = lstencil.begin(); it != lstencil.end(); it++) {
            stencil[iCell].Primary[count] =  *it;
            count++;
        }

        lstencil.clear();
	}

	// 1b) Select secondary neighbours

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

    	// Second layer

    	for (i = 0; i < stencil[iCell].nPrimary; ++i) {

    		c = stencil[iCell].Primary[i];

    		for (j = 0; j < stencil[c].nPrimary; ++j) {

				already_counted = false;

    			for (k = 0; k < stencil[iCell].nPrimary; ++k) {

    				if (stencil[c].Primary[j] == iCell || stencil[c].Primary[j] == stencil[iCell].Primary[k])
    					already_counted = true;

    			}

    			if (!already_counted) {
    				lstencil.insert(stencil[c].Primary[j]);
    			}
    		}
    	}

        stencil[iCell].nSecondary = lstencil.size();
        stencil[iCell].Secondary = new int[stencil[iCell].nSecondary];

        count = 0;

        for (std::set<int>::iterator it = lstencil.begin(); it != lstencil.end(); it++) {
            stencil[iCell].Secondary[count] =  *it;
            count++;
        }

        lstencil.clear();
	}

	// 1c) Select tertiary neighbours

	// For WENO3 we do not require tertiary neighbors, so I am commenting out this part to save storage space.
	// I will add it if required for say WENO5

	/*
	for (iCell = 0; iCell < tria.nCell; ++iCell) {

		count = 0;

    	for (i = 0; i < stencil[iCell].nSecondary; ++i) {

    		c = stencil[iCell].Secondary[i];

    		for (int j = 0; j < stencil[c].nPrimary; ++j) {

				already_counted = false;

    			for (k = 0; k < stencil[iCell].nSecondary; ++k) {

    				if (stencil[c].Primary[j] == iCell || stencil[c].Primary[j] == stencil[iCell].Secondary[k])
    					already_counted = true;

    			}

    			for (k = 0; k < stencil[iCell].nPrimary; ++k) {

    				if (stencil[c].Primary[j] == iCell || stencil[c].Primary[j] == stencil[iCell].Primary[k])
    					already_counted = true;

    			}

    			if (!already_counted) {
    				lstencil.insert(stencil[c].Primary[j]);
    			}
    		}
    	}

        stencil[iCell].nTertiary = lstencil.size();
        stencil[iCell].Tertiary = new int[stencil[iCell].nTertiary];


        count = 0;

        for (std::set<int>::iterator it = lstencil.begin(); it != lstencil.end(); it++) {
            stencil[iCell].Tertiary[count] =  *it;
            count++;
        }

    	lstencil.clear();
	}
	*/

	// 1d) Select stencils (using sectorial search)

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

		x0 = tria.cell[iCell].xc[0]; y0 = tria.cell[iCell].xc[1]; // Centroid of the triangle

		// a) Central stencil

		count = 0;

		// Add all the primary neighbours

		for (int i = 0; i < stencil[iCell].nPrimary; ++i) {
			int_store[count] = stencil[iCell].Primary[i];
			count++;
		}

		// Add the direct neighbours of primary neighbours

    	for (i = 0; i < stencil[iCell].nPrimary; ++i) {

    		c = stencil[iCell].Primary[i];

    		for (j = 0; j < GeometryInfo::faces_per_cell; ++j) {

				already_counted = false;

    			for (k = 0; k < stencil[iCell].nPrimary; ++k) {

    				if (tria.cell[c].nghbr[j] == iCell ||
    					tria.cell[c].nghbr[j] == stencil[iCell].Primary[k] ||
    					tria.cell[c].nghbr[j] == -1) {

    					already_counted = true;
    				}
    			}

    			if (!already_counted) {
    				int_store[count] = tria.cell[c].nghbr[j];
    				count++;
    			}
    		}
    	}

    	stencil[iCell].nCentral = count;
		stencil[iCell].Central = new int[stencil[iCell].nCentral];

		for (i = 0; i < stencil[iCell].nCentral; ++i) {
			stencil[iCell].Central[i] = int_store[i];
		}

		// b) Sectorial stencil 1

		count = 0;

		x1 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[0]; y1 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[1];
		x2 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[0]; y2 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[1];

		// Search among the primary neighbours

		for (i = 0; i < stencil[iCell].nPrimary; ++i) {

			x = tria.cell[stencil[iCell].Primary[i]].xc[0]; y = tria.cell[stencil[iCell].Primary[i]].xc[1];

			if ( point_lies_in_sector(x0,y0,x1,y1,x2,y2,x,y) ) {

				int_store[count] = stencil[iCell].Primary[i];
				count++;
			}
		}

		// Search among secondary neighbours

		for (int i = 0; i < stencil[iCell].nSecondary; ++i) {

			x = tria.cell[stencil[iCell].Secondary[i]].xc[0]; y = tria.cell[stencil[iCell].Secondary[i]].xc[1];

			if ( point_lies_in_sector(x0,y0,x1,y1,x2,y2,x,y) ) {

				int_store[count] = stencil[iCell].Secondary[i];
				count++;
			}
		}

		stencil[iCell].nSector1 = count;
		stencil[iCell].Sector1 = new int[stencil[iCell].nSector1];

		for (int i = 0; i < stencil[iCell].nSector1; ++i) {
			stencil[iCell].Sector1[i] = int_store[i];
		}

		// c) Sectorial stencil 2

		count = 0;

		x1 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[0]; y1 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[1];
		x2 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[0]; y2 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[1];

		// Search among the primary neighbours

		for (i = 0; i < stencil[iCell].nPrimary; ++i) {

			x = tria.cell[stencil[iCell].Primary[i]].xc[0]; y = tria.cell[stencil[iCell].Primary[i]].xc[1];

			if ( point_lies_in_sector(x0,y0,x1,y1,x2,y2,x,y) ) {

				int_store[count] = stencil[iCell].Primary[i];
				count++;
			}
		}

		// Search among secondary neighbours

		for (int i = 0; i < stencil[iCell].nSecondary; ++i) {

			x = tria.cell[stencil[iCell].Secondary[i]].xc[0]; y = tria.cell[stencil[iCell].Secondary[i]].xc[1];

			if ( point_lies_in_sector(x0,y0,x1,y1,x2,y2,x,y) ) {

				int_store[count] = stencil[iCell].Secondary[i];
				count++;
			}
		}

		stencil[iCell].nSector2 = count;
		stencil[iCell].Sector2 = new int[stencil[iCell].nSector2];

		for (int i = 0; i < stencil[iCell].nSector2; ++i) {
			stencil[iCell].Sector2[i] = int_store[i];
		}

		// d) Sectorial stencil 3

		count = 0;

		x1 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[0]; y1 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[1];
		x2 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[0]; y2 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[1];

		// Search among the primary neighbours

		for (i = 0; i < stencil[iCell].nPrimary; ++i) {

			x = tria.cell[stencil[iCell].Primary[i]].xc[0]; y = tria.cell[stencil[iCell].Primary[i]].xc[1];

			if ( point_lies_in_sector(x0,y0,x1,y1,x2,y2,x,y) ) {

				int_store[count] = stencil[iCell].Primary[i];
				count++;
			}
		}

		// Search among secondary neighbours

		for (int i = 0; i < stencil[iCell].nSecondary; ++i) {

			x = tria.cell[stencil[iCell].Secondary[i]].xc[0]; y = tria.cell[stencil[iCell].Secondary[i]].xc[1];

			if ( point_lies_in_sector(x0,y0,x1,y1,x2,y2,x,y) ) {

				int_store[count] = stencil[iCell].Secondary[i];
				count++;
			}
		}

		stencil[iCell].nSector3 = count;
		stencil[iCell].Sector3  = new int[stencil[iCell].nSector3];

		for (int i = 0; i < stencil[iCell].nSector3; ++i) {
			stencil[iCell].Sector3[i] = int_store[i];
		}

	}

    delete[] int_store;

    // To visualize stencils for a cell look_at, un-comment this part

  

    double *u = new double[tria.nCell];
    int look_at = 3750;

    for (iCell = 0; iCell < tria.nCell; ++iCell) {
    	if (iCell == look_at)
    		u[iCell] = 1.0;
    	else
    		u[iCell] = 0.0;
    }

    for (int l  = 0; l < stencil[look_at].nSector1; ++l) {
    	u[stencil[look_at].Sector1[l]] = 2.0;
    }

    for (int l  = 0; l < stencil[look_at].nSector2; ++l) {
    	u[stencil[look_at].Sector2[l]] = 3.0;
    }

    for (int l  = 0; l < stencil[look_at].nSector3; ++l) {
    	u[stencil[look_at].Sector3[l]] = 4.0;
    }

    tria.plot_triangulation(u);

    delete[] u;

	std::cout << " Done.";
}

//----------------------------------------------------------------------------
// Compute the QR decomposition matrices for each cell and each of its
// stencil and store them. These matrices are completely geometry dependent
// and hence they need to be computed only once for static meshes.
//----------------------------------------------------------------------------

void HyPE::compute_reconstruction_matrices() {

	std::cout << "\n    >> Pre-computing stencil matrices ...";

	int i, j, c, iCell;
	double x0,y0,x1,y1,x2,y2;
	double X0,X1,X2,Y0,Y1,Y2;
	double xi0,xi1,xi2,eta0,eta1,eta2;
	const int min_cells_in_sector = 7;  // For third order scheme 5 is sufficient, but I am putting 7 considering the
                                        // well-conditioning of the matrix formed, as suggested in Dumber, Kaser et. al.
	Matrix A_central, A_S1, A_S2, A_S3, A_S4;

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

		stencil[iCell].ntSectors = 0;

    	// Vertices of target cell

    	X0 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[0]; Y0 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[1];
    	X1 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[0]; Y1 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[1];
    	X2 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[0]; Y2 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[1];

    	// a) -----------------------------  Central stencil -----------------------------

    	A_central.reinit(stencil[iCell].nCentral, nDOF-1);

    	for (i = 0; i < stencil[iCell].nCentral; ++i) {

    		// Vertices of cell in stencil (in physical coordinates)

    		c = stencil[iCell].Central[i];

        	x0 = tria.vrtx[tria.cell[c].vrtx[0]].x[0]; y0 = tria.vrtx[tria.cell[c].vrtx[0]].x[1];
        	x1 = tria.vrtx[tria.cell[c].vrtx[1]].x[0]; y1 = tria.vrtx[tria.cell[c].vrtx[1]].x[1];
        	x2 = tria.vrtx[tria.cell[c].vrtx[2]].x[0]; y2 = tria.vrtx[tria.cell[c].vrtx[2]].x[1];

        	// Convert to reference coordinates

        	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x0,y0,xi0,eta0);
        	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x1,y1,xi1,eta1);
        	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x2,y2,xi2,eta2);

        	for (j = 1; j < nDOF; ++j) {
        		A_central(i,j-1) = basis_average_over_triangle(xi0,eta0,xi1,eta1,xi2,eta2,j);
        	}
    	}

    	stencil[iCell].QR_central.initialize(A_central);

    	if (!(stencil[iCell].QR_central.is_full_rank())) {
    		std::cerr << "Central stencil matrix for cell " << iCell << " is rank deficient." << std::endl;
    		std::cerr << "I cannot proceed further. Getting out." << std::endl;
    		std::exit(EXIT_FAILURE);
    	}

    	// b) -----------------------------  Sector 1 -----------------------------

    	if (stencil[iCell].nSector1 < min_cells_in_sector) {
    		stencil[iCell].is_admissible_S1 = false;
    	}

    	else {

    		A_S1.reinit(stencil[iCell].nSector1, nDOF-1);

        	for (i = 0; i < stencil[iCell].nSector1; ++i) {

        		// Vertices of cell in stencil (in physical coordinates)

        		c = stencil[iCell].Sector1[i];

            	x0 = tria.vrtx[tria.cell[c].vrtx[0]].x[0]; y0 = tria.vrtx[tria.cell[c].vrtx[0]].x[1];
            	x1 = tria.vrtx[tria.cell[c].vrtx[1]].x[0]; y1 = tria.vrtx[tria.cell[c].vrtx[1]].x[1];
            	x2 = tria.vrtx[tria.cell[c].vrtx[2]].x[0]; y2 = tria.vrtx[tria.cell[c].vrtx[2]].x[1];

            	// Convert to reference coordinates

            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x0,y0,xi0,eta0);
            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x1,y1,xi1,eta1);
            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x2,y2,xi2,eta2);

            	for (j = 1; j < nDOF; ++j) {
            		A_S1(i,j-1) = basis_average_over_triangle(xi0,eta0,xi1,eta1,xi2,eta2,j);
            	}
        	}

        	stencil[iCell].QR_S1.initialize(A_S1);

        	if (stencil[iCell].QR_S1.is_full_rank()) {
        		stencil[iCell].is_admissible_S1 = true;
        		stencil[iCell].ntSectors++;
        	}

        	else {
        		stencil[iCell].is_admissible_S1 = false;
        	}
    	}

    	// c) -----------------------------  Sector 2 -----------------------------

    	if (stencil[iCell].nSector2 < min_cells_in_sector) {
    		stencil[iCell].is_admissible_S2 = false;
    	}

    	else {

    		A_S2.reinit(stencil[iCell].nSector2, nDOF-1);

        	for (i = 0; i < stencil[iCell].nSector2; ++i) {

        		// Vertices of cell in stencil (in physical coordinates)

        		c = stencil[iCell].Sector2[i];

            	x0 = tria.vrtx[tria.cell[c].vrtx[0]].x[0]; y0 = tria.vrtx[tria.cell[c].vrtx[0]].x[1];
            	x1 = tria.vrtx[tria.cell[c].vrtx[1]].x[0]; y1 = tria.vrtx[tria.cell[c].vrtx[1]].x[1];
            	x2 = tria.vrtx[tria.cell[c].vrtx[2]].x[0]; y2 = tria.vrtx[tria.cell[c].vrtx[2]].x[1];

            	// Convert to reference coordinates

            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x0,y0,xi0,eta0);
            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x1,y1,xi1,eta1);
            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x2,y2,xi2,eta2);

            	for (j = 1; j < nDOF; ++j) {
            		A_S2(i,j-1) = basis_average_over_triangle(xi0,eta0,xi1,eta1,xi2,eta2,j);
            	}
        	}

        	stencil[iCell].QR_S2.initialize(A_S2);

        	if (stencil[iCell].QR_S2.is_full_rank()) {
        		stencil[iCell].is_admissible_S2 = true;
        		stencil[iCell].ntSectors++;
        	}

        	else {
        		stencil[iCell].is_admissible_S2 = false;
        	}
    	}

    	// d) -----------------------------  Sector 3 -----------------------------

    	if (stencil[iCell].nSector3 < min_cells_in_sector) {
    		stencil[iCell].is_admissible_S3 = false;
    	}

    	else {

    		A_S3.reinit(stencil[iCell].nSector3, nDOF-1);

        	for (i = 0; i < stencil[iCell].nSector3; ++i) {

        		// Vertices of cell in stencil (in physical coordinates)

        		c = stencil[iCell].Sector3[i];

            	x0 = tria.vrtx[tria.cell[c].vrtx[0]].x[0]; y0 = tria.vrtx[tria.cell[c].vrtx[0]].x[1];
            	x1 = tria.vrtx[tria.cell[c].vrtx[1]].x[0]; y1 = tria.vrtx[tria.cell[c].vrtx[1]].x[1];
            	x2 = tria.vrtx[tria.cell[c].vrtx[2]].x[0]; y2 = tria.vrtx[tria.cell[c].vrtx[2]].x[1];

            	// Convert to reference coordinates

            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x0,y0,xi0,eta0);
            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x1,y1,xi1,eta1);
            	physical_to_reference(X0,Y0,X1,Y1,X2,Y2,x2,y2,xi2,eta2);

            	for (j = 1; j < nDOF; ++j) {
            		A_S3(i,j-1) = basis_average_over_triangle(xi0,eta0,xi1,eta1,xi2,eta2,j);
            	}
        	}

        	stencil[iCell].QR_S3.initialize(A_S3);

        	if (stencil[iCell].QR_S3.is_full_rank()) {
        		stencil[iCell].is_admissible_S3 = true;
        		stencil[iCell].ntSectors++;
        	}

        	else {
        		stencil[iCell].is_admissible_S3 = false;
        	}
    	}

    	// Assign linear weights for each stencil

    	if (stencil[iCell].ntSectors == 0) { // Make sure at least one sectorial stencil is present.
    		std::cerr << "\nThe cell " << iCell << " has no sectorial stencils" << std::endl;
    		std::cerr << "I cannot proceed further. Getting out." << std::endl;
    		std::exit(EXIT_FAILURE);
    	}

	} // End of cell loop

	std::cout << " Done.\n";
}

//----------------------------------------------------------------------------
// Reconstruct the solution to get polynomial coefficients
//----------------------------------------------------------------------------

double smoothness_indicator(const Vector& u) {


	double IS = 2.5*u[0]*u[0] + 3.0*u[0]*u[1] + 4.5*u[1]*u[1] + 2.0*u[1]*u[2] + 2.0*u[1]*u[3]-
			    4.0*u[1]*u[4] + 96.0*u[2]*u[2]  + 38.0*u[2]*u[4]  + 212.0*u[4]*u[4] +
				(1./3.)*(2.0*u[0]*u[2] + 10.0*u[0]*u[3] - 4.0*u[0]*u[4] + 248.0*u[2]*u[3] + 320.0*u[3]*u[3] + 614.0*u[3]*u[4]);

	//IS = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3] + u[4]*u[4];

	return IS;
}

void HyPE::reconstruct() {

	int i, k, ck, iCell, iVar;

	const double small_num = 1.0e-12;

	Vector rhs(60);
	Vector sol_Central(nDOF-1);
	Vector sol_S1(nDOF-1);
	Vector sol_S2(nDOF-1);
	Vector sol_S3(nDOF-1);

	double w_Central, w_S1, w_S2, w_S3;
	double IS_Central, IS_S1, IS_S2, IS_S3;
	double total_wt;

	const double central_wt = 0.99;
	const double stencil_wt = (1.0-central_wt)/3.0;

	Array1D<double> Qbnd(nVar);
	Array1D<double> Vbnd(nVar);

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

		for (iVar = 0; iVar < nVar; ++iVar) {

			total_wt = 0.0;

			duh(iCell,iVar) = 0.0;

			// Loop through all the cells in the stencil to form right hand side

			// a) Central stencil

			for (k = 0; k < stencil[iCell].nCentral; ++k) {
				ck = stencil[iCell].Central[k];
				rhs[k] = uh(ck,iVar,0) - uh(iCell,iVar,0);
			}

			stencil[iCell].QR_central.solve(rhs, sol_Central);

			IS_Central = smoothness_indicator(sol_Central);
			w_Central = central_wt/( pow4(IS_Central+small_num) );
			total_wt += w_Central;

			// b) Sector 1

			if (stencil[iCell].is_admissible_S1) {

				for (k = 0; k < stencil[iCell].nSector1; ++k) {
					ck = stencil[iCell].Sector1[k];
					rhs[k] = uh(ck,iVar,0) - uh(iCell,iVar,0);
				}

				stencil[iCell].QR_S1.solve(rhs, sol_S1);
			}

			else {

				for (i = 0; i < 5; ++i) {
					sol_S1[i] = sol_Central[i];
				}
			}

			IS_S1 = smoothness_indicator(sol_S1);
			w_S1 = stencil_wt/( pow4(IS_S1+small_num) );
			total_wt += w_S1;

			// c) Sector 2

			if (stencil[iCell].is_admissible_S2) {

				for (k = 0; k < stencil[iCell].nSector2; ++k) {
					ck = stencil[iCell].Sector2[k];
					rhs[k] = uh(ck,iVar,0) - uh(iCell,iVar,0);
				}

				stencil[iCell].QR_S2.solve(rhs, sol_S2);
			}

			else {
				for (i = 0; i < 5; ++i) {
					sol_S2[i] = sol_Central[i];
				}
			}


			IS_S2 = smoothness_indicator(sol_S2);
			w_S2 = stencil_wt/( pow4(IS_S2+small_num) );
			total_wt += w_S2;

			// d) Sector 3

			if (stencil[iCell].is_admissible_S3) {

				for (k = 0; k < stencil[iCell].nSector3; ++k) {
					ck = stencil[iCell].Sector3[k];
					rhs[k] = uh(ck,iVar,0) - uh(iCell,iVar,0);
				}

				stencil[iCell].QR_S3.solve(rhs, sol_S3);
			}

			else {
				for (i = 0; i < 5; ++i) {
					sol_S3[i] = sol_Central[i];
				}
			}

			IS_S3 = smoothness_indicator(sol_S3);
			w_S3 = stencil_wt/( pow4(IS_S3+small_num) );
			total_wt += w_S3;

			// Normalize the weights

			w_Central = w_Central/total_wt;
			w_S1 = w_S1/total_wt;
			w_S2 = w_S2/total_wt;
			w_S3 = w_S3/total_wt;

			// Store in the solution variable

			for (i = 1; i < nDOF; ++i) {
				uh(iCell,iVar,i) = w_Central*sol_Central[i-1] + w_S1*sol_S1[i-1] +  w_S2*sol_S2[i-1] + w_S3*sol_S3[i-1];
			}

		} // Variable loop

		lim[iCell].is_corrupt = false;

		for (int f = 0; f < 3; ++f ) {

			for (int q = 0; q < NGP_face; ++q) {

				for (int iVar = 0; iVar < nVar; ++iVar) {
					Qbnd[iVar] = 0.0;
				}

				for (int iVar = 0; iVar < nVar; ++iVar) {

					for (int i = 0; i < nDOF; ++i) {
						Qbnd[iVar] += uh(iCell,iVar,i)*phiGP_Face(f,q,i);
					}
				}

				PDECons2Prim(Qbnd, Vbnd);

				if (Vbnd[0] < rho_floor || Vbnd[3] < prs_floor) {
					lim[iCell].is_corrupt  = true;
				}

			}
		}

		if (lim[iCell].is_corrupt) {
			for (iVar = 0; iVar < nVar; ++iVar) {
				for (i = 1; i < nDOF; ++i) {
					uh(iCell,iVar,i) = 0.0;
				}
			}
		}

	} // Cell loop

}

//----------------------------------------------------------------------------
// Preserve Density/Pressure positivity in cell.
// Reference: Self-adjusting, positivity preserving high order schemes for
// hydrodynamics and magneto-hydrodynamics. JCP 2012, D. Balsara
//----------------------------------------------------------------------------

// Solve pressure positivity equation


double solve_pressure_positivity_equation(double a2, double a1, double a0) {

	// Evaluate discriminant

	const double small_num = 1.0e-12;
	double dis = a1*a1 - 4.0*a2*a0;

	// Find the two roots of the equation

	double root1, root2;

	if (std::abs(a2) <= small_num) { // Equal roots

		if ( std::abs(a1) > small_num)
			root1 = -a0/a1;
		else
			root1 = 2.0;

		root2 = 2.0;
	}

	else if (dis >= 0.0) { // Two real roots
		root1 = 0.5*(-a1 + std::sqrt(dis))/a2;
        root2 = 0.5*(-a1 - std::sqrt(dis))/a2;
	}

	else { // Complex roots, set arbitrary values
		root1 = 2.0;
		root2 = 2.0;
	}

	// Now obtain root by polling the one or two roots.

	double x_root = 2.0;

	if ( (0.0 <= root1) && (root1 <= 1.0) )
		x_root  = std::min(x_root, root1);

	if ( (0.0 <= root2) && (root2 <= 1.0) )
	  x_root = std::min (x_root, root2);

	if ( x_root > 1.0)
		x_root = 0.0;

	return x_root;
}

void HyPE::preserve_positivity() {

	int iCell,iNb,Nb,iFace,iNode,iDOF,iVar;
	double irho,rho,rhou,rhov,u,v,E,u_nb,v_nb,nx,ny,temp;
	double del_rhou, del_rhov;
    double rho_min, rho_max, prs_min,e_min_ext, flatten_max;
	double rho_min_ext, rho_max_ext, prs_min_ext;
	double rho_node, prs_node, x_root,reduce_modes;
	double a0,a1,a2;
	const double kappa      = 0.4;
	const double expand_fac = 0.4;
	Array1D<double> Q_Node(nVar), V_Node(nVar);

	// 1) Find the divergence of velocity, pressure and speed of sound in each cell

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

		rho = uh(iCell,0,0);
		irho = 1.0/rho;
		u = irho*uh(iCell,1,0); v = irho*uh(iCell,2,0);

		lim[iCell].prs = (GAMMA -1.0)*( uh(iCell,3,0) - 0.5*rho*(u*u + v*v) );
		lim[iCell].c = std::sqrt(lim[iCell].prs*GAMMA*irho);

		lim[iCell].div_v = 0.0;

		for (iNb = 0; iNb < GeometryInfo::faces_per_cell; ++iNb) {

			if (tria.cell[iCell].nghbr[iNb] != -1) {
				Nb = tria.cell[iCell].nghbr[iNb];

                u_nb = uh(Nb,1,0)/uh(Nb,0,0);
                v_nb = uh(Nb,2,0)/uh(Nb,0,0);

                iFace = tria.cell[iCell].face[iNb];

                nx = tria.face[iFace].nv[0]; ny = tria.face[iFace].nv[1];

                if (tria.face[iFace].cell[1] == iCell) {
                	nx = -nx; ny = -ny;
                }

                lim[iCell].div_v += tria.face[iFace].area*( (u_nb-u)*nx + (v-v_nb)*ny );

			}

		}

		lim[iCell].div_v = lim[iCell].div_v/tria.cell[iCell].vol;

	}  // End of cell loop

	// 2) Find the minimum speed of sound in each cell and the flattener variable

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

		lim[iCell].cmin = lim[iCell].c;

        for (iNb = 0; iNb < stencil[iCell].nPrimary; ++iNb) {

        	Nb = stencil[iCell].Primary[iNb];

            if (lim[Nb].c < lim[iCell].cmin) {
            	lim[iCell].cmin = lim[Nb].c;
            }

        }

        temp = -(lim[iCell].div_v/kappa*lim[iCell].cmin + 1.0);

        lim[iCell].flatten = std::min(1.0, std::max (0.0, temp));
	}

	// 3) Limit the solution in each cell

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

        reduce_modes = 1.0;
        lim[iCell].is_corrupt = false;

        // Find the minimum and maximum of density, pressure and flattener in the cell neighbourhood

        rho  = uh(iCell,0,0);
        rhou = uh(iCell,1,0);
        rhov = uh(iCell,2,0);
        E    = uh(iCell,3,0);
        flatten_max = lim[iCell].flatten;
        rho_min     = uh(iCell,0,0);
        rho_max     = uh(iCell,0,0);
        prs_min     = lim[iCell].prs;
        //prs_max     = lim[iCell].prs;

        for (iNb = 0; iNb < stencil[iCell].nPrimary; ++iNb) {

        	Nb = stencil[iCell].Primary[iNb];

        	// Maximum flattener in the neighbourhood

            if (lim[Nb].flatten > flatten_max)
            	flatten_max = lim[Nb].flatten;

            // Minimum density in the neighbourhood

            if (uh(Nb,0,0) < rho_min)
            	rho_min = uh(Nb,0,0);

            // Maximum density in the neighbourhood

            if (uh(Nb,0,0) > rho_max)
            	rho_max = uh(Nb,0,0);

        	// Minimum pressure in the neighbourhood

            if (lim[Nb].prs < prs_min)
            	prs_min = lim[Nb].prs;

            /*
            // Maximum pressure in the neighbourhood

            if (lim[Nb].prs > prs_max)
            	prs_max = lim[Nb].prs;
            */
        }

        // Extend ranges

        rho_min_ext = rho_min*(1.0 - expand_fac + expand_fac*flatten_max);
        rho_max_ext = rho_max*(1.0 + expand_fac - expand_fac*flatten_max);

        prs_min_ext = prs_min*(1.0 - expand_fac + expand_fac*flatten_max);
        //prs_max_ext = prs_max*(1.0 + expand_fac - expand_fac*flatten_max);

        if (rho_min_ext < rho_floor)
            rho_min_ext = rho_min;

        if (prs_min_ext < prs_floor)
            prs_min_ext = prs_min;

        // Bring density within range

        for (iNode = 0; iNode < N_Node; ++iNode) {

        	rho_node = 0.0;

			for (iDOF = 0; iDOF < nDOF; ++iDOF)
				 rho_node += uh(iCell,0,iDOF)*phiNode(iNode,iDOF);

            if (rho_node > rho_max_ext*1.000001) {
                lim[iCell].is_corrupt = true;
                temp = (rho_max_ext-rho)/(rho_node-rho);
                reduce_modes = std::min(reduce_modes, temp);
            }

            if (rho_node < rho_min_ext*0.999999) {
                lim[iCell].is_corrupt = true;
                temp = (rho_min_ext-rho)/(rho_node-rho);
                reduce_modes = std::min(reduce_modes, temp);
            }

        }

        if (reduce_modes < 0.0)
        	reduce_modes = 0.0;

        if (reduce_modes > 1.0)
        	reduce_modes = 1.0;

        for (iVar = 0; iVar < nVar; ++iVar) {
        	for (iDOF = 1; iDOF < nDOF; ++iDOF) {
        		uh(iCell,iVar,iDOF) = reduce_modes*uh(iCell,iVar,iDOF);
        	}
        }

        // Now check for pressure positivity

        reduce_modes = 1.0;
        e_min_ext = prs_min_ext/(GAMMA-1.0);

        for (iNode = 0; iNode < N_Node; ++iNode) {

        	for (iVar = 0; iVar < nVar; ++iVar)
        		Q_Node[iVar] = 0.0;

        	for (iVar = 0; iVar < nVar; ++iVar)
        		for (iDOF = 0; iDOF < nDOF; ++iDOF)
        			Q_Node[iVar] += uh(iCell,iVar,iDOF)*phiNode(iNode,iDOF);

        	PDECons2Prim(Q_Node, V_Node);

        	prs_node = V_Node[3];

        	if (prs_node < prs_min_ext) {

        		// Solve the pressure positivity equation

        		lim[iCell].is_corrupt = true;

        		temp = rhou*rhou + rhov*rhov;
        		del_rhou = Q_Node[1] - rhou;
        		del_rhov = Q_Node[2] - rhov;

        		a0 = 2.0*rho*E - temp - 2.0*rho*e_min_ext;
        		a1 = 2.0*rho*(Q_Node[3]-E) + 2.0*E*(Q_Node[0]-rho) - 2.0*(del_rhou*rhou + del_rhov*rhov) - 2.0*e_min_ext*(Q_Node[0]-rho);
        		a2 = 2.0*(Q_Node[0] - rho)*(Q_Node[3] - E) - (del_rhou*del_rhou + del_rhov*del_rhov);

        		x_root = solve_pressure_positivity_equation(a2, a1, a0);

        		reduce_modes = std::min(reduce_modes,0.95*x_root);

        	}
        }

        if (reduce_modes < 0.0)
        	reduce_modes = 0.0;

        if (reduce_modes > 1.0)
        	reduce_modes = 1.0;

        for (iVar = 0; iVar < nVar; ++iVar) {
        	for (iDOF = 1; iDOF < nDOF; ++iDOF) {
        		uh(iCell,iVar,iDOF) = reduce_modes*uh(iCell,iVar,iDOF);
        	}
        }

	} // End of cell loop


}


//----------------------------------------------------------------------------
// Find upwind flux on each face of the mesh and add its contribution to rhs
//----------------------------------------------------------------------------

void HyPE::compute_rhs() {

	int iFace, q, iDOF, iVar, LCell, RCell,lfl, lfr;

	Array2D<double> QL_GP(NGP_face,nVar);
	Array2D<double> QR_GP(NGP_face,nVar);
	Array2D<double> F_GP(NGP_face,nVar);
	Array1D<double> QL(nVar), QR(nVar), Flux(nVar);

	Array1D<double> Vb(nVar);

	// Free stream conditions

	Vb[0] = 8.0;
	Vb[1] = 8.25;
	Vb[2] = 0.0;
	Vb[3] = 116.5;

	double s, smax = 0.0;

	// Loop over all the faces of the triangulation

	for (iFace = 0; iFace < tria.nFace; ++iFace) {

		LCell = tria.face[iFace].cell[0];
		lfl = tria.face[iFace].lf[0];

		RCell = tria.face[iFace].cell[1];
		lfr = tria.face[iFace].lf[1];

		//-----------------------------------------------------------------------

		// Left state

		for (q = 0; q < NGP_face; ++ q) {

			for (iVar = 0; iVar < nVar; ++ iVar) {
				QL_GP(q,iVar) = 0.0;
			}


			for (iVar = 0; iVar < nVar; ++iVar) {

				for (iDOF = 0; iDOF < nDOF; ++iDOF) {
					QL_GP(q,iVar) += uh(LCell,iVar,iDOF)*phiGP_Face(lfl,q,iDOF);
				}

			}
		}

		//-----------------------------------------------------------------------

		if (tria.face[iFace].at_boundary) {

			// Boundary state

			for (q = 0; q < NGP_face; ++ q) {

				for (iVar = 0; iVar < nVar; ++ iVar) {
					QL(iVar) = QL_GP(q,iVar);
				}

				get_right_state(QL, Vb, tria.face[iFace].bcond, tria.face[iFace].nv[0], tria.face[iFace].nv[1], QR);

				for (iVar = 0; iVar < nVar; ++ iVar)
					QR_GP(q,iVar) = QR(iVar);

			}
		}

		else { // Internal faces

			// Right state

			for (q = 0; q < NGP_face; ++ q) {

				for (iVar = 0; iVar < nVar; ++ iVar)
					QR_GP(q,iVar) = 0.0;

				for (iVar = 0; iVar < nVar; ++iVar)
					for (iDOF = 0; iDOF < nDOF; ++iDOF)
						QR_GP(q,iVar) += uh(RCell,iVar,iDOF)*phiGP_Face(lfr,NGP_face-q-1,iDOF);
			}
		}

		// Solve Riemann problem at each quadrature point

		for (q = 0; q < NGP_face; ++q) {

			for (iVar = 0; iVar < nVar; ++iVar) {
				QL[iVar] = QL_GP(q,iVar);
				QR[iVar] = QR_GP(q,iVar);
			}

			s = LLFRiemannSolver(QL, QR, tria.face[iFace].nv[0], tria.face[iFace].nv[1], Flux); if (s>smax) smax = s;

			for (iVar = 0; iVar < nVar; ++iVar) {
				F_GP(q,iVar) = Flux[iVar];
			}
		}

		//-----------------------------------------------------------------------

		// Add the contribution to rhs

        for (q = 0; q < NGP_face; ++q) { // Subtract -> because the normal is outward facing
            for (iVar = 0; iVar < nVar; ++ iVar) {
            	duh(LCell,iVar) -= (tria.face[iFace].area*wGP_Face(lfl,q)*F_GP(q,iVar))/tria.cell[LCell].vol;
            }
        }

        if (! (tria.face[iFace].at_boundary) ) {

            for (q = 0; q < NGP_face; ++q) { // Add -> because the normal is inward facing
                for (iVar = 0; iVar < nVar; ++ iVar) {
                	duh(RCell,iVar) += (tria.face[iFace].area*wGP_Face(lfr,NGP_face-q-1)*F_GP(NGP_face-q-1,iVar))/tria.cell[RCell].vol;
                }
            }

        }
	}

	// Set time step

	if (rk_stage == 1) {
		dt = 0.5*CFL*tria.heffi_min/smax;
	}
}

/*
//  First order - for debugging

void HyPE::compute_rhs() {

	int iFace, iVar, LCell, RCell;

	Array1D<double> QL(nVar), QR(nVar), Flux(nVar);

	Array1D<double> Vb(nVar);

	// Free stream conditions

	Vb[0] = 1.4;
	Vb[1] = 3.0;
	Vb[2] = 0.0;
	Vb[3] = 1.0;

	double s, smax = 0.0;

	// Loop over all the faces of the triangulation

	for (iFace = 0; iFace < tria.nFace; ++iFace) {

		LCell = tria.face[iFace].cell[0];
		RCell = tria.face[iFace].cell[1];

		//-----------------------------------------------------------------------

		// Left state

		for (iVar = 0; iVar < nVar; ++iVar) {
			QL[iVar] = uh(LCell,iVar,0);
		}

		//-----------------------------------------------------------------------

		if (tria.face[iFace].at_boundary) {

			// Boundary state

			get_right_state(QL, Vb, tria.face[iFace].bcond, tria.face[iFace].nv[0], tria.face[iFace].nv[1], QR);
		}

		else { // Internal faces

			// Right state

			for (iVar = 0; iVar < nVar; ++iVar) {
				QR[iVar] = uh(RCell,iVar,0);
			}
		}

		// Solve Riemann problem at each quadrature point

		s = LLFRiemannSolver(QL, QR, tria.face[iFace].nv[0], tria.face[iFace].nv[1], Flux); if (s>smax) smax = s;

		//-----------------------------------------------------------------------

		// Add the contribution to rhs

		for (iVar = 0; iVar < nVar; ++iVar) {
			duh(LCell,iVar) -= tria.face[iFace].area*Flux[iVar]/tria.cell[LCell].vol;
		}


        if (! (tria.face[iFace].at_boundary) ) {

        	for (iVar = 0; iVar < nVar; ++iVar) {
        		duh(RCell,iVar) += tria.face[iFace].area*Flux[iVar]/tria.cell[RCell].vol;
        	}
        }

	} // End of face loop


	// Set time step

	if (rk_stage == 1) {
		dt = 0.5*CFL*tria.heffi_min/smax;
	}
}
*/

//----------------------------------------------------------------------------
// Update one-step of the solution using SSPRK(3,3) method
//----------------------------------------------------------------------------

void HyPE::step_SSPRK33() {

	int iCell, iVar;

	// Stage 1

	rk_stage = 1;
	reconstruct();
	//preserve_positivity();
	compute_rhs();

	for (iCell = 0; iCell < tria.nCell; ++iCell) {
		for (iVar = 0; iVar < nVar; ++iVar) {
            uh0(iCell,iVar)  = uh(iCell,iVar,0);
            uh(iCell,iVar,0) = uh0(iCell,iVar) + dt*duh(iCell,iVar);
		}
	}

	// Stage 2

	rk_stage = 2;
	reconstruct();
	//preserve_positivity();
	compute_rhs();

	for (iCell = 0; iCell < tria.nCell; ++iCell) {
		for (iVar = 0; iVar < nVar; ++iVar) {
            uh(iCell,iVar,0) =  0.25*(3.0*uh0(iCell,iVar) + uh(iCell,iVar,0) + dt*duh(iCell,iVar));
		}
	}

	// Stage 3

	rk_stage = 3;
	reconstruct();
	//preserve_positivity();
	compute_rhs();

	for (iCell = 0; iCell < tria.nCell; ++iCell) {
		for (iVar = 0; iVar < nVar; ++iVar) {
            uh(iCell,iVar,0) =  (1./3.)*(uh0(iCell,iVar) + 2.0*uh(iCell,iVar,0) + 2.0*dt*duh(iCell,iVar));
		}
	}

}

//----------------------------------------------------------------------------
// Plot solution in vtk format
//----------------------------------------------------------------------------

void HyPE::plot() {

	const unsigned int digits = 8;

	const std::string filename = "sol-" + int_to_string (timestep, digits) + ".vtk";

    std::ofstream vtk;
    vtk.open (filename);
    vtk.flags( std::ios::dec | std::ios::scientific );
    vtk.precision(6);
    Array1D<double> Q(nVar), V(nVar);

    if ( !( vtk.is_open() ) ) {
        std::cerr << "Error. Unable to open file: " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::cout << "Writing data to " << filename << std::endl;

    // Obtain primitive variables at mesh vertices

    for (int v = 0; v < tria.nVrtx; ++v) {

		for (int iVar = 0; iVar < nVar; ++iVar) {
			Q[iVar] = 0.0;
		}

    	for (int c = 0; c < tria.vrtx[v].nc; ++c) {

    		for (int iVar = 0; iVar < nVar; ++iVar) {
    			Q[iVar] += uh(tria.vrtx[v].cell[c],iVar,0)/static_cast<double>(tria.vrtx[v].nc);
    		}
    	}

    	PDECons2Prim(Q, V);

		for (int iVar = 0; iVar < nVar; ++iVar) {
			wh(v,iVar) = V[iVar];
		}
    }

    vtk << "# vtk DataFile Version 3.0" << "\n";
    vtk << "Unstructured Mesh" << "\n";
    vtk << "ASCII" << "\n";
    vtk << "\nDATASET UNSTRUCTURED_GRID" << "\n";
    vtk << "\nFIELD FieldData 1" << "\n";
    vtk << "TIME 1 1 double" << "\n";
    vtk << time << "\n";

    // Point data

    vtk << "POINTS " << tria.nVrtx << " double" << "\n";

    for (int v = 0; v < tria.nVrtx ; ++v)
        vtk << tria.vrtx[v].x[0] << " " << tria.vrtx[v].x[1] << " " << 0.0 << "\n";

    // Cell connectivity

    vtk << "CELLS " << tria.nCell << " " << 4*tria.nCell  << "\n";

    for (int c = 0; c < tria.nCell; ++c) {

        vtk << 3 << " ";

        for (int local_v = 0; local_v < GeometryInfo::vertices_per_cell; ++local_v )
            vtk << tria.cell[c].vrtx[local_v] << " ";

        vtk << "\n";
    }

    // Cell type - 5 for triangles

    vtk << "CELL_TYPES " << tria.nCell << "\n";

    for (int c = 0; c < tria.nCell; ++c)
        vtk << 5 << "\n";

    vtk << "POINT_DATA " << tria.nVrtx << "\n";
    vtk << "SCALARS Density double 1" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";

    for (int v = 0; v < tria.nVrtx; ++v)
    	vtk << wh(v,0) << std::endl;

    vtk << "\n";

    vtk << "SCALARS Pressure double" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";
    for (int v = 0; v < tria.nVrtx; ++v)
        vtk << wh(v,3) << "\n";

    vtk << "\n";

    vtk << "VECTORS Velocity double" << "\n";
    for (int v = 0; v < tria.nVrtx; ++v)
        vtk << wh(v,1) << " " << wh(v,2) << " " << 0.0 << "\n";

    vtk << "CELL_DATA " << tria.nCell << "\n";
    vtk << "SCALARS Tcells double 1" << "\n";
    vtk << "LOOKUP_TABLE default" << "\n";

    for (int c = 0; c < tria.nCell; ++c) {
    	if (lim[c].is_corrupt)
    		vtk << 1.0 << std::endl;
    	else
    		vtk << 0.0 << std::endl;
    }

    vtk.close();
}

//----------------------------------------------------------------------------
// Find the maximum error in the solution compared to the exact solution
//----------------------------------------------------------------------------

double HyPE::max_error() const {

	double max_norm = 0.0;
	double error;
	int q, iCell, iVar;
	double X0, Y0, X1, Y1, X2, Y2, x, y;
	const int NGP_Vol = 6; // STRANG5, order 6, degree of precision 4 (for initializng the solution)

	double xiGP[]  = {0.816847572980459,0.091576213509771,0.091576213509771,0.108103018168070,0.445948490915965,0.445948490915965};
	double etaGP[] = {0.091576213509771,0.816847572980459,0.091576213509771,0.445948490915965,0.108103018168070,0.445948490915965};
	double wGP[]   = {0.109951743655322,0.109951743655322,0.109951743655322,0.223381589678011,0.223381589678011,0.223381589678011};

	double q0[nVar];
	double qexact[nVar];

	for (iCell = 0; iCell < tria.nCell; ++iCell) {

    	// Vertices of the cell

    	X0 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[0]; Y0 = tria.vrtx[tria.cell[iCell].vrtx[0]].x[1];
    	X1 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[0]; Y1 = tria.vrtx[tria.cell[iCell].vrtx[1]].x[1];
    	X2 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[0]; Y2 = tria.vrtx[tria.cell[iCell].vrtx[2]].x[1];

		for (iVar = 0; iVar < nVar; ++iVar)
			qexact[iVar] = 0.0;

    	for (q = 0; q < NGP_Vol; ++ q) {

    		reference_to_physical(X0, Y0, X1, Y1, X2, Y2, xiGP[q], etaGP[q], x, y);

    		exact_solution(x,y,tend,q0);

    		for (iVar = 0; iVar < nVar; ++iVar) {
    			qexact[iVar] += wGP[q]*q0[iVar];
    		}
    	}

    	error = std::abs(qexact[0] - uh(iCell,0,0));

    	if (error > max_norm)
    		max_norm = error;
	}

	return max_norm;
}

//----------------------------------------------------------------------------
// Put everything together and run the problem
//----------------------------------------------------------------------------

void HyPE::run() {

	clock_t start, end;
	double cpu_time_used;

	select_stencils();
	compute_reconstruction_matrices();

	// This is just to set the initial time step

	reconstruct();
	//preserve_positivity();
	compute_rhs();

	start = clock();
	
	while (time < tend) {

		if (time+dt > tend)
			dt = tend-time;

		if (timestep % write_interval == 0)
			plot();

		step_SSPRK33();
		printf ("t = %4.3e\n", time);
		time+= dt;
		timestep++;
	}

	printf ("t = %4.3e\n", time);

	printf ("Max-Error = %4.7e\n", max_error());

	plot();

	end = clock();
	cpu_time_used = ((double) (end - start))/CLOCKS_PER_SEC;

	std::cout << "CPU time used (s) = " << cpu_time_used << std::endl;
	std::cout << "Finished successfully. Good bye." << std::endl;
	std::cout << "------------------------------------------------------------------------" << std::endl;
}

