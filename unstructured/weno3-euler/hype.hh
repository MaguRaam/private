/*
 * hype.hh
 *      Author: sunder
 */

#ifndef HYPE_HH_
#define HYPE_HH_

#include "la.hh"
#include "array.hh"
#include "mappings.hh"
#include "utilities.hh"
#include "quadrature.hh"
#include "triangulation.hh"

/* Structure storing all the information about stencils */

struct Stencil {

    int nPrimary;          // Number of primary neighbours (which have vertex common with the target cell)
    int *Primary;          // List of primary neighbours

    int nSecondary;        // Number of secondary neighbours (which have vertex common with the primary neighbours)
    int *Secondary;        // List of secondary neighbours

    int nTertiary;         // Number of tertiary neighbours (which have vertex common with the secondary neighbours)
    int *Tertiary;         // List of tertiary neighbours

    int nCentral;          // Number of cells in the central stencil
    int *Central;          // List of all the cells in the central stencil
    QRdcmp QR_central;     // Least-squares solution for central stencil

    int nSector1;          // Number of cells in the first sector (formed by centroid, vertex 1 and vertex 2)
    int *Sector1;          // List of cells in the first sector
    bool is_admissible_S1; // Whether Sector 1 stencil is admissible
    QRdcmp QR_S1;          // Least-squares solution for first sectorial stencil

    int nSector2;          // Number of cells in the second sector (formed by centroid, vertex 2 and vertex 3)
    int *Sector2;          // List of cells in the second sector
    bool is_admissible_S2; // Whether Sector 2 stencil is admissible
    QRdcmp QR_S2;          // Least-squares solution for first sectorial stencil

    int nSector3;          // Number of cells in the third sector (formed by centroid, vertex 3 and vertex 1)
    int *Sector3;          // List of cells in the third sector
    bool is_admissible_S3; // Whether Sector 3 stencil is admissible
    QRdcmp QR_S3;          // Least-squares solution for first sectorial stencil

    int ntSectors;         // Total number admissible sectors present

    ~Stencil() {
        delete[] Primary;
        delete[] Secondary;
        //delete[] Tertiary;
        delete[] Central;
        delete[] Sector1;
        delete[] Sector2;
        delete[] Sector3;
    }

};

/* Structure storing data related to the limiter */

struct Limiter {
	double c, cmin;    // Speed of sound in a cell and its minimum in the neighbourhood
	double prs;        // Pressure in the cell
	double div_v;      // Divergence of velocity the cell
	double flatten;    // Flattener in the cell
	bool is_corrupt;   // Whether the cell is corrupt
};

/* Main class for solving the Hyperbolic PDE */

class HyPE {

	/* Data related to the solution */

	int N;                // Degree of approximation
	int nDOF;             // Number of degrees of freedom
	Stencil* stencil;     // Stencils and associated data
	Triangulation tria;   // Mesh object
	Array3D<double> uh;   // Conserved variables along with their spatial gradients
	Array2D<double> wh;   // Primitive variables at mesh vertices
	Array2D<double> uh0;  // Copy of solution vector for using in time stepping
	Array2D<double> duh;  // RHS for updating the solution
	Array1D<Limiter> lim; // Array storing data related to the limiter

	/* Time related data */

	double CFL;         // CFL constant < 1
	double time;        // Current time of the simulation
	double dt;          // Time step size
	int timestep;       // Time step
	int write_interval; // Number of times steps after which to write output files
	int rk_stage;       // Stage of the RK method
	double tend;        // Final time of the simulation

	/* Quadrature related data */

	int NGP_face;
	Array2D<double> xiGP_Face;
	Array2D<double> etaGP_Face;
	Array2D<double> wGP_Face;
	Array3D<double> phiGP_Face;
	int N_Node;
	Array1D<double> xiNode;
	Array1D<double> etaNode;
	Array2D<double> phiNode;

	/* Backbone of the class */

	void select_stencils();
	void compute_reconstruction_matrices();
	void reconstruct();
	void compute_rhs();
	void step_SSPRK33();
	void plot();

	/* Additional monitoring */

	double max_error() const;
	void preserve_positivity();

public:

	/* Constructor and destructor */

	HyPE();
	~HyPE();

	void run();
};

/* Initial condition function and exact solutions */

void icond(double,double,Array1D<double>&);
void exact_solution(double,double,double,double*);

/* Boundary conditions */

void get_right_state(const Array1D<double>&, const Array1D<double>&, int, double, double, Array1D<double>&);


/* PDE related functions */

void PDECons2Prim(const Array1D<double>&, Array1D<double>&);
void PDEPrim2Cons(const Array1D<double>&, Array1D<double>&);
double PDEFlux(const Array1D<double>&, double,double, Array1D<double>&);
double LLFRiemannSolver(const Array1D<double>&, const Array1D<double>&, double, double, Array1D<double>&);

#endif /* HYPE_HH_ */
