/***********************************************************
** Nanyang Technological University, Singapore 637459.
** 2015.
*
**	Raman Monte Carlo simulation of  
**	photon distribution in cuboid  in ANSI Standard C.
****
**	Starting Date:	06/2013.
**	Current Date:	04/2015.
*
**	Vijitha Periyasamy, B.E.
**	Manojit Pramanik, Ph.D.
**	Biomedical Imaging Laboratory,
**  School of Chemical and Biomedical Engineering,
**	Nanyang Technological University, Singapore 637459.
*
*	This program was based on:
**	(1)	MCML-Monte Carlo modeling of photon transport  
**  in multi-layered tissues, L.-H. Wang, S. L. Jacques,
**	and L.-Q. Zheng Computer Methods and programs in 
**	Biomedicine, 47, 131-146 (1995)
*
**  (2) Monte Carlo simulation of light transport in turbid 
**  medium with embedded object - spherical, cylindrical, 
**  ellipsoidal, or cuboidal object embedded within 
**  multilayered tissues, V. Periyasamy and M. Pramanik, 
**  Journal of Biomedical Optics 19(4), 045003 (2014).
*
**	(3) Experimentally validated Raman Monte Carlo simulation  
**	for a cuboid object to obtain Raman Spectroscopic signature
**  for hidden material, V. Periyasamy, S. Sil, G. Dhal, F. Ariese,
**  S. Umapathy, and M. Pramanik, Journal of Raman Spectroscopy, (2015)
*
****
*	General Naming Conventions:
*	Preprocessor names: all capital letters, 
*		e.g. #define PREPROCESSORS
*	Globals: first letter of each word is capital, no 
*		underscores, 
*		e.g. short GlobalVar;
*	Dummy variables:  first letter of each word is capital,
*		and words are connected by underscores, 
*		e.g. void NiceFunction(char Dummy_Var);
*	Local variables:  all lower cases, words are connected 
*		by underscores,
*		e.g. short local_var;
*	Function names or data types:  same as Globals.
*
****
*	Dimension of length: cm.
*
****/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#define PI 3.1415926
#define WEIGHT 1E-4		/* Critical weight for roulette. */
#define CHANCE 0.1		/* Chance of roulette survival. */
#define STRLEN 256		/* String length. */

#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)
#define MIN(x,y) ((x<=y)? x:y)

/****************** Stuctures *****************************/

/****
*	Structure used to describe a photon packet.
****/
typedef struct {
	double x, y, z;		/* Cartesian coordinates.[cm] */
	double ux, uy, uz;	/* directional cosines of a photon. */
	double w;			/* weight. */
	Boolean dead;		/* 1 if photon is terminated. */
	Boolean inMat;		/** 1 if the photon is in material. **/
	double s;			/* current step size. [cm]. */
	double sleft;		/* step size left. dimensionless [-]. */
	Boolean isRaman;	    /**	1 if photon is converted to Raman photon **/
	short ramanLayer;	/**	Layer in which photon is converted to Raman photon **/
	Boolean justRaman;	/** To perform isotropic scattering if the photon is just converted to Raman **/
} PhotonStruct;

/****

*	Structure used to describe the geometry and optical
*	properties of a layer.
*	
*	cos_crit0 and cos_crit1 are the cosines of the 
*	critical angle of total internal reflection for the
*	upper boundary and lower boundary respectively.
*	They are set to zero if no total internal reflection
*	exists.
*	They are used for computation speed.
**	Includeed Raman probability for the layers
****/
typedef struct {
	double z0, z1;		/* z coordinates of a layer. [cm] */
	double n;			/* refractive index of a layer. */
	double mua;			/* absorption coefficient. [1/cm] */
	double mus;			/* scattering coefficient. [1/cm] */
	double g;		    /* anisotropy. */
	double ramProb;		/** Raman probability. **/

	double cos_crit0,	cos_crit1;	
} LayerStruct;


/****
**	Structure used to describe the geometry of box.
******/
typedef struct {
	double lt;				/** Length along x-axis **/
	double bt;				/** Breadth along y-axis **/
	double ht;				/** Heigth along z-axis **/
	double tn;				/** Thickness of container **/
	double XMin,YMin,ZMin;	/** Minmum cordinates of box **/
	double XMax,YMax,ZMax;	/** Maxmum cordinates of box **/
} BoxStruct;

/****
*	Input parameters for each independent run.
*
*	z and r are for the cylindrical coordinate system. [cm]
*	a is for the angle alpha between the photon exiting 
*	direction and the surface normal. [radian]
*
*	The grid line separations in z, r, and alpha
*	directions are dz, dr, and da respectively.  The numbers 
*	of grid lines in z, r, and alpha directions are
*	nz, nr, and na respectively.
*
**	x, y and z are for cartesian coordinate system.
**	The grid line separations in x,y and z
**	directions are dx, dy, and dz respectively.  The numbers 
**	of grid lines in x, y, and z directions are
**	nx, ny, and nz respectively.
*	
**	The member boxspecs stores the dimension of box.
*	
*	The member layerspecs will point to an array of 
*	structures which store parameters of each layer. 
*	This array has (2 + 2) elements. One element is 
*	for a container and the other for material.
*	The layers 0 and 3 are for top ambient 
*	medium and the bottom ambient medium respectively.
****/
typedef struct {
	char	 out_fname[STRLEN];	/* output file name. */
	char	 out_fformat;		/* output file format. */
	/* 'A' for ASCII, */
	/* 'B' for binary. */
	long	 num_photons; 		/* to be traced. */
	double Wth; 				/* play roulette if photon */
	/* weight < Wth.*/

	double dz;					/* z grid separation.[cm] */ 
	double dr;					/* r grid separation.[cm] */
	double da;					/* alpha grid separation. */
	/* [radian] */
	short nz;					/* array range 0..nz-1. */
	short nr;					/* array range 0..nr-1. */
	short na;					/* array range 0..na-1. */

	double dx;					/** x grid separation.[cm] **/ 
	double dy;					/** y grid separation.[cm] **/

	short nx;					/** array range 0..nx-1. **/
	short ny;					/** array range 0..ny-1. **/

	short	num_layers;			/* number of layers. */
	
	LayerStruct * layerspecs;	/* layer parameters. */	
	BoxStruct * boxspecs;		/** box dimensions. **/	
	Boolean ramanScatter;		/** 1 is to have isotropic scatter of Raman **/
} InputStruct;

/****
*	Structures for scoring physical quantities. 
*	z and r represent z and r coordinates of the 
*	cylindrical coordinate system. [cm]
**	x,y and z represent coordinates of the 
**	cartesian coordinate system. [cm]
*	a is the angle alpha between the photon exiting 
*	direction and the normal to the surfaces. [radian]
****/
typedef struct {
	double    Rsp;		/* specular reflectance. [-] */
	double ** Rd_ra;	/* 2D distribution of diffuse */
	/* reflectance. [1/(cm2 sr)] */

	double ** A_rz;		/* 2D probability density in turbid */
	/* media over r & z. [1/cm3] */
	double ** Tt_ra;	    /* 2D distribution of transmittance in ambient layer*/
	double ** RdCon_xy;		/** 2D distribution of diffuse reflectance from container **/
	double ** TtCon_xy;		/** 2D distribution of transmittance from container in ambient medium **/
	double ** TtCon_yz;		/** 2D distribution of transmittance from container in constant x surface **/
	double ** TtCon_xz;		/** 2D distribution of transmittance from container in constant y surface **/	
	double ** RdMat_xy;		/** 2D distribution of diffuse reflectance from container **/
	double ** TtMat_xy;		/** 2D distribution of transmittance from material in ambient surface **/
	double ** TtMat_yz;		/** 2D distribution of transmittance from material in constant x surface **/
	double ** TtMat_xz;		/** 2D distribution of transmittance from material in constant y surface **/
	long RamanPhotons;		/** Track number of Raman photons **/
} OutStruct;

/***********************************************************
*	Routine prototypes for dynamic memory allocation and 
*	release of arrays and matrices.
*	Modified from Numerical Recipes in C.
****/
double  *AllocVector(short, short);
double **AllocMatrix(short, short,short, short);
void 	 FreeVector(double *, short, short);
void 	 FreeMatrix(double **, short, short, short, short);
void 	 nrerror(char *);
