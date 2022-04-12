/***********************************************************
** Nanyang Technological University, Singapore 637459.
** 2015.
*
*	Launch, move, and record photon weight.
*
**	main program for Raman Monte Carlo simulation of photon
**	distribution in Cuboid.
*
**	This code modified from the original MCML code
**	developed by Dr Lihong Wang et al.
****/

#include "mcml.h"

#define STANDARDTEST 0
/* testing program using fixed rnd seed. */

#define PARTIALREFLECTION 0     
/* 1=split photon, 0=statistical reflection. */

#define COSZERO (1.0-1.0E-12)	
/* cosine of about 1e-6 rad. */

#define COS90D  1.0E-6		
/* cosine of about 1.57 - 1e-6 rad. */


/***********************************************************
*	A random number generator from Numerical Recipes in C.
****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

float ran3(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
			inext=0;
			inextp=31;
			*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/***********************************************************
*	Generate a random number between 0 and 1.  Take a 
*	number as seed the first time entering the function.  
*	The seed is limited to 1<<15.  
*	We found that when idum is too large, ran3 may return 
*	numbers beyond 0 and 1.
****/
double RandomNum(void)
{
	static Boolean first_time=1;
	static int idum;	/* seed for ran3. */

	if(first_time) {
#if STANDARDTEST /* Use fixed seed to test the program. */
		idum = - 1;
#else
		idum = -(int)time(NULL)%(1<<15);
		/* use 16-bit integer as the seed. */
#endif
		ran3(&idum);
		first_time = 0;
		idum = 1;
	}

	return( (double)ran3(&idum) );
}

/***********************************************************
*	Compute the specular reflection. 
*
*	If the first layer is a turbid medium, use the Fresnel
*	reflection from the boundary of the first layer as the 
*	specular reflectance.
*
*	If the first layer is glass, multiple reflections in
*	the first layer is considered to get the specular
*	reflectance.
*
*	The subroutine assumes the Layerspecs array is correctly 
*	initialized.
****/
double Rspecular(LayerStruct * Layerspecs_Ptr)
{
	double r1, r2;
	/* direct reflections from the 1st and 2nd layers. */
	double temp;

	temp =(Layerspecs_Ptr[0].n - Layerspecs_Ptr[1].n)
		/(Layerspecs_Ptr[0].n + Layerspecs_Ptr[1].n);
	r1 = temp*temp;

	if((Layerspecs_Ptr[1].mua == 0.0) 
		&& (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
			temp = (Layerspecs_Ptr[1].n - Layerspecs_Ptr[2].n)
				/(Layerspecs_Ptr[1].n + Layerspecs_Ptr[2].n);
			r2 = temp*temp;
			r1 = r1 + (1-r1)*(1-r1)*r2/(1-r1*r2);
	}

	return (r1);	
}

/***********************************************************
**	Initialize a photon packet. With pointers for Raman
****/
void LaunchPhoton(double Rspecular,
	LayerStruct  * Layerspecs_Ptr,
	PhotonStruct * Photon_Ptr)
{

	double radRand, radius, angRand, aziAng; 
	FILE * fp;

	Photon_Ptr->w	 	= 1.0 - Rspecular;	
	Photon_Ptr->dead 	= 0;
	Photon_Ptr->inMat	= 0;
	Photon_Ptr->isRaman = 0;
	Photon_Ptr->s		= 0;
	Photon_Ptr->sleft	= 0;
	Photon_Ptr->ramanLayer = 0;
	Photon_Ptr->justRaman = 0;

	Photon_Ptr->x 	= 0.0;	
	Photon_Ptr->y	= 0.0;	
	Photon_Ptr->z	= 0.0;	

	Photon_Ptr->ux	= 0.0;	
	Photon_Ptr->uy	= 0.0;	
	Photon_Ptr->uz	= sqrt(1-Photon_Ptr->ux*Photon_Ptr->ux-Photon_Ptr->uy*Photon_Ptr->uy);	

	if((Layerspecs_Ptr[1].mua == 0.0) 
		&& (Layerspecs_Ptr[1].mus == 0.0))  { /* glass layer. */
			Photon_Ptr->inMat 	= 1;
			Photon_Ptr->z	= Layerspecs_Ptr[2].z0;	
	}
}

/***********************************************************
*	Choose (sample) a new theta angle for photon propagation
*	according to the anisotropy.
*
*	If anisotropy g is 0, then
*		cos(theta) = 2*rand-1.
*	otherwise
*		sample according to the Henyey-Greenstein function.
*
*	Returns the cosine of the polar deflection angle theta.
****/
double SpinTheta(double g, PhotonStruct* Photon_Ptr)
{
	double cost;

	if(g == 0.0 || Photon_Ptr->justRaman==1) {
		cost = 2*RandomNum() -1;
		Photon_Ptr->justRaman=0;
	}
	else {
		double temp = (1-g*g)/(1-g+2*g*RandomNum());
		cost = (1+g*g - temp*temp)/(2*g);
		if(cost < -1) cost = -1;
		else if(cost > 1) cost = 1;
	}
	return(cost);
}


/***********************************************************
*	Choose a new direction for photon propagation by 
*	sampling the polar deflection angle theta and the 
*	azimuthal angle psi.
*
*	Note:
*  	theta: 0 - pi so sin(theta) is always positive 
*  	feel free to use sqrt() for cos(theta).
* 
*  	psi:   0 - 2pi 
*  	for 0-pi  sin(psi) is + 
*  	for pi-2pi sin(psi) is - 
****/
void Spin(double g,
	PhotonStruct * Photon_Ptr)
{
	double cost, sint;	/* cosine and sine of the */
	/* polar deflection angle theta. */
	double cosp, sinp;	/* cosine and sine of the */
	/* azimuthal angle psi. */
	double ux = Photon_Ptr->ux;
	double uy = Photon_Ptr->uy;
	double uz = Photon_Ptr->uz;
	double psi;

	cost = SpinTheta(g,Photon_Ptr);
	sint = sqrt(1.0 - cost*cost);	
	/* sqrt() is faster than sin(). */

	psi = 2.0*PI*RandomNum(); /* spin psi 0-2pi. */
	cosp = cos(psi);
	if(psi<PI)
		sinp = sqrt(1.0 - cosp*cosp);	
	/* sqrt() is faster than sin(). */
	else
		sinp = - sqrt(1.0 - cosp*cosp);	

	if(fabs(uz) > COSZERO)  { 	/* normal incident. */
		Photon_Ptr->ux = sint*cosp;
		Photon_Ptr->uy = sint*sinp;
		Photon_Ptr->uz = cost*SIGN(uz);	
		/* SIGN() is faster than division. */
	}
	else  {		/* regular incident. */
		double temp = sqrt(1.0 - uz*uz);
		Photon_Ptr->ux = sint*(ux*uz*cosp - uy*sinp)
			/temp + ux*cost;
		Photon_Ptr->uy = sint*(uy*uz*cosp + ux*sinp)
			/temp + uy*cost;
		Photon_Ptr->uz = -sint*cosp*temp + uz*cost;
	}
}

/***********************************************************
*	Move the photon s away in the current layer of medium.  
****/
void Hop(PhotonStruct *	Photon_Ptr)
{
	double s = Photon_Ptr->s;

	Photon_Ptr->x += s*Photon_Ptr->ux;
	Photon_Ptr->y += s*Photon_Ptr->uy;
	Photon_Ptr->z += s*Photon_Ptr->uz;
}			

/***********************************************************
*	If uz != 0, return the photon step size in glass, 
*	Otherwise, return 0.
*
*	The step size is the distance between the current 
*	position and the boundary in the photon direction.
*
*	Make sure uz !=0 before calling this function.
****/
void StepSizeInGlass(PhotonStruct *  Photon_Ptr,
	InputStruct  *  In_Ptr)
{
	double dl_b;	/* step size to boundary. */
	short  layer = Photon_Ptr->inMat;
	double uz = Photon_Ptr->uz;

	/* Stepsize to the boundary. */	
	if(uz>0.0)
		dl_b = (In_Ptr->layerspecs[layer].z1 - Photon_Ptr->z)
		/uz;
	else if(uz<0.0)
		dl_b = (In_Ptr->layerspecs[layer].z0 - Photon_Ptr->z)
		/uz;
	else
		dl_b = 0.0;

	Photon_Ptr->s = dl_b;
}

/***********************************************************
*	Pick a step size for a photon packet when it is in 
*	tissue.
*	If the member sleft is zero, make a new step size 
*	with: -log(rnd)/(mua+mus).
*	Otherwise, pick up the leftover in sleft.
*
*	Layer is the index to layer.
*	In_Ptr is the input parameters.
****/
void StepSizeInTissue(PhotonStruct * Photon_Ptr,
	InputStruct  * In_Ptr)
{
	short  layer = Photon_Ptr->inMat;
	double mua = In_Ptr->layerspecs[layer+1].mua;
	double mus = In_Ptr->layerspecs[layer+1].mus;

	if(Photon_Ptr->sleft == 0.0) {  /* make a new step. */
		double rnd;
		do rnd = RandomNum(); 
		while( rnd <= 0.0 );    /* avoid zero. */
		Photon_Ptr->s = -log(rnd)/(mua+mus);
	}
	else {	/* take the leftover. */
		Photon_Ptr->s = Photon_Ptr->sleft/(mua+mus);
		Photon_Ptr->sleft = 0.0;
	}
}

/***********************************************************
*	Check if the step will hit the boundary.
*	Return 1 if hit boundary.
*	Return 0 otherwise.
*
* 	If the projected step hits the boundary, the members
*	s and sleft of Photon_Ptr are updated.
****/
Boolean HitBoundary(PhotonStruct *  Photon_Ptr,
	InputStruct  *  In_Ptr)
{
	double dl_b;  /* length to boundary. */
	short  layer = Photon_Ptr->inMat + 1;
	double uz = Photon_Ptr->uz;
	Boolean hit;

	/* Distance to the boundary. */
	if(uz>0.0)
		dl_b = (In_Ptr->layerspecs[layer].z1 
		- Photon_Ptr->z)/uz;	/* dl_b>0. */
	else if(uz<0.0)
		dl_b = (In_Ptr->layerspecs[layer].z0 
		- Photon_Ptr->z)/uz;	/* dl_b>0. */

	if(uz != 0.0 && Photon_Ptr->s > dl_b) {
		/* not horizontal & crossing. */
		double mut = In_Ptr->layerspecs[layer].mua 
			+ In_Ptr->layerspecs[layer].mus;

		Photon_Ptr->sleft = (Photon_Ptr->s - dl_b)*mut;
		Photon_Ptr->s    = dl_b;
		hit = 1;
	}
	else
		hit = 0;

	return(hit);
}

/***********************************************************
**	Check if the given distance moved by a photon 
**	lies on the surface of box.
**	Return 1 if hits surface.
**	Return 0 otherwise.
*****/

Boolean checkOnPlane(double bXMin,double bYMin,double bZMin,double bXMax,double bYMax,double bZMax,
	double x,double y,double z,double ux,double uy,double uz, double dist) 
{
	double newX,newY,newZ;
	newX = (x + ux*dist);	/*Find the point after the photon travels the distance */
	newY = (y + uy*dist);
	newZ = (z + uz*dist);
	if((bXMin <= newX && newX <= bXMax)&& /* Check if the point is on box surface */
		(bYMin <= newY && newY <= bYMax)&&
		(bZMin <= fabs(newZ) && newZ <= bZMax)){
			return(1);
	}
	else return(0);
}
/***********************************************************
**	Check if the step will hit the surface of box.
**	Return 1 if hits boundary.
**	Return 0 otherwise.
**	The minimum distance to boundary is saved.
****/

Boolean hitBox(int boxNum, InputStruct  *  In_Ptr, 
	PhotonStruct *  Photon_Ptr, double * minDistBox) 
{
	double bXMin, bYMin, bZMin, bXMax, bYMax, bZMax;
	double x, y, z, ux, uy, uz;
	double distX1, distX2, distY1, distY2, distZ1, distZ2;
	double minDistX = 9999, minDistY = 9999, minDistZ = 9999;
	float newX, newY, newZ;
	double minDist, tempMin;
	double ERR = 1E-9, HIGHVAL = 99999;

	/** Get the box dimensions **/
	bXMin = In_Ptr->boxspecs[boxNum].XMin;
	bYMin = In_Ptr->boxspecs[boxNum].YMin;
	bZMin = In_Ptr->boxspecs[boxNum].ZMin;
	bXMax = In_Ptr->boxspecs[boxNum].XMax;
	bYMax = In_Ptr->boxspecs[boxNum].YMax;
	bZMax = In_Ptr->boxspecs[boxNum].ZMax;

	x = Photon_Ptr->x;
	y = Photon_Ptr->y;
	z = Photon_Ptr->z;

	ux = Photon_Ptr->ux;
	uy = Photon_Ptr->uy;
	uz = Photon_Ptr->uz;

	/**	Find the distance of photon with respect to 6 directions. Discard the negative distance	**/
	if(ux != 0) {
		distX1 = (bXMin-x)/ux;
		if(distX1 < ERR || !checkOnPlane(bXMin,bYMin,bZMin,bXMax,bYMax,bZMax,x,y,z,ux,uy,uz,distX1)){
			distX1 = HIGHVAL;		/** Setting a high value because point doesn't lie on surface with the distance **/
		}
		distX2 = (bXMax-x)/ux;
		if(distX2 < ERR || !checkOnPlane(bXMin,bYMin,bZMin,bXMax,bYMax,bZMax,x,y,z,ux,uy,uz,distX2)){
			distX2 = HIGHVAL;
		}
		minDistX = MIN(distX1,distX2);
	}

	if(uy != 0) {
		distY1 = (bYMin-y)/uy;
		if(distY1 < ERR || !checkOnPlane(bXMin,bYMin,bZMin,bXMax,bYMax,bZMax,x,y,z,ux,uy,uz,distY1)){
			distY1 = HIGHVAL;
		}
		distY2 = (bYMax-y)/uy;
		if(distY2 < ERR || !checkOnPlane(bXMin,bYMin,bZMin,bXMax,bYMax,bZMax,x,y,z,ux,uy,uz,distY2)){
			distY2 = HIGHVAL;
		}
		minDistY = MIN(distY1,distY2);
	}

	if (uz !=0) {
		distZ1 = (bZMin-z)/uz;
		if(distZ1 < ERR || !checkOnPlane(bXMin,bYMin,bZMin,bXMax,bYMax,bZMax,x,y,z,ux,uy,uz,distZ1)){
			distZ1 = HIGHVAL;
		}
		distZ2 = (bZMax-z)/uz;
		if(distZ2 < ERR || !checkOnPlane(bXMin,bYMin,bZMin,bXMax,bYMax,bZMax,x,y,z,ux,uy,uz,distZ2)){
			distZ2 = HIGHVAL;
		}
		minDistZ = MIN(distZ1,distZ2);
	}

	if (minDistX == HIGHVAL && minDistY == HIGHVAL && minDistZ == HIGHVAL) {
		return 0; /** No hit when there is no valid distance. **/
	}

	minDist = MIN(MIN(minDistX,minDistY),minDistZ); /** Minimum of the distance in x,y and z directions **/

	if (minDist <= Photon_Ptr->s) {
		*minDistBox = minDist;
		return(1);		/** Hit if minimum distance is less than step-size **/
	}
	return(0);
}

/***********************************************************
**	Check if the step will hit the surface of box1 or box2.
**	Return 1 if hits boundary.
**	Return 0 otherwise.
*
** 	If the projected step hits the boundary, the members
**	s and sleft of Photon_Ptr are updated.
****/
Boolean HitBoundaryBox(PhotonStruct *  Photon_Ptr,
	InputStruct  *  In_Ptr)
{
	Boolean hitBox1 = 0,hitBox2 = 0;
	double minDistBox1, minDistBox2, hitDist;
	double mut;

	hitBox2 = hitBox(1,In_Ptr,Photon_Ptr,&minDistBox2);		/** Check if photon hits inner box **/
	if (Photon_Ptr->inMat == 0) {
		hitBox1 = hitBox(0,In_Ptr,Photon_Ptr,&minDistBox1);	/** Check if photon hits outer box and if it is in outer box **/
	}

	/**	Get minimum of the hit distance from outer and inner box when there is a hit **/
	if(hitBox1 && hitBox2) {
		hitDist = MIN(minDistBox1,minDistBox2);
	} 
	else if(hitBox1 == 1) {
		hitDist = minDistBox1;
	}
	else if(hitBox2==1) {
		hitDist = minDistBox2;
	}
	else {
		return(0);
	}

	/**	Update the photon pointer when there is a hit. **/
	mut = In_Ptr->layerspecs[Photon_Ptr->inMat+1].mua + In_Ptr->layerspecs[Photon_Ptr->inMat+1].mus;

	Photon_Ptr->sleft = (Photon_Ptr->s - hitDist)*mut;
	Photon_Ptr->s    = hitDist;
	return(1);
}


/***********************************************************
*	Drop photon weight inside the tissue (not glass).
*
*  The photon is assumed not dead. 
*
*	The weight drop is dw = w*mua/(mua+mus).
*
*	The dropped weight is assigned to the absorption array 
*	elements.
****/
void Drop(InputStruct  *	In_Ptr, 
	PhotonStruct *	Photon_Ptr,
	OutStruct *		Out_Ptr)
{
	double dwa;		/* absorbed weight.*/
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	double izd, ird, ixd, iyd;	/* LW 5/20/98. To avoid out of short range.*/
	long  iz, ir, ix, iy;	/* index to z & r. */
	short  layer = Photon_Ptr->inMat+1;
	double mua, mus;		

	/* compute array indices. */
	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd > In_Ptr->nz-1) iz = In_Ptr->nz-1;
	else iz = izd;

	ixd = x/In_Ptr->dx;

	if(ixd > 0) ixd = ixd + 0.5;
	else ixd = ixd - 0.5;
	if(fabs(ixd) > (In_Ptr->nx-1)/2) ix = SIGN(ixd)*(In_Ptr->nx-1)/2;
	else ix = ixd;
	ix = (In_Ptr->nx-1)/2 + ix;

	iyd = y/In_Ptr->dy;
	if(iyd > 0) iyd = iyd + 0.5;
	else iyd = iyd - 0.5;
	if(fabs(iyd) > (In_Ptr->ny-1)/2) iy = SIGN(iyd)*(In_Ptr->ny-1)/2;
	else iy = iyd;
	iy = (In_Ptr->ny-1)/2 + iy;

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird > In_Ptr->nr-1) ir = In_Ptr->nr-1;
	else ir = ird;

	/* update photon weight. */
	mua = In_Ptr->layerspecs[layer].mua;
	mus = In_Ptr->layerspecs[layer].mus;
	dwa = Photon_Ptr->w * mua/(mua+mus);
	Photon_Ptr->w -= dwa;

	/* assign dwa to the absorption array element. */
	Out_Ptr->A_rz[ir][iz] += dwa;
}

/***********************************************************
*	The photon weight is small, and the photon packet tries 
*	to survive a roulette.
****/
void Roulette(PhotonStruct * Photon_Ptr)
{
	if(Photon_Ptr->w == 0.0)	
		Photon_Ptr->dead = 1;
	else if(RandomNum() < CHANCE) /* survived the roulette.*/
		Photon_Ptr->w /= CHANCE;
	else 
		Photon_Ptr->dead = 1;
}

/***********************************************************
*	Compute the Fresnel reflectance.
*
*	Make sure that the cosine of the incident angle a1
*	is positive, and the case when the angle is greater 
*	than the critical angle is ruled out.
*
* 	Avoid trigonometric function operations as much as
*	possible, because they are computation-intensive.
****/
double RFresnel(double n1,	/* incident refractive index.*/
	double n2,	/* transmit refractive index.*/
	double ca1,	/* cosine of the incident */
	/* angle. 0<a1<90 degrees. */
	double * ca2_Ptr)  /* pointer to the */
	/* cosine of the transmission */
	/* angle. a2>0. */
{
	double r;

	if(n1==n2) {			  	/** matched boundary. **/
		*ca2_Ptr = ca1;
		r = 0.0;
	}
	else if(ca1>COSZERO) {	/** normal incident. **/
		*ca2_Ptr = ca1;
		r = (n2-n1)/(n2+n1);
		r *= r;
	}
	else if(ca1<COS90D)  {	/** very slant. **/
		*ca2_Ptr = 0.0;
		r = 1.0;
	}
	else  {			  		/** general. **/
		double sa1, sa2;	
		/* sine of the incident and transmission angles. */
		double ca2;

		sa1 = sqrt(1-ca1*ca1);
		sa2 = n1*sa1/n2;
		if(sa2>=1.0) {	
			/* double check for total internal reflection. */
			*ca2_Ptr = 0.0;
			r = 1.0;
		}
		else  {
			double cap, cam;	/* cosines of the sum ap or */
			/* difference am of the two */
			/* angles. ap = a1+a2 */
			/* am = a1 - a2. */
			double sap, sam;	/* sines. */

			*ca2_Ptr = ca2 = sqrt(1-sa2*sa2);

			cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
			cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
			sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
			sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
			r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
			/* rearranged for speed. */
		}
	}
	return(r);
}

/***********************************************************
*	Record the photon weight exiting the first layer(uz<0), 
*	no matter whether the layer is glass or not, to the 
*	reflection array.
*
*	Update the photon weight as well.
****/
void RecordR(double			Refl,	/* reflectance. */
	InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ir, ia;	/* index to r & angle. */
	double ird, iad;	/* LW 5/20/98. To avoid out of short range.*/

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird > In_Ptr->nr-1) ir = In_Ptr->nr-1;
	else ir = ird;

	iad = acos(-Photon_Ptr->uz)/In_Ptr->da;
	if(iad > In_Ptr->na-1) ia = In_Ptr->na-1;
	else ia = iad;

	/* assign photon to the reflection array element. */
	Out_Ptr->Rd_ra[ir][ia] += Photon_Ptr->w*(1.0-Refl);

	Photon_Ptr->w *= Refl;
}
/***********************************************************
**	Record the photon which became Raman in conatiner exiting from launch surface(z=0), 
*	no matter whether the layer is glass or not, to the 
*	reflection array. 
****/
void RecordRCon_xy(InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ix, iy;	
	double ixd, iyd;	

	ixd = x/In_Ptr->dx;
	if(ixd > 0) ixd = ixd + 0.5;
	else ixd = ixd - 0.5;
	if(fabs(ixd) > (In_Ptr->nx-1)/2) ix = SIGN(ixd)*(In_Ptr->nx-1)/2;
	else ix = ixd;

	ix = (In_Ptr->nx-1)/2 + ix;

	iyd = y/In_Ptr->dy;
	if(iyd > 0) iyd = iyd + 0.5;
	else iyd = iyd - 0.5;
	if(fabs(iyd) > (In_Ptr->ny-1)/2) iy = SIGN(iyd)*(In_Ptr->ny-1)/2;
	else iy = iyd;

	iy = (In_Ptr->ny-1)/2 + iy;

	/* assign photon to the reflection array element. */
	Out_Ptr->RdCon_xy[ix][iy] += Photon_Ptr->w;
}

/***********************************************************
**	Record the photon which became Raman in conatiner exiting from surface where x=XMAX, 
****/
void RecordTCon_yz(InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double y = Photon_Ptr->y;
	double z = Photon_Ptr->z;
	short  iy, iz;	
	double iyd, izd;

	iyd = y/In_Ptr->dy;
	if(iyd > 0) iyd = iyd + 0.5;
	else iyd = iyd - 0.5;
	if(fabs(iyd) > (In_Ptr->ny-1)/2) iy = SIGN(iyd)*(In_Ptr->ny-1)/2;
	else iy = iyd;

	iy = (In_Ptr->ny-1)/2 + iy;

	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd > In_Ptr->nz-1) iz = In_Ptr->nz-1;
	else iz = izd;

	/* assign photon to the reflection array element. */
	Out_Ptr->TtCon_yz[iy][iz] += Photon_Ptr->w;
}

/***********************************************************
**	Record the photon which became Raman in conatiner exiting from surface where y=YMax, 
*****/
void RecordTCon_xz(InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ix, iz;	
	double ixd, izd;	

	ixd = x/In_Ptr->dx;
	if(ixd > 0) ixd = ixd + 0.5;
	else ixd = ixd - 0.5;
	if(fabs(ixd) > (In_Ptr->nx-1)/2) ix = SIGN(ixd)*(In_Ptr->nx-1)/2;
	else ix = ixd;

	ix = (In_Ptr->nx-1)/2 + ix;

	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd > In_Ptr->nz-1) iz = In_Ptr->nz-1;
	else iz = izd;

	/* assign photon to the reflection array element. */
	Out_Ptr->TtCon_xz[ix][iz] += Photon_Ptr->w;
}

/***********************************************************
**	Record the photon which became Raman in conatiner exiting from ambient surface(z=ht), 
*	no matter whether the layer is glass or not, to the 
*	reflection array.
****/
void RecordTCon_xy(InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ix, iy;	
	double ixd, iyd;	

	ixd = x/In_Ptr->dx;
	if(ixd > 0) ixd = ixd + 0.5;
	else ixd = ixd - 0.5;
	if(fabs(ixd) > (In_Ptr->nx-1)/2) ix = SIGN(ixd)*(In_Ptr->nx-1)/2;
	else ix = ixd;

	ix = (In_Ptr->nx-1)/2 + ix;

	iyd = y/In_Ptr->dy;
	if(iyd > 0) iyd = iyd + 0.5;
	else iyd = iyd - 0.5;
	if(fabs(iyd) > (In_Ptr->ny-1)/2) iy = SIGN(iyd)*(In_Ptr->ny-1)/2;
	else iy = iyd;

	iy = (In_Ptr->ny-1)/2 + iy;

	/* assign photon to the reflection array element. */
	Out_Ptr->TtCon_xy[ix][iy] += Photon_Ptr->w;
}
/***********************************************************
**	Record the photon which became Raman in material exiting from launch surface(z=0), 
*	no matter whether the layer is glass or not, to the 
*	reflection array.
****/
void RecordRMat_xy(InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ix, iy;	
	double ixd, iyd;	

	ixd = x/In_Ptr->dx;
	if(ixd > 0) ixd = ixd + 0.5;
	else ixd = ixd - 0.5;
	if(fabs(ixd) > (In_Ptr->nx-1)/2) ix = SIGN(ixd)*(In_Ptr->nx-1)/2;
	else ix = ixd;

	ix = (In_Ptr->nx-1)/2 + ix;

	iyd = y/In_Ptr->dy;
	if(iyd > 0) iyd = iyd + 0.5;
	else iyd = iyd - 0.5;
	if(fabs(iyd) > (In_Ptr->ny-1)/2) iy = SIGN(iyd)*(In_Ptr->ny-1)/2;
	else iy = iyd;

	iy = (In_Ptr->ny-1)/2 + iy;

	/* assign photon to the reflection array element. */
	Out_Ptr->RdMat_xy[ix][iy] += Photon_Ptr->w;
}

/***********************************************************
**	Record the photon which became Raman in material exiting from surface where x=XMAX, 
****/
void RecordTMat_yz(InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double y = Photon_Ptr->y;
	double z = Photon_Ptr->z;
	short  iy, iz;	
	double iyd, izd;

	iyd = y/In_Ptr->dy;
	if(iyd > 0) iyd = iyd + 0.5;
	else iyd = iyd - 0.5;
	if(fabs(iyd) > (In_Ptr->ny-1)/2) iy = SIGN(iyd)*(In_Ptr->ny-1)/2;
	else iy = iyd;

	iy = (In_Ptr->ny-1)/2 + iy;

	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd > In_Ptr->nz-1) iz = In_Ptr->nz-1;
	else iz = izd;

	/* assign photon to the reflection array element. */
	Out_Ptr->TtMat_yz[iy][iz] += Photon_Ptr->w;
}

/***********************************************************
**	Record the photon which became Raman in material exiting from surface where y=YMax, 
*****/
void RecordTMat_xz(InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double z = Photon_Ptr->z;
	short  ix, iz;	
	double ixd, izd;

	ixd = x/In_Ptr->dx;
	if(ixd > 0) ixd = ixd + 0.5;
	else ixd = ixd - 0.5;
	if(fabs(ixd) > (In_Ptr->nx-1)/2) ix = SIGN(ixd)*(In_Ptr->nx-1)/2;
	else ix = ixd;

	ix = (In_Ptr->nx-1)/2 + ix;

	izd = Photon_Ptr->z/In_Ptr->dz;
	if(izd > In_Ptr->nz-1) iz = In_Ptr->nz-1;
	else iz = izd;

	/* assign photon to the reflection array element. */
	Out_Ptr->TtMat_xz[ix][iz] += Photon_Ptr->w;
}

/***********************************************************
**	Record the photon which became Raman in material exiting from launch surface(z=0), 
*	no matter whether the layer is glass or not, to the 
*	reflection array.
****/
void RecordTMat_xy(InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ix, iy;	
	double ixd, iyd;	

	ixd = x/In_Ptr->dx;
	if(ixd > 0) ixd = ixd + 0.5;
	else ixd = ixd - 0.5;
	if(fabs(ixd) > (In_Ptr->nx-1)/2) ix = SIGN(ixd)*(In_Ptr->nx-1)/2;
	else ix = ixd;

	ix = (In_Ptr->nx-1)/2 + ix;

	iyd = y/In_Ptr->dy;
	if(iyd > 0) iyd = iyd + 0.5;
	else iyd = iyd - 0.5;
	if(fabs(iyd) > (In_Ptr->ny-1)/2) iy = SIGN(iyd)*(In_Ptr->ny-1)/2;
	else iy = iyd;

	iy = (In_Ptr->ny-1)/2 + iy;

	/* assign photon to the reflection array element. */
	Out_Ptr->TtMat_xy[ix][iy] += Photon_Ptr->w;
}

/***********************************************************
*	Record the photon weight exiting the last layer(uz>0), 
*	no matter whether the layer is glass or not, to the 
*	transmittance array.
*
*	Update the photon weight as well.
****/
void RecordT(double 		Refl,
	InputStruct  *	In_Ptr,
	PhotonStruct *	Photon_Ptr,
	OutStruct *	Out_Ptr)
{
	double x = Photon_Ptr->x;
	double y = Photon_Ptr->y;
	short  ir, ia;	/* index to r & angle. */
	double ird, iad;	/* LW 5/20/98. To avoid out of short range.*/

	ird = sqrt(x*x+y*y)/In_Ptr->dr;
	if(ird > In_Ptr->nr-1) ir = In_Ptr->nr-1;
	else ir = ird;

	iad = acos(Photon_Ptr->uz)/In_Ptr->da; /* LW 1/12/2000. Removed -. */
	if(iad > In_Ptr->na-1) ia = In_Ptr->na-1;
	else ia = iad;

	/* assign photon to the transmittance array element. */
	Out_Ptr->Tt_ra[ir][ia] += Photon_Ptr->w*(1.0-Refl);

	Photon_Ptr->w *= Refl;
}

/***********************************************************
**	Record the photon weight exiting the cuboid surfaces 
*	Update the photon weight as well.
****/
void RecordRT(InputStruct  *      In_Ptr,
	PhotonStruct *       Photon_Ptr,
	OutStruct *   Out_Ptr)
{     
	double SMALL = 1.0E-9;

	float x = (float)Photon_Ptr->x;
	float y = (float)Photon_Ptr->y;
	float z = (float)Photon_Ptr->z;


	if(Photon_Ptr->ramanLayer == 1) {
		if (fabs(z) < SMALL){
			RecordRCon_xy(In_Ptr,Photon_Ptr,Out_Ptr);
			RecordR(0.0,In_Ptr,Photon_Ptr,Out_Ptr);
		}
		else if (fabs(x-In_Ptr->boxspecs[0].XMax) < SMALL) {
			RecordTCon_yz(In_Ptr,Photon_Ptr,Out_Ptr);
		}
		else if (fabs(y-In_Ptr->boxspecs[0].YMax) < SMALL) {
			RecordTCon_xz(In_Ptr,Photon_Ptr,Out_Ptr);
		}
		else if (fabs(z-In_Ptr->boxspecs[0].ZMax) < SMALL) {
			RecordTCon_xy(In_Ptr,Photon_Ptr,Out_Ptr);
			RecordT(0.0,In_Ptr,Photon_Ptr,Out_Ptr);
		}
		else {
		}
	}
	else {
		if (fabs(z) < SMALL){
			RecordRMat_xy(In_Ptr,Photon_Ptr,Out_Ptr);
			RecordR(0.0,In_Ptr,Photon_Ptr,Out_Ptr);
		}
		else if (fabs(x-In_Ptr->boxspecs[0].XMax) < SMALL) {
			RecordTMat_yz(In_Ptr,Photon_Ptr,Out_Ptr);
		}
		else if (fabs(y-In_Ptr->boxspecs[0].YMax) < SMALL) {
			RecordTMat_xz(In_Ptr,Photon_Ptr,Out_Ptr);
		}
		else if (fabs(z-In_Ptr->boxspecs[0].ZMax) < SMALL) {
			RecordTMat_xy(In_Ptr,Photon_Ptr,Out_Ptr);
			RecordT(0.0,In_Ptr,Photon_Ptr,Out_Ptr);
		}
		else {
		}
	}
}

/***********************************************************
*	Decide whether the photon will be transmitted or 
*	reflected on the upper boundary (uz<0) of the current 
*	layer.
*
*	If "layer" is the first layer, the photon packet will 
*	be partially transmitted and partially reflected if 
*	PARTIALREFLECTION is set to 1,
*	or the photon packet will be either transmitted or 
*	reflected determined statistically if PARTIALREFLECTION 
*	is set to 0.
*
*	Record the transmitted photon weight as reflection.  
*
*	If the "layer" is not the first layer and the photon 
*	packet is transmitted, move the photon to "layer-1".
*
*	Update the photon parmameters.
****/
void CrossUpOrNot(InputStruct  *	In_Ptr, 
	PhotonStruct *	Photon_Ptr,
	OutStruct *		Out_Ptr)
{
	double uz = Photon_Ptr->uz; /* z directional cosine. */
	double uz1;	/* cosines of transmission alpha. always positive. */
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->inMat + 1;
	double ni = In_Ptr->layerspecs[layer].n;
	double nt = In_Ptr->layerspecs[layer-1].n;

	/* Get r. */
	if( - uz <= In_Ptr->layerspecs[layer].cos_crit0) 
		r = 1.0;		      /* total internal reflection. */
	else r = RFresnel(ni, nt, -uz, &uz1);

#if PARTIALREFLECTION
	if(layer == 1 && r<1.0) {	/* partially transmitted. */
		Photon_Ptr->uz = -uz1;	/* transmitted photon. */
		RecordR(r, In_Ptr, Photon_Ptr, Out_Ptr);
		Photon_Ptr->uz = -uz;	/* reflected photon. */
	}		
	else if(RandomNum() > r) {/* transmitted to layer-1. */
		Photon_Ptr->layer--;
		Photon_Ptr->ux *= ni/nt;
		Photon_Ptr->uy *= ni/nt;
		Photon_Ptr->uz = -uz1;
	}
	else			      		/* reflected. */
		Photon_Ptr->uz = -uz;
#else
	if(RandomNum() > r) {		/* transmitted to layer-1. */
		if(layer == 1)  {
			Photon_Ptr->uz = -uz1;
			RecordR(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
			Photon_Ptr->dead = 1;
		}
		else {
			Photon_Ptr->inMat = 1-Photon_Ptr->inMat;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz = -uz1;
		}
	}
	else 						/* reflected. */
		Photon_Ptr->uz = -uz;
#endif
}

/***********************************************************
*	Decide whether the photon will be transmitted  or be 
*	reflected on the bottom boundary (uz>0) of the current 
*	layer.
*
*	If the photon is transmitted, move the photon to 
*	"layer+1". If "layer" is the last layer, record the 
*	transmitted weight as transmittance. See comments for 
*	CrossUpOrNot.
*
*	Update the photon parmameters.
****/
void CrossDnOrNot(InputStruct  *	In_Ptr, 
	PhotonStruct *	Photon_Ptr,
	OutStruct *		Out_Ptr)
{
	double uz = Photon_Ptr->uz; /* z directional cosine. */
	double uz1;	/* cosines of transmission alpha. */
	double r = 0.0;	/* reflectance */
	short  layer = Photon_Ptr->inMat;
	double ni = In_Ptr->layerspecs[layer].n;
	double nt = In_Ptr->layerspecs[layer+1].n;

	/* Get r. */
	if( uz <= In_Ptr->layerspecs[layer].cos_crit1) 
		r = 1.0;		/* total internal reflection. */
	else r = RFresnel(ni, nt, uz, &uz1);

#if PARTIALREFLECTION	
	if(layer == In_Ptr->num_layers && r<1.0) {
		Photon_Ptr->uz = uz1;
		RecordT(r, In_Ptr, Photon_Ptr, Out_Ptr);
		Photon_Ptr->uz = -uz;
	}
	else if(RandomNum() > r) {/* transmitted to layer+1. */
		Photon_Ptr->layer++;
		Photon_Ptr->ux *= ni/nt;
		Photon_Ptr->uy *= ni/nt;
		Photon_Ptr->uz = uz1;
	}
	else 						/* reflected. */
		Photon_Ptr->uz = -uz;
#else
	if(RandomNum() > r) {		/* transmitted to layer+1. */
		if(layer == In_Ptr->num_layers) {
			Photon_Ptr->uz = uz1;
			RecordT(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
			Photon_Ptr->dead = 1;
		}
		else {
			Photon_Ptr->inMat = 1-Photon_Ptr->inMat;
			Photon_Ptr->ux *= ni/nt;
			Photon_Ptr->uy *= ni/nt;
			Photon_Ptr->uz = uz1;
		}
	}
	else 						/* reflected. */
		Photon_Ptr->uz = -uz;
#endif
}

/***********************************************************
****/
void CrossOrNot(InputStruct  *	In_Ptr, 
	PhotonStruct *	Photon_Ptr,
	OutStruct    *	Out_Ptr)
{
	if(Photon_Ptr->uz < 0.0)
		CrossUpOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
	else
		CrossDnOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
}



/***********************************************************
**	Check if the photon hits the surface of the box and escapes.
****/
void CrossOrNotBox(InputStruct  * In_Ptr,
	PhotonStruct *       Photon_Ptr,
	OutStruct    *       Out_Ptr)
{
	double SMALL = 1.0E-9;
	float x = (float)Photon_Ptr->x;
	float y = (float)Photon_Ptr->y;
	float z = (float)Photon_Ptr->z;

	if (Photon_Ptr->inMat == 1) {     /** If photon is in material change to container **/
		Photon_Ptr->inMat = 0;
	}
	else {
		/** Check if photon lies on outer surface of the container and if so record it and kill **/
		if (fabs(z) < SMALL || fabs(z-In_Ptr->boxspecs[0].ht) < SMALL ||
			fabs(x-In_Ptr->boxspecs[0].lt/2) < SMALL ||
			fabs(x+In_Ptr->boxspecs[0].lt/2) < SMALL ||
			fabs(y-In_Ptr->boxspecs[0].bt/2) < SMALL ||
			fabs(y+In_Ptr->boxspecs[0].bt/2) < SMALL)  {

				if (Photon_Ptr->x > In_Ptr->boxspecs[0].XMax) {
					Photon_Ptr->x = In_Ptr->boxspecs[0].XMax;
				}
				if (Photon_Ptr->x < In_Ptr->boxspecs[0].XMin) {
					Photon_Ptr->x = In_Ptr->boxspecs[0].XMin;
				}

				if (Photon_Ptr->y > In_Ptr->boxspecs[0].YMax) {
					Photon_Ptr->y = In_Ptr->boxspecs[0].YMax;
				}
				if (Photon_Ptr->y < In_Ptr->boxspecs[0].YMin) {
					Photon_Ptr->y = In_Ptr->boxspecs[0].YMin;
				}

				if (Photon_Ptr->z > In_Ptr->boxspecs[0].ZMax) {
					Photon_Ptr->z = In_Ptr->boxspecs[0].ZMax;
				}
				if (Photon_Ptr->z < 0) {
					Photon_Ptr->z = 0;
				}

				x = (float)Photon_Ptr->x;
				y = (float)Photon_Ptr->y;
				z = (float)Photon_Ptr->z;

				Photon_Ptr->dead = 1;
				if(Photon_Ptr->isRaman == 1 ){
					RecordRT(In_Ptr,Photon_Ptr,Out_Ptr);
				}     
				else {
					RecordR(0.0,In_Ptr, Photon_Ptr, Out_Ptr);
				}
		}
		else {
			Photon_Ptr->inMat = 1;
		}
	}
}


/***********************************************************
*	Move the photon packet in glass layer.
*	Horizontal photons are killed because they will
*	never interact with tissue again.
****/
void HopInGlass(InputStruct  * In_Ptr,
	PhotonStruct * Photon_Ptr,
	OutStruct    * Out_Ptr)
{
	double dl;     /* step size. 1/cm */

	if(Photon_Ptr->uz == 0.0) { 
		/* horizontal photon in glass is killed. */
		Photon_Ptr->dead = 1;
	}
	else {
		StepSizeInGlass(Photon_Ptr, In_Ptr);
		Hop(Photon_Ptr);
		CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
	}
}

/***********************************************************
*	Set a step size, move the photon, drop some weight, 
*	choose a new photon direction for propagation.  
*
*	When a step size is long enough for the photon to 
*	hit an interface, this step is divided into two steps. 
*	First, move the photon to the boundary free of 
*	absorption or scattering, then decide whether the 
*	photon is reflected or transmitted.
*	Then move the photon in the current or transmission 
*	medium with the unfinished stepsize to interaction 
*	site.  If the unfinished stepsize is still too long, 
*	repeat the above process.  
****/
void HopDropSpinInTissue(InputStruct  *  In_Ptr,
	PhotonStruct *  Photon_Ptr,
	OutStruct    *  Out_Ptr)
{
	/** Recheck the boundary conditions **/
	if (Photon_Ptr->x > In_Ptr->boxspecs[0].XMax || Photon_Ptr->x < In_Ptr->boxspecs[0].XMin ||
		Photon_Ptr->y > In_Ptr->boxspecs[0].YMax || Photon_Ptr->y < In_Ptr->boxspecs[0].YMin ||
		Photon_Ptr->z > In_Ptr->boxspecs[0].ZMax || Photon_Ptr->z < In_Ptr->boxspecs[0].ZMin) {

			if (Photon_Ptr->x > In_Ptr->boxspecs[0].XMax) {
				Photon_Ptr->x = In_Ptr->boxspecs[0].XMax;
			}
			if (Photon_Ptr->x < In_Ptr->boxspecs[0].XMin) {
				Photon_Ptr->x = In_Ptr->boxspecs[0].XMin;
			}

			if (Photon_Ptr->y > In_Ptr->boxspecs[0].YMax) {
				Photon_Ptr->y = In_Ptr->boxspecs[0].YMax;
			}
			if (Photon_Ptr->y < In_Ptr->boxspecs[0].YMin) {
				Photon_Ptr->y = In_Ptr->boxspecs[0].YMin;
			}

			if (Photon_Ptr->z > In_Ptr->boxspecs[0].ZMax) {
				Photon_Ptr->z = In_Ptr->boxspecs[0].ZMax;
			}
			if (Photon_Ptr->z < In_Ptr->boxspecs[0].ZMin) {
				Photon_Ptr->z = In_Ptr->boxspecs[0].ZMin;
			}

			Photon_Ptr->dead = 1;
			if(Photon_Ptr->isRaman == 1 ){
				RecordRT(In_Ptr,Photon_Ptr,Out_Ptr);
			}     
			else {
				RecordR(0.0,In_Ptr, Photon_Ptr, Out_Ptr);
			}

	} else {
		StepSizeInTissue(Photon_Ptr, In_Ptr);

		if(HitBoundaryBox(Photon_Ptr, In_Ptr)) {
			Hop(Photon_Ptr);     /* move to boundary plane. */
			CrossOrNotBox(In_Ptr, Photon_Ptr, Out_Ptr);
		}
		else {
			Hop(Photon_Ptr);
			Drop(In_Ptr, Photon_Ptr, Out_Ptr);
			Spin(In_Ptr->layerspecs[Photon_Ptr->inMat+1].g,
				Photon_Ptr);
		}
	}

}

/***********************************************************
**	Check if the random number is greater than Raman probability.
**  If it is so convert the photon to Raman photon and record 
* * the layer where the laser photon becomes Raman photon.  
****/

void checkForRaman(InputStruct  *  In_Ptr,
	PhotonStruct *  Photon_Ptr,
	OutStruct    *  Out_Ptr)
{
	double randNum = RandomNum();
	double ramanProb = In_Ptr->layerspecs[Photon_Ptr->inMat+1].ramProb;
	if ( randNum < ramanProb && Photon_Ptr->isRaman == 0){
		Photon_Ptr->isRaman = 1;
		if (In_Ptr->ramanScatter == 1){
			Photon_Ptr->justRaman = 1;
		}
		Out_Ptr->RamanPhotons = Out_Ptr->RamanPhotons + 1;
		Photon_Ptr->ramanLayer = Photon_Ptr->inMat + 1;
	}
}

/***********************************************************
****/
void HopDropSpin(InputStruct  *  In_Ptr,
	PhotonStruct *  Photon_Ptr,
	OutStruct    *  Out_Ptr)
{
	short layer = Photon_Ptr->inMat + 1;

	if((In_Ptr->layerspecs[layer].mua == 0.0) 
		&& (In_Ptr->layerspecs[layer].mus == 0.0)) 
		/* glass layer. */
		HopInGlass(In_Ptr, Photon_Ptr, Out_Ptr);
	else {
		HopDropSpinInTissue(In_Ptr, Photon_Ptr, Out_Ptr);
		if (Photon_Ptr->isRaman !=1) {
			checkForRaman(In_Ptr,Photon_Ptr, Out_Ptr);
		}
	}
	if( Photon_Ptr->w < In_Ptr->Wth && !Photon_Ptr->dead) 
		Roulette(Photon_Ptr);
}
