/***********************************************************
** Nanyang Technological University, Singapore 637459.
** 2015.
*
**	main program for Raman Monte Carlo simulation of photon
**	distribution in Cuboid.
*
**	Input/output of data.
*
**	This code modified from the original MCML code
**	developed by Dr Lihong Wang et al.
****/

#include "mcml.h"

/***********************************************************
*	Structure used to check against duplicated file names.
****/
struct NameList {
	char name[STRLEN];
	struct NameList * next;
};

typedef struct NameList NameNode;
typedef NameNode * NameLink;


/***********************************************************
*	Center a string according to the column width.
****/
char * CenterStr(short  Wid,
	char * InStr,
	char * OutStr)
{
	size_t nspaces;	/* number of spaces to be filled */
	/* before InStr. */

	nspaces = (Wid - strlen(InStr))/2;
	if(nspaces<0) nspaces = 0;

	strcpy(OutStr, "");
	while(nspaces--)  strcat(OutStr, " ");

	strcat(OutStr, InStr);

	return(OutStr);
}

/**************************************************************************
*	Center a string according to the column width.
****/
#define COLWIDTH 80
void
	CtrPuts(char *InStr)
{
	short       nspaces;		/* number of spaces to be left-filled. */
	char        outstr[STRLEN];

	nspaces = (COLWIDTH - strlen(InStr)) / 2;
	if (nspaces < 0)
		nspaces = 0;

	strcpy(outstr, "");
	while (nspaces--)
		strcat(outstr, " ");

	strcat(outstr, InStr);

	puts(outstr);
}

/***********************************************************
*	Print some messages before starting simulation.  
*	e.g. author, address, program version, year.
****/
#define COLWIDTH 80
void ShowVersion(char *version)
{
	char str[STRLEN];

	CtrPuts(" ");
	CtrPuts("\Raman Monte Carlo Simulation of Light Transport in Cuboid \n");
	CtrPuts(version);
	CtrPuts(" ");

	CtrPuts("Manojit Pramanik, Ph.D.");
	CtrPuts("Biomedical Imaging Laboratory (BILab)");
	CtrPuts("School of Chemical and Biomedical Engineering");
	CtrPuts("Nanyang Technological University (NTU), Singapore");
	CtrPuts("Email: manojit@ntu.edu.sg");
	CtrPuts("Website:http://www.ntu.edu.sg/home/manojit/ ");

	puts("\n");	
}
#undef COLWIDTH

/***********************************************************
*	Get a filename and open it for reading, retry until 
*	the file can be opened.  '.' terminates the program.
*      
*	If Fname != NULL, try Fname first.
****/
FILE *GetFile(char *Fname)
{
	FILE * file=NULL;
	Boolean firsttime=1;

	do {
		if(firsttime && Fname[0]!='\0') { 
			/* use the filename from command line */
			firsttime = 0;
		}
		else {
			printf("Input filename(or . to exit):");
			scanf("%s", Fname);
			firsttime = 0;
		}

		if(strlen(Fname) == 1 && Fname[0] == '.') 
			exit(1);			/* exit if no filename entered. */

		file = fopen(Fname, "r");
	}  while(file == NULL);

	return(file);
}

/***********************************************************
*	Kill the ith char (counting from 0), push the following 
*	chars forward by one.
****/
void KillChar(size_t i, char * Str)
{
	size_t sl = strlen(Str);

	for(;i<sl;i++) Str[i] = Str[i+1];
}

/***********************************************************
*	Eliminate the chars in a string which are not printing 
*	chars or spaces.
*
*	Spaces include ' ', '\f', '\t' etc.
*
*	Return 1 if no nonprinting chars found, otherwise 
*	return 0.
****/
Boolean CheckChar(char * Str)
{
	Boolean found = 0;	/* found bad char. */
	size_t sl = strlen(Str);
	size_t i=0;

	while(i<sl) 
		if (Str[i]<0 || Str[i]>255)
			nrerror("Non-ASCII file\n");
		else if(isprint(Str[i]) || isspace(Str[i])) 
			i++;
		else {
			found = 1;
			KillChar(i, Str);
			sl--;
		}

		return(found);	
}

/***********************************************************
*	Return 1 if this line is a comment line in which the 
*	first non-space character is "#".
*
*	Also return 1 if this line is space line.
****/
Boolean CommentLine(char *Buf)
{
	size_t spn, cspn;

	spn = strspn(Buf, " \t");	
	/* length spanned by space or tab chars. */

	cspn = strcspn(Buf, "#\n");
	/* length before the 1st # or return. */

	if(spn == cspn) 	/* comment line or space line. */
		return(1);
	else				/* the line has data. */	 
		return(0);		
}

/***********************************************************
*	Skip space or comment lines and return a data line only.
****/
char * FindDataLine(FILE *File_Ptr)
{
	static char buf[STRLEN];	/* LW 1/11/2000. Added static. */

	buf[0] = '\0';
	do {	/* skip space or comment lines. */
		if(fgets(buf, 255, File_Ptr) == NULL)  {
			printf("Incomplete data\n");
			buf[0]='\0';
			break;
		}
		else
			CheckChar(buf);
	} while(CommentLine(buf));

	return(buf);
}

/***********************************************************
*	Skip file version, then read number of runs.
****/
short ReadNumRuns(FILE* File_Ptr)
{
	char buf[STRLEN];
	short n=0;

	FindDataLine(File_Ptr); /* skip file version. */

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') nrerror("Reading number of runs\n");
	sscanf(buf, "%hd",&n);	        
	return(n);
}


/***********************************************************
*	Read the file name and the file format.
*
*	The file format can be either A for ASCII or B for
*	binary.
****/
void ReadFnameFormat(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in file name and format. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading file name and format.\n");
	sscanf(buf, "%s %c", 
		In_Ptr->out_fname, &(In_Ptr->out_fformat) );
	if(toupper(In_Ptr->out_fformat) != 'B') 
		In_Ptr->out_fformat = 'A';
}


/***********************************************************
*	Read the number of photons.
****/
void ReadNumPhotons(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in number of photons. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading number of photons.\n");
	sscanf(buf, "%ld", &In_Ptr->num_photons);
	if(In_Ptr->num_photons<=0) 
		nrerror("Nonpositive number of photons.\n");
}


/***********************************************************
*	Read the members dz and dr.
****/
void ReadDzDr(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in dz, dr. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') nrerror("Reading dz, dr.\n");
	sscanf(buf, "%lf%lf", &In_Ptr->dz, &In_Ptr->dr);
	if(In_Ptr->dz<=0) nrerror("Nonpositive dz.\n");
	if(In_Ptr->dr<=0) nrerror("Nonpositive dr.\n");
}


/***********************************************************
*	Read the members nz, nr, na.
****/
void ReadNzNrNa(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in number of dz, dr, da. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading number of dz, dr, da's.\n");
	sscanf(buf, "%hd%hd%hd", 
		&In_Ptr->nz, &In_Ptr->nr, &In_Ptr->na);
	if(In_Ptr->nz<=0) 
		nrerror("Nonpositive number of dz's.\n");
	if(In_Ptr->nr<=0) 
		nrerror("Nonpositive number of dr's.\n");
	if(In_Ptr->na<=0) 
		nrerror("Nonpositive number of da's.\n");
	In_Ptr->da = 0.5*PI/In_Ptr->na;
}

/***********************************************************
**	Read the members dx and dy.
****/
void ReadDxDy(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in dz, dr. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') nrerror("Reading dx, dy.\n");
	sscanf(buf, "%lf%lf", &In_Ptr->dx, &In_Ptr->dy);
	if(In_Ptr->dx<=0) nrerror("Nonpositive dx.\n");
	if(In_Ptr->dy<=0) nrerror("Nonpositive dy.\n");
}


/***********************************************************
**	Read the members nx, ny.
****/
void ReadNxNy(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in number of dz, dr, da. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading number of dx, dy's.\n");
	sscanf(buf, "%hd%hd", 
		&In_Ptr->nx, &In_Ptr->ny);
	if(In_Ptr->nx<=0) 
		nrerror("Nonpositive number of dx's.\n");
	if(In_Ptr->ny<=0) 
		nrerror("Nonpositive number of dy's.\n");
}

/***********************************************************
**	Read isoRaman, which =1 for Isotropic scatter of Raman photon
****/
void ReadRecordRaman(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Read isoRaman's.\n");
	sscanf(buf, "%d", 
		&In_Ptr->ramanScatter);
	if(In_Ptr->ramanScatter!=0 && In_Ptr->ramanScatter!=1 ) 
		nrerror("Invalid isotropic scatter of Raman\n");
}

/***********************************************************
*	Read the number of layers.
****/
void ReadNumLayers(FILE *File_Ptr, InputStruct *In_Ptr)
{
	char buf[STRLEN];

	/** read in number of layers. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading number of layers.\n");
	sscanf(buf, "%hd", &In_Ptr->num_layers);
	if(In_Ptr->num_layers<=0) 
		nrerror("Nonpositive number of layers.\n");
}


/***********************************************************
*	Read the refractive index n of the ambient.
****/
void ReadAmbient(FILE *File_Ptr, 
	LayerStruct * Layer_Ptr,
	char *side)
{
	char buf[STRLEN], msg[STRLEN];
	double n;

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') {
		sprintf(msg, "Rading n of %s ambient.\n", side);
		nrerror(msg);
	}

	sscanf(buf, "%lf", &n );
	if(n<=0) nrerror("Wrong n.\n");
	Layer_Ptr->n = n;
}


/***********************************************************
*	Read the parameters of one layer.
*
*	Return 1 if error detected.
*	Return 0 otherwise.
*
*	*Z_Ptr is the z coordinate of the current layer, which
*	is used to convert thickness of layer to z coordinates
*	of the two boundaries of the layer.
****/
Boolean ReadOneLayer(FILE *File_Ptr, 
	LayerStruct * Layer_Ptr,
	double *Z_Ptr)
{
	char buf[STRLEN], msg[STRLEN];
	double d, n, mua, mus, g;	/* d is thickness. */

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') return(1);	/* error. */

	sscanf(buf, "%lf%lf%lf%lf%lf", &n, &mua, &mus, &g, &d);
	if(d<0 || n<=0 || mua<0 || mus<0 || g<0 || g>1) 
		return(1);			/* error. */

	Layer_Ptr->n	= n;
	Layer_Ptr->mua = mua;	
	Layer_Ptr->mus = mus;	
	Layer_Ptr->g   = g;
	Layer_Ptr->z0	= *Z_Ptr;
	*Z_Ptr += d;
	Layer_Ptr->z1	= *Z_Ptr;

	return(0);
}

/***********************************************************
*	Read the parameters of one layer at a time.
****/
void ReadLayerSpecs(FILE *File_Ptr, 
	short Num_Layers,
	LayerStruct ** Layerspecs_PP)
{
	char msg[STRLEN];
	short i=0;
	double z = 0.0;	/* z coordinate of the current layer. */

	/* Allocate an array for the layer parameters. */
	/* layer 0 and layer Num_Layers + 1 are for ambient. */
	*Layerspecs_PP = (LayerStruct *)
		malloc((unsigned) (Num_Layers+2)*sizeof(LayerStruct));
	if (!(*Layerspecs_PP)) 
		nrerror("allocation failure in ReadLayerSpecs()");

	ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "top"); 
	for(i=1; i<=Num_Layers; i++)  
		if(ReadOneLayer(File_Ptr, &((*Layerspecs_PP)[i]), &z)) {
			sprintf(msg, "Error reading %hd of %hd layers\n", 
				i, Num_Layers);
			nrerror(msg);
		}
		ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "bottom"); 
}

/***********************************************************
*	Compute the critical angles for total internal
*	reflection according to the relative refractive index
*	of the layer.
*	All layers are processed.
****/
void CriticalAngle( short Num_Layers, 
	LayerStruct ** Layerspecs_PP)
{
	short i=0;
	double n1, n2;

	for(i=1; i<=Num_Layers; i++)  {
		n1 = (*Layerspecs_PP)[i].n;
		n2 = (*Layerspecs_PP)[i-1].n;
		(*Layerspecs_PP)[i].cos_crit0 = n1>n2 ? 
			sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;

		n2 = (*Layerspecs_PP)[i+1].n;
		(*Layerspecs_PP)[i].cos_crit1 = n1>n2 ? 
			sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
	}
}


/***********************************************************
**	Read the parameters of one layer of cuboid.
*
*	Return 1 if error detected.
*	Return 0 otherwise.
*
**	Raman probability is read from input file
*****/
Boolean ReadOneLayerBox(FILE *File_Ptr, 
	LayerStruct * Layer_Ptr)
{
	char buf[STRLEN], msg[STRLEN];
	double n, mua, mus, g, ramProb;

	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') return(1);	/* error. */

	sscanf(buf, "%lf%lf%lf%lf%lf", &n, &mua, &mus, &g, &ramProb);
	if(n<=0 || mua<0 || mus<0 || g<0 || g>1) 
		return(1);			/* error. */

	Layer_Ptr->n	= n;
	Layer_Ptr->mua = mua;	
	Layer_Ptr->mus = mus;	
	Layer_Ptr->g   = g;
	Layer_Ptr->ramProb=ramProb;

	return(0);
}

/***********************************************************
**	Read the length,breadth,height and thickness of box.
****/
void ReadBoxDim(FILE *File_Ptr, BoxStruct * Box_Ptr)
{
	char buf[STRLEN];
	double lt,bt,ht,tn;
	/** read in lt,bt,ht and tn of box. **/
	strcpy(buf, FindDataLine(File_Ptr));
	if(buf[0]=='\0') 
		nrerror("Reading lt,bt,ht and tn of box.\n");
	sscanf(buf, "%lf%lf%lf%lf", &lt,&bt,&ht,&tn);
	if(lt<=0 || bt<=0 || ht<=0 ||tn<=0) 
		nrerror("Nonpositive dimension of the box.\n");
	Box_Ptr->lt=lt;
	Box_Ptr->bt=bt;
	Box_Ptr->ht=ht;
	Box_Ptr->tn=tn;
}
/***********************************************************
**	Function to call the function to read length,breadth,height and thickness of box.
****/
void ReadBoxParm(FILE *File_Ptr, BoxStruct ** Boxspecs_PP)
{
	*Boxspecs_PP = (BoxStruct *)
		malloc((unsigned) (2)*sizeof(BoxStruct));

	ReadBoxDim(File_Ptr,&((*Boxspecs_PP[0])));
}

/***********************************************************
**	Read the parameters of box.
****/
void ReadLayerSpecsBox(FILE *File_Ptr, 
	short Num_Layers,
	BoxStruct ** Boxspecs_PP,
	LayerStruct ** Layerspecs_PP)
{
	char msg[STRLEN];
	short i=0;
	double z = 0.0;	/* z coordinate of the current layer. */


	/* Allocate an array for the layer parameters. */
	/* layer 0 and layer Num_Layers + 1 are for ambient. */

	*Layerspecs_PP = (LayerStruct *)

		malloc((unsigned) (Num_Layers+2)*sizeof(LayerStruct));
	if (!(*Layerspecs_PP)) 
		nrerror("allocation failure in ReadLayerSpecs()");

	ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "top"); 
	for(i=1; i<=Num_Layers; i++)  
		if(ReadOneLayerBox(File_Ptr, &((*Layerspecs_PP)[i]))) {
			sprintf(msg, "Error reading %hd of %hd layers\n", 
				i, Num_Layers);
			nrerror(msg);
		}
		ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "bottom"); 
}


/***********************************************************
*	Read in the input parameters for one run.
****/
void ReadParm(FILE* File_Ptr, InputStruct * In_Ptr)
{
	char buf[STRLEN];

	In_Ptr->Wth = WEIGHT;

	ReadFnameFormat(File_Ptr, In_Ptr);
	ReadNumPhotons(File_Ptr, In_Ptr);
	ReadDzDr(File_Ptr, In_Ptr);
	ReadNzNrNa(File_Ptr, In_Ptr);
	ReadDxDy(File_Ptr, In_Ptr);
	ReadNxNy(File_Ptr, In_Ptr);
	ReadRecordRaman(File_Ptr, In_Ptr);
	In_Ptr->num_layers=2;	/** Container is layer 1 and aterial is layer 2 **/
	ReadBoxParm(File_Ptr, &In_Ptr->boxspecs);
	ReadLayerSpecsBox(File_Ptr, In_Ptr->num_layers,&In_Ptr->boxspecs, &In_Ptr->layerspecs);
	CriticalAngle(In_Ptr->num_layers, &In_Ptr->layerspecs);
}

/***********************************************************
*	Return 1, if the name in the name list.
*	Return 0, otherwise.
****/
Boolean NameInList(char *Name, NameLink List)
{
	while (List != NULL) {
		if(strcmp(Name, List->name) == 0) 
			return(1);
		List = List->next;
	};
	return(0);
}

/***********************************************************
*	Add the name to the name list.
****/
void AddNameToList(char *Name, NameLink * List_Ptr)
{
	NameLink list = *List_Ptr;

	if(list == NULL) {	/* first node. */
		*List_Ptr = list = (NameLink)malloc(sizeof(NameNode));
		strcpy(list->name, Name);
		list->next = NULL;
	}
	else {				/* subsequent nodes. */
		/* Move to the last node. */
		while(list->next != NULL)
			list = list->next;

		/* Append a node to the list. */
		list->next = (NameLink)malloc(sizeof(NameNode));
		list = list->next;
		strcpy(list->name, Name);
		list->next = NULL;
	}
}

/***********************************************************
*	Check against duplicated file names.
*
*	A linked list is set up to store the file names used
*	in this input data file.
****/
Boolean FnameTaken(char *fname, NameLink * List_Ptr)
{
	if(NameInList(fname, *List_Ptr))
		return(1);
	else {
		AddNameToList(fname, List_Ptr);
		return(0);
	}
}

/***********************************************************
*	Free each node in the file name list.
****/
void FreeFnameList(NameLink List)
{
	NameLink next;

	while(List != NULL) {
		next = List->next;
		free(List);
		List = next;
	}
}

/***********************************************************
*	Check the input parameters for each run.
****/
void CheckParm(FILE* File_Ptr, InputStruct * In_Ptr)
{
	short i_run;
	short num_runs;	/* number of independent runs. */
	NameLink head = NULL;
	Boolean name_taken;/* output files share the same */
	/* file name.*/
	char msg[STRLEN];

	num_runs = ReadNumRuns(File_Ptr);
	for(i_run=1; i_run<=num_runs; i_run++)  {
		printf("Checking input data for run %hd\n", i_run);
		ReadParm(File_Ptr, In_Ptr);

		name_taken = FnameTaken(In_Ptr->out_fname, &head);
		if(name_taken) 
			sprintf(msg, "file name %s duplicated.\n", 
			In_Ptr->out_fname);

		free(In_Ptr->layerspecs);
		free(In_Ptr->boxspecs);
		if(name_taken) nrerror(msg);
	}
	FreeFnameList(head);
	rewind(File_Ptr);
}


/***********************************************************
*	Allocate the arrays in OutStruct for one run, and 
*	array elements are automatically initialized to zeros.
****/
void InitOutputData(InputStruct In_Parm, 
	OutStruct * Out_Ptr)
{
	short nz = In_Parm.nz;
	short nr = In_Parm.nr;
	short na = In_Parm.na;
	short nx = In_Parm.nx;
	short ny = In_Parm.ny;
	short nl = In_Parm.num_layers;	
	/* remember to use nl+2 because of 2 for ambient. */

	if(nz<=0 || nr<=0 || na<=0 || nl<=0) 
		nrerror("Wrong grid parameters.\n");

	/* Init pure numbers. */
	Out_Ptr->Rsp = 0.0;

	/* Allocate the arrays and the matrices. */

	Out_Ptr->Rd_ra = AllocMatrix(0,nr-1,0,na-1);
	Out_Ptr->Tt_ra = AllocMatrix(0,nr-1,0,na-1);
	Out_Ptr->RdCon_xy = AllocMatrix(0,nx-1,0,ny-1);
	Out_Ptr->TtCon_xy = AllocMatrix(0,nx-1,0,ny-1);
	Out_Ptr->TtCon_yz = AllocMatrix(0,ny-1,0,nz-1);
	Out_Ptr->TtCon_xz = AllocMatrix(0,nx-1,0,nz-1);
	Out_Ptr->RdMat_xy = AllocMatrix(0,nx-1,0,ny-1);
	Out_Ptr->TtMat_xy = AllocMatrix(0,nx-1,0,ny-1);
	Out_Ptr->TtMat_yz = AllocMatrix(0,ny-1,0,nz-1);
	Out_Ptr->TtMat_xz = AllocMatrix(0,nx-1,0,nz-1);
	Out_Ptr->RamanPhotons=0;

	Out_Ptr->A_rz  = AllocMatrix(0,nr-1,0,nz-1);

	/**	Calculate the box dimensions **/
	In_Parm.boxspecs[0].XMin=-In_Parm.boxspecs[0].lt/2;
	In_Parm.boxspecs[0].YMin=-In_Parm.boxspecs[0].bt/2;
	In_Parm.boxspecs[0].ZMin=0;
	In_Parm.boxspecs[0].XMax=In_Parm.boxspecs[0].lt/2;
	In_Parm.boxspecs[0].YMax=In_Parm.boxspecs[0].bt/2;
	In_Parm.boxspecs[0].ZMax=In_Parm.boxspecs[0].ht;

	In_Parm.boxspecs[1].XMin=-In_Parm.boxspecs[0].lt/2+In_Parm.boxspecs[0].tn;
	In_Parm.boxspecs[1].YMin=-In_Parm.boxspecs[0].bt/2+In_Parm.boxspecs[0].tn;
	In_Parm.boxspecs[1].ZMin=In_Parm.boxspecs[0].tn;
	In_Parm.boxspecs[1].XMax=In_Parm.boxspecs[0].lt/2-In_Parm.boxspecs[0].tn;
	In_Parm.boxspecs[1].YMax=In_Parm.boxspecs[0].bt/2-In_Parm.boxspecs[0].tn;
	In_Parm.boxspecs[1].ZMax=In_Parm.boxspecs[0].ht-In_Parm.boxspecs[0].tn;
}

/***********************************************************
*	Undo what InitOutputData did.
*  i.e. free the data allocations.
****/
void FreeData(InputStruct In_Parm, OutStruct * Out_Ptr)
{
	short nz = In_Parm.nz;
	short nr = In_Parm.nr;
	short na = In_Parm.na;
	short nx = In_Parm.nx;
	short ny = In_Parm.ny;
	short nl = In_Parm.num_layers;	
	/* remember to use nl+2 because of 2 for ambient. */

	free(In_Parm.layerspecs);

	FreeMatrix(Out_Ptr->Rd_ra, 0,nr-1,0,na-1);
	FreeMatrix(Out_Ptr->RdCon_xy, 0,nx-1,0,ny-1);
	FreeMatrix(Out_Ptr->TtCon_xy, 0,nx-1,0,ny-1);
	FreeMatrix(Out_Ptr->TtCon_yz, 0,ny-1,0,nz-1);
	FreeMatrix(Out_Ptr->TtCon_xz, 0,nx-1,0,nz-1);
	FreeMatrix(Out_Ptr->RdMat_xy, 0,nx-1,0,ny-1);
	FreeMatrix(Out_Ptr->TtMat_xy, 0,nx-1,0,ny-1);
	FreeMatrix(Out_Ptr->TtMat_yz, 0,ny-1,0,nz-1);
	FreeMatrix(Out_Ptr->TtMat_xz, 0,nx-1,0,nz-1);

	FreeMatrix(Out_Ptr->A_rz, 0, nr-1, 0,nz-1);

	FreeMatrix(Out_Ptr->Tt_ra, 0,nr-1,0,na-1);
}

/***********************************************************
*	Scale Rd and Tt properly.
*
*	"a" stands for angle alpha.
****
*	Scale Rd(r,a) and Tt(r,a) by
*      (area perpendicular to photon direction)
*		x(solid angle)x(No. of photons).
*	or
*		[2*PI*r*dr*cos(a)]x[2*PI*sin(a)*da]x[No. of photons]
*	or
*		[2*PI*PI*dr*da*r*sin(2a)]x[No. of photons]
****
*	Scale Rd(r) and Tt(r) by
*		(area on the surface)x(No. of photons).
****
*	Scale Rd(a) and Tt(a) by
*		(solid angle)x(No. of photons).
****/
void ScaleRdTt(InputStruct In_Parm, OutStruct *	Out_Ptr)
{
	short nr = In_Parm.nr;
	short na = In_Parm.na;
	double dr = In_Parm.dr;
	double da = In_Parm.da;
	short nx = In_Parm.nx;
	short ny = In_Parm.ny;
	double dx = In_Parm.dx;
	double dy = In_Parm.dy;
	short nz = In_Parm.nz;
	double dz = In_Parm.dz;
	short ir,ia,ix,iy,iz;
	double scale1,scale2,scale3;

	scale1 = 4.0*PI*PI*dr*sin(da/2)*dr*In_Parm.num_photons;
	/* The factor (ir+0.5)*sin(2a) to be added. */

	for(ir=0; ir<nr; ir++)  
		for(ia=0; ia<na; ia++) {
			scale2 = 1.0/((ir+0.5)*sin(2.0*(ia+0.5)*da)*scale1);
			Out_Ptr->Rd_ra[ir][ia] *= scale2;
			Out_Ptr->Tt_ra[ir][ia] *= scale2;
		}

		scale3 = (1.0)/(In_Parm.num_photons);
		for(ix=0; ix<nx; ix++)  
			for(iy=0; iy<ny; iy++) {
				Out_Ptr->RdCon_xy[ix][iy] *= scale3;
				Out_Ptr->TtCon_xy[ix][iy] *= scale3;
				Out_Ptr->RdMat_xy[ix][iy] *= scale3;
				Out_Ptr->TtMat_xy[ix][iy] *= scale3;
			}
			scale3 = (1.0)/(In_Parm.num_photons);
			for(iy=0; iy<ny; iy++)  
				for(iz=0; iz<nz; iz++) {
					Out_Ptr->TtCon_yz[iy][iz] *= scale3;
					Out_Ptr->TtMat_yz[iy][iz] *= scale3;
				}
				scale3 = (1.0)/(In_Parm.num_photons);
				for(ix=0; ix<nx; ix++)  
					for(iz=0; iz<nz; iz++) {
						Out_Ptr->TtCon_xz[ix][iz] *= scale3;
						Out_Ptr->TtMat_xz[ix][iz] *= scale3;
					}
}

/***********************************************************
*	Scale absorption arrays properly.
****/
void ScaleA(InputStruct In_Parm, OutStruct * Out_Ptr)
{
	short nz = In_Parm.nz;
	short nr = In_Parm.nr;
	double dz = In_Parm.dz;
	double dr = In_Parm.dr;
	short nl = In_Parm.num_layers;
	short iz,ir;
	short il;
	double scale1;

	/* Scale A_rz. */
	scale1 = 2.0*PI*dr*dr*dz*In_Parm.num_photons;	
	/* volume is 2*pi*(ir+0.5)*dr*dr*dz.*/ 
	/* ir+0.5 to be added. */
	for(iz=0; iz<nz; iz++) 
		for(ir=0; ir<nr; ir++) 
			Out_Ptr->A_rz[ir][iz] /= (ir+0.5)*scale1;
}

/***********************************************************
*	Sum and scale results of current run.
****/
void SumScaleResult(InputStruct In_Parm, 
	OutStruct * Out_Ptr)
{
	ScaleRdTt(In_Parm, Out_Ptr);
	ScaleA(In_Parm, Out_Ptr);
}

/***********************************************************
*	Write the version number as the first string in the 
*	file.
*	Use chars only so that they can be read as either 
*	ASCII or binary.
****/
void WriteVersion(FILE *file, char *Version)
{
	fprintf(file, 
		"%s \t# Version number of the file format.\n\n", 
		Version);
	fprintf(file, "####\n# Data categories include: \n");
	fprintf(file, "# Rd_ra,RdCon_xy,TtCon_xy,TtCon_yz,TtCon_xz,\nRdMat_xy,TtMat_xy,TtMat_yz,TtMat_xz,Tt_ra \n####\n\n");
}

/***********************************************************
*	Write the input parameters to the file.
****/
void WriteInParm(FILE *file, InputStruct In_Parm)
{
	short i;

	fprintf(file, 
		"InParm \t\t\t# Input parameters. cm is used.\n");

	fprintf(file, 
		"%s \tA\t\t# output file name, ASCII.\n", 
		In_Parm.out_fname);
	fprintf(file, 
		"%ld \t\t\t# No. of photons\n", In_Parm.num_photons);

	fprintf(file, 
		"%G\t%G\t\t# dz, dr [cm]\n", In_Parm.dz,In_Parm.dr);
	fprintf(file, "%hd\t%hd\t%hd\t# No. of dz, dr, da.\n", 
		In_Parm.nz, In_Parm.nr, In_Parm.na);
	fprintf(file, 
		"%G\t%G\t\t# dx, dy [cm]\n", In_Parm.dx,In_Parm.dy);
	fprintf(file, "%hd\t%hd\t\t# No. of dx, dy.\n", 
		In_Parm.nx, In_Parm.ny);
	fprintf(file, "%d\t\t\t# Isotropic scatter of Raman.\n\n", 
		In_Parm.ramanScatter);
	fprintf(file, 
		"%G\t%G\t%G\t%G\t\t# lt, bt, ht, tn of box[cm]\n", In_Parm.boxspecs[0].lt,In_Parm.boxspecs[0].bt,In_Parm.boxspecs[0].ht,In_Parm.boxspecs[0].tn);

	fprintf(file, 
		"#n\tmua\tmus\tg\tramProb\t# One line for each layer\n"); 
	fprintf(file, 
		"%G\t\t\t\t\t# n for medium above\n", 
		In_Parm.layerspecs[0].n); 
	for(i=1; i<=In_Parm.num_layers; i++)  {
		LayerStruct s;
		s = In_Parm.layerspecs[i];
		fprintf(file, "%G\t%G\t%G\t%G\t%G\t# layer %hd\n",
			s.n, s.mua, s.mus, s.g, s.ramProb, i);
	}
	fprintf(file, "%G\t\t\t\t\t# n for medium below\n\n", 
		In_Parm.layerspecs[i].n); 
}


/***********************************************************
*	5 numbers each line.
****/
void WriteRd_ra(FILE * file, 
	short Nr,
	short Na,
	OutStruct Out_Parm)
{
	short ir, ia;

	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Rd[r][angle]. [1/(cm2sr)].",
		"# Rd[0][0], [0][1],..[0][na-1]",
		"# Rd[1][0], [1][1],..[1][na-1]",
		"# ...",
		"# Rd[nr-1][0], [nr-1][1],..[nr-1][na-1]",
		"Rd_ra");

	for(ir=0;ir<Nr;ir++)
		for(ia=0;ia<Na;ia++)  {
			fprintf(file, "%12.4E \n", Out_Parm.Rd_ra[ir][ia]);
		}

		fprintf(file, "\n");
}

/***********************************************************
**  WriteRd_xy .
****/
void WriteRT_xyz(FILE * file, 
	short Nx,
	short Ny,
	short Nz,
	OutStruct Out_Parm)
{
	short ix, iy,iz;

	fprintf(file,"\nTotal Raman photons = %ld\n",Out_Parm.RamanPhotons);
	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Rd[x][y]. [1/(cm2sr)].",
		"# Rd[0][0], [0][1],..[0][nx-1]",
		"# Rd[1][0], [1][1],..[1][nx-1]",
		"# ...",
		"# Rd[ny-1][0], [ny-1][1],..[nx-1][ny-1]",
		"Rd_xy");
	fprintf(file,"\n%d\t%d\t\t#nx,ny\n",Nx,Ny);
	fprintf(file,"RdCon\n");
	for(ix=0;ix<Nx;ix++) {
		for(iy=0;iy<Ny;iy++)  {
			fprintf(file, "%12.4E\n", Out_Parm.RdCon_xy[ix][iy]);
		}
	}
	fprintf(file, "\n");
	fprintf(file,"RdMat\n");
	for(ix=0;ix<Nx;ix++) {
		for(iy=0;iy<Ny;iy++)  {
			fprintf(file, "%12.4E\n", Out_Parm.RdMat_xy[ix][iy]);
		}
	}
	fprintf(file, "\n");
	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Tt[x][y]. [1/(cm2sr)].",
		"# Tt[0][0], [0][1],..[0][nx-1]",
		"# Tt[1][0], [1][1],..[1][nx-1]",
		"# ...",
		"# Tt[ny-1][0], [ny-1][1],..[nx-1][ny-1]",
		"Tt_xy");
	fprintf(file,"\n%d\t%d\t\t#nx,ny\n",Nx,Ny);
	fprintf(file,"TtCon_xy\n");
	for(ix=0;ix<Nx;ix++){
		for(iy=0;iy<Ny;iy++)  {
			fprintf(file, "%12.4E\n", Out_Parm.TtCon_xy[ix][iy]);
		}
	}
	fprintf(file, "\n");
	fprintf(file,"TtMat_xy\n");
	for(ix=0;ix<Nx;ix++) {
		for(iy=0;iy<Ny;iy++)  {
			fprintf(file, "%12.4E\n", Out_Parm.TtMat_xy[ix][iy]);
		}
	}
	fprintf(file, "\n");
	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Tt[y][z]. [1/(cm2sr)].",
		"# Tt[0][0], [0][1],..[0][ny-1]",
		"# Tt[1][0], [1][1],..[1][ny-1]",
		"# ...",
		"# Tt[nz-1][0], [nz-1][1],..[ny-1][nz-1]",
		"Tt_yz");
	fprintf(file,"\n%d\t%d\t\t#ny,nz\n",Ny,Nz);
	fprintf(file, "TtCon_yz\n");
	for(iy=0;iy<Ny;iy++){
		for(iz=0;iz<Nz;iz++)  {
			fprintf(file, "%12.4E\n", Out_Parm.TtCon_yz[iy][iz]);
		}
	}
	fprintf(file, "TtMat_yz\n");
	fprintf(file, "\n");
	for(iy=0;iy<Ny;iy++) {
		for(iz=0;iz<Nz;iz++)  {
			fprintf(file, "%12.4E\n", Out_Parm.TtMat_yz[iy][iz]);
		}
	}
	fprintf(file, "\n");
	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Tt[x][z]. [1/(cm2sr)].",
		"# Tt[0][0], [0][1],..[0][nz-1]",
		"# Tt[1][0], [1][1],..[1][nz-1]",
		"# ...",
		"# Tt[nz-1][0], [nz-1][1],..[nx-1][nz-1]",
		"Tt_xz");
	fprintf(file,"\n%d\t%d\t\t#nx,nz\n",Nx,Nz);
	fprintf(file, "TtCon_xz\n");
	for(ix=0;ix<Nx;ix++) {
		for(iz=0;iz<Nz;iz++)  {
			fprintf(file, "%12.4E\n", Out_Parm.TtCon_xz[ix][iz]);
		}
	}
	fprintf(file, "\n");
	fprintf(file, "TtMat_xz\n");
	for(ix=0;ix<Nx;ix++) {
		for(iz=0;iz<Nz;iz++)  {
			fprintf(file, "%12.4E\n", Out_Parm.TtMat_xz[ix][iz]);
		}
	}
	fprintf(file, "\n");
}



/***********************************************************
*	5 numbers each line.
****/
void WriteTt_ra(FILE * file, 
	short Nr,
	short Na,
	OutStruct Out_Parm)
{
	short ir, ia;

	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
		"# Tt[r][angle]. [1/(cm2sr)].",
		"# Tt[0][0], [0][1],..[0][na-1]",
		"# Tt[1][0], [1][1],..[1][na-1]",
		"# ...",
		"# Tt[nr-1][0], [nr-1][1],..[nr-1][na-1]",
		"Tt_ra");

	for(ir=0;ir<Nr;ir++)
		for(ia=0;ia<Na;ia++)  {
			fprintf(file, "%12.4E ", Out_Parm.Tt_ra[ir][ia]);
			if( (ir*Na + ia + 1)%5 == 0) fprintf(file, "\n");
		}

		fprintf(file, "\n");
}

/***********************************************************
*	5 numbers each line.
****/
void WriteA_rz(FILE * file, 
	short Nr,
	short Nz,
	OutStruct Out_Parm)
{
	short iz, ir;

	fprintf(file, 
		"%s\n%s\n%s\n%s\n%s\n%s\n", /* flag. */
		"# A[r][z]. [1/cm3]",
		"# A[0][0], [0][1],..[0][nz-1]",
		"# A[1][0], [1][1],..[1][nz-1]",
		"# ...",
		"# A[nr-1][0], [nr-1][1],..[nr-1][nz-1]",
		"A_rz");

	for(ir=0;ir<Nr;ir++)
		for(iz=0;iz<Nz;iz++)  {
			fprintf(file, "%12.4E ", Out_Parm.A_rz[ir][iz]);
			if( (ir*Nz + iz + 1)%5 == 0) fprintf(file, "\n");
		}

		fprintf(file, "\n");
}

/***********************************************************
****/
void WriteResult(InputStruct In_Parm, 
	OutStruct Out_Parm,
	char * TimeReport)
{
	FILE *file;

	file = fopen(In_Parm.out_fname, "w");
	if(file == NULL) nrerror("Cannot open file to write.\n");

	if(toupper(In_Parm.out_fformat) == 'A') 
		WriteVersion(file, "A1");
	else 
		WriteVersion(file, "B1");

	fprintf(file, "# %s", TimeReport);
	fprintf(file, "\n");

	WriteInParm(file, In_Parm);

	/* 2D arrays. */
	//WriteA_rz(file, In_Parm.nr, In_Parm.nz, Out_Parm);
	//WriteRd_ra(file, In_Parm.nr, In_Parm.na, Out_Parm);
	WriteRT_xyz(file, In_Parm.nx, In_Parm.ny, In_Parm.nz, Out_Parm);	/** Reflectance or transmittance is 6 planes **/
	//WriteTt_ra(file, In_Parm.nr, In_Parm.na, Out_Parm);

	fclose(file);
}
