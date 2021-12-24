// Type: real
typedef double real;

// Type: real
typedef long integer;

// Type: array
struct realArray
{
	long n;
	double *v;
};

void makeRealArray ( realArray *a, const integer n )
{
	a->v = (double *) malloc(n * sizeof(double));
}

void freeRealArray ( realArray *a )
{
	free( a->v );
}

// Type: vector
struct Nvector
{
	double *x, *y, *z;
};

void makeVector ( struct Nvector *v, const integer n )
{
	//v->n = n;
	//v->nMax = n;
	v->x = (double *) malloc(n * sizeof(double));
	v->y = (double *) malloc(n * sizeof(double));
	v->z = (double *) malloc(n * sizeof(double));
	//zeroVector ( v, n );
}

void freeVector ( struct Nvector *v )
{
	free( v->x );
	free( v->y );
	free( v->z );
}

/*
void EsperaEnter()  // Definicao da funcao "EsperaEnter"
{
    int tecla;
    printf("Pressione ENTER\n");
    do
    {
        tecla = getch();
        
    } while(tecla != 13); // 13 e' o codigo ASCII do ENTER
}
*/
void printVTK (std::ofstream &file, Nvector r, Nvector f, Nvector df, realArray pnd, realArray pndNew,
	Nvector grad, Nvector grad1, Nvector gradCPM, realArray divU, realArray divMPS, realArray div1,
	realArray d2f, realArray lapMPS, realArray lapCPM, realArray lap2, realArray lapLDD, int nPart)
{
	//int flag = 1;
	
	file << "# vtk DataFile Version 2.0";
	file << "\nvtk_file";
	file << "\nASCII";
	file << "\nDATASET POLYDATA";
	
	file << "\nPOINTS " << nPart << " float\n";	
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << r.x[i] << " " << r.y[i] << " " << r.z[i];
		//fprintf (output, "\t%lf %lf %lf", r.x[i], r.y[i], r.z[i]);
	
	file << "\nPOINT_DATA " << nPart;
	//fprintf(output, "\nPOINT_DATA %d", nPart);
	
	file << "\nVECTORS Function float\n";
	//fprintf(output, "\nVECTORS Function float\n");
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << f.x[i] << " " << f.y[i] << " " << f.z[i];
		//fprintf (output, "\t%lf %lf %lf", f.x[i], f.y[i], f.z[i]);
		
	file << "\nVECTORS Gradient float\n";
	//fprintf(output, "\nVECTORS Gradient float\n");
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << df.x[i] << " " << df.y[i] << " " << df.z[i];
		//fprintf (output, "\t%lf %lf %lf", df.x[i], df.y[i], df.z[i]);  
	
	file << "\nSCALARS PND float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << pnd.v[i];
	
	file << "\nSCALARS PNDnew float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << pndNew.v[i];
	
	file << "\nVECTORS GradMPS float\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << grad.x[i] << " " << grad.y[i] << " " << grad.z[i];
        
        file << "\nVECTORS Grad1stOrd float\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << grad1.x[i] << " " << grad1.y[i] << " " << grad1.z[i];
	
	file << "\nVECTORS GradCPM float\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << gradCPM.x[i] << " " << gradCPM.y[i] << " " << gradCPM.z[i];
	
	file << "\nSCALARS Divergence float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << divU.v[i];
		
	file << "\nSCALARS DivMPS float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << divMPS.v[i];
	
	file << "\nSCALARS Div1stOrd float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << div1.v[i];
		
	file << "\nSCALARS Laplacian float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << d2f.v[i];
	
	file << "\nSCALARS LapMPS float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << lapMPS.v[i];
	
	file << "\nSCALARS LapCPM float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << lapCPM.v[i];
	
	file << "\nSCALARS Lap2 float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << lap2.v[i];

	file << "\nSCALARS LapLDD float\nLOOKUP_TABLE default\n";
	
	for (int i = 0; i < nPart; i++)
		file << "\t" << lapLDD.v[i];
}

void printVTK2 (FILE *fp, Nvector r, realArray f, Nvector u, Nvector df, realArray pnd,
        realArray pndNew, realArray bc, Nvector normal, Nvector grad, Nvector grad1, 
        Nvector gradCPM, realArray divf, realArray divMPS, realArray div1, realArray div1New, realArray d2f, 
        realArray lapMPS, realArray lapMPS1, realArray lapCPM, realArray lap2, realArray lapLDD, int nPart)
{
	//int flag = 1;
        
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "vtk_file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET POLYDATA\n");
	
	fprintf(fp,"POINTS %d float\n",nPart);
	
        for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f %f %f", r.x[i], r.y[i], r.z[i]);
	
	fprintf(fp, "\nPOINT_DATA %d", nPart);
	
        fprintf(fp, "\nSCALARS Function float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", f.v[i]);
        
	fprintf(fp, "\nVECTORS U float\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f %f %f", u.x[i], u.y[i], u.z[i]);
		
	fprintf(fp, "\nVECTORS Gradient float\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f %f %f", df.x[i], df.y[i], df.z[i]);  
	
	fprintf(fp, "\nSCALARS PND float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", pnd.v[i]);
	
	fprintf(fp, "\nSCALARS PNDnew float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", pndNew.v[i]);
	
        fprintf(fp, "\nSCALARS BC float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", bc.v[i]);
        
        fprintf(fp, "\nVECTORS Normal float\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f %f %f", normal.x[i], normal.y[i], normal.z[i]);
        
	fprintf(fp, "\nVECTORS GradMPS float\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f %f %f", grad.x[i], grad.y[i], grad.z[i]);
        
        fprintf(fp, "\nVECTORS Grad1stOrd float\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f %f %f", grad1.x[i], grad1.y[i], grad1.z[i]);
	
	fprintf(fp, "\nVECTORS GradCPM float\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f %f %f", gradCPM.x[i], gradCPM.y[i], gradCPM.z[i]); 
	
	fprintf(fp, "\nSCALARS Divergence float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", divf.v[i]);
	
	fprintf(fp, "\nSCALARS DivMPS float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", divMPS.v[i]);	
	
	fprintf(fp, "\nSCALARS Div1stOrd float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", div1.v[i]);
        
        fprintf(fp, "\nSCALARS Div1stOrdNew float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", div1New.v[i]);
		
	fprintf(fp, "\nSCALARS Laplacian float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", d2f.v[i]);
		
	fprintf(fp, "\nSCALARS LapMPS float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", lapMPS.v[i]);
        
        fprintf(fp, "\nSCALARS LapMPS1stOrd float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", lapMPS1.v[i]);
		
	fprintf(fp, "\nSCALARS LapCPM float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", lapCPM.v[i]);
		
	fprintf(fp, "\nSCALARS Lap2 float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", lap2.v[i]);

	fprintf(fp, "\nSCALARS LapLDD float\nLOOKUP_TABLE default\n");
	
	for (int i = 0; i < nPart; i++)
		fprintf (fp, "\t%f", lapLDD.v[i]);
}