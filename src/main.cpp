/* 
 * File:   main.cpp
 * Author: Rubens
 *
 * Created on 29 de Outubro de 2015, 19:51
 */

// Gradiente and Laplacian operators

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
//#include <conio.h>
//#include <dos.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "operators.h"
#include "invMatrix.h"

using namespace std;

int nPart;
double eps = 1.0e-8;

double weight(double r, double re, int wijType) {
	switch (wijType) {
		case 0:
			return re/r - 1.0;
		case 1:
			return re/r + r/re - 2.0;
		case 2:
			return re/r - r/re;
		case 3:
			return (1.0-r/re)*(1.0-r/re)*(1.0-r/re);
		case 4:
			return (1.0-r/re)*(1.0-r/re);
		default:
			return re/r - 1.0;
	}
}

void initVals(	struct Nvector *r, struct Nvector *rn, struct Nvector *u, struct Nvector *df, struct Nvector *gradMPS, struct Nvector *gradCPM, struct Nvector *gradMPSadd, struct Nvector *gradMPS1, struct Nvector *normal,
				struct realArray *f, struct realArray *divU, struct realArray *d2f, struct realArray *pndS, struct realArray *pndL, struct realArray *pndNew, struct realArray *lambda, struct realArray *lambda2, 
				struct realArray *lambda3, struct realArray *divMPS, struct realArray *divMPS1, struct realArray *lapMPS, struct realArray *lapMPS1, struct realArray *lapCPM, struct realArray *lap2, struct realArray *lapLDD,
				struct realArray *bc, struct realArray *divMPS1New, double dp, int x, int y, int order, int irreg) {
	int n = 0;
	for (int j = 0; j <= y; j++) {
	for (int i = 0; i <= x; i++, n++) {
		// Particles position
		r->x[n] = i * dp;
		r->y[n] = j * dp;
		rn->x[n] = r->x[n];
		rn->y[n] = r->y[n];
		double e = rand() % 50;
		if (irreg == 0) e = 0;
		rn->x[n] += (e/500.0) * dp;
		rn->y[n] += (e/500.0) * dp;
		r->z[n] = rn->z[n] = 0.0;

		// Avoid position y = zero
		if(fabs(rn->y[n])<eps && irreg == 1) rn->y[n] += 0.01 * dp;

		double flagX = 0.0;
		double flagY = 1.0;
		// Scalar Function
		f->v[n] = flagX * pow(rn->x[n], order) + flagY * pow(rn->y[n],order);

		// Analytical gradient
		df->x[n] = flagX * double(order) * pow(rn->x[n], order-1);
		df->y[n] = flagY * double(order) * pow(rn->y[n], order-1);

		// Analytical Laplacian
		if (order == 1) {
			d2f->v[n] = 0.0;
		}
		else {
			d2f->v[n] =	flagX * double(order) * (double(order)-1.0) * pow(rn->x[n], order-2) +
						flagY * double(order) * (double(order)-1.0) * pow(rn->y[n], order-2);
		}

		flagX = 1.0;
		flagY = 1.0;
		// Vectorial Function
		u->x[n] = flagX * pow(rn->x[n], order);
		u->y[n] = flagY * pow(rn->y[n], order);

		// Analytical divergence
		//divf->v[n] = double(order)*pow(rn->x[n],order-1) + double(order)*pow(rn->y[n],order-1);
		divU->v[n] = 	flagX * double(order) * pow(rn->x[n], order-1) + 
						flagY * double(order) * pow(rn->y[n], order-1);
		
		// Set variables to ZERO
		pndS->v[n] = pndL->v[n] = pndNew->v[n] = lambda->v[n] = lambda2->v[n] = lambda3->v[n] = 
		divMPS->v[n] = divMPS1->v[n] = lapMPS->v[n] = lapMPS1->v[n] = lapCPM->v[n] =
		lap2->v[n] = lapLDD->v[n] = bc->v[n] = divMPS1New->v[n] = 0.0;

		gradMPS->x[n] = gradMPS->y[n] = gradMPS->z[n] = gradCPM->x[n] = gradCPM->y[n] = gradCPM->z[n] = 
		gradMPSadd->x[n] = gradMPSadd->y[n] = gradMPSadd->z[n] = gradMPS1->x[n] = gradMPS1->y[n] = gradMPS1->z[n] =
		normal->x[n] = normal->y[n] = normal->z[n] = 0.0;
	}}
}

void calcPnd(struct realArray *pndS, struct realArray *pndL, double &pndS0, double &pndL0, struct Nvector r, double reS2, double reL2, int wijType) {
	for (int i = 0; i < nPart; i++) {
		for (int j = 0; j < nPart; j++) {
			if (i != j) {
				double dx, dy, dist, dist2, reS, reL;
				dx = r.x[j] - r.x[i];
				dy = r.y[j] - r.y[i];
				dist2 = dx*dx + dy*dy;
				if (dist2<=eps) {
					printf("dist = 0 !!! i = %d, j =%d", i, j);
					//EsperaEnter();
				}
				else if (dist2 < reL2) {
					dist = sqrt(dist2);
					reL = sqrt(reL2);
					pndL->v[i] += weight(dist, reL, wijType);
					if (dist2 < reS2) {
						reS = sqrt(reS2);
						pndS->v[i] += weight(dist, reS, wijType);
					}
				}
			}
		}

		if (pndS->v[i] > pndS0)
			pndS0 = pndS->v[i];
		if (pndL->v[i] > pndL0)
			pndL0 = pndL->v[i];
	}
}

void pndCorrection(struct realArray *pndNew, struct Nvector r, double reS2, int wijType) {
	for (int i=0; i<nPart; i++) {
		double dx, dy, dist, dist2, reS, beta0, beta1, beta2, wij;

		double** Aux = new double* [3];
		for (int k = 0; k < 3; k++)
			Aux[k] = new double [3];

		double** Coef = new double* [3];
		for (int k = 0; k < 3; k++)
			Coef[k] = new double [3];

		for (int ii = 0; ii < 3; ii++) {
		for (int jj = 0; jj < 3; jj++) {
			Aux[ii][jj] = 0.0;
			Coef[ii][jj] = 0.0;
		}}

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);
				Aux[0][0] += wij*1;     Aux[0][1] += wij*dx;    Aux[0][2] += wij*dy;
				Aux[1][0] += wij*dx;    Aux[1][1] += wij*dx*dx; Aux[1][2] += wij*dx*dy;
				Aux[2][0] += wij*dy;    Aux[2][1] += wij*dx*dy; Aux[2][2] += wij*dy*dy;
			}
		}}

		MatrixInversion(Aux, 3, Coef);

		beta0 = Coef[0][0];
		beta1 = Coef[1][0];
		beta2 = Coef[2][0];
		double valor = 0.0;
		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);
				valor = (beta0 + beta1*dx + beta2*dy) * wij;
				printf("v = %f ", valor);
				pndNew->v[i] += (beta0 + beta1*dx + beta2*dy) * wij;
			}
		}}
		//double valor = pndNew->v[i];
		printf("\n\n");
	}
}

void calcLambda(struct realArray *lambda1, struct realArray *lambda2, struct realArray *lambda3, 
	double &lambda0, struct Nvector r, double reS2, double reL2, int wijType) {
	for (int i = 0; i < nPart; i++) {
		double dx, dy, dist, dist2, reS, reL, wij, wr, wr2, wr3, wtS, wtL;
		wr = wr2 = wr3 = wtS = wtL = 0.0;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < eps) {
				printf("dist = 0 !!! i = %d, j =%d", i, j);
				//EsperaEnter();
			}
			else if (dist2 < reL2) {
				dist = sqrt(dist2);
				reL = sqrt(reL2);
				wij = weight(dist, reL, wijType);
				wr2 += wij * dist2;
				wr3 += wij * dist2 * dist;
				wtL += wij;
			}
			else if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);
				wr  += wij * dist;
				wtS += wij;
			}
		}}
		lambda1->v[i] = wr2/wtL;
		lambda2->v[i] = wr/wtS;
		lambda3->v[i] = wr3;

		if (lambda1->v[i] > lambda0)
			lambda0 = lambda1->v[i];
	}
}

void calcGradient(struct Nvector *grad, struct Nvector *gradAdd, struct Nvector *grad1, struct Nvector r, 
	struct realArray f, struct realArray pnd, double pnd0, int dim, double reS2, int wijType) {
	double C1 = dim / pnd0;
	for (int i = 0; i < nPart; i++) {
		double dx, dy, dist, dist2, reS, wij;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);

				// Subtraction
				grad->x[i] += ( f.v[j] - f.v[i] ) * dx * wij / dist2;
				grad->y[i] += ( f.v[j] - f.v[i] ) * dy * wij / dist2;

				//grad->x[i] += ( pnd.v[i] / pnd.v[j] * f.v[j] - pnd.v[j] / pnd.v[i] * f.v[i] ) * dx * wij / dist2;
				//grad->y[i] += ( pnd.v[i] / pnd.v[j] * f.v[j] - pnd.v[j] / pnd.v[i] * f.v[i] ) * dy * wij / dist2;

				// Addition
				gradAdd->x[i] += ( f.v[j] + f.v[i] ) * dx * wij / dist2;
				gradAdd->y[i] += ( f.v[j] + f.v[i] ) * dy * wij / dist2;

				// Addition
				//gradAdd->x[i] += ( pnd.v[i] / pnd.v[j] * f.v[j] + pnd.v[j] / pnd.v[i] * f.v[i] ) * dx * wij / dist2;
				//gradAdd->y[i] += ( pnd.v[i] / pnd.v[j] * f.v[j] + pnd.v[j] / pnd.v[i] * f.v[i] ) * dy * wij / dist2;
			}
		}}

		grad->x[i] *= C1;
		grad->y[i] *= C1;
		gradAdd->x[i] *= C1;
		gradAdd->y[i] *= C1;
		grad1->x[i] = grad->x[i];
		grad1->y[i] = grad->y[i];
		//grad1->x[i] = gradAdd->x[i];
		//grad1->y[i] = gradAdd->y[i];
	}
}

void calcDivergence(struct realArray *div, struct Nvector r, struct Nvector u,
	double pnd0, int dim, double reS2, int wijType) {
	double C1 = dim / pnd0;
	for (int i = 0; i < nPart; i++) {
		double dx, dy, dist, dist2, reS, wij;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);
				div->v[i] += ( (u.x[j] - u.x[i]) * dx + (u.y[j] - u.y[i]) * dy ) * wij / dist2;
			}
		}}

		div->v[i] *= C1;
	}
}

void calcLaplacian(struct realArray *lap, struct Nvector r, struct realArray f,
	struct realArray lambda, double pnd0, double lambda0, int dim, double reL2, int wijType)
{
	double C1 = 2.0 * dim / ( pnd0 * lambda0 );
	for (int i = 0; i < nPart; i++) {
		double dx, dy, dist, dist2, reL, wij;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reL2) {
				dist = sqrt(dist2);
				reL = sqrt(reL2);
				wij = weight(dist, reL, wijType);

				//lap->v[i] += 2.0 * dim * ( f.v[j] - f.v[i] ) * wij / ( pnd0 * lambda.v[i] );
				lap->v[i] += ( f.v[j] - f.v[i] ) * wij;
				//lap->v[i] += 2.0 * dim * ( f.x[j] - f.x[i] ) * wij / ( pnd.v[i] * lambdaA );
				/// Laplaciano n = 3
				//lap->v[i] += 2.0 * dim * ( f.x[j] - f.x[i] ) * dist * wij / lambda3.v[i];
				// Gotoh Khayyer
				//lap->v[i] += 3.0 * ( f.x[j] - f.x[i] ) * reL / ( pnd.v[i] * dist2 * dist );
			}
		}}
		lap->v[i] *= C1;
	}
}

void gradCorrection(struct Nvector *grad, struct Nvector r, struct realArray pnd, double pnd0, 
	int dim, double dp, double reS2, int wijType) {
	double C1 = dim / pnd0;
	for (int i=0; i<nPart; i++) {
		double xxij, yyij, xyij, detC, gradx, grady, wij, dx, dy, dist, dist2, reS;// dx0, dy0;
		double corrGrad[3];
		corrGrad[0]=corrGrad[1]=corrGrad[2]=0.0;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);

				corrGrad[0] += ( dx * dx / dist2 ) * wij;
				corrGrad[1] += ( dx * dy / dist2 ) * wij;
				corrGrad[2] += ( dy * dy / dist2 ) * wij;
			}
		}}
		corrGrad[0] *= C1;
		corrGrad[1] *= C1;
		corrGrad[2] *= C1;

		//detC = ( corrGrad[0] * corrGrad[2] - corrGrad[1] * corrGrad[1] ) * dim / pnd.v[i];
		detC = ( corrGrad[0] * corrGrad[2] - corrGrad[1] * corrGrad[1] );// * dim / pnd0;
		if(fabs(detC) < eps) fprintf(stderr, "detC[i=%d] = %.4e\n", i, detC);
		if(fabs(detC) >= eps) {
			xxij =   corrGrad[2]/detC;
			yyij =   corrGrad[0]/detC;
			xyij = - corrGrad[1]/detC;
			//fprintf(stderr, "xy = %f, %f, %f\n", xij2, xyij, yij2);
			gradx = grad->x[i] * xxij + grad->y[i] * xyij;
			grady = grad->x[i] * xyij + grad->y[i] * yyij;
			grad->x[i] = gradx;
			grad->y[i] = grady;
			//fprintf(stderr, "OK");
		}
	}
}

void divCorrection(struct realArray *div1, struct Nvector u, struct Nvector r, double dp, double reS2, int wijType) {
	for (int i=0; i<nPart; i++) {
		double xxij, yyij, xyij, yxij, detC, wij, dx, dy, dist, dist2, reS;// dx0, dy0;
		double corrDiv[4], Div[4];//
		corrDiv[0]=corrDiv[1]=corrDiv[2]=corrDiv[3]=0.0;
		Div[0]=Div[1]=Div[2]=Div[3]=0.0;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);

				corrDiv[0] += ( dx * dx / dist2 ) * wij;
				corrDiv[1] += ( dx * dy / dist2 ) * wij;
				corrDiv[2] += ( dy * dx / dist2 ) * wij;
				corrDiv[3] += ( dy * dy / dist2 ) * wij;

				Div[0] += ( dx * (u.x[j] - u.x[i]) / dist2) * wij;
				Div[1] += ( dx * (u.y[j] - u.y[i]) / dist2) * wij;
				Div[2] += ( dy * (u.x[j] - u.x[i]) / dist2) * wij;
				Div[3] += ( dy * (u.y[j] - u.y[i]) / dist2) * wij;
			}
		}}

		detC = corrDiv[0] * corrDiv[3] - corrDiv[2] * corrDiv[1];
		if(fabs(detC) < eps) fprintf(stderr, "detC[i=%d] = %.4e\n", i, detC);
		if(fabs(detC) >= eps) {
			xxij =   corrDiv[3] / detC;
			xyij = - corrDiv[1] / detC;
			yxij = - corrDiv[2] / detC;
			yyij =   corrDiv[0] / detC;
			div1->v[i] = xxij * Div[0] + xyij * Div[2] + yxij * Div[1] + yyij * Div[3];
			//fprintf(stderr, "OK");
		}
	}
}

void divCorrectionNew(struct realArray *div1, struct Nvector u, struct Nvector r, double dp, double reS2, int wijType) {
	for (int i=0; i<nPart; i++) {
		
		double C[4]; // Corrective matrix
		/* | C[0] C[1] | */
		/* | C[2] C[3] | */
		C[0]=C[1]=C[2]=C[3]=0.0;

		double dx, dy, dist, dist2, reS, wij, detC;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);

				C[0] += (dx * dx / dist2) * wij;
				C[1] += (dx * dy / dist2) * wij;
				C[2] += (dy * dx / dist2) * wij;
				C[3] += (dy * dy / dist2) * wij;
			}
		}}

		detC = C[0] * C[3] - C[2] * C[1];
		if(fabs(detC) < eps) fprintf(stderr, "detC[i=%d] = %.4e\n", i, detC);

		double invC[4]; // Inverse of corrective matrix
		if(fabs(detC) >= eps) {
			invC[0] =   C[3] / detC;
			invC[1] = - C[1] / detC;
			invC[2] = - C[2] / detC;
			invC[3] =   C[0] / detC;

			// Divergence correction
			double Crij[2]; // Corrective direction vector
			Crij[0]=Crij[1]=0.0;

			for (int j = 0; j < nPart; j++) {
			if (i != j) {
				dx = r.x[j] - r.x[i];
				dy = r.y[j] - r.y[i];
				dist2 = dx*dx + dy*dy;
				if (dist2 < reS2) {
					dist = sqrt(dist2);
					reS = sqrt(reS2);
					wij = weight(dist, reS, wijType);

					Crij[0] = (invC[0] * dx / dist) + (invC[1] * dy / dist);
					Crij[1] = (invC[2] * dx / dist) + (invC[3] * dy / dist);

					div1->v[i] += ( ((u.x[j] - u.x[i]) / dist) * Crij[0] + ((u.y[j] - u.y[i]) / dist) * Crij[1] ) * wij;
				}
			}}
		}
	}
}

void lapCorrection(struct realArray *lap1, struct realArray lap, struct Nvector grad, struct Nvector r,
	struct realArray f, struct realArray lambda, double pnd0, double lambda0, int dim, double reL2, int wijType) {
	double C1 = 2.0 * dim / (pnd0 * lambda0);
	for (int i = 0; i < nPart; i++) {

		double dx, dy, dist, dist2, reL, wij, Li[2];
		Li[0] = Li[1] = 0;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reL2) {
				dist = sqrt(dist2);
				reL = sqrt(reL2);
				wij = weight(dist, reL, wijType);

				//Li[0] += 2.0 * dim * dx * wij / ( pnd0 * lambda.v[i] );
				//Li[1] += 2.0 * dim * dy * wij / ( pnd0 * lambda.v[i] );
				Li[0] += dx * wij;
				Li[1] += dy * wij;
			}
		}}

		Li[0] *= C1;
		Li[1] *= C1;
		lap1->v[i] = lap.v[i] - (Li[0]*grad.x[i] + Li[1]*grad.y[i]);
	}
}

void lapHigherOrder(struct realArray *lap2, struct Nvector r, struct realArray f, struct realArray lambda, 
	double **Aux, double **Coef, double pnd0, double lambda0, int dim, double dp, double reL2, int wijType) {
	for (int i = 0; i < nPart; i++) {

		//int nSize = 3;
		int nSize = 5;
		for (int ii = 0; ii < nSize; ii++) {
		for (int jj = 0; jj < nSize; jj++) {
			Aux[ii][jj] = 0.0;
			Coef[ii][jj] = 0.0;
		}}

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			double dx = r.x[j] - r.x[i];
			double dy = r.y[j] - r.y[i];
			double dist2 = dx*dx + dy*dy;
			if (dist2 < reL2) {
				double dist = sqrt(dist2);
				double reL = sqrt(reL2);
				//dx /= dist;
				//dy /= dist;
				double wij = weight(dist, reL, wijType);
				//double wij2 = wij * wij;
				//double wij2 = pow(1/dist,3) * pow(1/dist,3);
				//double dx2 = dx*dx;			double dy2 = dy*dy;
				//double dx3 = dx*dx*dx;		double dy3 = dy*dy*dy;
				//double dx4 = dx*dx*dx*dx;	double dy4 = dy*dy*dy*dy;

				//Aux[0][0] += wij*dx4;		Aux[0][1] += wij*dx3*dy;	Aux[0][2] += wij*dx2*dy2;
				//Aux[1][0] += wij*dx3*dy;	Aux[1][1] += wij*dx2*dy2;	Aux[1][2] += wij*dx*dy3;
				//Aux[2][0] += wij*dx2*dy2;	Aux[2][1] += wij*dx*dy3;	Aux[2][2] += wij*dy4;

				double P1 = dx / dist; double P2 = dy / dist; double P3 = dx * dx / (dp * dist);
				double P4 = dy * dy / (dp * dist); double P5 = dx * dy / (dp * dist);

				Aux[0][0] += wij*P1*P1; Aux[0][1] += wij*P1*P2; Aux[0][2] += wij*P1*P3; Aux[0][3] += wij*P1*P4; Aux[0][4] += wij*P1*P5;
				Aux[1][0] += wij*P2*P1; Aux[1][1] += wij*P2*P2; Aux[1][2] += wij*P2*P3; Aux[1][3] += wij*P2*P4; Aux[1][4] += wij*P2*P5;
				Aux[2][0] += wij*P3*P1; Aux[2][1] += wij*P3*P2; Aux[2][2] += wij*P3*P3; Aux[2][3] += wij*P3*P4; Aux[2][4] += wij*P3*P5;
				Aux[3][0] += wij*P4*P1; Aux[3][1] += wij*P4*P2; Aux[3][2] += wij*P4*P3; Aux[3][3] += wij*P4*P4; Aux[3][4] += wij*P4*P5;
				Aux[4][0] += wij*P5*P1; Aux[4][1] += wij*P5*P2; Aux[4][2] += wij*P5*P3; Aux[4][3] += wij*P5*P4; Aux[4][4] += wij*P5*P5;
			}
		}}

		for (int ii = 0; ii < nSize; ii++) {
		for (int jj = 0; jj < nSize; jj++) {
			Aux[ii][jj] /= pnd0;
		}}
		
		// Inverse
		MatrixInversion(Aux, nSize, Coef);

		//for (int ii = 0; ii < 5; ii++)
		//    for (int jj = 0; jj < 5; jj++)
		//    {
		//        printf ("C %d %d = %f ", ii, jj, Coef[ii][jj]);
		//    }
		//    printf ("\n");

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			double dx = (r.x[j] - r.x[i]);
			double dy = (r.y[j] - r.y[i]);
			double dist2 = dx*dx + dy*dy;
			if (dist2 < reL2) {
				double dist = sqrt(dist2);
				double reL = sqrt(reL2);
				//dx /= dist;
				//dy /= dist;
				double wij = weight(dist, reL, wijType);
				double P1 = dx / dist; double P2 = dy / dist; double P3 = dx * dx / (dp * dist);
				double P4 = dy * dy / (dp * dist); double P5 = dx * dy / (dp * dist);

				// Laplacian
				//lap2->v[i] += wij * 2.0 * (f.v[j] - f.v[i]) * ( Coef[0][0] * dx * dx + Coef[0][1] * dx * dy + Coef[0][2] * dy * dy +
				//Coef[2][0] * dx * dx + Coef[2][1] * dx * dy + Coef[2][2] * dy * dy ) / dist2;
				lap2->v[i] += wij * (f.v[j] - f.v[i]) * (	(Coef[2][0] + Coef[3][0]) * P1 + 
															(Coef[2][1] + Coef[3][1]) * P2 + 
															(Coef[2][2] + Coef[3][2]) * P3 +
															(Coef[2][3] + Coef[3][3]) * P4 +
															(Coef[2][4] + Coef[3][4]) * P5 ) / (dp * dist);
			}
		}}

		lap2->v[i] *= 2.0 / pnd0;
	}

}

void lap2ndLDD(struct realArray *lapLDD, struct Nvector r, struct realArray f, struct realArray lambda, 
	double **Aux, double **Coef, double pnd0, double lambda0, int dim, double dp, double reL2, int wijType) {
	for (int i = 0; i < nPart; i++) {

		int nSize = 2;
		for (int ii = 0; ii < nSize; ii++) {
		for (int jj = 0; jj < nSize; jj++) {
			Aux[ii][jj] = 0.0;
			Coef[ii][jj] = 0.0;
		}}
		
		double Oi[2];
		Oi[0] = Oi[1] = 0;
		
		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			double dx = r.x[j] - r.x[i];
			double dy = r.y[j] - r.y[i];
			double dist2 = dx * dx + dy * dy;
			if (dist2 < reL2) {
				double dist = sqrt(dist2);
				double reL = sqrt(reL2);
				double wij = weight(dist, reL, wijType);
				double dx2 = dx * dx;	double dy2 = dy * dy;	double dxy = dx * dy;

				Aux[0][0] += wij * dx2 / (dp * dist);	Aux[0][1] += wij * dxy / (dp * dist);
				Aux[1][0] += wij * dxy / (dp * dist);	Aux[1][1] += wij * dy2 / (dp * dist);

				Oi[0] += wij * dx / (dp * dist);	Oi[1] += wij * dy / (dp * dist);
			}
		}}

		// Inverse
		MatrixInversion(Aux, nSize, Coef);

		double wr2xBO = 0.0;
		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			double dx = (r.x[j] - r.x[i]);
			double dy = (r.y[j] - r.y[i]);
			double dist2 = dx*dx + dy*dy;
			if (dist2 < reL2) {
				double dist = sqrt(dist2);
				double reL = sqrt(reL2);
				double wij = weight(dist, reL, wijType);
				
				// Laplacian
				double xBO = 1 - dx * (Coef[0][0] * Oi[0] + Coef[0][1] * Oi[1]) - dy * (Coef[1][0] * Oi[0] + Coef[1][1] * Oi[1]);
				lapLDD->v[i] += wij * (f.v[j] - f.v[i]) * xBO;
				wr2xBO += wij * dist2 * xBO;
			}
		}}

		lapLDD->v[i] *= 2.0 * dim / wr2xBO;
	}

}

void calcGradLapCPM(struct Nvector *grad, struct realArray *lap, struct Nvector r, struct realArray f, struct realArray lambda, 
	double **Aux, double **Coef, double pnd0, double lambda0, int dim, double dp, double reS2, double reL2, int wijType) {
	for (int i = 0; i < nPart; i++) {

		for (int ii = 0; ii < 5; ii++) {
		for (int jj = 0; jj < 5; jj++) {
			Aux[ii][jj] = 0.0;
			Coef[ii][jj] = 0.0;
		}}

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			double dx = r.x[j] - r.x[i];
			double dy = r.y[j] - r.y[i];
			double dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				double dist = sqrt(dist2);
				double reS = sqrt(reS2);
				dx /= dp; // / dist;
				dy /= dp; // / dist;
				double wij2 = weight(dist, reS, wijType) * weight(dist, reS, wijType);
				double dx2 = dx*dx;			double dy2 = dy*dy;
				double dx3 = dx*dx*dx;		double dy3 = dy*dy*dy;
				double dx4 = dx*dx*dx*dx;	double dy4 = dy*dy*dy*dy;

				Aux[0][0] += wij2*dx2;			Aux[0][1] += wij2*dx*dy;		Aux[0][2] += 0.5*wij2*dx3;		Aux[0][3] += wij2*dx2*dy;		Aux[0][4] += 0.5*wij2*dx*dy2;
				Aux[1][0] += wij2*dx*dy;		Aux[1][1] += wij2*dy2;			Aux[1][2] += 0.5*wij2*dx2*dy;	Aux[1][3] += wij2*dx*dy2;		Aux[1][4] += 0.5*wij2*dy3;
				Aux[2][0] += 0.5*wij2*dx3;		Aux[2][1] += 0.5*wij2*dx2*dy;	Aux[2][2] += 0.25*wij2*dx4;		Aux[2][3] += 0.5*wij2*dx3*dy;	Aux[2][4] += 0.25*wij2*dx2*dy2;
				Aux[3][0] += wij2*dx2*dy;		Aux[3][1] += wij2*dx*dy2;		Aux[3][2] += 0.5*wij2*dx3*dy;	Aux[3][3] += wij2*dx2*dy2;		Aux[3][4] += 0.5*wij2*dx*dy3;
				Aux[4][0] += 0.5*wij2*dx*dy2;	Aux[4][1] += 0.5*wij2*dy3;		Aux[4][2] += 0.25*wij2*dx2*dy2;	Aux[4][3] += 0.5*wij2*dx*dy3;	Aux[4][4] += 0.25*wij2*dy4;
			}
		}}

		// Inverse
		MatrixInversion(Aux, 5, Coef);

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			double dx = r.x[j] - r.x[i];
			double dy = r.y[j] - r.y[i];
			double dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				double dist = sqrt(dist2);
				double reS = sqrt(reS2);
				dx /= dp; // / dist;
				dy /= dp; // / dist;
				double wij2 = weight(dist, reS, wijType) * weight(dist, reS, wijType);
				//double wij2 = pow(1/dist,3) * pow(1/dist,3);
				double dx2 = dx*dx;	double dy2 = dy*dy;	double dxy = dx*dy;

				// Gradient
				grad->x[i] += wij2*(Coef[0][0]*dx + Coef[0][1]*dy + 0.5*Coef[0][2]*dx2 + Coef[0][3]*dxy + 0.5*Coef[0][4]*dy2) * (f.v[j] - f.v[i]) / dp; // / dist;//reS;// / dist;
				grad->y[i] += wij2*(Coef[1][0]*dx + Coef[1][1]*dy + 0.5*Coef[1][2]*dx2 + Coef[1][3]*dxy + 0.5*Coef[1][4]*dy2) * (f.v[j] - f.v[i]) / dp; // / dist;//reS;// / dist;

				// Laplacian
				lap->v[i] += wij2*((Coef[2][0]+Coef[4][0])*dx + (Coef[2][1]+Coef[4][1])*dy + 0.5*(Coef[2][2]+Coef[4][2])*dx2 + 
				(Coef[2][3]+Coef[4][3])*dxy + 0.5*(Coef[2][4]+Coef[4][4])*dy2) * (f.v[j] - f.v[i]) / (dp*dp); // / dist2;
			}
		}}
	}
}

void freeSurafce(struct realArray *bc, struct realArray pnd, double &pnd0, double beta) {
	for (int i = 0; i < nPart; i++) {
		if (pnd.v[i]<beta*pnd0)
			bc->v[i]=1;
		else
		bc->v[i]=0;
	}
}

void calcNormal(struct Nvector *normal, struct realArray bc, struct Nvector r, struct realArray pnd,
	double pnd0, int dim, double reS2, double dPart, int wijType) {
	/*
	// 1st way - neighboor positions
	double neighDist = dPart*1.2;

	for (int i = 0; i < nPart; i++) {
	if (bc.v[i]==1) {
		double dx, dy, dz, dist, orient[3], vec[3];
		orient[0]=0.0; orient[1]=0.0; orient[2]=0.0;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			// Only solid particles or dummies
			//if(prop[mm].type!=solid||i==j||m!=mm)
			//	continue;
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dz = 0;
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				vec[0]=-dx;
				vec[1]=-dy;
				vec[2]=-dz;
				if(fabs(vec[0])<neighDist&&fabs(vec[1])<neighDist&&fabs(vec[2])<neighDist) {
					orient[0]+=vec[0];
					orient[1]+=vec[1];
					orient[2]+=vec[2];
				}
			}
		}}
		/// Normalizing normals
		double normalMod;
		normalMod=sqrt(orient[0]*orient[0]+orient[1]*orient[1]+orient[2]*orient[2]);
		if(normalMod>=eps) {
			normal->x[i]=orient[0]/normalMod;
			normal->y[i]=orient[1]/normalMod;
			normal->z[i]=orient[2]/normalMod;
		}
		else {
			fprintf(stderr, "\nWarning: The norm of a normal vector is 0 \n\n");
			normal->x[i]=0.0;
			normal->y[i]=0.0;
			normal->z[i]=0.0;
		}
		if((isnan(normal->x[i])!=0) || (isnan(normal->y[i])!=0) || (isnan(normal->z[i])!=0)) {
			fprintf(stderr, "\nWarning: normal is NaN 0 ! \n\n");
			normal->x[i]=0.0;
			normal->y[i]=0.0;
			normal->z[i]=0.0;
		}
	}}
	*/
	// 2nd way - gradient of PDN
	for (int i = 0; i < nPart; i++) {
	if (bc.v[i]==1) {
		double dx, dy, dist, dist2, reS, wij, orient[3];
		orient[0]=0.0; orient[1]=0.0; orient[2]=0.0;

		for (int j = 0; j < nPart; j++) {
		if (i != j) {
			dx = r.x[j] - r.x[i];
			dy = r.y[j] - r.y[i];
			dist2 = dx*dx + dy*dy;
			if (dist2 < reS2) {
				dist = sqrt(dist2);
				reS = sqrt(reS2);
				wij = weight(dist, reS, wijType);

				// Subtraction
				orient[0] -= dim * ( pnd.v[j] - pnd.v[i] ) * dx * wij / ( pnd0 * dist2 );
				orient[1] -= dim * ( pnd.v[j] - pnd.v[i] ) * dy * wij / ( pnd0 * dist2 );

				// Addition
				//normal->x[i] += dim * ( pnd.v[j] + pnd.v[i] ) * dx * wij / ( pnd.v[i] * dist2 );
				//normal->y[i] += dim * ( pnd.v[j] + pnd.v[i] ) * dy * wij / ( pnd.v[i] * dist2 );
			}
		}}
		
		/// Normalizing normals
		double normalMod;
		normalMod=sqrt(orient[0]*orient[0]+orient[1]*orient[1]+orient[2]*orient[2]);
		if(normalMod>=eps) {
			normal->x[i]=orient[0]/normalMod;
			normal->y[i]=orient[1]/normalMod;
			normal->z[i]=orient[2]/normalMod;
		}
		else {
			fprintf(stderr, "\nWarning: The norm of a normal vector is 0 \n\n");
			normal->x[i]=0.0;
			normal->y[i]=0.0;
			normal->z[i]=0.0;
		}
		if((isnan(normal->x[i])!=0) || (isnan(normal->y[i])!=0) || (isnan(normal->z[i])!=0)) {
			fprintf(stderr, "\nWarning: normal is NaN 0 ! \n\n");
			normal->x[i]=0.0;
			normal->y[i]=0.0;
			normal->z[i]=0.0;
		}
	}}
}


int main(int argc, char** argv) 
{
	double t0, tf, tempo, dp, rS, rL, pndSZero, pndLZero, lambdaZero;
	int x, y, dim, wijType, order;
	struct timeval tempo_inicio,tempo_fim;

	t0 = tf = tempo = pndSZero = pndLZero = lambdaZero = 0.0;
	dim = 2;
	// irregular arrangements of points (no = 0, yes = 1)
	int irreg = 1;

	printf("Digite o numero de particulas em X e Y (Ex: 20 20): ");
	scanf("%d %d", &x, &y);
	printf("Digite a distancia entre particulas dp (m): ");
	scanf("%lf", &dp);
	printf("Digite o raio small e large (Ex: 2.1 3.1): ");
	scanf("%lf %lf", &rS, &rL);
	printf("Escolha a função peso\n0: re/r-1\n1: re/r+r/re-2\n2: re/r-r/re\n3: (1-r/re)^3\n4: (1-r/re)^2\n");
	scanf("%d", &wijType);
	printf("Escolha a ordem da função exemplo (1, 2, 3, 4): ");
	scanf("%d", &order);

	char fname[256];
	stringstream strs;
	strs << order;
	string temp_str = strs.str();
	char const* pchar = temp_str.c_str(); //dont use cast
	char str1[100] = "f";
	char str2[100] = "_rand.vtk";
	strcat(str1, pchar);
	strcat(str1, str2);
	FILE *fp;
	fp = fopen(str1, "w");

	nPart = (x+1)*(y+1);
	double reS = rS * dp;
	double reL = rL * dp;
	double reS2 = reS * reS;
	double reL2 = reL * reL;
	printf("\n reS: %.3e reS^2: %.3e", reS, reS2);
	printf("\n reL: %.3e reL^2: %.3e", reL, reL2);
	double beta = 0.85;
	eps *= dp;

	struct realArray bc, f, divU, divMPS, divMPS1, divMPS1New, pndSmall, pndLarge, pndNew, d2f;
	struct realArray lambda, lambda2, lambda3, lapMPS, lapMPS1, lapCPM, lap2, lapLDD;
	struct Nvector r, rn, u, un, normal, df, gradMPS, gradMPSadd, gradMPS1, gradCPM;

	makeRealArray (&bc, nPart);
	makeRealArray (&f, nPart);
	makeRealArray (&divU, nPart);
	makeRealArray (&divMPS, nPart);
	makeRealArray (&divMPS1, nPart);
	makeRealArray (&divMPS1New, nPart);
	makeRealArray (&pndSmall, nPart);
	makeRealArray (&pndLarge, nPart);
	makeRealArray (&pndNew, nPart);
	makeRealArray (&lambda, nPart);
	makeRealArray (&lambda2, nPart);
	makeRealArray (&lambda3, nPart);
	makeRealArray (&d2f, nPart);
	makeRealArray (&lapMPS, nPart);
	makeRealArray (&lapMPS1, nPart);
	makeRealArray (&lapCPM, nPart);
	makeRealArray (&lap2, nPart);
	makeRealArray (&lapLDD, nPart);

	makeVector (&r, nPart);
	makeVector (&rn, nPart);
	makeVector (&u, nPart);
	makeVector (&un, nPart);
	makeVector (&df, nPart);
	makeVector (&normal, nPart);
	makeVector (&gradMPS, nPart);
	makeVector (&gradMPSadd, nPart);
	makeVector (&gradMPS1, nPart);
	makeVector (&gradCPM, nPart);

	gettimeofday(&tempo_inicio,NULL);

	// Initial values
	initVals(	&r, &rn, &u, &df, &gradMPS, &gradCPM, &gradMPSadd, &gradMPS1, &normal,
				&f, &divU, &d2f, &pndSmall, &pndLarge, &pndNew, &lambda, &lambda2, &lambda3, 
				&divMPS, &divMPS1, &lapMPS, &lapMPS1, &lapCPM, &lap2, &lapLDD, &bc, &divMPS1New,
				dp, x, y, order, irreg);

	// Particle Number Density - PND
	calcPnd(&pndSmall, &pndLarge, pndSZero, pndLZero, rn, reS2, reL2, wijType);

	// Corrected PND
	//pndCorrection (&pndNew, r, reS2, wijType );

	// Lambda
	// Numerical
	calcLambda(&lambda, &lambda2, &lambda3, lambdaZero, rn, reS2, reL2, wijType);
	printf("\n lambda Numeric: %.4e ", lambdaZero);

	// Analytical
	double ra=1.0/1.7724539;			// pi*ra^2=disa^2 //
	double lambdaA = (	pow(reL2, 2)/12.0-reL * pow(ra, 3)/3.0 +
						pow(ra,4)/4.0)/(reL2/2.0-reL * ra + pow(ra, 2)/2.0);
	printf("\n lambda Analytic: %.4e ", lambdaA);

	// Boundary condition
	freeSurafce(&bc, pndSmall, pndSZero, beta);

	// Normal vector
	calcNormal(&normal, bc, rn, pndSmall, pndSZero, dim, reS2, dp, wijType);

	/// MPS operators
	// Gradient
	calcGradient(&gradMPS, &gradMPSadd, &gradMPS1, rn, f, pndSmall, pndSZero, dim, reS2, wijType);

	// Divergence
	calcDivergence(&divMPS, rn, u, pndSZero, dim, reS2, wijType);

	// Laplacian
	calcLaplacian(&lapMPS, rn, f, lambda, pndLZero, lambdaZero, dim, reL2, wijType);

	// MPS 1st Order
	// 2011 - Enhancement of stability and accuracy of the moving particle semi-implicit method
	// Corrected Gradient - subtraction
	gradCorrection(&gradMPS1, rn, pndSmall, pndSZero, dim, dp, reS2, wijType);
	//gradCorrection(&gradMPS, r, pndSmall, pndSZero, dim, dp, reS2, wijType);

	// Corrected Gradient - addition
	//gradCorrection(&gradPs, r, pndSmall, dim, reS2, wijType);

	// Corrected Divergence
	divCorrection(&divMPS1, u, rn, dp, reS2, wijType);

	// Corrected Divergence New
	// 2018 - An accurate and stable multiphase moving particle semi-implicit 
	// method based on a corrective matrix for all particle interaction models
	divCorrectionNew(&divMPS1New, u, rn, dp, reS2, wijType);

	/// Corrected Laplacian
	// 2018 - An accurate and stable multiphase moving particle semi-implicit 
	// method based on a corrective matrix for all particle interaction models
	lapCorrection(&lapMPS1, lapMPS, gradMPS1, rn, f, lambda, pndLZero, lambdaZero, dim, reL2, wijType);

	/// Laplacian MPS high order
	// 2019 - The truncation and stabilization error in multiphase moving 
	// particle semi-implicit method based on corrective matrix: Which is dominant ?
	//int nSize = 3;
	int nSize = 5;
	double** Aux = new double* [nSize];
	for (int i = 0; i < nSize; i++)
		Aux[i] = new double [nSize];
	double** Coef = new double* [nSize];
	for (int i = 0; i < nSize; i++)
		Coef[i] = new double [nSize];

	lapHigherOrder(&lap2, rn, f, lambda, Aux, Coef, pndLZero, lambdaZero, dim, dp, reL2, wijType);

	/// Laplacian MPS from Lagrangian differencing dynamics (LDD)
	// 2018 - A class of renormalised meshless Laplacians for boundary value problems
	double** Aux2 = new double* [2];
	for (int i = 0; i < 2; i++)
		Aux2[i] = new double [2];
	double** Coef2 = new double* [2];
	for (int i = 0; i < 2; i++)
		Coef2[i] = new double [2];
	
	lap2ndLDD(&lapLDD, rn, f, lambda, Aux2, Coef2, pndLZero, lambdaZero, dim, dp, reL2, wijType);

	//Aux[0][0] = 2; Aux[0][1] = 4; Aux[0][2] = 5; Aux[0][3] = 6; Aux[0][4] = 3;
	//Aux[1][0] = 2; Aux[1][1] = 4; Aux[1][2] = 8; Aux[1][3] = 6; Aux[1][4] = 9;
	//Aux[2][0] = 4; Aux[2][1] = 8; Aux[2][2] = 7; Aux[2][3] = 9; Aux[2][4] = 6;
	//Aux[3][0] = 4; Aux[3][1] = 1; Aux[3][2] = 5; Aux[3][3] = 9; Aux[3][4] = 6;
	//Aux[4][0] = 3; Aux[4][1] = 5; Aux[4][2] = 98; Aux[4][3] = 4; Aux[4][4] = 1;

	//MatrixInversion(Aux, 5, Coef);

	//for (int i = 0; i < 5; i++)
	//{
	//        for (int j = 0; j < 5; j++)
	//                printf ("%f ", Coef[i][j]);
	//        printf ("\n");
	//}

	/// CPM
	// 2012 - A new particle method for simulation of incompressible free surface flow problems
	double** Aux3 = new double* [5];
	for (int i = 0; i < 5; i++)
		Aux3[i] = new double [5];

	double** Coef3 = new double* [5];
	for (int i = 0; i < 5; i++)
		Coef3[i] = new double [5];

	calcGradLapCPM(&gradCPM, &lapCPM, rn, f, lambda, Aux3, Coef3, pndLZero, lambdaZero, dim, dp, reS2, reL2, wijType);

	gettimeofday(&tempo_fim,NULL);
	tf = (double)tempo_fim.tv_usec + ((double)tempo_fim.tv_sec * (1000000.0));
	t0 = (double)tempo_inicio.tv_usec + ((double)tempo_inicio.tv_sec * (1000000.0));
	tempo = (tf - t0) / 1000;

	printf("\nTempo gasto em milissegundos %.3f\n",tempo);

	//std::ostringstream name;
	//name << "case_x2.vtk" ;
	//std::string filename=name.str();
	//std::ofstream file(filename.c_str(), std::fstream::out);

	//printVTK (file, r, f, df, pndSmall, pndNew, gradP, gradCPM, divf, divMPS, div1, d2f, lapMPS, lapCPM, lap2, nPart);

	//char palavra[20];
	//printf("Escreva uma palavra qualquer: ");
	//scanf("%s", palavra);

	printVTK2(	fp, rn, f, u, df, pndSmall, pndNew, bc, normal, gradMPS, gradMPS1, gradCPM, divU, 
				divMPS, divMPS1, divMPS1New, d2f, lapMPS, lapMPS1, lapCPM, lap2, lapLDD, nPart);

	fclose(fp);
	printf("Dados gravados com sucesso!");

	freeRealArray (&bc);
	freeRealArray (&f);
	freeRealArray(&divU);
	freeRealArray(&divMPS);
	freeRealArray(&divMPS1);
	freeRealArray(&divMPS1New);
	freeRealArray(&pndSmall);
	freeRealArray(&pndLarge);
	freeRealArray(&lambda);
	freeRealArray(&lambda3);
	freeRealArray(&d2f);
	freeRealArray(&lapMPS);
	freeRealArray(&lapMPS1);
	freeRealArray(&lapCPM);
	freeRealArray(&lap2);
	freeRealArray(&lapLDD);
	//freeVector (&r0ij);
	freeVector (&r);
	freeVector (&rn);
	freeVector (&u);
	freeVector (&un);
	freeVector (&normal);
	freeVector (&df);
	freeVector (&gradMPS);
	freeVector (&gradMPSadd);
	freeVector (&gradMPS1);
	freeVector (&gradCPM);

	delete[] Aux;
	delete[] Aux3;
	delete[] Coef;
	delete[] Coef3;

	//EsperaEnter();
	return 0;
}