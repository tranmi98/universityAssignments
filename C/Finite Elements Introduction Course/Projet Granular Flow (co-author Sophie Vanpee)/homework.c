// Auteurs : Minh-Phuong Tran et Sophie Vanpee

#include"fem.h"
#include <float.h>
femCouetteProblem * theProblemg;
//===========================================================================
//Grains
//============================================================================
#ifndef NOCONTACTITERATE

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)
{
	int n = myGrains->n;
	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *r = myGrains->r;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double *dvBoundary = myGrains->dvBoundary;
	double *dvContacts = myGrains->dvContacts;
	double rIn = myGrains->radiusIn;
	double rOut = myGrains->radiusOut;
	double zeta = 0.0;

	if (iter != 0) {
		// Collision entre 2 billes ?
		double distx = 0;
		double disty = 0;
		double ndist = 0;
		double gam = 0;
		double normx = 0;
		double normy = 0;
		double vn = 0;

		double dv = 0;
		double DeltaVix = 0;
		double DeltaViy = 0;
		double DeltaVjx = 0;
		double DeltaVjy = 0;
		int p = 0;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {

				distx = x[j] - x[i];
				disty = y[j] - y[i];
				ndist = sqrt(pow(distx, 2) + pow(disty, 2));
				gam = ndist - (r[i] + r[j]);
				normx = distx / ndist;
				normy = disty / ndist;
				vn = (vx[i] * normx + vy[i] * normy) - (vx[j] * normx + vy[j] * normy);

				dv = fmax(0.0, vn + dvContacts[p] - (gam / dt)) - dvContacts[p];
				DeltaVix = -dv * normx*(m[j] / (m[i] + m[j]));
				DeltaViy = -dv * normy*(m[j] / (m[i] + m[j]));
				DeltaVjx = dv * normx*(m[i] / (m[i] + m[j]));
				DeltaVjy = dv * normy*(m[i] / (m[i] + m[j]));
				dvContacts[p] += dv;
				vx[i] += DeltaVix;
				vy[i] += DeltaViy;
				vx[j] += DeltaVjx;
				vy[j] += DeltaVjy;
				zeta = fmax(zeta, fabs(dv));
				p++;

			}
		}

		// Collision avec le mur ?
		double gamOut = 0;
		double gamIn = 0;
		double dvout = 0;
		double dvin = 0;
		double dist = 0;

		for (int i = 0; i < n; i++) {
			dist = sqrt(x[i] * x[i] + y[i] * y[i]);
			gamOut = rOut - dist - r[i];
			gamIn = -rIn + dist - r[i];
			normx = x[i] / dist;
			normy = y[i] / dist;
			vn = (vx[i] * normx + vy[i] * normy);

			dvout = fmax(0.0, vn + dvBoundary[i] - (gamOut / dt));
			dvin = fmax(0.0, -vn - dvBoundary[i] - (gamIn / dt));
			dv = dvout - dvin - dvBoundary[i];
			DeltaVix = -dv * normx;
			DeltaViy = -dv * normy;
			dvBoundary[i] += dv;
			vx[i] += DeltaVix;
			vy[i] += DeltaViy;
			zeta = fmax(zeta, fabs(dv));


		}

	}
	else {

		for (int i = 0; i < n; i++) {
			dvBoundary[i] = 0;
			dvContacts[i] = 0;
		}
		for (int i = n; i < (n*(n - 1) / 2); i++) {
			dvContacts[i] = 0;
		}

	}



	return zeta;

}

#endif
#ifndef NOUPDATE

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femCouetteProblem *theProblem)
{
	int n = myGrains->n;
	int i, iter = 0;
	double zeta;
	double *x = myGrains->x;
	double *y = myGrains->y;
	double *m = myGrains->m;
	double *vy = myGrains->vy;
	double *vx = myGrains->vx;
	double *u = myGrains->u;
	double *v = myGrains->v;
	double gamma = myGrains->gamma;
	double gx = myGrains->gravity[0];
	double gy = myGrains->gravity[1];
	double *X = theProblem->mesh->X;
	double *Y = theProblem->mesh->Y;
	int *elem = theProblem->mesh->elem;
	double phi[3];
	int *triangles = femTriangleBilles(x, y, myGrains->r, n);
	for (int i = 0; i < n; i++) {
		
		u[i] = theProblem->soluce[elem[triangles[i]]]; 
		v[i] = theProblem->soluce2[elem[triangles[i]]];
	}

	free(triangles);
	// 
	// -1- Calcul des nouvelles vitesses des grains sur base de la gravité et de la trainee
	//
	double DeltaVx = 0;
	double DeltaVy = 0;
	for (int i = 0; i < myGrains->n; i++) {
		DeltaVx = gx * dt - (gamma*dt*(vx[i] - u[i]) / m[i]);
		DeltaVy = gy * dt - (gamma*dt*(vy[i] - v[i]) / m[i]);
		vx[i] += DeltaVx;
		vy[i] += DeltaVy;

	}


	// -2- Correction des vitesses pour tenir compte des contacts        
	//       
	do {
		zeta = femGrainsContactIterate(myGrains, dt, iter);
		iter++;
	} while ((zeta > tol / dt && iter < iterMax) || iter == 1);
	printf("iterations = %4d : error = %14.7e \n", iter - 1, zeta);

	//  
	// -3- Calcul des nouvelles positions sans penetrations de points entre eux
	//
	for (i = 0; i < n; ++i) {
		x[i] += vx[i] * dt;
		y[i] += vy[i] * dt;
	}
}


#endif
int funx(const void *a, const void *b) {
	double Xa = theProblemg->mesh->X[*(int *)a];
	double Xb = theProblemg->mesh->X[*(int *)b];
	if (Xa > Xb) {
		return 1;
	}
	else if (Xa < Xb) {
		return -1;
	}
	else {
		return 0;
	}
}
int funy(const void *a, const void *b) {
	double Ya = theProblemg->mesh->Y[*(int *)a];
	double Yb = theProblemg->mesh->Y[*(int *)b];
	if (Ya > Yb) {
		return 1;
	}
	else if (Ya < Yb) {
		return -1;
	}
	else {
		return 0;
	}
}
void femCouetteRenumber(femCouetteProblem *theProblem, femRenumType renumType)
{
	int i;

	int n = (theProblem->mesh->nNode);
	theProblemg = theProblem;
	int *tab;

	switch (renumType) {
	case FEM_NO:
		for (i = 0; i < theProblem->mesh->nNode; i++)
			theProblem->number[i] = i;
		break;
		
	case FEM_XNUM:
		tab = (int *)malloc(n * sizeof(int));
		for (i = 0; i < n; i++) {
			tab[i] = i;
		}

		qsort(tab, n, sizeof(int), funx);

		for (i = 0; i < n; i++) {
			theProblem->number[tab[i]] = i;
		}
		free(tab);
		break;

	case FEM_YNUM:

		tab = (int *)malloc(n * sizeof(int));
		for (i = 0; i < n; i++) {
			tab[i] = i;
		}

		qsort(tab, n, sizeof(int), funy);

		for (i = 0; i < n; i++) {
			theProblem->number[tab[i]] = i;
		}
		free(tab);
		break;


		//
		// end
		//

	default: Error("Unexpected renumbering option");
	}
}

int femCouetteComputeBand(femCouetteProblem *theProblem)
{
	femMesh *theMesh = theProblem->mesh;
	int nloc = theMesh->nLocalNode;
	int m = theMesh->nElem;
	int max1, min1, i, j;
	int myBand = 0;
	for (i = 0; i < m; i++) {
		max1 = theProblem->number[theMesh->elem[i*nloc]];
		min1 = max1;
		for (j = 0; j < nloc; j++) {
			max1 = fmax(theProblem->number[theMesh->elem[i*nloc + j]], max1);
			min1 = fmin(theProblem->number[theMesh->elem[i*nloc + j]], min1);
		}
		myBand = fmax(myBand, (max1 - min1));
	}

	return myBand + 1;
}

