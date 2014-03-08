#include <Math/QuasiRandom.h>
#include <Math/Random.h>
#include <TMatrixD.h>
#include <stdio.h>
#include <cmath>
#include <set>
#include <cassert>

#include "OutputManager.hh"
#include "PathUtils.hh"
#include "strutils.hh"
#include <TGraph.h>

using namespace ROOT::Math;
using namespace std;

#define N_DIM 6
#define N_COEFFS ((N_DIM+1)*(N_DIM+2))/2
#define N_SEARCHPTS 200000
#define N_SEARCHBALL 1000

double searchpts[N_SEARCHPTS][N_DIM];	//< potential search points
double best_pts[N_COEFFS][N_DIM];		//< best points ever found
double ptvals[N_SEARCHPTS][N_COEFFS];	//< coefficient values at each point
QuasiRandomNiederreiter rNied(N_DIM);	//< random generator
unsigned int pt_idx[N_COEFFS];			//< selected indices of points to use
set<unsigned int> pt_idx_set;			//< set of used point indices, to avoid duplication
double best_val = 0;					//< optimum value found so far
double best_ever = 0;					//< best value ever found

TMatrixD M(N_COEFFS,N_COEFFS);			//< matrix of fit vectors
TMatrixD* Mm[N_COEFFS];					//< minors for fast row-swapping determinant calculation
double DetMm[N_COEFFS];					//< determinants of minors

double rsearch = pow(pow(2.,N_DIM)/N_SEARCHPTS,1./N_DIM)/2.;

// random search point in sphere
void sphere_point(double* v) {
	while(true) {
		rNied.Next(v);
		double s = 0;
		for(unsigned int n=0; n<N_DIM; n++) {
			v[n] = (v[n]-0.5)*2;
			s += v[n]*v[n];
		}
		if(s<1) break;
	}
}

/// initialize matrix minors array
void init_Mm() {
	for(unsigned int n=0; n<N_COEFFS; n++)
		Mm[n] = new TMatrixD(N_COEFFS-1,N_COEFFS-1);
}

/// calculate determinants of minors
void calcMinorDet(unsigned int n_coeff) {
	for(unsigned int n=0; n<N_COEFFS; n++) {
		unsigned int ii=0;
		for(unsigned int i=0; i<N_COEFFS; i++) {
			if(i==n_coeff) continue;
			unsigned int jj=0;
			for(unsigned int j=0; j<N_COEFFS; j++) {
				if(j==n) continue;
				(*Mm[n])[ii][jj++] = M[i][j];
			}
			ii++;
		}
		DetMm[n] = Mm[n]->Determinant();
	}
}

/// select new vector to try; return new Det^2
inline double set_vector(unsigned int n_coeff, unsigned int n_pt) {
	pt_idx_set.erase(pt_idx[n_coeff]);
	pt_idx_set.insert(n_pt);
	pt_idx[n_coeff] = n_pt;
	double D = 0;
	double s = 1;
	for(unsigned int n=0; n<N_COEFFS; n++) {
		M[n_coeff][n] = ptvals[n_pt][n];
		D += s*M[n_coeff][n]*DetMm[n];
		s *= -1;
	}
	return D*D;
}

void new_searchpts(bool use_prev) {
	assert(N_COEFFS + N_SEARCHBALL*N_COEFFS*2 < N_SEARCHPTS);
	
	printf("Selecting new search points...\n");
	
	best_val = 0;
	unsigned int i=0;
	
	// search around previous best points
	double searchball[N_DIM];
	if(use_prev) {
		for(unsigned int bp=0; bp<N_COEFFS; bp++) {
			for(unsigned int n=0; n<N_DIM; n++)
					searchpts[i][n] = best_pts[bp][n];
			i++;
		}
		for(unsigned int bp=0; bp<N_COEFFS; bp++) {
			for(unsigned int sb=0; sb<N_SEARCHBALL; sb++) {
				double s = 0;
				sphere_point(searchball);
				for(unsigned int n=0; n<N_DIM; n++) {
					searchpts[i][n] = best_pts[bp][n]+searchball[n];
					s += searchpts[i][n]*searchpts[i][n];
				}
				if(s<1) i++;
			}
		}
	}
	
	// add new points
	while(i<N_SEARCHPTS)
		sphere_point(searchpts[i++]);
	
	printf("Coefficient values at search points...\n");
	// coefficient values at search point
	for(i=0; i<N_SEARCHPTS; i++) {
		unsigned int m=0;
		ptvals[i][m++] = 1.;
		for(unsigned int j=0; j<N_DIM; j++) {
			ptvals[i][m++] = searchpts[i][j];
			for(unsigned int k=j; k<N_DIM; k++)
				ptvals[i][m++] = searchpts[i][j]*searchpts[i][k];
		}
	}
	
	for(unsigned int n=0; n<N_COEFFS; n++)
		set_vector(n,n);
}

/// scan through all possible points for given coefficient vector
bool search_coeff(unsigned int n_coeff) {
	unsigned int c0 = pt_idx[n_coeff];
	unsigned int c1 = c0;
	calcMinorDet(n_coeff);
	for(unsigned int i=0; i<N_SEARCHPTS; i++) {
		if(pt_idx_set.count(i)) continue;
		double v = set_vector(n_coeff,i);
		if(v>best_val) {
			best_val = v;
			c1 = i;
		}
	}
	set_vector(n_coeff,c1);
	//printf("%i %i->%i %.5f\n",n_coeff,c0,c1,best_val);
	return c1 != c0;
}


unsigned int scan_pass() {
	unsigned int nchanged = 0;
	for(unsigned int n=0; n<N_COEFFS; n++)
		nchanged += search_coeff(n);
	return nchanged;
}

void display_best() {
	printf("Current best = %.5f\n",best_ever);
	for(unsigned int i=0; i<N_COEFFS; i++) {
		for(unsigned int j=0; j<N_DIM; j++)
			printf("\t%.3f",best_pts[i][j]);
		printf("\n");
	}
}

void full_scan() {
	printf("Performing scan...\n");
	unsigned int nch;
	while( (nch=scan_pass()) ) { }
	if(best_val>best_ever) {
		best_ever = best_val;
		for(unsigned int i=0; i<N_COEFFS; i++)
			for(unsigned int j=0; j<N_DIM; j++)
				best_pts[i][j] = searchpts[pt_idx[i]][j];
		display_best();
	} else {
		printf("Scan best %f\n",best_val);
	}
}







int main(int argc, char *argv[]) {
	
	init_Mm();
	new_searchpts(false);
	
	unsigned int n_iterations = 20;
	for(unsigned int i=0; i<n_iterations; i++) {
		full_scan();
		new_searchpts(true);
	}
	
	display_best();
	
	OutputManager OM("OptQuadFitPts_"+itos(N_DIM),getEnvSafe("UCNA_ANA_PLOTS")+"/test/OptQuadFitPts");
	for(unsigned int a1 = 0; a1<N_DIM; a1++) {
		for(unsigned int a2 = a1+1; a2<N_DIM; a2++) {
			TGraph g(N_COEFFS);
			for(unsigned int i=0; i<N_COEFFS; i++) g.SetPoint(i,best_pts[i][a1],best_pts[i][a2]);
			g.SetMinimum(-1);
			g.SetMaximum(1);
			g.Draw("AP");
			g.GetXaxis()->SetLimits(-1,1);
			g.SetMarkerStyle(33);
			g.SetMarkerSize(0.75);
			g.Draw("AP");
			OM.printCanvas("fit_pts_"+itos(N_DIM)+"_"+itos(a1)+itos(a2));
		}
	}
	
	for(unsigned int i=0; i<N_COEFFS; i++) {
		Stringmap m;
		m.insert("fit_pt",vtos(&best_pts[i][0],&best_pts[i][N_DIM]));
		m.insert("n",itos(i));
		OM.qOut.insert("fit_pt",m);
	}
	
	Stringmap m;
	m.insert("n_dim",N_DIM);
	m.insert("n_coeffs",N_COEFFS);
	m.insert("best_val",best_ever);
	m.insert("n_searchpts",N_SEARCHPTS);
	m.insert("n_searchball",N_SEARCHBALL);
	m.insert("rsearch",rsearch);
	m.insert("n_iterations",n_iterations);
	OM.qOut.insert("opt_search",m);
	
	OM.write();
	
	return 0;
}
