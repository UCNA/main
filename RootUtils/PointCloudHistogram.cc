#include "PointCloudHistogram.hh"
#include <cassert>
#include <Math/QuasiRandom.h>
#include <Math/Random.h>

//-------------------------

void KDTreeSet::finalize() {
	printf("Building %lu-point, %u-dimensional kd-tree...\n", fData[0].size(), ndim);
	T = new TKDTree<int,float>(fData[0].size(), ndim, 1);
	for(unsigned int j=0; j<ndim; j++) T->SetData(j, &fData[j][0]);
	T->Build();
	T->SetOwner(kTRUE);
}

void KDTreeSet::fillPointRange(unsigned int npts, const float* xlo, const float* xhi, const float* dens) {
	
	assert(!T);

	ROOT::Math::QuasiRandomNiederreiter r(ndim);
	ROOT::Math::RandomMT r0;
	std::vector<double> x(ndim);
	
	// generate bin points
	printf("Generating %u point %u-dimensional cloud...\n", npts, ndim);
	unsigned int i = 0;
	while(i<npts) {
		if(false)
			r.Next(&x[0]);
		else
			r0.RndmArray(ndim,&x[0]);
		if(dens) {
			unsigned int j = 0;
			for(; j<ndim; j++) if(x[j]>dens[j]) break;
			if(j==ndim) continue;
		}
		for(unsigned int j=0; j<ndim; j++) fData[j].push_back(xlo[j] + (xhi[j]-xlo[j])*x[j]);
		i++;
	}
}

//---------------------------

void PointCloudHistogram::Fill(const float* x, float v) {
	int idx;
	float dist;
	myTree->T->FindNearestNeighbors(x, 1, &idx, &dist);
	std::map<int, float>::iterator it = bins.find(idx);
	if(it != bins.end()) it->second += v;
	else bins.insert(std::pair<int,float>(idx,v));
}

void PointCloudHistogram::project(const float* v, TGraph& g) const {
	unsigned int i=0;
	for(std::map<int,float>::const_iterator it = bins.begin(); it != bins.end(); it++) {
		double s = 0;
		for(unsigned int j=0; j<myTree->ndim; j++) s += myTree->fData[j][it->first] * v[j];
		g.SetPoint(i++, s, it->second);
	}
	g.Sort();
}

void PointCloudHistogram::project(const float* v, TH1& h) const {
	for(std::map<int,float>::const_iterator it = bins.begin(); it != bins.end(); it++) {
		double s = 0;
		for(unsigned int j=0; j<myTree->ndim; j++) s += myTree->fData[j][it->first] * v[j];
		h.Fill(s,it->second);
	}
}

void PointCloudHistogram::getPoint(unsigned int i, float* x) const {
	for(unsigned int j=0; j<myTree->ndim; j++) x[j] = myTree->fData[j][i];
}
