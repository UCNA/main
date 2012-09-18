#include "GraphUtils.hh"
#include "strutils.hh"
#include <cassert>
#include <math.h> 
#include <TROOT.h>
#include <TDirectory.h>
#include <cfloat>

Stringmap histoToStringmap(const TH1* h) {
	assert(h);
	Stringmap m;
	m.insert("nbins",h->GetNbinsX());
	m.insert("name",h->GetName());
	m.insert("title",h->GetTitle());
	std::vector<float> binEdges;
	std::vector<float> binConts;
	std::vector<float> binErrs;
	for(int i=0; i<=h->GetNbinsX()+1; i++) {
		binConts.push_back(h->GetBinContent(i));
		binErrs.push_back(h->GetBinError(i));
		if(i<=h->GetNbinsX())
			binEdges.push_back(h->GetBinLowEdge(i+1));
	}
	m.insert("binEdges",vtos(binEdges));
	m.insert("binErrs",vtos(binErrs));
	m.insert("binConts",vtos(binConts));
	return m;
}

TH1F* stringmapToTH1F(const Stringmap& m) {
	std::string hName = m.getDefault("name","hFoo");
	std::string hTitle = m.getDefault("name","hFoo");
	unsigned int nBins = (unsigned int)(m.getDefault("nbins",0));
	assert(nBins >= 1);
	std::vector<Float_t> binEdges = sToFloats(m.getDefault("binEdges",""));
	std::vector<float> binConts = sToFloats(m.getDefault("binConts",""));
	std::vector<float> binErrs = sToFloats(m.getDefault("binErrs",""));
	assert(binEdges.size()==nBins+1);
	assert(binConts.size()==nBins+2);
	assert(binErrs.size()==nBins+2);
	
	// TODO does this work right?
	TH1F* h = new TH1F(hName.c_str(),hTitle.c_str(),nBins,&binEdges[0]);
	for(unsigned int i=0; i<=nBins+1; i++) {
		h->SetBinContent(i,binConts[i]);
		h->SetBinError(i,binErrs[i]);
	}
	return h;
}

Stringmap graphToStringmap(const TGraph& g) {
	Stringmap m;
	m.insert("npts",g.GetN());
	std::vector<float> xs;
	std::vector<float> ys;
	double x,y;
	for(int i=0; i<g.GetN(); i++) {
		g.GetPoint(i,x,y);
		xs.push_back(x);
		ys.push_back(y);
	}
	m.insert("x",vtos(xs));
	m.insert("y",vtos(ys));
	return m;	
}

TGraphErrors* TH1toTGraph(const TH1& h) {
	TGraphErrors* g = new TGraphErrors(h.GetNbinsX()-2);
	for(int i=0; i<h.GetNbinsX()-2; i++) {
		g->SetPoint(i,h.GetBinCenter(i+1),h.GetBinContent(i+1));
		g->SetPointError(i,0.0,h.GetBinError(i+1));
	}
	return g;
}


TH1F* cumulativeHist(const TH1F& h, bool normalize) {
	TH1F* c = new TH1F(h);
	int n = h.GetNbinsX()-2;
	float ecum2 = 0;
	c->SetBinContent(0,0);
	c->SetBinError(0,0);
	for(int i=1; i<=n+1; i++) {
		c->SetBinContent(i,c->GetBinContent(i-1)+h.GetBinContent(i));
		ecum2 += h.GetBinError(i);
		c->SetBinError(i,sqrt(ecum2));
	}
	if(normalize)
		c->Scale(1.0/c->GetBinContent(n));
	return c;
}

TGraph* invertGraph(const TGraph* g) {
	assert(g);
	TGraph* gi = new TGraph(g->GetN());
	double x,y;
	for(int i=0; i<g->GetN(); i++) {
		g->GetPoint(i,x,y);
		gi->SetPoint(i,y,x);
	}
	return gi;
}

TGraph* combine_graphs(const std::vector<TGraph*> gs) {
	unsigned int npts = 0;
	for(unsigned int n=0; n<gs.size(); n++)
		npts += gs[n]->GetN();
	TGraph* g = new TGraph(npts);
	npts = 0;
	double x,y;
	for(unsigned int n=0; n<gs.size(); n++) {
		for(int n2 = 0; n2 < gs[n]->GetN(); n2++) {
			gs[n]->GetPoint(n2,x,y);
			g->SetPoint(npts++,x,y);
		}
	}
	return g;
}

TGraphErrors* merge_plots(const std::vector<TGraphErrors*>& pin, const std::vector<int>& toffset) {
	printf("Merging %i graphs...\n",(int)pin.size());
	unsigned int npts = 0;
	for(unsigned int n=0; n<pin.size(); n++)
		npts += pin[n]->GetN();
	TGraphErrors* tg = new TGraphErrors(npts);
	npts = 0;
	double x,y;
	for(unsigned int n=0; n<pin.size(); n++) {
		for(int n2 = 0; n2 < pin[n]->GetN(); n2++) {
			pin[n]->GetPoint(n2,x,y);
			tg->SetPoint(npts,(x+toffset[n])/3600.0,y);
			tg->SetPointError(npts,pin[n]->GetErrorX(n2)/3600.0,pin[n]->GetErrorY(n2));
			++npts;
		}
	}
	tg->GetXaxis()->SetTitle("Time [Hours]");
	return tg;
}

void drawTogether(std::vector<TGraphErrors*>& gs, float ymin, float ymax, TCanvas* C, const char* outname, const char* graphTitle) {
	if(!gs.size())
		return;
	for(unsigned int t=0; t<gs.size(); t++)
		gs[t]->SetLineColor(t+1);
	gs[0]->SetMinimum(ymin);
	gs[0]->SetMaximum(ymax);
	gs[0]->SetTitle(graphTitle);
	gs[0]->Draw("AP");
	for(unsigned int i=1; i<gs.size(); i++)
		gs[i]->Draw("P");
	C->Print(outname);
	
}

TGraph* matchHistoShapes(const TH1F& h1, const TH1F& h2) {
	TH1F* c1 = cumulativeHist(h1,true);
	TH1F* c2 = cumulativeHist(h2,true);
	TGraph* c2g = TH1toTGraph(*c2);
	delete(c2);
	TGraph* c2i = invertGraph(c2g);
	delete(c2g);
	int n = h1.GetNbinsX()-2;
	TGraph* T = new TGraph(n);
	for(int i=1; i<=n; i++)
		T->SetPoint(i-1,c1->GetBinCenter(i),c2i->Eval(c1->GetBinContent(i)));
	delete(c1);
	delete(c2i);
	return T;
}

void scale(TGraphErrors& tg, float s) {
	double x,y;
	for(int i=0; i<tg.GetN(); i++) {
		tg.GetPoint(i,x,y);
		tg.SetPoint(i,x,s*y);
		tg.SetPointError(i,tg.GetErrorX(i),s*tg.GetErrorY(i));
	}
}

TGraph* derivative(TGraph& g) {
	g.Sort();
	TGraph* d = new TGraph(g.GetN()-1);
	double x1,y1,x2,y2;
	g.GetPoint(0,x1,y1);
	for(int i=0; i<g.GetN()-1; i++) {
		g.GetPoint(i+1,x2,y2);
		d->SetPoint(i,0.5*(x1+x2),(y2-y1)/(x2-x1));
		x1 = x2;
		y1 = y2;
	}
	return d;
}

void transformAxis(TGraph& g, TGraph& T, bool useJacobean) {
	double x,y,j;
	j = 1.0;
	TGraph* d = NULL;
	if(useJacobean)
		d = derivative(T);
	for(int i=0; i<g.GetN(); i++) {
		g.GetPoint(i,x,y);
		if(useJacobean)
			j=d->Eval(x);
		g.SetPoint(i,T.Eval(x),j*y);
	}
	if(useJacobean)
		delete(d);
}

TGraphErrors* interpolate(TGraphErrors& tg, float dx) {
	std::vector<float> xnew;
	std::vector<float> ynew;
	std::vector<float> dynew;
	double x0,x1,y,dy0,dy1;
	
	// sort input points by x value
	tg.Sort();
	
	// interpolate each interval of the original graph
	for(int i=0; i<tg.GetN()-1; i++) {
		tg.GetPoint(i,x0,y);
		tg.GetPoint(i+1,x1,y);
		dy0 = tg.GetErrorY(i);
		dy1 = tg.GetErrorY(i+1);
		// determine number of points for this interval
		int ninterp = (x1-x0>dx)?int((x1-x0)/dx):1;
		for(int n=0; n<ninterp; n++) {
			float l = float(n)/float(ninterp);
			xnew.push_back(x0+(x1-x0)*l);
			ynew.push_back(tg.Eval(xnew.back()));
			dynew.push_back(sqrt(ninterp)*((1-l)*dy0+l*dy1));
		}
	}
	
	// fill interpolated output graph
	TGraphErrors* gout = new TGraphErrors(xnew.size());
	for(unsigned int i=0; i<xnew.size(); i++) {
		gout->SetPoint(i,xnew[i],ynew[i]);
		gout->SetPointError(i,0,dynew[i]);
	}
	return gout;
}

double invCDF(TH1* h, double p) {
	assert(h);
	unsigned int nbins = h->GetNbinsX()-2;
	if(p<=0) return 0;
	if(p>=1) return h->GetBinLowEdge(nbins+1);
	Double_t* cdf = h->GetIntegral();
	unsigned int mybin = std::upper_bound(cdf,cdf+nbins,p)-cdf;
	assert(mybin>0);
	assert(mybin<=nbins);
	double l = (p-cdf[mybin-1])/(cdf[mybin]-cdf[mybin-1]);
	return h->GetBinLowEdge(mybin)*(1.0-l)+h->GetBinLowEdge(mybin+1)*l;
}

void fixNaNs(TH1* h) {
	unsigned int nb = h->GetNbinsX();
	for(unsigned int i=0; i<=nb+1; i++) {
		if(!(h->GetBinContent(i)==h->GetBinContent(i))) {
			printf("NaN found in bin %i/%i!\n",i,nb);
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}
}

std::vector<TH2F*> sliceTH3(const TH3& h3, AxisDirection d) {
	assert(d<=Z_DIRECTION);
	
	TAxis* Ax1 = d==X_DIRECTION?h3.GetYaxis():h3.GetXaxis();
	TAxis* Ax2 = d==Z_DIRECTION?h3.GetYaxis():h3.GetZaxis();
	TAxis* Ax3 = d==X_DIRECTION?h3.GetXaxis():d==Y_DIRECTION?h3.GetYaxis():h3.GetZaxis();
	const unsigned int n1 = Ax1->GetNbins();
	const unsigned int n2 = Ax2->GetNbins();
	const unsigned int n3 = Ax3->GetNbins();
	
	std::vector<TH2F*> h2s;
	for(unsigned int z = 0; z <= n3+1; z++) {
		TH2F* h2 = new TH2F((std::string(h3.GetName())+"_"+itos(z)).c_str(),
							h3.GetTitle(),
							n1,
							Ax1->GetBinLowEdge(1),
							Ax1->GetBinUpEdge(n1),
							n2,
							Ax2->GetBinLowEdge(1),
							Ax2->GetBinUpEdge(n2));
		if(h3.GetSumw2()) h2->Sumw2();
		h2->GetXaxis()->SetTitle(Ax1->GetTitle());
		h2->GetYaxis()->SetTitle(Ax2->GetTitle());
		for(unsigned int x=0; x <= n1+1; x++) {
			for(unsigned int y=0; y <= n2+1; y++) {
				if(d==X_DIRECTION) {
					h2->SetBinContent(x,y,h3.GetBinContent(z,x,y));
					h2->SetBinError(x,y,h3.GetBinError(z,x,y));
				} else if(d==Y_DIRECTION) {
					h2->SetBinContent(x,y,h3.GetBinContent(x,z,y));
					h2->SetBinError(x,y,h3.GetBinError(x,z,y));
				} else if(d==Z_DIRECTION) {
					h2->SetBinContent(x,y,h3.GetBinContent(x,y,z));
					h2->SetBinError(x,y,h3.GetBinError(x,y,z));
				}
			}
		}
		h2s.push_back(h2);
	}
	return h2s;
}

std::vector<TH1F*> sliceTH2(const TH2& h2, AxisDirection d, bool includeOverflow) {
	assert(d==X_DIRECTION || d==Y_DIRECTION);
	std::vector<TH1F*> h1s;
	const unsigned int nx = h2.GetNbinsX();
	const unsigned int ny = h2.GetNbinsY();
	const TAxis* axs = d==X_DIRECTION?h2.GetYaxis():h2.GetXaxis();
	const unsigned int nz = d==X_DIRECTION?nx:ny;
	const unsigned int nn = d==X_DIRECTION?ny:nx;
	
	for(unsigned int z = 0; z <= nz+1; z++) {
		if(!includeOverflow && (z==0 || z==nz+1)) continue;
		TH1F* h1 = new TH1F((std::string(h2.GetName())+"_"+itos(z)).c_str(),
							h2.GetTitle(),
							nn,
							axs->GetBinLowEdge(1),
							axs->GetBinUpEdge(nn));
		if(h2.GetSumw2()) h1->Sumw2();
		h1->GetXaxis()->SetTitle(axs->GetTitle());
		for(unsigned int n=0; n <= nn+1; n++) {
			if(d==X_DIRECTION) {
				h1->SetBinContent(n,h2.GetBinContent(z,n));
				h1->SetBinError(n,h2.GetBinError(z,n));
			} else { 
				h1->SetBinContent(n,h2.GetBinContent(n,z));
				h1->SetBinError(n,h2.GetBinError(n,z));
			}
		}
		h1s.push_back(h1);
	}
	return h1s;	
}

std::vector<unsigned int> equipartition(const std::vector<float>& elems, unsigned int n) {
	std::vector<float> cumlist;
	for(unsigned int i=0; i<elems.size(); i++)
		cumlist.push_back(i?cumlist[i-1]+elems[i]:elems[i]);
	
	std::vector<unsigned int> part;
	part.push_back(0);
	for(unsigned int i=1; i<n; i++) {
		double x0 = cumlist.back()*float(i)/float(n);
		unsigned int p = (unsigned int)(std::upper_bound(cumlist.begin(),cumlist.end(),x0)-cumlist.begin()-1);
		if(p != part.back()) part.push_back(p);
	}
	part.push_back(elems.size());
	return part;
}

TH1* projectTH2(const TH2& h, double nb, double cx, double cy) {
	TAxis* Ax = h.GetXaxis();
	TAxis* Ay = h.GetYaxis();
	double x0 = Ax->GetXmin();
	double x1 = Ax->GetXmax();
	double y0 = Ay->GetXmin();
	double y1 = Ay->GetXmax();
	TH1D* hOut = new TH1D((h.GetName()+std::string("_Projected")).c_str(),"Projected Histogram",nb,
						  cx*(cx>0?x0:x1)+cy*(cy>0?y0:y1),cx*(cx>0?x1:x0)+cy*(cy>0?y1:y0));
	for(int bx=1; bx<=Ax->GetNbins(); bx++)
		for(int by=1; by<=Ay->GetNbins(); by++)
			hOut->Fill(cx*Ax->GetBinCenter(bx)+cy*Ay->GetBinCenter(by),h.GetBinContent(bx,by));
	return hOut;
}

TH1* histsep(const TH1& h1, const TH1& h2) {
	int nb = h1.GetNbinsX();
	assert(nb==h2.GetNbinsX());
	TH1* hDiv = (TH1*)h1.Clone("hDivision");
	hDiv->SetBinContent(0,0);
	hDiv->SetBinContent(nb+1,0);
	for(int b=1; b<=nb; b++)
		hDiv->SetBinContent(b,hDiv->GetBinContent(b-1)+h2.GetBinContent(b));
	double c = 0;
	for(int b=nb; b>=1; b--) {
		c += (b==nb?0:h1.GetBinContent(b+1));
		hDiv->SetBinContent(b,hDiv->GetBinContent(b)+c);
		hDiv->SetBinError(b,0);
	}
	return hDiv;
}

void histoverlap(const TH1& h1, const TH1& h2, double& o, double& xdiv) {
	int nb = h1.GetNbinsX();
	assert(nb==h2.GetNbinsX());
	double* csum = new double[nb+2];
	csum[0] = csum[nb+1] = 0;
	for(int b=1; b<=nb; b++)
		csum[b] = csum[b-1]+h2.GetBinContent(b);
	double c = 0;
	int bmn = nb;
	for(int b=nb; b>=1; b--) {
		c += (b==nb?0:h1.GetBinContent(b+1));
		csum[b] += c;
		if(csum[b] <= csum[bmn]) bmn=b;
	}
	o = csum[bmn];
	xdiv = h1.GetBinLowEdge(bmn+1);
	delete[] csum;
}

