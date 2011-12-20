#include "KurieFitter.hh"

float_err kuriePlotter(TH1* spectrum, float endpoint, TGraphErrors** tgout, float targetEP, float fitStart, float fitEnd) {
	
	assert(spectrum);
	printf("Endpoint %.3f (%i points)\n",endpoint,spectrum->GetNbinsX());
	fflush(stdout);
	
	// set up parameters
	float x,y,y0,dy;
	TGraphErrors* tg_uninterp = new TGraphErrors(spectrum->GetNbinsX());
	
	// convert histogram to estimated true energy scale; get fit parameters estimate; make Kurie plot
	float tgmaxy = 0;
	float tgmaxx = 0;
	for(int i=0; i<spectrum->GetNbinsX(); i++) {
		
		x = spectrum->GetBinCenter(i+1)*targetEP/endpoint;
		y0 = spectrum->GetBinContent(i+1);
		
		if(x<0 || y0 < 0) {
			tg_uninterp->SetPoint(i,0,0);
			tg_uninterp->SetPointError(i,0,0);
			continue;
		}
		
		y = sqrt(y0/(sqrt(x*x+2*m_e*x)*(x+m_e)));
		dy = y/(2*y0)*spectrum->GetBinError(i+1);
		
		tg_uninterp->SetPoint(i,x,y);
		tg_uninterp->SetPointError(i,0,dy);
		if(fitStart <= x && x <= fitEnd && y>tgmaxy) {
			tgmaxy = y;
			tgmaxx = x;
		}
	}
	float slope = tgmaxy/(tgmaxx-endpoint);
	TGraphErrors* tg = interpolate(*tg_uninterp,2.5);
	delete(tg_uninterp);
	if(!slope || !(slope==slope)) {
		printf("\n**** Kurie Endpoint slope estimation failed (%g)!\n",slope);
		if(tgout)
			*tgout = tg;
		else
			delete(tg);
		return float_err(0,0);
	}
	
	// fit for endpoint
	TF1 lf("linfit","[0]*(x-[1])",fitStart,fitEnd);
	lf.SetParameter(1,endpoint);
	lf.SetParameter(0,tgmaxy/(tgmaxx-endpoint));
	int fiterr = tg->Fit(&lf,"QR");
	if(fiterr && fiterr != 4) {
		lf.SetParameter(1,endpoint);
		lf.SetParameter(0,slope);
		printf("\n**** Kurie Endpoint fit failed (err=%i)!\n",fiterr);
		tg->Fit(&lf,"VR");
		if(tgout)
			*tgout = tg;
		else
			delete(tg);
		return float_err(0,0);
	}
	
	assert(lf.GetParameter(1));
	
	if(tgout) {
		// normalize to endpoint
		scale(*tg,1.0/lf.Eval(500.0));
		lf.SetParameter(0,lf.GetParameter(0)/lf.Eval(500.0));
		fiterr = tg->Fit(&lf,"QR");
		if(fiterr && fiterr != 4) {
			printf("\n**** Kurie Endpoint re-fit failed (err=%i)!\n",fiterr);
			return float_err(0,0);
		}		
		tg->SetTitle("Kurie Plot");
		tg->SetMarkerColor(4);
		tg->SetMarkerStyle(21);
		tg->SetMinimum(0.0);
		tg->SetMaximum(4.0);
		*tgout = tg;
	} else
		delete(tg);
	
	assert(lf.GetParameter(1));
	
	return float_err(lf.GetParameter(1)*endpoint/targetEP,lf.GetParError(1)*endpoint/targetEP);
	
	/*
	float a,b,da,db;
	a = lf.GetParameter(0);
	b = lf.GetParameter(1);
	da = lf.GetParError(0);
	db = lf.GetParError(1);
	
	float eperr = sqrt(b*b*da*da+a*a*db*db)/(b*b)*(endpoint/neutronBetaEp);
	endpoint = -a/b*(endpoint/neutronBetaEp);
	printf(" => %.3f~%.3f\n",endpoint,eperr);
	 
	 return float_err(endpoint,eperr);
	*/
	
	
}

float_err kurieIterator(TH1* spectrum, float iguess, TGraphErrors** tgout, float targetEP,
						float fitStart, float fitEnd, unsigned int nmax, float etol) {
	
	float_err oldGuess = 0;
	float_err newGuess = iguess;
	unsigned int ntries = 0;
	
	while(fabs(newGuess.x-oldGuess.x) > etol && ntries < nmax) {
		oldGuess = newGuess;
		newGuess = kuriePlotter(spectrum,newGuess.x,NULL,targetEP,fitStart,fitEnd);
		if(newGuess.x < iguess/2.0 || newGuess.x > 2.0*iguess || !(newGuess.x == newGuess.x)) {
			if(newGuess.x == newGuess.x)
				printf("*** ERROR *** Crazy new guess %f far from old %f\n",newGuess.x,iguess);
			else
				printf("*** ERROR *** Kurie fit returned NaN!\n");
			newGuess = 0;
			break;
		}
		newGuess = 0.5*(oldGuess+newGuess);	// step half-way for smoother convergence
		++ntries;
	}
	
	if(ntries == nmax)
		printf("\n********** Warning: FAILED CONVERGENCE **************\n\n");
	if(newGuess.x) {
		float_err finalAnswer = 0.5*(newGuess+kuriePlotter(spectrum,newGuess.x,tgout,targetEP,fitStart,fitEnd));
		printf("----------------> %.3f +/- %.3f\n",finalAnswer.x,finalAnswer.err);
		return finalAnswer;
	}
	
	return float_err(0,0);		
}

