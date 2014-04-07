//2-19-02 J. Yuan: magnet field for UCNA experiment
////////////////////////////////////////////////////

#include "bmField.hh"
#include <cmath>
#include <fstream>

extern "C" {
	      void b_field_(double *bx, double *by, double *bz, 
			     const double *x, const double *y, const double *z);
	    }

bmField::bmField(const TString filename): A(0.796), B(0.204), C(192.0), D(0.08491), addAFP(false), rmax2(20*20*cm2), fieldScale(1.0), RobbyField(true)  {

  if(!RobbyField) {LoadFieldMap(filename);}
}

void bmField::LoadFieldMap(const TString filename) {
	
	Bpoints.clear();
	Zpoints.clear();
	
	if(filename==""){
		//default field profile
		addPoint(-3.0*m,0.6*tesla);
		addPoint(-2.2*m,0.6*tesla);
		addPoint(-1.5*m,1.0*tesla);
		addPoint(1.5*m,1.0*tesla);
		addPoint(2.2*m,0.6*tesla);
		addPoint(3.0*m,0.6*tesla);
	} else {
		// load profile from file
		ifstream fin;
		fin.open(filename.Data());
		if(!fin) {
			G4cout<<"Can not open "<<filename.Data()<<G4endl;
			exit(1);
		} 
		TString stmp;
		stmp.ReadLine(fin);  //skip the first title line
		
		TString stmp1, stmp2;
		while(fin){
			fin>>stmp1>>stmp2;
			if(stmp1.IsNull()) break;
			else {
				G4double z = atof(stmp1.Data())*m;
				G4double B = atof(stmp2.Data())*tesla;
				addPoint(z,B);
				cout<<z/m<<" "<<B/tesla<<endl;
			}
		}
		fin.close();
	}
}

void addAFPFringeField(const G4double Point[3], G4double *Bfield) {
	double z0 = 100.-Point[0]/cm;								// z distance from dipole center
	double l = sqrt(Point[1]*Point[1]+Point[2]*Point[2])/cm;	// perpendicular distance from dipole center
	double r = sqrt(z0*z0+l*l);									// total distance from dipole center
	const double m = 1.2e3*tesla;								// dipole strength
	
	double r3 = r*r*r;
	double b0 = 3*m*(z0/r)/r3;
	Bfield[0] += b0*(z0/r)-m/r3;
	Bfield[1] += b0*(Point[1]/cm/r);
	Bfield[2] += b0*(Point[2]/cm/r);
}

G4double bmField::calculate_BZ(const double X, const double Y, const double Z) const
{
  double P = sqrt(X*X + Y*Y);
  double ARGZ, CH, TH, D0BZ, D2BZ, D4BZ, BZ;
  ARGZ = D*(abs(Z)-C);
  CH = cosh(ARGZ)*cosh(ARGZ);
  TH = tanh(ARGZ);
  D0BZ = A - B*TH;
  D2BZ = 2.*A*D*D*TH/CH;
  D4BZ = 8.*A*D*D*D*D*(TH*TH*TH/CH - 2.*TH/(CH*CH));
  
  BZ = D0BZ - 0.25*D2BZ*P*P + D4BZ*P*P*P*P/32.;
  return BZ;
}

double bmField::calculate_BP(const double X, const double Y, const double Z) const
{
  double P = sqrt(X*X + Y*Y);
  double ARGZ, CH, TH, DBZ, D3BZ, D5BZ, BP;
  
  ARGZ = D*(abs(Z) - C);

  CH   = cosh(ARGZ)*cosh(ARGZ);
  TH   = tanh(ARGZ)*tanh(ARGZ);

  DBZ  = -B*D/CH;

  D3BZ = 2*B*D*D*D/CH*(1./CH- 2.*TH);
     
  D5BZ = 8.*B*D*D*D*D*D/CH*(6.*TH/CH - (TH - 2./CH) * (2.*TH + 1./CH));

  BP   = -DBZ*P/2. + D3BZ*P*P*P/16. + D5BZ*P*P*P*P*P/C;
  return BP;
}


void bmField::GetFieldValue(const G4double Point[3], G4double *Bfield) const {
  
  if (RobbyField) 
    {

      if(Point[0]*Point[0]+Point[1]*Point[1]>rmax2 || !fieldScale) {
	Bfield[0] = Bfield[1] = Bfield[2] = 0;
	return;
      }
      
      /*
      double x = Point[0]/cm;
      double y = Point[1]/cm;
      double z = Point[2]/cm;
      double P = sqrt(x*x + y*y);
      double psi = asin(y/P);
      
      G4double Bz = calculate_BZ(x, y, z)*tesla;
      G4double Bp = calculate_BP(x, y, z)*tesla;
      //G4double Bz = calculate_BZ(Point[0]/cm, Point[1]/cm, Point[2]/cm);
      //G4double Bp = calculate_BP(Point[0]/cm, Point[1]/cm, Point[2]/cm);
      G4double Bx = 0.0*tesla, By = 0.0*tesla;

      if (x==0.0 && y==0.0)
	{
	  Bx=0.0;
	  By=0.0;
	}
  
      else
	{
	  Bx=Bp*cos(psi);
	  By=Bp*sin(psi);
	  if (x<0.0){Bx=-1.0*abs(Bx);}
	  if (y<0.0){By=-1.0*abs(By);}
	}

      if (z<0.0)
	{
	  Bx=-1.0*Bx;
	  By=-1.0*By;
	}
      Bfield[0] = Bx; 
      Bfield[1] = By;
      Bfield[2] = Bz;
      */
      
      G4double p0 = Point[0]/cm, p1 = Point[1]/cm, p2 = Point[2]/cm;
      double b0 = 0., b1=0., b2=0.;
      //cout << Bfield[0] << " " << Bfield[1] << " " << Bfield[2] << endl;
      b_field_(&b0, &b1, &b2, &p0, &p1, &p2);
      Bfield[0] = b0*tesla;
      Bfield[1] = b1*tesla;
      Bfield[2] = b2*tesla;
    }
      
      
  else 
    {
      G4double z=Point[2];	// point z
      unsigned int zindex = int(lower_bound(Zpoints.begin(), Zpoints.end(), z)-Zpoints.begin());	// location in points list
      
      // no field defined outside experimental volume
      if(zindex==0 || zindex>=Zpoints.size() || Point[0]*Point[0]+Point[1]*Point[1]>rmax2 || !fieldScale) {
	Bfield[0] = Bfield[1] = Bfield[2] = 0;
	return;
      }
	
      // interpolate between defined regions
      G4double base = 0.5*(Bpoints[zindex-1]+Bpoints[zindex]);// midpoint value
      G4double amp = 0.5*(Bpoints[zindex-1]-Bpoints[zindex]);	// variation amplitude between ends
      G4double dz = Zpoints[zindex]-Zpoints[zindex-1];		// z distance between ends
      G4double l = (z-Zpoints[zindex-1])/dz;					// fractional distance between ends
      
      Bfield[2] = base*fieldScale;
      if(amp) {
	Bfield[2] += amp*cos(l*M_PI)*fieldScale; // interpolate B_z component with cosine
	// B_r component to obey Maxwell equation grad dot B = dB_z/dz + 1/r d(r B_r)/dr = 0
	G4double Brtemp = amp*M_PI*sin(l*M_PI)/(2*dz)*fieldScale;
	Bfield[0]=Point[0]*Brtemp;
	Bfield[1]=Point[1]*Brtemp;
      } else {
	Bfield[0]=Bfield[1]=0.0;
      }
    }	
  if(addAFP) addAFPFringeField(Point,Bfield);
}
