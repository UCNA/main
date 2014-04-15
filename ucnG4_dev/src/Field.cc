#include <cmath>
#include <fstream>

#include "Field.hh"
#include "FieldMessenger.hh"

#include <TString.h>

#include <G4SystemOfUnits.hh>

Field::Field(const G4String& filename): myMessenger(new FieldMessenger(this)), rmax2(20*20*cm2), fieldScale(1.0), afp_m(0.) {
	LoadFieldMap(filename);
}

void Field::LoadFieldMap(const G4String& filename) {
	
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
		fin.open(filename.c_str());
		if(!fin) {
			G4cout << "Can not open " << filename << G4endl;
			exit(1);
		}
		G4cout << "Loading field profile from " << filename << "...\n";
		
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

void Field::addAFPFringeField(const G4double Point[3], G4double *Bfield) const {
	double z0 = Point[0]-(-280.*cm);							// z distance from dipole center TODO: what is actual distance??
	double l = sqrt(Point[1]*Point[1]+Point[2]*Point[2]);		// perpendicular distance from dipole center
	double r = sqrt(z0*z0+l*l);									// total distance from dipole center
	
	double r3 = r*r*r;
	double mu_0_4pi_m = afp_m * m*m*m * tesla * 1e-7;			// afp_m [A*m^2] * (mu_0/4pi = 1e-7 T*m/A)
	double x = 3 * mu_0_4pi_m * z0/(r3*r*r);
	
	Bfield[0] += x*z0 - mu_0_4pi_m/r3;
	Bfield[1] += x*Point[1];
	Bfield[2] += x*Point[2];
}

void Field::GetFieldValue(const G4double Point[3], G4double *Bfield) const {
	
	G4double z=Point[2];	// point z
	unsigned int zindex = int(lower_bound(Zpoints.begin(), Zpoints.end(), z)-Zpoints.begin());	// location in points list
	
	if(zindex==0 || zindex>=Zpoints.size() || Point[0]*Point[0]+Point[1]*Point[1]>rmax2 || !fieldScale) {
		// no field defined outside experimental volume
		Bfield[0] = Bfield[1] = Bfield[2] = 0;
	} else {
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
	
	if(afp_m) addAFPFringeField(Point,Bfield);
}
