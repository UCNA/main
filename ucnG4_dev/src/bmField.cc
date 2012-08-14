//2-19-02 J. Yuan: magnet field for UCNA experiment
////////////////////////////////////////////////////

#include "bmField.hh"
#include <cmath>
#include <fstream>

bmField::bmField(const TString filename): addAFP(false), rmax2(20*20*cm2), fieldScale(1.0) {
	LoadFieldMap(filename);
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

void bmField::GetFieldValue(const G4double Point[3], G4double *Bfield) const {
	
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
	
	if(addAFP) addAFPFringeField(Point,Bfield);
}
