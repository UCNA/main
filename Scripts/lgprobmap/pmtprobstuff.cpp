int sign(double value) {
    return ( value > 0) - (value < 0);
}


void twoDvec2txt(const char* filename,double number1[],double number2[], double value[],int length){     
 
  ofstream myfile;
  myfile.open (filename);
 	for(int i=0; i<length; i++) {
  myfile << number1[i] <<" "<< number2[i] << " " << value[i]<<"\n"; 
  }
  myfile.close();
}

