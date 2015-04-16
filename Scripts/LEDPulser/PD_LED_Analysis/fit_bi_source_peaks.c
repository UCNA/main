// Quick script to fit bi calibration source peak to a gaussian

#include <iostream>
using namespace std;

//int main (int argc, char **argv)
int fit_bi_source_peaks(int argc, char **argv)
{
  TString filename = argv[1];
  TChain h1("h1");
  
  if (not h1.Add(filename))
    {
      printf("Could not find root file full%s.root\n", filename);
      exit(1);
    }
  
  cout << "Found " << filename << endl;
  
  h9998.Draw();

}
