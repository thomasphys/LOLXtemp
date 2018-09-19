#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "TGraph.h"
#include "TMath.h"

namespace LOLXReadData{
    extern TGraph* PMMA_nREAL;
    extern int state_PMMA_nREAL;
    extern TGraph* PMMA_nIMAG;
    extern int state_PMMA_nIMAG;
    extern TGraph* Filter_Trans;
    extern int state_Filter_Trans;
    extern TGraph* eStop_Col;
    extern int state_eStop_Col;
    extern TGraph* eStop_Rad;
    extern int state_eStop_Rad;
    extern TGraph* eStop_Total;
    extern int state_eStop_Total;
    extern TGraph* eStop_DensCor;
    extern int state_eStop_DensCor;
    extern TGraph* SiPM_eff;
    extern int state_SiPM_eff;
    extern TGraph* xenon_Index;
    extern int state_xenon_Index;
    
    int ReadData(char* filename,TGraph* &graph,std::string target,int unit);
    double GetValue(int state,TGraph* graph,double eV);
    
    double GetPMMA_nREAL(double eV);
    double GetPMMA_nIMAG(double eV);
    double GetFilter_Trans(double eV);
    double GeteStop_Collision(double eV);
    double GeteStop_Radiation(double eV);
    double GeteStop_Total(double eV);
    double GeteStop_DensityCorrection(double eV);
    double GetXenon_n(double eV);
    double GetXenon_Scintillation(double eV);
    double GetXenon_Absorption(double eV);
    double GetXenon_Rayleigh(double eV);
    double GetXenon_IndexofRefraction(double eV);
    double GetXenon_n_stitch(double eV);
    double GetSiPM_Efficiency(double eV);
    
}
