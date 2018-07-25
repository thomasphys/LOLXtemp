#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TGraph.h"

namespace LOLXData{
    extern TGraph* PMMA_nREAL;
    extern TGraph* PMMA_nIMAG;
    
    void GetPMMA_nREAL(double eV);
    void GetPMMA_nIMAG(double eV);
    
    
}
