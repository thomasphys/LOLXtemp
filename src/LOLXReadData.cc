#include "LOLXReadData.hh"

namespace DMaterials {
    
    TGraph* PMMA_nREAL=NULL;
    TGraph* PMMA_nIMAG=NULL;
    
    void GetPMMA_nREAL(double eV){
        if(!PMMA_nREAL){
            std::ifstream in(Form("%s/PMMAOpticalProperties.txt",std::getenv("LOLXSTLDIR")));
            // Check if object is valid
            if(!in){
                std::cerr << "Cannot open the File : "<<fileName<<std::endl;
                return false;
            }
            
            std::string str;
            std::string target("PMMA nREAL");
            // Read the next line from File untill it reaches the end.
            while (std::getline(in, str))
            {
                if(str == target)break;
            }
            
            std::string extrapolation;
            std::getline(in, extrapolation);
            std::string wavelengths;
            std::getline(in, wavelengths);
            std::string values;
            std::getline(in,values);
            
            std::replace(wavelengths.begin(), wavelengths.end(), ',', ' ');
            std::replace(values.begin(), values.end(), ',', ' ');
            
            vector<double> wavelengths_array;
            stringstream ss(wavelengths);
            double temp;
            while (ss >> temp) wavelengths_array.push_back(temp);
            
            vector<double> value_array;
            ss = stringstream(values);
            while (ss >> temp) value_array.push_back(temp);
            
            assert(wavelength_array.size() == value_array.size());
            
            PMMA_nREAL = new TGraph();
            
            for(int i=wavelength_array.size()-1; i>-1; --i){
                PMMA_nREAL->SetPoint(i,1240./wavelength_array[i],value_array[i]);
            }
            //Close The File
            in.close();
        }
        return PMMA_nREAL->Eval(eV);
    }
    
    void GetPMMA_nIMAG(double eV){
        if(!PMMA_nIMAG){
                std::ifstream in(Form("%s/PMMAOpticalProperties.txt",std::getenv("LOLXSTLDIR")));
                // Check if object is valid
                if(!in){
                    std::cerr << "Cannot open the File : "<<fileName<<std::endl;
                    return false;
                }
                
                std::string str;
                std::string target("PMMA nIMAG");
                // Read the next line from File untill it reaches the end.
                while (std::getline(in, str))
                {
                    if(str == target)break;
                }
                
                std::string extrapolation;
                std::getline(in, extrapolation);
                std::string wavelengths;
                std::getline(in, wavelengths);
                std::string values;
                std::getline(in,values);
                
                std::replace(wavelengths.begin(), wavelengths.end(), ',', ' ');
                std::replace(values.begin(), values.end(), ',', ' ');
                
                vector<double> wavelengths_array;
                stringstream ss(wavelengths);
                double temp;
                while (ss >> temp) wavelengths_array.push_back(temp);
                
                vector<double> value_array;
                ss = stringstream(values);
                while (ss >> temp) value_array.push_back(temp);
                
                assert(wavelength_array.size() == value_array.size());
            
                PMMA_nIMAG = new TGraph();
            
                for(int i=wavelength_array.size()-1; i>-1; --i){
                    PMMA_nIMAG->SetPoint(i,1240./wavelength_array[i],value_array[i]);
                }
                //Close The File
                in.close();
        }
        
        return PMMA_nIMAG->Eval(eV);
    }
    
}
