#include "LOLXReadData.hh"

namespace LOLXReadData{
    
    TGraph* PMMA_nREAL=NULL;
    int state_PMMA_nREAL=0;
    TGraph* PMMA_nIMAG=NULL;
    int state_PMMA_nIMAG=0;
    TGraph* Filter_Trans=NULL;
    int state_Filter_Trans=0;
    TGraph* eStop_Col=NULL;
    int state_eStop_Col=0;
    TGraph* eStop_Rad=NULL;
    int state_eStop_Rad=0;
    TGraph* eStop_Total=NULL;
    int state_eStop_Total=0;
    TGraph* eStop_DensCor=NULL;
    int state_eStop_DensCor=0;
    TGraph* SiPM_eff = NULL;
    int state_SiPM_eff = 0;
    TGraph* xenon_Index = NULL;
    int state_xenon_Index = 0;
    
    int ReadData(char* filename,TGraph* &graph,std::string target,int unit){
        std::ifstream in(Form("%s/%s",std::getenv("LOLXSTLDIR"),filename));
        // Check if object is valid
        if(!in){
            std::cerr << "Cannot open the File : "<<Form("%s/%s",std::getenv("LOLXSTLDIR"),filename)<<std::endl;
            return false;
        }
        //printf("target = %s\n",target.data());
        std::string str;
        // Read the next line from File untill it reaches the end.
        while (std::getline(in, str))
        {
            if(str == target)break;
        }
        
        int state = 0;
        std::string LogX("LogX");
        std::string LOGX("LOGX");
        std::string logx("logx");
        std::string Logx("Logx");
        std::string logX("logX");
        std::string LogY("LogY");
        std::string LOGY("LOGY");
        std::string logy("logy");
        std::string logY("logY");
        std::string Logy("Logy");
        std::string LogXY("LogXY");
        std::string LOGXY("LOGXY");
        std::string Logxy("Logxy");
        std::string logxy("logxy");
        std::string logXY("logXY");
        
        std::string reference;
        std::getline(in, reference);
        std::string extrapolation;
        std::getline(in, extrapolation);
        std::string inputunit;
        std::getline(in, inputunit);
        std::string x;
        std::getline(in, x);
        std::string y;
        std::getline(in,y);
        
        if(extrapolation == LogX) state = 1;
        else if(extrapolation == LOGX) state = 1;
        else if(extrapolation == logx) state = 1;
        else if(extrapolation == Logx) state = 1;
        else if(extrapolation == logX) state = 1;
        else if(extrapolation == LogY) state = 2;
        else if(extrapolation == LOGY) state = 2;
        else if(extrapolation == logy) state = 2;
        else if(extrapolation == Logy) state = 2;
        else if(extrapolation == LogY) state = 2;
        else if(extrapolation == LogXY) state = 3;
        else if(extrapolation == LOGXY) state = 3;
        else if(extrapolation == logxy) state = 3;
        else if(extrapolation == Logxy) state = 3;
        else if(extrapolation == logXY) state = 3;
        
        std::replace(x.begin(),x.end(),',',' ');
        std::replace(y.begin(),y.end(),',',' ');
        
        std::vector<double> x_array;
        std::stringstream ss(x);
        double temp;
        while (ss >> temp){
            if(state == 1 || state == 3){
                if(unit == 1) x_array.push_back(1240./TMath::Log10(temp));
                else x_array.push_back(TMath::Log10(temp));
            }else{
                if(unit == 1) x_array.push_back(1240./temp);
                else x_array.push_back(temp);
            }
            //printf("x %d = %f\n",(int)x_array.size(),x_array[x_array.size()-1]);
        }
        std::vector<double> y_array;
        ss = std::stringstream(y);
        while (ss >> temp){
            if(state == 2 || state == 3){
                y_array.push_back(TMath::Log10(temp));
            }else{
                y_array.push_back(temp);
            }
            //printf("y %d = %f\n",(int)y_array.size(),y_array[y_array.size()-1]);
        }
        assert(x_array.size() == y_array.size());
        
        graph = new TGraph();
        
        if(unit==1){
            for(int i=x_array.size()-1; i>-1; --i){
                //printf("%d %d %f %f\n",i,(x_array.size()-1)-i,x_array[i],y_array[i]);
                graph->SetPoint((x_array.size()-1)-i,x_array[i],y_array[i]);
                //printf("%d %f %f\n",i,x_array[i],graph->Eval(x_array[i]-0.01));
            }
        }else{
            for(int i=0; i<x_array.size(); ++i){
                //printf("%d %f %f\n",i,x_array[i],y_array[i]);
                graph->SetPoint(i,x_array[i],y_array[i]);
            }
        }
        //Close The File
        in.close();
        
        return state;
    }
    
    double GetValue(int state, TGraph* graph,double eV){
        if(state==0) return graph->Eval(eV);
        if(state==1) return graph->Eval(TMath::Log10(eV));
        if(state==2) return TMath::Power(10.,graph->Eval(eV));
        if(state==3) return TMath::Power(10.,graph->Eval(TMath::Log10(eV)));
    }
    
    double GetPMMA_nREAL(double eV){
        if(!PMMA_nREAL) state_PMMA_nREAL = ReadData("PMMAOpticalProperties.txt",PMMA_nREAL,std::string("PMMA nREAL"),1);
        return GetValue(state_PMMA_nREAL,PMMA_nREAL,eV);
    }
    
    double GetPMMA_nIMAG(double eV){
        if(!PMMA_nIMAG) state_PMMA_nIMAG = ReadData("PMMAOpticalProperties.txt",PMMA_nIMAG,std::string("PMMA nIMAG"),1);
        return GetValue(state_PMMA_nIMAG,PMMA_nIMAG,eV);
    }
    
    double GetFilter_Trans(double eV){
        if(!Filter_Trans) state_Filter_Trans = ReadData("FilterOpticalProperties.txt",Filter_Trans,std::string("Filter Transmission"),1);
        return GetValue(state_Filter_Trans,Filter_Trans,eV);
    }
    
    double GeteStop_Collision(double eV){
        if(!eStop_Col) state_eStop_Col = ReadData("XenonData.txt",eStop_Col,std::string("eStop Collision"),0);
        return GetValue(state_eStop_Col,eStop_Col,eV);
    }
    
    double GeteStop_Radiation(double eV){
        if(!eStop_Rad) state_eStop_Rad = ReadData("XenonData.txt",eStop_Rad,std::string("eStop Radiation"),0);
        return GetValue(state_eStop_Rad,eStop_Rad,eV);
    }
    
    double GeteStop_Total(double eV){
        if(!eStop_Total) state_eStop_Total = ReadData("XenonData.txt",eStop_Total,std::string("eStop Total"),0);
        return GetValue(state_eStop_Total,eStop_Total,eV);
    }
    
    double GeteStop_DensityCorrection(double eV){
        if(!eStop_DensCor) state_eStop_DensCor = ReadData("XenonData.txt",eStop_DensCor,std::string("eStop DensityCorrection"),0);
        return GetValue(state_eStop_DensCor,eStop_DensCor,eV);
    }
    
    double GetXenon_n(double eV){
        //from https://arxiv.org/pdf/1502.04213.pdf  AS Feb 21, 2018
        double WL = 1240./eV;
        if(WL>700) return 1.365;
        if(WL<147.) return TMath::Sqrt(1.5+ 0.38*147.*147./(147.*147.-146.9*146.9)+ 0.009*147.*147./(147.*147.-827*827));
        return TMath::Sqrt(1.5+ 0.38*WL*WL/(WL*WL-146.9*146.9)+ 0.009*WL*WL/(WL*WL-827*827));
    }
    
    double GetXenon_Scintillation(double eV){
        double WL = 1240./eV;
        return TMath::Gaus(WL, 178., 5., true);
    }
    
    double GetXenon_Absorption(double eV){
        //fix
        return 200.;
    }
    
    double GetXenon_Rayleigh(double eV){
        //fix
        return 35.;
    }
    
    double GetXenon_IndexofRefraction(double eV){
        if(!xenon_Index) state_xenon_Index = ReadData("XenonData.txt",xenon_Index,std::string("Index of Refraction"),0);
        return GetValue(state_xenon_Index,xenon_Index,eV);
    }
    
    double GetXenon_n_stitch(double eV){
        if(eV<8.) return GetXenon_n(eV);
        if(eV>11.0) return 1.0;
        else return TMath::Max(1.0,GetXenon_IndexofRefraction(eV)*GetXenon_n(8.0)/GetXenon_IndexofRefraction(8.0));
    }
    
    double GetSiPM_Efficiency(double eV){
        if(!SiPM_eff) state_eStop_DensCor = ReadData("SiPMEfficiency.txt",SiPM_eff,std::string("SiPM Efficiency"),1);
        //printf("e = %f v = %f\n",eV,SiPM_eff->Eval(eV));
        return GetValue(state_SiPM_eff,SiPM_eff,eV);
    }
    
}
