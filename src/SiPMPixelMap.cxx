#include "SiPMPixelMap.hh"

SiPMPixelMap::SiPMPixelMap();
SiPMPixelMap::~SiPMPixelMap();

bool timesort (std::pair<double,double> i,std::pair<double,double> j) { return (i.first<j.first); }

int SiPMPixelMap::getpixel(TVector3 pos)
{
    //GetPosition on SiPM
    pos = pos - center;
    TVector3 vert = normal*pos;
    pos = pos - vert;
    pos = pos.Rotate(angle,normal);
    double x = pos.X()+halfx;
    double y = pos.Y()+halfy;
    
    int xind = TMath::Max(0,TMath::Min(x/xpitch,nxpix-1));
    int yind = TMath::Max(0,TMath::Min(y/ypitch,nypix-1));
    
    return nxpix*yind+xind;
}

void SiPMPixelMap::AddPhotonHit(double time,double wavelength,TVector3 pos)
{
    int pixelind = getpixel(pos);
    int i=0;
    while(i<pixelnumer.size() && pixelnumer[i] != pixelind)++i;
    
    if(i<pixelnumer.size()){
        hittimes[i].push_back(std::pair<double,double>(time,wavelength));
    }else{
        pixelnumer.push_back(pixelind);
        hittimes.push_back(std::vector<std::pair<double,double> >(1,std::pair<double,double>(time,wavelength)));
    }
    ++totalcounts;
}

void SiPMPixelMap::TimeOrderHits()
{
    if(totalcounts>1){
        hittimes_ordered = hittimes;
    }
    
    for(int i=0; i<hittimes.size().size(); ++i){
        hittimes_ordered.push_back(std::sort(hittimes[i].begin(),hittimes[i].end(),timesort));
    }
}

double SiPMPixelMap::GetCharge(double deltaT){
    if(deltaT<0.0) return rand->Gaus(SPE,sigma);
    return rand->Gaus(SPE(1-TMath::Exp(deltaT/recoverytime)),sigma);
}

void SiPMPixelMap::FillOverlaps(){
    for(int i=0; i<hittimes_ordered.size(); ++i){
        time_charge.push_back(std::pair<double,double>(hittimes_ordered[i][0].first,GetCharge(-1.)));
        for(int j=1; j<hittimes_ordered[i].size(); ++j){
            time_charge.push_back(std::pair<double,double>(hittimes_ordered[i][0].first,GetCharge(hittimes_ordered[i][j]-hittimes_ordered[i][j-1])));
        }
    }
}

double SiPMPixelMap::GetSiPMCurrent(double t){
    double value = 0.;
    for(int i=0; i<time_charge.size(); ++i){
        value += time_charge[i].second*pulsefunc->Eval(t-time_charge[i]);
    }
    return value;
}

void SiPMPixelMap::ClearSiPM()
{
    pixelnumber.clear();
    hittimes.clear();
    hittimes_ordered.clear();
    time_charge.clear();
    totalcounts=0;
}
