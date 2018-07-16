#ifndef __DS_Event__
#define __DS_Event__

#include <vector>
#include <string>
#include <utility>
#include <iostream>

#include <TObject.h>
#include "MPPC.h" 

		class Event : public TObject
		{
			public:
			Event() : TObject() { Init(); }
      			Event(const Event &rhs) : TObject(rhs) { Init(); CopyObj(rhs); }
      			virtual ~Event() { Destroy(); }
      			virtual Event &operator=(const Event &rhs) {TObject::operator=(rhs); CopyObj(rhs); return *this; }

      			/// Clone this object. This overrides the default ROOT TObject::Clone,
      			//  as that version corrupts old QT and Pulse objects.
      			virtual TObject* Clone(const char* opt = "") const { (void)opt; return new Event(*this); }
      			virtual void Clear(Option_t* = "") { Destroy(); Init(); }	
			
			void SetSensorNum(int nsens){PhotoSensors = std::vector<MPPC*>(nsens,new MPPC());}
			
			void AddMPPC(MPPC* _sensor){PhotoSensors.push_back(_sensor);}
            MPPC* GetMPPC(int i){return PhotoSensors[i];}
            int GetMPPCCount(){return PhotoSensors.size();}
			void SetEventID(long int _eventid){eventid = _eventid;}

			void SetMCSintillation(int _MCSintillation){MCSintillation = _MCSintillation;}
			int GetMCSintillation(){return MCSintillation;}
			void SetMCCherenkov(int _MCCherenkov){MCCherenkov = _MCCherenkov;}
                        int GetMCCherenkov(){return MCCherenkov;}

			void AddCherenkovWavelength(double _wl){CherenkovWavelength.push_back(_wl);}
                        void AddElectronBeta(double _beta){ElectronBeta.push_back(_beta);}
			int GetCherenkovSpectrumCount(){return (int)CherenkovWavelength.size();}
			double GetCherenkovWavelength(int i){CherenkovWavelength[i];}
                        double GetElectronBeta(int i){ElectronBeta[i];}

            ClassDef(Event,1);

			private:

			void Init()
			{
				eventseed = 0;
			}

			void Destroy()
			{
				PhotoSensors.resize(0);
				ElectronBeta.resize(0);
				CherenkovWavelength.resize(0);
			}

			virtual void CopyObj(const Event &rhs){
				eventid = rhs.eventid;
				eventseed = rhs.eventseed;
				PhotoSensors = rhs.PhotoSensors;
			}
			
			long int eventid;
			long int eventseed; //event random seed for MC regeneration.
			std::vector<MPPC*> PhotoSensors; /// this is where the different hit sensors are stored.
			int MCSintillation;
			int MCCherenkov;
			std::vector<double> CherenkovWavelength;
			std::vector<double> ElectronBeta;
			
		};

#endif 
				
