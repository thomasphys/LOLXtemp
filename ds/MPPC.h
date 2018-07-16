#ifndef __DS_MPPC__
#define __DS_MPPC__

#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <TObject.h>

		class MPPC : public TObject
		{
			public:

			enum mppcmodeindex
			{
				HAMAMATSU,
				MODEL_COUNT	
			};

			MPPC() : TObject() { Init(); }
			MPPC(int mppcid) : TObject() {Init(); id = mppcid;}
      			MPPC(const MPPC &rhs) : TObject(rhs) { Init(); CopyObj(rhs); }
      			virtual ~MPPC() { Destroy(); }
      			virtual MPPC &operator=(const MPPC &rhs) {TObject::operator=(rhs); CopyObj(rhs); return *this; }

      			/// Clone this object. This overrides the default ROOT TObject::Clone,
      			//  as that version corrupts old QT and Pulse objects.
      			virtual TObject* Clone(const char* opt = "") const { (void)opt; return new MPPC(*this); }
      			virtual void Clear(Option_t* = "") { Destroy(); Init(); }	
			
			void SetID(int _id){id = _id;}
			unsigned int GetID(){return id;}
			void AddSensorHit(double time){hittimes.push_back(time);}
			unsigned int GetHitCount(){return (unsigned int)hittimes.size();}
			float GetHitTime(int n){if(n<0 || n>hittimes.size()) return -1.0; else return hittimes[n];}
            
            ClassDef(MPPC,1);

			private:

			void Init()
			{
				id = 0;
				mppcmodel = 0;
				hittimes = std::vector<float>(0,0.0);
			}
			void Destroy()
			{
				hittimes.resize(0);
			}

			virtual void CopyObj(const MPPC &rhs){
				mppcmodel = rhs.mppcmodel;
				id = rhs.id;
				hittimes = rhs.hittimes;
			}

			unsigned short mppcmodel;
			unsigned int id;
			std::vector<float> hittimes; /// photon hit times. 
		};

#endif
