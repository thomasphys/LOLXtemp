// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "Event.h"
#include "MPPC.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_MPPC(void *p = 0);
   static void *newArray_MPPC(Long_t size, void *p);
   static void delete_MPPC(void *p);
   static void deleteArray_MPPC(void *p);
   static void destruct_MPPC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MPPC*)
   {
      ::MPPC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MPPC >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MPPC", ::MPPC::Class_Version(), "MPPC.h", 10,
                  typeid(::MPPC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MPPC::Dictionary, isa_proxy, 4,
                  sizeof(::MPPC) );
      instance.SetNew(&new_MPPC);
      instance.SetNewArray(&newArray_MPPC);
      instance.SetDelete(&delete_MPPC);
      instance.SetDeleteArray(&deleteArray_MPPC);
      instance.SetDestructor(&destruct_MPPC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MPPC*)
   {
      return GenerateInitInstanceLocal((::MPPC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MPPC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Event(void *p = 0);
   static void *newArray_Event(Long_t size, void *p);
   static void delete_Event(void *p);
   static void deleteArray_Event(void *p);
   static void destruct_Event(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Event*)
   {
      ::Event *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Event >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Event", ::Event::Class_Version(), "Event.h", 12,
                  typeid(::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Event::Dictionary, isa_proxy, 4,
                  sizeof(::Event) );
      instance.SetNew(&new_Event);
      instance.SetNewArray(&newArray_Event);
      instance.SetDelete(&delete_Event);
      instance.SetDeleteArray(&deleteArray_Event);
      instance.SetDestructor(&destruct_Event);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Event*)
   {
      return GenerateInitInstanceLocal((::Event*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Event*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr MPPC::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MPPC::Class_Name()
{
   return "MPPC";
}

//______________________________________________________________________________
const char *MPPC::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MPPC*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MPPC::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MPPC*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MPPC::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MPPC*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MPPC::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MPPC*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Event::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Event::Class_Name()
{
   return "Event";
}

//______________________________________________________________________________
const char *Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Event*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Event*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Event*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Event*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void MPPC::Streamer(TBuffer &R__b)
{
   // Stream an object of class MPPC.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MPPC::Class(),this);
   } else {
      R__b.WriteClassBuffer(MPPC::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MPPC(void *p) {
      return  p ? new(p) ::MPPC : new ::MPPC;
   }
   static void *newArray_MPPC(Long_t nElements, void *p) {
      return p ? new(p) ::MPPC[nElements] : new ::MPPC[nElements];
   }
   // Wrapper around operator delete
   static void delete_MPPC(void *p) {
      delete ((::MPPC*)p);
   }
   static void deleteArray_MPPC(void *p) {
      delete [] ((::MPPC*)p);
   }
   static void destruct_MPPC(void *p) {
      typedef ::MPPC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MPPC

//______________________________________________________________________________
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(Event::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Event(void *p) {
      return  p ? new(p) ::Event : new ::Event;
   }
   static void *newArray_Event(Long_t nElements, void *p) {
      return p ? new(p) ::Event[nElements] : new ::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_Event(void *p) {
      delete ((::Event*)p);
   }
   static void deleteArray_Event(void *p) {
      delete [] ((::Event*)p);
   }
   static void destruct_Event(void *p) {
      typedef ::Event current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Event

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 216,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 216,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEMPPCmUgR_Dictionary();
   static void vectorlEMPPCmUgR_TClassManip(TClass*);
   static void *new_vectorlEMPPCmUgR(void *p = 0);
   static void *newArray_vectorlEMPPCmUgR(Long_t size, void *p);
   static void delete_vectorlEMPPCmUgR(void *p);
   static void deleteArray_vectorlEMPPCmUgR(void *p);
   static void destruct_vectorlEMPPCmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<MPPC*>*)
   {
      vector<MPPC*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<MPPC*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<MPPC*>", -2, "vector", 216,
                  typeid(vector<MPPC*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEMPPCmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<MPPC*>) );
      instance.SetNew(&new_vectorlEMPPCmUgR);
      instance.SetNewArray(&newArray_vectorlEMPPCmUgR);
      instance.SetDelete(&delete_vectorlEMPPCmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEMPPCmUgR);
      instance.SetDestructor(&destruct_vectorlEMPPCmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<MPPC*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<MPPC*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEMPPCmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<MPPC*>*)0x0)->GetClass();
      vectorlEMPPCmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEMPPCmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEMPPCmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MPPC*> : new vector<MPPC*>;
   }
   static void *newArray_vectorlEMPPCmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MPPC*>[nElements] : new vector<MPPC*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEMPPCmUgR(void *p) {
      delete ((vector<MPPC*>*)p);
   }
   static void deleteArray_vectorlEMPPCmUgR(void *p) {
      delete [] ((vector<MPPC*>*)p);
   }
   static void destruct_vectorlEMPPCmUgR(void *p) {
      typedef vector<MPPC*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<MPPC*>

namespace {
  void TriggerDictionaryInitialization_Dict_Impl() {
    static const char* headers[] = {
"Event.h",
"MPPC.h",
0
    };
    static const char* includePaths[] = {
"/phys/root/rootbuild/include",
"/phys/root/rootbuild/include",
"/home/tmcelroy/lolxsim/ds/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$MPPC.h")))  __attribute__((annotate("$clingAutoload$Event.h")))  MPPC;
class __attribute__((annotate("$clingAutoload$Event.h")))  Event;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "Event.h"
#include "MPPC.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"Event", payloadCode, "@",
"MPPC", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Dict() {
  TriggerDictionaryInitialization_Dict_Impl();
}
