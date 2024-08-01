// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MyDict
#define R__NO_DEPRECATION

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

// Header files passed as explicit arguments
#include "TPseudoRestFrame.h"
#include "RooTauLeptonInvisible.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace RooFit {
   namespace ROOTDict {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *RooFit_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("RooFit", 0 /*version*/, "RooPrintable.h", 64,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &RooFit_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *RooFit_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static TClass *TPseudoRestFrame_Dictionary();
   static void TPseudoRestFrame_TClassManip(TClass*);
   static void *new_TPseudoRestFrame(void *p = 0);
   static void *newArray_TPseudoRestFrame(Long_t size, void *p);
   static void delete_TPseudoRestFrame(void *p);
   static void deleteArray_TPseudoRestFrame(void *p);
   static void destruct_TPseudoRestFrame(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPseudoRestFrame*)
   {
      ::TPseudoRestFrame *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TPseudoRestFrame));
      static ::ROOT::TGenericClassInfo 
         instance("TPseudoRestFrame", "TPseudoRestFrame.h", 23,
                  typeid(::TPseudoRestFrame), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TPseudoRestFrame_Dictionary, isa_proxy, 4,
                  sizeof(::TPseudoRestFrame) );
      instance.SetNew(&new_TPseudoRestFrame);
      instance.SetNewArray(&newArray_TPseudoRestFrame);
      instance.SetDelete(&delete_TPseudoRestFrame);
      instance.SetDeleteArray(&deleteArray_TPseudoRestFrame);
      instance.SetDestructor(&destruct_TPseudoRestFrame);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPseudoRestFrame*)
   {
      return GenerateInitInstanceLocal((::TPseudoRestFrame*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TPseudoRestFrame*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TPseudoRestFrame_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TPseudoRestFrame*)0x0)->GetClass();
      TPseudoRestFrame_TClassManip(theClass);
   return theClass;
   }

   static void TPseudoRestFrame_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","TPseudoRestFrame.h");
   }

} // end of namespace ROOT

namespace ROOT {
   static void delete_RooTauLeptonInvisible(void *p);
   static void deleteArray_RooTauLeptonInvisible(void *p);
   static void destruct_RooTauLeptonInvisible(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooTauLeptonInvisible*)
   {
      ::RooTauLeptonInvisible *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooTauLeptonInvisible >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooTauLeptonInvisible", ::RooTauLeptonInvisible::Class_Version(), "RooTauLeptonInvisible.h", 20,
                  typeid(::RooTauLeptonInvisible), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooTauLeptonInvisible::Dictionary, isa_proxy, 4,
                  sizeof(::RooTauLeptonInvisible) );
      instance.SetDelete(&delete_RooTauLeptonInvisible);
      instance.SetDeleteArray(&deleteArray_RooTauLeptonInvisible);
      instance.SetDestructor(&destruct_RooTauLeptonInvisible);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooTauLeptonInvisible*)
   {
      return GenerateInitInstanceLocal((::RooTauLeptonInvisible*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooTauLeptonInvisible*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooTauLeptonInvisible::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooTauLeptonInvisible::Class_Name()
{
   return "RooTauLeptonInvisible";
}

//______________________________________________________________________________
const char *RooTauLeptonInvisible::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooTauLeptonInvisible*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooTauLeptonInvisible::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooTauLeptonInvisible*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooTauLeptonInvisible::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooTauLeptonInvisible*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooTauLeptonInvisible::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooTauLeptonInvisible*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPseudoRestFrame(void *p) {
      return  p ? new(p) ::TPseudoRestFrame : new ::TPseudoRestFrame;
   }
   static void *newArray_TPseudoRestFrame(Long_t nElements, void *p) {
      return p ? new(p) ::TPseudoRestFrame[nElements] : new ::TPseudoRestFrame[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPseudoRestFrame(void *p) {
      delete ((::TPseudoRestFrame*)p);
   }
   static void deleteArray_TPseudoRestFrame(void *p) {
      delete [] ((::TPseudoRestFrame*)p);
   }
   static void destruct_TPseudoRestFrame(void *p) {
      typedef ::TPseudoRestFrame current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TPseudoRestFrame

//______________________________________________________________________________
void RooTauLeptonInvisible::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooTauLeptonInvisible.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RooTauLeptonInvisible::Class(),this);
   } else {
      R__b.WriteClassBuffer(RooTauLeptonInvisible::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooTauLeptonInvisible(void *p) {
      delete ((::RooTauLeptonInvisible*)p);
   }
   static void deleteArray_RooTauLeptonInvisible(void *p) {
      delete [] ((::RooTauLeptonInvisible*)p);
   }
   static void destruct_RooTauLeptonInvisible(void *p) {
      typedef ::RooTauLeptonInvisible current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RooTauLeptonInvisible

namespace {
  void TriggerDictionaryInitialization_MyDict_Impl() {
    static const char* headers[] = {
"TPseudoRestFrame.h",
"RooTauLeptonInvisible.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/belle.cern.ch/el7/externals/v01-10-02/Linux_x86_64/common/root/include/",
"/gpfs/home/belle2/johancol/tau3x1/RooFitClass/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MyDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate(R"ATTRDUMP(file_name@@@TPseudoRestFrame.h)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TPseudoRestFrame.h")))  TPseudoRestFrame;
class __attribute__((annotate(R"ATTRDUMP(file_name@@@RooTauLeptonInvisible.h)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RooTauLeptonInvisible.h")))  RooTauLeptonInvisible;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MyDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "TPseudoRestFrame.h"
#include "RooTauLeptonInvisible.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"", payloadCode, "@",
"RooTauLeptonInvisible", payloadCode, "@",
"RooTauLeptonInvisible::fgIsA", payloadCode, "@",
"TPseudoRestFrame", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MyDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MyDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MyDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MyDict() {
  TriggerDictionaryInitialization_MyDict_Impl();
}
