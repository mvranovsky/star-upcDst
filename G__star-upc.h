/********************************************************************
* G__star-upc.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error G__star-upc.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtableG__starmIupc();
extern void G__cpp_setup_inheritanceG__starmIupc();
extern void G__cpp_setup_typetableG__starmIupc();
extern void G__cpp_setup_memvarG__starmIupc();
extern void G__cpp_setup_globalG__starmIupc();
extern void G__cpp_setup_memfuncG__starmIupc();
extern void G__cpp_setup_funcG__starmIupc();
extern void G__set_cpp_environmentG__starmIupc();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"
#include "StUPCTofHit.h"
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"
#include "StUPCRpsCluster.h"
#include "StPicoPhysicalHelix.h"
#include "StPicoHelix.h"
#include "StUPCV0.h"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__G__starmIupcLN_TClass;
extern G__linked_taginfo G__G__starmIupcLN_TBuffer;
extern G__linked_taginfo G__G__starmIupcLN_TMemberInspector;
extern G__linked_taginfo G__G__starmIupcLN_TObject;
extern G__linked_taginfo G__G__starmIupcLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__G__starmIupcLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__G__starmIupcLN_TClonesArray;
extern G__linked_taginfo G__G__starmIupcLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__G__starmIupcLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__G__starmIupcLN_StUPCTrack;
extern G__linked_taginfo G__G__starmIupcLN_StUPCBemcCluster;
extern G__linked_taginfo G__G__starmIupcLN_StUPCTofHit;
extern G__linked_taginfo G__G__starmIupcLN_TIterator;
extern G__linked_taginfo G__G__starmIupcLN_StUPCVertex;
extern G__linked_taginfo G__G__starmIupcLN_TParticle;
extern G__linked_taginfo G__G__starmIupcLN_TArrayI;
extern G__linked_taginfo G__G__starmIupcLN_StUPCEvent;
extern G__linked_taginfo G__G__starmIupcLN_TLorentzVector;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTBaselEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTBaselEdoublegR;
extern G__linked_taginfo G__G__starmIupcLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__G__starmIupcLN_TElementActionTlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TElementPosActionTlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTRow_constlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTRowlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTDiag_constlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTColumn_constlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTFlat_constlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTSub_constlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTSparseRow_constlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTSparseDiag_constlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTColumnlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTDiaglEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTFlatlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTSublEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTSparseRowlEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TMatrixTSparseDiaglEfloatgR;
extern G__linked_taginfo G__G__starmIupcLN_TVector3;
extern G__linked_taginfo G__G__starmIupcLN_StUPCTrackcLcLFlag;
extern G__linked_taginfo G__G__starmIupcLN_StUPCTrackcLcLPart;
extern G__linked_taginfo G__G__starmIupcLN_StUPCRpsCluster;
extern G__linked_taginfo G__G__starmIupcLN_StUPCRpsTrack;
extern G__linked_taginfo G__G__starmIupcLN_StUPCRpsTrackPoint;
extern G__linked_taginfo G__G__starmIupcLN_StRPEvent;
extern G__linked_taginfo G__G__starmIupcLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__G__starmIupcLN_StUPCRpsTrackcLcLStRpsTrackType;
extern G__linked_taginfo G__G__starmIupcLN_StUPCRpsTrackcLcLStRpsAngles;
extern G__linked_taginfo G__G__starmIupcLN_StUPCRpsTrackPointcLcLStUPCRpsTrackPointQuality;
extern G__linked_taginfo G__G__starmIupcLN_StPicoHelix;
extern G__linked_taginfo G__G__starmIupcLN_pairlEdoublecOdoublegR;
extern G__linked_taginfo G__G__starmIupcLN_StPicoPhysicalHelix;
extern G__linked_taginfo G__G__starmIupcLN_StUPCV0;

/* STUB derived class for protected member access */
