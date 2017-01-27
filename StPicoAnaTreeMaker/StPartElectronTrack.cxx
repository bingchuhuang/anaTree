#include "StPartElectronTrack.h"
#include "StMessMgr.h"
#include "TVector2.h"
#include "TMath.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StBTofUtil/tofPathLength.hh"
#include "PhysicalConstants.h"
#include <climits>

ClassImp(StPartElectronTrack)

   //----------------------------------------------------------------------------------
StPartElectronTrack::StPartElectronTrack() : mId(-1), mPMom(0., 0., 0.), mGMom(0.,0.,0.),mDca(32768),
   mDcaXY(32768),mDcaZ(32768),
   mNHitsFit(0), mNHitsDedx(0), 
   mNSigmaElectron(32768),
   mDnDx(0),mDnDxError(0),
   mBeta(0), mLocalY(32768),
   mBTOWADC0(0), mBTOWE0(0), mBTOWE(0),
   mBEMCDistZ(32768), mBEMCDistPhi(32768), mBSMDNEta(0), mBSMDNPhi(0),
   mBTOWId(0),  mEmcTrgId(-1)
{

}

/////////////////////////////////////////////////////////////////////////////////////////
// t - the global track.  p - the associated primary track from the first primary vertex
/////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------------
   StPartElectronTrack::StPartElectronTrack(StPicoDst *picoDst, StPicoTrack* t, Int_t idx)
: mId(idx), mPMom(0., 0., 0.), mGMom(0.,0.,0.),
   mNHitsFit(0), mNHitsDedx(0), 
   mNSigmaElectron(32768),
   mBeta(0), mLocalY(32768)
{
   mPMom = t->pMom();
   int q      = t->charge();
   //mDedx      = (t->dEdx()*1000. > 65536) ? 65536 : (UShort_t)(TMath::Nint(t->dEdx()*1000.));
   mNHitsFit  = t->nHitsFit()*q;
   //mNHitsMax  = t->nHitsMax();
   mNHitsDedx = (UChar_t)(t->nHitsDedx());
   mNSigmaElectron = (fabs(t->nSigmaElectron() * 1000.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaElectron() * 1000.));


   StThreeVectorF vertexPos = picoDst->event()->primaryVertex();
   StPhysicalHelixD helix = t->helix(picoDst->event()->bField());
   mGMom = t->gMom(vertexPos,picoDst->event()->bField());
   StThreeVectorF dcaPoint = helix.at(helix.pathLength(vertexPos.x(), vertexPos.y()));
   float dcaZ = (dcaPoint.z() - vertexPos.z())*10000.;
   float dcaXY = (helix.geometricSignedDistance(vertexPos.x(),vertexPos.y()))*10000.;
   mDcaZ = dcaZ>32768?32768:(Short_t)dcaZ;
   mDcaXY = dcaXY>32768?32768:(Short_t)dcaXY;

   mDnDx = t->dNdx();
   mDnDxError = t->dNdxError();

   double thePath = helix.pathLength(vertexPos);
   StThreeVectorF dcaPos = helix.at(thePath);

   double dca = (dcaPos-vertexPos).mag();
   bool isHft = t->isHFTTrack();
   if(isHft) dca *= -1;
   mDca = fabs(dca*10000.)>32768? 32768: (Short_t)(dca*10000.);

   int index2TofPid = t->bTofPidTraitsIndex();
   mLocalY = 32768;
   mBeta = 0;
   if (index2TofPid>=0){
      StPicoBTofPidTraits *tofPid = picoDst->btofPidTraits(index2TofPid);
      //mTofMatchFlag = tofPid->btofMatchFlag();
      Float_t mom = mPMom.mag();
      mLocalY = tofPid->btofYLocal()*1000;
      //mLocalZ = tofPid->btofZLocal()*1000;
      Float_t beta = tofPid->btofBeta();
      if(beta<1e-4||beta>=(USHRT_MAX-1)/20000){
         Float_t tof = tofPid->btof();
         StThreeVectorF btofHitPos = tofPid->btofHitPos();
         float L = tofPathLength(&vertexPos, &btofHitPos, helix.curvature()); 
         beta = L/(tof*(c_light/1.0e9));
      }
      mBeta = (UShort_t)(beta*20000);
   }
   

   int index2EmcPid = t->bemcPidTraitsIndex();
   if (index2EmcPid>=0){
      StPicoBEmcPidTraits *emcPid = picoDst->bemcPidTraits(index2EmcPid);
      mBTOWADC0 = emcPid->bemcAdc0();
      mBTOWE0 = emcPid->bemcE0()*1000;
      mBTOWE = emcPid->bemcE()*1000;
      mBEMCDistZ = emcPid->bemcZDist()*1000;
      mBEMCDistPhi = emcPid->bemcPhiDist()*10000;
      mBSMDNEta = emcPid->bemcSmdNEta();
      mBSMDNPhi = emcPid->bemcSmdNPhi();
      mBTOWDistEta = emcPid->btowEtaDist()*10000;
      mBTOWDistPhi = emcPid->btowPhiDist()*10000;

      mBTOWId = emcPid->btowId();
      mEmcTrgId = -1;
   }else{
      mBTOWADC0 = 0;
      mBTOWE0 = 0;
      mBTOWE = 0;
      mBEMCDistZ = 32768;
      mBEMCDistPhi = 32768;
      mBSMDNEta = 0;
      mBSMDNPhi = 0;
      mBTOWDistEta = 32768;
      mBTOWDistPhi = 32768;

      mBTOWId = 0;
      mEmcTrgId = -1;
   }
}

//----------------------------------------------------------------------------------
StPartElectronTrack::~StPartElectronTrack()
{
   /* noop */
}
//----------------------------------------------------------------------------------
void StPartElectronTrack::Print(const Char_t *option) const
{
   if (strcmp(option, "tpc") == 0 || strcmp(option, "") == 0)
   {
      LOG_INFO << "id=" << id() <<" charge = "<<charge()
         << endm;
      LOG_INFO << "pMom=" << pMom() << endm;
      LOG_INFO << "gpt =" << gPt() << " gEta = "<<gEta()<<" gPhi = "<<gPhi()<<endm;
      LOG_INFO << " nHitsFit = " << nHitsFit() << " nHitsdEdx = " << nHitsDedx() << endm;
      LOG_INFO << " nSigma E = " << nSigmaElectron() << endm;
      LOG_INFO << " beta = "<<beta()<<endm;
   }
}
