#include "StHadronTrack.h"
#include "StMessMgr.h"
#include "TVector2.h"
#include "TMath.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoDstMaker/StPicoEmcPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StBTofUtil/tofPathLength.hh"
#include "PhysicalConstants.h"
#include <climits>

ClassImp(StHadronTrack)

//----------------------------------------------------------------------------------
StHadronTrack::StHadronTrack() : mId(0), mDca(0), mPMom(0., 0., 0.),  mGMom(0., 0., 0.), 
   mNHitsFit(0), mNHitsDedx(0), mIsHftTrack(0), mNSigmaPion(32768), mNSigmaElectron(32768),   
    mBeta(0), mLocalY(-999)
{

}

/////////////////////////////////////////////////////////////////////////////////////////
// t - the global track.  p - the associated primary track from the first primary vertex
/////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------------
StHadronTrack::StHadronTrack(StPicoDst *picoDst, StPicoTrack* t)
   : mId(0), mDca(0), mPMom(0., 0., 0.),  mGMom(0., 0., 0.), 
   mNHitsFit(0), mNHitsDedx(0), mIsHftTrack(0), mNSigmaPion(32768), mNSigmaElectron(32768),   
    mBeta(0), mLocalY(-999)
{
      mId        = (UShort_t)t->id();
      mPMom      = t->pMom();
      int q      = t->charge();
      mNHitsFit  = t->nHitsFit()*q;
      mNHitsDedx = (UChar_t)(t->nHitsDedx());
      mNSigmaPion = (fabs(t->nSigmaPion() * 100.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaPion() * 100.));
      mNSigmaElectron = (fabs(t->nSigmaElectron() * 100.) > 32768) ? 32768 : (Short_t)(TMath::Nint(t->nSigmaElectron() * 100.));

	   StThreeVectorF vertexPos = picoDst->event()->primaryVertex();
      StPhysicalHelixD helix = t->helix();
      double thePath = helix.pathLength(vertexPos);
      StThreeVectorF dcaPos = helix.at(thePath);
      mDca = (dcaPos-vertexPos).mag()*1000;

      int index2TofPid = t->bTofPidTraitsIndex();
      if (index2TofPid>=0){
        StPicoBTofPidTraits *tofPid = picoDst->btofPidTraits(index2TofPid);
        //mTofMatchFlag = tofPid->btofMatchFlag();
        Float_t mom = mPMom.mag();
        mLocalY = tofPid->btofYLocal();
        //mLocalZ = tofPid->btofZLocal();
        Float_t beta = tofPid->btofBeta();
        if(beta<1e-4||beta>=(USHRT_MAX-1)/20000){
           Float_t tof = tofPid->btof();
           StThreeVectorF btofHitPos = tofPid->btofHitPos();
           float L = tofPathLength(&vertexPos, &btofHitPos, helix.curvature()); 
           beta = L/(tof*(c_light/1.0e9));
        }
        mBeta = (UShort_t)(beta*20000);
      }
}

//----------------------------------------------------------------------------------
StHadronTrack::~StHadronTrack()
{
   /* noop */
}
//----------------------------------------------------------------------------------
void StHadronTrack::Print(const Char_t *option) const
{
      LOG_INFO << "id=" << id() << endm;
      LOG_INFO << "pMom=" << pMom() << endm;
      LOG_INFO << " nHitsFit = " << nHitsFit() << " nHitsdEdx = " << nHitsDedx() << endm;
}


