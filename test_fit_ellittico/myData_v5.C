#define myData_cxx
#include "myData.h"
#include <TH2.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <vector>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <fstream>
#include <TGraph.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TF1.h>

#include "TROOT.h"
#include "TRandom.h"
static constexpr Int_t kMaxcl = 11414; 
static constexpr Int_t kMaxgr = 2482; 
static constexpr Int_t kMaxmt = 121; 
static constexpr Int_t kMaxim = 11414; 
static constexpr Int_t kMaxfr = 320; 
#include "TGraph.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TCanvas.h"
#include "TEllipse.h"

#include <cmath>
#include <iostream>


using namespace std;


const double pi=TMath::Pi();

//// SETTINGS
double x_pix_size=0; //um
double y_pix_size=0; //um
double len_view_x=0; //um
double len_view_y=0; //um
double thr_ldust_br =0;
double thr_ldust_area =0;
double fid_cut_par=0;
double maxcl=0;
const int npol=8;
int cut_nofit=1;
int cut_goodzone=1;
int cut_ncl=30;
double cut_minor=0.100;
int cut_isolated=-1;
int cut_npol=8;
int cut_mtrk=-1;
double cut_bar_l=0.03;
double cut_reg_u=2.3;
double cut_bar_u=1;
double cut_reg_l=0;
int cut_view=-1;
int channel=0;
//////////////

//Double_t xcorr[8] = {0.1813,0.0065,0.1731,0.3163,0.3770,0.3360,-0.1816,-1.148};
//Double_t ycorr[8] = {0.9882,0.5992,0.4274,0.2614,-0.2331,-0.6426,-0.6527,-0.4959};

Double_t major_mean[8] = {0.1894,0.1880,0.18621,0.1832,0.1806,0.1811,0.1849,0.1885};
Double_t phi_elliptic[8] = {0.264,0.610,0.942,1.310,1.753,2.216,2.627,2.923};
Double_t xy_gap[8] = {0.001634,0.000952,0.000263,-0.000476,-0.001146,-0.001490,-0.000592,0.0006743};

///// TREE VARIABLES /////////////
Bool_t eLargeDust=false;
Bool_t eGoodZone=false;
Bool_t eShadow=false;
Bool_t eClustTypeDust=false;
Int_t eChannel=0;
Int_t eGrainID=0;
Int_t ePolID=0;
Int_t eHeaderID=0;
Int_t eViewID=0;
Int_t eBfcPolID=0;
Int_t eBfcGap=0;
Int_t eNcl=0;
Int_t eNgr=0;
Int_t eNclFr=0;
Int_t ePuls=0;
Int_t eBfcID=0;
Int_t eNumIm=0;
Int_t eImCols=0;
Int_t eImRows=0;
Int_t eFlag=0;
Double_t eImArea=0;
Double_t eGrainPol=0;
Double_t eGrainx=0;
Double_t eGrainy=0;
Double_t ePreGrainx=0;
Double_t ePreGrainy=0;
Double_t eGrainz=0;
Double_t eClx=0;
Double_t eCly=0;
Double_t eClz=0;
Double_t eGrainrRms=0;
Double_t eGrainxRms=0;
Double_t eGrainyRms=0;
Double_t eGrainzRms=0;
Double_t eGrainMinRms=0;
Double_t eGrainMajRms=0;
Double_t eGrainPhiRms=0;
Double_t eGrainMin=0;
Double_t eGrainMaj=0;
Double_t eGrainEll=0;
Double_t eGrainPhi=0;
Double_t eGrainTheta=0;
Double_t eEllPrjX=0;
Double_t eEllPrjY=0;
Double_t eClustx=0;
Double_t eClusty=0;
Double_t eClustz=0;
Double_t eClustMin=0;
Double_t eClustMaj=0;
Double_t eClustEll=0;
Double_t eClustPhi=0;
Double_t eClustMaxPeak=0;
Double_t eClustMeanBkg=0;
Double_t eVolume=0;
Double_t eBfcVolume=0;
Double_t eFrBfVolume=0;
Double_t eArea=0;
Double_t eAreaRms=0;
Double_t eGrainArea=0;
Double_t eBfcArea=0;
Double_t eFrBfArea=0;
Double_t eBfcBorderFrame=0;
Double_t eXView=0;
Double_t eYView=0;
Double_t eZ=0;
Double_t eZlen=0;
Int_t eChain=0;
Int_t eSameFrame=0;
Int_t eNFrame=0;
Int_t eFrame=0;
Int_t eSetFrame=0;
Int_t eDeltaFrame=0;
Double_t eClDist=0;
Int_t eSetID=0;
Int_t eSetNCopy=0;
Int_t eSetNStatic=0;
Double_t eSetGapZ=0;
Double_t eSetGapZ2=0;
Double_t eSetGap=0;
Double_t eSetPath=0;
Double_t eSetMeanPath=0;
Double_t eSetRmsPath=0;
Double_t eSetMaxPath=0;
Double_t eSetMaxDist=0;
Double_t eSetNpeaks=0;
Double_t eSetNpeaksMax=0;
Double_t eSetNpeaksNcopy=0;
Double_t eSetNpeaksDist=0;
Double_t eSetNpeaksPhi=0;
Double_t eSetNpeaksDVol=0;
Double_t eSetNpeaksDNpx=0;
Double_t eSetNpeaksDBri=0;
Double_t eSetNpeaksMaxPhiAmp=0;
Double_t eSetMaxPol1=0;
Double_t eSetMaxPol2=0;
Double_t eSetMaxPixDist=0;
Double_t eSetMaxBar=0;
Double_t eSetBrAmp=0;
Double_t eSetBkgAmp=0;
Double_t eSetPeakAmp=0;
Double_t eSetMaxAmp=0;
Double_t eSetRmsAmp=0;
Double_t eSetMeanAmp=0;
Double_t eSetBrMaxPol=0;
Double_t eSetBrMinPol=0;
Double_t eSetXRms=0;
Double_t eSetYRms=0;
Double_t eSetXBar=0;
Double_t eSetYBar=0;
Double_t eSetXMaxBar=0;
Double_t eSetYMaxBar=0;
Double_t eSetXMinBar=0;
Double_t eSetYMinBar=0;
Double_t eSetPhiRms=0;
Double_t eSetPhiMean=0;
Double_t eSetChi2=0;
Double_t eMinDistGrain=0;
Double_t eSetPhi=0;
Double_t eSetPhiBar=0;
Double_t eSetPhiMaxAmp=0;
Double_t eSetTheBar=0;
Double_t eSetPhiFit=0;
Double_t eSetMeanBright=0;
Double_t eIsolated=0;
Double_t eSetNpxRms=0;
Double_t eSetVolRms=0;
Double_t eSetVolRatio=0;
Double_t eBfcMeanBkg=0;
Double_t eBfcSigPeak=0;
Double_t eBfcWeigth=0;
Double_t eSclMeanBkg=0;
Double_t eSclSigPeak=0;
Double_t eGrMeanBkg=0;
Double_t eGrSigPeak=0;
Double_t eGrSclMeanBkg=0;
Double_t eGrSclSigPeak=0;
Double_t eGrWeigth=0;
Double_t eDeltaPhi=0;
/////// Microtracks /////////
Int_t eMTrk=0;
Int_t eMTID=-1;
Int_t eMTGr=-1;
Int_t eMTNfr=-1;
Double_t eMTPhi=-10;
Double_t eMTThe=0;
Double_t eMTLen=-1;
Double_t eMTChi2=-1;
Double_t eMTVolDif=-1;
////////////////////////////////

Double_t eGrFitMin=0;
Double_t eGrFitMaj=0;
Double_t eGrFitEll=0;
Double_t eGrFitPhi=0;


//// OTHER DEFINITIONS /////////
Int_t viewID=0;
Int_t hID=0;
Int_t sum_vid=0;
Int_t sum_aid=0;
Double_t xlast=0.;
Double_t ylast=0.;
Double_t zlast=0.;

Double_t xmin=0.;
Double_t ymin=0.;
Double_t zmin=0.;

Double_t minor_good=0.;
Double_t ell_good=0.;

Double_t tmp_gr_dist=0;
Double_t tmp_gap=-1;
Int_t tmp_frame=-1;
Int_t n_ifr=1;
Int_t tmp_ifr=-1;
Int_t tmp_same_frame=-1;

Double_t xydist=0;
Double_t hdist=0;
Double_t tmp_rdist=-1;
Double_t tmp_rdist2=-1;
Double_t tmp_rdist_ldust=-1;
Double_t tmp_mtrk_dist=0;
//Double_t tmp_mtrk_dist[kMaxmt]={};
Double_t clust_tracks_dist=1000;
Double_t rdist_min=100;
Double_t rdist=0.;
Double_t phiang2=0.;
Double_t npeak_vol_dif=0.;
Double_t npeak_npx_dif=0.;
Double_t npeak_bri_dif=0.;
Double_t mtrk_rdist[4]={};
Double_t dist_mtrk=0.;
Double_t bfcl_rdist=0.;
Double_t bfcl_rdist2=0.;
Double_t bfcl_rdist3=0.;
Double_t bfcl_rdist_pol=0.;
Double_t bfcl_rdist_ldust=0.;
Double_t rdist3=0.;
Double_t xdist=0.;
Double_t ydist=0.;
Double_t zdist=0.;
Int_t NoClust_first=-1;
Int_t xbin=-1;
Int_t ybin=-1;
////////////////////////
const Int_t totXbin=80;
const Int_t totYbin=80;
Float_t new_phi_bar=0;
Int_t ngr_cell[totXbin][totYbin];
Float_t phi_projX[totXbin][totYbin];
Float_t phi_projY[totXbin][totYbin];
Float_t phi_projX_adi[totXbin][totYbin];
Float_t phi_projY_adi[totXbin][totYbin];
Float_t radius=0;

Double_t xgrain=0;
Double_t ygrain=0;
Double_t zgrain=0;

Int_t min_area_x=0;
Int_t min_area_y=0;
Int_t max_area_x=0;
Int_t max_area_y=0;

Double_t scale_result[4]={};
Double_t gr_scale_result[2]={};
Double_t im_info[2]={};
Double_t fit_info[3]={};

Int_t tmp_gr=-1;
Int_t index_pol=0;
Int_t index_2pol=0;

Double_t pre_step_x=0;
Double_t pre_step_y=0;
Double_t first_step_x=0;
Double_t first_step_y=0;
Int_t match=0;
Double_t scale_ref=0;

const int nsides=4;
Double_t area_cut=0;
Double_t area_scan=0;
Int_t iLarge_dust=0;
UInt_t minclust=6;
Long64_t index_clust=0;
Long64_t index_fr_clust=0;
Double_t fiducial_cut=0;


///////// ARRAYS ///////////////////
Double_t grain_Ox[kMaxgr]={};
Double_t grain_Oy[kMaxgr]={};
Double_t grain_rx[kMaxgr]={};
Double_t grain_ry[kMaxgr]={};
Double_t grain_area[kMaxgr]={};
Double_t bfcl_Area[kMaxgr]={};
Double_t bfcl_Vol[kMaxgr]={};
Double_t im_Cols[kMaxgr]={};
Double_t im_Rows[kMaxgr]={};
Double_t im_Area[kMaxgr]={};
Double_t x_gr_rms[kMaxgr]={};
Double_t y_gr_rms[kMaxgr]={};
Double_t z_gr_rms[kMaxgr]={};
Double_t min_gr_rms[kMaxgr]={};
Double_t maj_gr_rms[kMaxgr]={};
Double_t phi_gr_rms[kMaxgr]={};
Double_t area_gr_rms[kMaxgr]={};
Double_t xdust[kMaxgr]={};
Double_t ydust[kMaxgr]={};
Double_t zdust[kMaxgr]={};
Int_t gr_chain[kMaxgr]={};
Double_t dist_cl_gr[kMaxgr]={};
Int_t pol_id[kMaxgr]={};
Double_t set_chi2[kMaxgr]={};
Double_t min_dist_grain[kMaxgr]={};
   
bool grain[kMaxgr]={};
bool shadow[kMaxgr]={false};
bool ctd[kMaxgr]={false};
bool ldust[kMaxgr]={false};
bool bfc_border[kMaxgr]={false};
Int_t same_grain[kMaxgr]={};
Int_t gr_copy[kMaxgr]={};
Int_t gr_same_pos[kMaxgr]={};
Double_t gr_npeaks[kMaxgr]={};
Int_t gr_180[kMaxgr]={};
   
Int_t bfc_gap[kMaxgr]={};
Int_t gr_npeaks_npol[kMaxgr]={};
Double_t gr_min_maj[kMaxgr]={};
Float_t  gr_npeaks_max_phi[kMaxgr]={};
Double_t gr_mean_dist_bar[kMaxgr]={};
Double_t gr_mean_bkg_amp[kMaxgr]={};
Double_t gr_max_br_amp[kMaxgr]={};
Double_t gr_max_peak_amp[kMaxgr]={};
Double_t gr_max_br_pol[kMaxgr]={};
Double_t gr_min_br_pol[kMaxgr]={};
Double_t twopeak_phi[kMaxgr]={};
Double_t twopeak_dist[kMaxgr]={};
Double_t twopeak_dvol[kMaxgr]={};
Double_t twopeak_dnpx[kMaxgr]={};
Double_t twopeak_dbri[kMaxgr]={};
Double_t gr_vol_ref[kMaxgr]={};
Double_t gr_vol_mean[kMaxgr]={};
Double_t gr_vol_mean_bar[kMaxgr]={};
Double_t gr_vol_rms_bar[kMaxgr]={};
Double_t gr_npx_rms_bar[kMaxgr]={};
Double_t gr_vol_rms[kMaxgr]={};
Double_t gr_npx_rms[kMaxgr]={};
Double_t gr_phi_rms[kMaxgr]={};
Double_t gr_x_rms[kMaxgr]={};
Double_t gr_y_rms[kMaxgr]={};
Double_t gr_z_rms[kMaxgr]={};
Double_t gr_x_mean[kMaxgr]={};
Double_t gr_y_mean[kMaxgr]={};
Double_t gr_z_mean[kMaxgr]={};
Double_t gr_mean_path[kMaxgr]={};
Double_t gr_rms_path[kMaxgr]={};
Double_t gr_max_path[kMaxgr]={};
Double_t gr_phi_mean[kMaxgr]={};
Double_t gr_phi_max[kMaxgr]={};
Double_t gr_x_maxbar[kMaxgr]={};
Double_t gr_y_maxbar[kMaxgr]={};
Double_t gr_x_minbar[kMaxgr]={};
Double_t gr_y_minbar[kMaxgr]={};
Double_t gr_npx_sum[kMaxgr]={};
Double_t gr_vol_sum[kMaxgr]={};
Double_t gr_ncl_sum[kMaxgr]={};
Double_t gr_path[kMaxgr]={};
Double_t gr_gap[kMaxgr]={};
Double_t chi2_pol[kMaxgr]={};
Double_t phi_set[kMaxgr]={};
Double_t phi_set_bar[kMaxgr]={};
Double_t the_set_bar[kMaxgr]={};
Double_t phi_set_fit[kMaxgr]={};
Double_t gr_max_dist[kMaxgr]={};
Double_t gr_max_dist_bar[kMaxgr]={};
Double_t gr_maxpol1[kMaxgr]={};
Double_t gr_maxpol2[kMaxgr]={};
Double_t gr_isolated[kMaxgr]={};
Double_t gr_gap_z[kMaxgr]={};
Double_t gr_gap_z2[kMaxgr]={};
Double_t q_line[kMaxgr]={};
Double_t m_line[kMaxgr]={};
Double_t gr_max_amp[kMaxgr]={};
Double_t gr_rms_amp[kMaxgr]={};
Double_t gr_mean_amp[kMaxgr]={};
Double_t maxpixdist[kMaxgr]={};
Int_t set_first[kMaxgr]={};
Int_t gr_fr[kMaxgr]={};
Int_t set_fr[kMaxgr]={};
Double_t mt_dif_vol[kMaxgr]={};
Double_t phi_mt[kMaxgr]={};
int frbf_ent[kMaxgr][npol]={};
double cut_area[kMaxgr][nsides]={};   
double cut_area2[kMaxgr][nsides]={};   
double ld_area[kMaxgr][nsides]={};   
Bool_t mtrk[kMaxmt]={};

Double_t xb_frbf[kMaxgr][npol]={};
Double_t yb_frbf[kMaxgr][npol]={};
Double_t vol_frbf[kMaxgr][npol]={};
Double_t npx_frbf[kMaxgr][npol]={};
Double_t cl_x_pos[kMaxgr][npol]={};
Double_t cl_y_pos[kMaxgr][npol]={};
Double_t cl_z_pos[kMaxgr][npol]={};
Double_t cl_min_bf[kMaxgr][npol]={};
Double_t cl_maj_bf[kMaxgr][npol]={};
Double_t cl_ell_bf[kMaxgr][npol]={};
Double_t cl_phi_bf[kMaxgr][npol]={};
Double_t cl_sig_peak[kMaxgr][npol]={};
Double_t cl_bkg_mean[kMaxgr][npol]={};
Int_t same_frame[kMaxgr][npol]={};
Int_t num_peaks[kMaxgr][npol]={};
Int_t num_tot_peaks[kMaxgr]={};  
Int_t ipol_gr[kMaxgr][npol]={};
/////////////////////////


DMRRun         *aRun;
DMRView        *view;
DMRImageCl     *image;

myData vv;

//
// "NUMERICALLY STABLE DIRECT LEAST SQUARES FITTING OF ELLIPSES"
// Radim Halir, Jan Flusser
// http://autotrace.sourceforge.net/WSCG98.pdf
//
// http://en.wikipedia.org/wiki/Ellipse
//
// An "algebraic distance" of a point "(x, y)" to a given conic:
//   F(x, y) = A * (x - X0)^2 + B * (x - X0) * (y - Y0) + C * (y - Y0)^2
//           + D * (x - X0) + E * (y - Y0) + F
//
// Ellipse-specific constraints:
//   F(x, y) = 0
//   B^2 - 4 * A * C < 0
//
// input parameter is a a pointer to a "TGraph" with at least 6 points
//
// returns a "TVectorD" ("empty" in case any problems encountered):
// ellipse[0] = "X0"
// ellipse[1] = "Y0"
// ellipse[2] = "A"
// ellipse[3] = "B"
// ellipse[4] = "C"
// ellipse[5] = "D"
// ellipse[6] = "E"
// ellipse[7] = "F"
//
TVectorD fit_ellipse(TGraph *g)
{
  TVectorD ellipse;

  if (!g) return ellipse; // just a precaution
  if (g->GetN() < 6) return ellipse; // just a precaution

  Int_t i;
  Double_t tmp;

  Int_t N = g->GetN();
  Double_t xmin, xmax, ymin, ymax, X0, Y0;
  g->ComputeRange(xmin, ymin, xmax, ymax);
#if 1 /* 0 or 1 */
  X0 = (xmax + xmin) / 2.0;
  Y0 = (ymax + ymin) / 2.0;
#else /* 0 or 1 */
  X0 = Y0 = 0.0;
#endif /* 0 or 1 */

  TMatrixD D1(N, 3); // quadratic part of the design matrix
  TMatrixD D2(N, 3); // linear part of the design matrix

  for (i = 0; i < N; i++) {
    Double_t x = (g->GetX())[i] - X0;
    Double_t y = (g->GetY())[i] - Y0;
    D1[i][0] = x * x;
    D1[i][1] = x * y;
    D1[i][2] = y * y;
    D2[i][0] = x;
    D2[i][1] = y;
    D2[i][2] = 1.0;
  }

  // quadratic part of the scatter matrix
  TMatrixD S1(TMatrixD::kAtA, D1);
  // combined part of the scatter matrix
  TMatrixD S2(D1, TMatrixD::kTransposeMult, D2);
  // linear part of the scatter matrix
  TMatrixD S3(TMatrixD::kAtA, D2);
  S3.Invert(&tmp); S3 *= -1.0;
  if (tmp == 0.0) {
    std::cout << "fit_ellipse : linear part of the scatter matrix is singular!" << std::endl;
    return ellipse;
  }
  // for getting a2 from a1
  TMatrixD T(S3, TMatrixD::kMultTranspose, S2);
  // reduced scatter matrix
  TMatrixD M(S2, TMatrixD::kMult, T); M += S1;
  // premultiply by inv(C1)
  for (i = 0; i < 3; i++) {
    tmp = M[0][i] / 2.0;
    M[0][i] = M[2][i] / 2.0;
    M[2][i] = tmp;
    M[1][i] *= -1.0;
  }
  // solve eigensystem
  TMatrixDEigen eig(M); // note: eigenvectors are not normalized
  const TMatrixD &evec = eig.GetEigenVectors();
  // const TVectorD &eval = eig.GetEigenValuesRe();
  if ((eig.GetEigenValuesIm()).Norm2Sqr() != 0.0) {
    std::cout << "fit_ellipse : eigenvalues have nonzero imaginary parts!" << std::endl;
    return ellipse;
  }
  // evaluate a’Ca (in order to find the eigenvector for min. pos. eigenvalue)
  for (i = 0; i < 3; i++) {
    tmp = 4.0 * evec[0][i] * evec[2][i] - evec[1][i] * evec[1][i];
    if (tmp > 0.0) break;
  }
  if (i > 2) {
    std::cout << "fit_ellipse : no min. pos. eigenvalue found!" << std::endl;
    // i = 2;
    return ellipse;
  }
  // eigenvector for min. pos. eigenvalue
  TVectorD a1(TMatrixDColumn_const(evec, i));
  tmp = a1.Norm2Sqr();
  if (tmp > 0.0) {
    a1 *= 1.0 / std::sqrt(tmp); // normalize this eigenvector
  } else {
    std::cout << "fit_ellipse : eigenvector for min. pos. eigenvalue is NULL!" << std::endl;
    return ellipse;
  }
  TVectorD a2(T * a1);

  // ellipse coefficients
  ellipse.ResizeTo(8);
  ellipse[0] = X0; // "X0"
  ellipse[1] = Y0; // "Y0"
  ellipse[2] = a1[0]; // "A"
  ellipse[3] = a1[1]; // "B"
  ellipse[4] = a1[2]; // "C"
  ellipse[5] = a2[0]; // "D"
  ellipse[6] = a2[1]; // "E"
  ellipse[7] = a2[2]; // "F"

  return ellipse;
}

//
// http://mathworld.wolfram.com/Ellipse.html
// http://mathworld.wolfram.com/QuadraticCurve.html
// http://mathworld.wolfram.com/ConicSection.html
//
// "Using the Ellipse to Fit and Enclose Data Points"
// Charles F. Van Loan
// http://www.cs.cornell.edu/cv/OtherPdf/Ellipse.pdf
//
// input parameter is a reference to a "TVectorD" which describes
// an ellipse according to the equation:
//   0 = A * (x - X0)^2 + B * (x - X0) * (y - Y0) + C * (y - Y0)^2
//     + D * (x - X0) + E * (y - Y0) + F
// conic[0] = "X0"
// conic[1] = "Y0"
// conic[2] = "A"
// conic[3] = "B"
// conic[4] = "C"
// conic[5] = "D"
// conic[6] = "E"
// conic[7] = "F"
//
// returns a "TVectorD" ("empty" in case any problems encountered):
// ellipse[0] = ellipse's "x" center ("x0")
// ellipse[1] = ellipse's "y" center ("y0")
// ellipse[2] = ellipse's "semimajor" axis along "x" ("a" > 0)
// ellipse[3] = ellipse's "semiminor" axis along "y" ("b" > 0)
// ellipse[4] = ellipse's axes rotation angle ("theta" = -45 ... 135 degrees)
//
TVectorD ConicToParametric(const TVectorD &conic)
{
  TVectorD ellipse;

  if (conic.GetNrows() != 8) {
    std::cout << "ConicToParametric : improper input vector length!" << std::endl;
    return ellipse;
  }

  Double_t a, b, theta;
  Double_t x0 = conic[0]; // = X0
  Double_t y0 = conic[1]; // = Y0

  // http://mathworld.wolfram.com/Ellipse.html
  Double_t A = conic[2];
  Double_t B = conic[3] / 2.0;
  Double_t C = conic[4];
  Double_t D = conic[5] / 2.0;
  Double_t F = conic[6] / 2.0;
  Double_t G = conic[7];

  Double_t J = B * B - A * C;
  Double_t Delta = A * F * F + C * D * D + J * G - 2.0 * B * D * F;
  Double_t I = - (A + C);

  // http://mathworld.wolfram.com/QuadraticCurve.html
  if (!( (Delta != 0.0) && (J < 0.0) && (I != 0.0) && (Delta / I < 0.0) )) {
    std::cout << "ConicToParametric : ellipse (real) specific constraints not met!" << std::endl;
    return ellipse;
  }

  x0 += (C * D - B * F) / J;
  y0 += (A * F - B * D) / J;

  Double_t tmp = std::sqrt((A - C) * (A - C) + 4.0 * B * B);
  a = std::sqrt(2.0 * Delta / J / (I + tmp));
  b = std::sqrt(2.0 * Delta / J / (I - tmp));

  theta = 0.0;
  if (B != 0.0) {
    tmp = (A - C) / 2.0 / B;
    theta = -45.0 * (std::atan(tmp) / TMath::PiOver2());
    if (tmp < 0.0) { theta -= 45.0; } else { theta += 45.0; }
    if (A > C) theta += 90.0;
  } else if (A > C) theta = 90.0;

  // try to keep "a" > "b"
  if (a < b) { tmp = a; a = b; b = tmp; theta -= 90.0; }
  // try to keep "theta" = -45 ... 135 degrees
  if (theta < -45.0) theta += 180.0;
  if (theta > 135.0) theta -= 180.0;

  // ellipse coefficients
  ellipse.ResizeTo(5);
  ellipse[0] = x0; // ellipse's "x" center
  ellipse[1] = y0; // ellipse's "y" center
  ellipse[2] = a; // ellipse's "semimajor" axis along "x"
  ellipse[3] = b; // ellipse's "semiminor" axis along "y"
  ellipse[4] = theta; // ellipse's axes rotation angle (in degrees)

  return ellipse;
}

//
// creates a test TGraph with an ellipse
//
TGraph *TestGraphDLSF(Bool_t randomize = kFALSE) {
  Int_t i;

  // define the test ellipse
  Double_t x0 = 4; // ellipse's "x" center
  Double_t y0 = 3; // ellipse's "y" center
  Double_t a = 2; // ellipse's "semimajor" axis along "x" (> 0)
  Double_t b = 1; // ellipse's "semiminor" axis along "y" (> 0)
  Double_t theta = 100; // ellipse's axes rotation angle (-45 ... 135 degrees)

  // gRandom->SetSeed(0);
  if (randomize) {
    x0 = 10.0 - 20.0 * gRandom->Rndm();
    y0 = 10.0 - 20.0 * gRandom->Rndm();
    a = 0.5 + 4.5 * gRandom->Rndm();
    b = 0.5 + 4.5 * gRandom->Rndm();
    theta = 180.0 - 360.0 * gRandom->Rndm();
  }

  const Int_t n = 100; // number of points
  Double_t x[n], y[n];
  Double_t dt = TMath::TwoPi() / Double_t(n);
  Double_t tmp;
  theta *= TMath::PiOver2() / 90.0; // degrees -> radians
  for (i = 0; i < n; i++) {
    x[i] = a * (std::cos(dt * Double_t(i)) + 0.1 * gRandom->Rndm() - 0.05);
    y[i] = b * (std::sin(dt * Double_t(i)) + 0.1 * gRandom->Rndm() - 0.05);
    // rotate the axes
    tmp = x[i];
    x[i] = x[i] * std::cos(theta) - y[i] * std::sin(theta);
    y[i] = y[i] * std::cos(theta) + tmp * std::sin(theta);
    // shift the center
    x[i] += x0;
    y[i] += y0;
  }

  // create the test TGraph
  TGraph *g = ((TGraph *)(gROOT->FindObject("g")));
  if (g) delete g;
  g = new TGraph(n, x, y);
  g->SetNameTitle("g", "test ellipse");

  return g;
}

//
// "ROOT Script" entry point (the same name as the "filename's base")
//
void fitEllipseTGraphDLSF(TGraph *g = ((TGraph *)0), double *fit_info=0)
{
  if (!g) g = TestGraphDLSF(kTRUE); // create a "random" ellipse

  // fit the TGraph
  TVectorD conic = fit_ellipse(g);
  TVectorD ellipse = ConicToParametric(conic);

#if 1 /* 0 or 1 */
  if ( ellipse.GetNrows() == 5 ) {
    std::cout << std::endl;
    std::cout << "x0 = " << ellipse[0] << std::endl;
    std::cout << "y0 = " << ellipse[1] << std::endl;
    std::cout << "a = " << ellipse[2] << std::endl;
    std::cout << "b = " << ellipse[3] << std::endl;
    std::cout << "theta = " << ellipse[4]*TMath::Pi()/180. << std::endl;
    std::cout << std::endl;

    fit_info[0]=ellipse[3]*2*27.52;
    fit_info[1]=ellipse[2]*2*27.52;
    fit_info[2]=ellipse[4]*TMath::Pi()/180.;
  }
#endif /* 0 or 1 */

#if 1 /* 0 or 1 */
  // draw everything
  TCanvas * c = ((TCanvas *)(gROOT->GetListOfCanvases()->FindObject("c")));
  if (c) { c->Clear(); } else { c = new TCanvas("c", "c"); }
  c->SetGrid(1, 1);
  g->Draw("A*");
  g->GetXaxis()->SetLimits(0,41);
  g->GetYaxis()->SetRangeUser(0,41);
  if ( ellipse.GetNrows() == 5 ) {
    TEllipse *e = new TEllipse(ellipse[0], ellipse[1], // "x0", "y0"
                               ellipse[2], ellipse[3], // "a", "b"
                               0, 360,
                               ellipse[4]); // "theta" (in degrees)
    e->SetFillStyle(0); // hollow
    e->SetLineColor(2);
    e->SetLineWidth(2);
    e->Draw();
  }
  c->Modified(); c->Update(); // make sure it's really drawn
#endif /* 0 or 1 */
  c->SaveAs("fitEll.root");
  return;
}

// end of file fitEllipseTGraphDLSF.cxx by Silesius Anonymus


TH2F *merged_histo(DMRRun *aRun, DMRView *view, int ihd, int igr, int *clset)
{

  const int dr=50;
  float pixX = aRun->GetHeader()->pixX;
  float pixY = aRun->GetHeader()->pixY;
  int nx     = aRun->GetHeader()->npixX;
  int ny     = aRun->GetHeader()->npixY;
  int x0,y0;

  int tmp_max=0;
  int tmp_min=255;

  int tmp_max2=-255;
  int tmp_min2=255;

  
  const int ncp=8;


  TH2F *hsum = new TH2F();
  hsum=0;
  TH2F *h[ncp];
  DMRImage imin; imin.SetImage( 2*dr, 2*dr );
  DMRImage imax; imax.SetImage( 2*dr, 2*dr );

  for(int i=0;i<ncp;i++){
    h[i]= new TH2F("","",2*dr,0,dr,2*dr,0,dr);
    if(clset[i]!=-1){
      DMRImageCl *im    = aRun->GetCLIM(ihd,clset[i],i,dr);
      if(im){
	h[i] = im->GetHist2();
	if(h[i]->GetMaximum()>tmp_max)tmp_max = h[i]->GetMaximum();
	if(h[i]->GetMinimum()<tmp_min)tmp_min = h[i]->GetMinimum();
      }
      delete im;
    }
  }
  
  for(int i=0;i<ncp;i++){

    if(clset[i]!=-1){
      
      DMRViewHeader   *hd = view->GetHD();
      DMRCluster      *cl = view->GetCL(clset[i]);
      DMRFrame        *frcl = view->GetFR(cl->ifr);       
      DMRFrame        *fr  = view->GetFR(frcl->iz,frcl->ipol);
      
      if(i==0){
	x0 = ((cl->x+hd->x)-fr->x)/pixX + nx/2;
	y0 = ((cl->y+hd->y)-fr->y)/pixY + ny/2;
      }
    
      DMRImageCl *im    = aRun->GetCLIMBFC(ihd,clset[i],i,dr,x0,y0);
      if(im){
	h[i] = im->GetHist2();
	//cout << i << " " << h[i]->GetMaximum() << endl;
	if(!hsum){
	  hsum = (TH2F*)(h[i]->Clone("hsum"));
	  hsum->SetTitle(Form("Sum over all polarizations"));
	  hsum->SetName(Form("hsum"));
	  imin.Max(*im);
	  imax.Max(*im);
	}
	else {
	  hsum->Add(h[i]);
	  imin.Min(*im);
	  imax.Max(*im);
	  //cout << "sum "<< hsum->GetMaximum() << endl;
	}
      }
      delete im;
      //delete hd;
      //delete cl;
      //delete frcl;
      //delete fr;
    }
    delete h[i];
  }
  hsum->Scale(1./8.);
  //delete view;
  return hsum;
  delete hsum;
}



void gr_spectrum(TH2F *h, double *gr_scale_result)
{
  int nbinsX;
  int nbinsY; 
  double cont=0;
  Double_t bkg_mean=0;
  Double_t sig_peak=0;

  TH1F *hspect = new TH1F("hspect","hs",255,0,255);

  nbinsX = h->GetNbinsX();
  nbinsY = h->GetNbinsY();
    
  for(int i=0;i<nbinsX;i++){
    for(int j=0;j<nbinsY;j++){
      cont = h->GetBinContent(i+1,j+1);
      if(cont>0)hspect->Fill(cont);
    }
  }
  
  if(hspect->GetEntries()!=0){
    hspect->Fit("gaus","Q");
    TF1 *gaus = hspect->GetFunction("gaus");
    bkg_mean = gaus->GetParameter(1);
    sig_peak = h->GetMaximum();
  }
  else {
    bkg_mean=-1;
    sig_peak=-1;
  }
  
  gr_scale_result[0]=bkg_mean;
  gr_scale_result[1]=sig_peak;
  
  //delete h;
  delete hspect;
}




void scaling(DMRRun *aRun, int ihd, int icl, int ipol, double *scale_result)
{

  double br[6] = {20,35,50,65,80,95};
  double r0[6] = {1.076,1.157,1.185,1.184,1.215,1.229};
  double r1[6] = {1.043,1.103,1.084,1.086,1.106,1.083};
  double r2[6] = {0.996,1.013,0.962,0.947,0.973,0.973};
  double r3[6] = {0.953,0.884,0.862,0.878,0.878,0.878};
  double r4[6] = {0.947,0.886,0.869,0.872,0.889,0.889};
  double r5[6] = {0.977,0.980,0.947,0.963,0.961,0.961};
  double r6[6] = {1.024,1.077,1.075,1.082,1.109,1.082};
  double r7[6] = {1.064,1.140,1.150,1.174,1.206,1.237};

  double brM[8][6];
 
  for(int j=0;j<6;j++){
    brM[0][j]=r0[j];
    brM[1][j]=r1[j];
    brM[2][j]=r2[j];
    brM[3][j]=r3[j];
    brM[4][j]=r4[j];
    brM[5][j]=r5[j];
    brM[6][j]=r6[j];
    brM[7][j]=r7[j];
  }

  
  int nbinsX;
  int nbinsY; 
  double cont=0;
  Double_t bkg_mean=0;
  Double_t sig_peak=0;
  double cont_S=0;
  Double_t bkg_mean_S=0;
  Double_t sig_peak_S=0;
  //float w[8] = {1,0.957,0.909,0.881,0.887,0.918,0.962,0.996};
  int dr=50;
  TH2F *h=0;
  TH2F *hs=0;
  TH1F *hspect = new TH1F("hspect","hs",255,0,255);
  TH1F *hspect_S = new TH1F("hspect_S","hs_S",255,0,255);
  hspect->Reset(0);
  hspect_S->Reset(0);
  DMRImageCl *im    = aRun->GetCLIM(ihd,icl,ipol,dr);
  
  if(im){
    h = im->GetHist2();
    hs = im->GetHist2();    

    for(int k=0;k<6;k++){
      if(k==0){
	if(hs->GetMaximum()<br[1])hs->Scale(1./brM[ipol][0]);
      }
      if(k>0 && k<5){
	if(hs->GetMaximum()>=br[k] && hs->GetMaximum()<br[k+1])hs->Scale(1./brM[ipol][k]);
      }
      if(k==5){
	if(hs->GetMaximum()>=br[5])hs->Scale(1./brM[ipol][5]);
      }
    }

    nbinsX = h->GetNbinsX();
    nbinsY = h->GetNbinsY();
    
    //if(icl==168) cout <<"hello1 "<< scale_result[0] << " " << scale_result[1] << " " << scale_result[2] << " " << scale_result[3] << endl;
    for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	cont = h->GetBinContent(i+1,j+1);
	if(cont>0)hspect->Fill(cont);
	cont_S = hs->GetBinContent(i+1,j+1);
	if(cont_S>0)hspect_S->Fill(cont_S);
      }
    }
    
    if(hspect->GetEntries()!=0){
      hspect->Fit("gaus","Q");
      TF1 *gaus = hspect->GetFunction("gaus");
      bkg_mean = gaus->GetParameter(1);
      sig_peak = h->GetMaximum();
    }
    else {
      bkg_mean=-1;
      sig_peak=-1;
    }
    if(hspect_S->GetEntries()!=0){
      hspect_S->Fit("gaus","Q");
      TF1 *gaus_S = hspect_S->GetFunction("gaus");
      bkg_mean_S = gaus_S->GetParameter(1);
      sig_peak_S = hs->GetMaximum();
    }
    else {
      bkg_mean_S=-1;
      sig_peak_S=-1;
    }
  }
  else{
    bkg_mean=-1;
    sig_peak=-1;
    bkg_mean_S=-1;
    sig_peak_S=-1;
  }

  scale_result[0]=bkg_mean;
  scale_result[1]=sig_peak;
  scale_result[2]=bkg_mean_S;
  scale_result[3]=sig_peak_S;

    
  
  delete im;  
  delete h;
  delete hspect;
  delete hs;
  delete hspect_S;
  
}

float maxpeak(DMRRun *aRun, int ihd, int icl, int ipol)
{
  int dr=15;
  TH2F *h=0;
  //cout << iv << " " << icl << " " << ipol << endl;
  DMRImageCl *im    = aRun->GetCLIM(ihd,icl,ipol,dr);
  if(im){
    h = im->GetHist2();
    float maxpeak = h->GetMaximum();
    return maxpeak;
  }
  else return -1;
}



void images_info(DMRRun *aRun, int ihd, int icl, int ipol, double *im_info)
{
  //cout << ihd << " " << icl <<" " << ipol << endl;
  int nbinsX=0;
  int nbinsY=0;
  double cont=0;
  Double_t bkg_mean=0;
  Double_t sig_peak=0;
  
  //TCanvas *c1 = new TCanvas();
  TH1F *hspect = new TH1F("hspect","hs",255,0,255);
  TF1 *gaus = new TF1("fit_gaus","gaus");
  gaus->SetParameter(1,23);
  gaus->SetParameter(2,1);

  int dr=15;
  TH2F *h = new TH2F();
  //cout << iv << " " << icl << " " << ipol << endl;
  DMRImageCl *im    = new DMRImageCl();
  im = aRun->GetCLIM(ihd,icl,ipol,dr);
  if(im){
    h = im->GetHist2();
    nbinsX = h->GetNbinsX();
    nbinsY = h->GetNbinsY();
    
    for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	cont = h->GetBinContent(i+1,j+1);
	if(cont>0)hspect->Fill(cont);
      }
    }
    
    
  if(hspect->GetEntries()!=0){
    double max_pos = hspect->GetBinCenter(hspect->GetMaximumBin());
    hspect->Fit(gaus,"Q","",max_pos-5,max_pos+5);
    //TF1 *gaus = hspect->GetFunction("gaus");
    bkg_mean = gaus->GetParameter(1);
    sig_peak = h->GetMaximum();
  }
  else {
    bkg_mean=-1;
    sig_peak=-1;
  }

  }
  else{
    bkg_mean=-1;
    sig_peak=-1;
  }

  im_info[0]=bkg_mean;
  im_info[1]=sig_peak;
  delete hspect;
  delete im;
  delete gaus;
  delete h;
  //delete c1;
}




int otsu_method(float *histogram, long int total_pixels) {
    double probability[256], mean[256];
    double max_between, between[256];
    int threshold;

    /*
    probability = class probability
    mean = class mean
    between = between class variance
    */

    for(int i = 0; i < 256; i++) {
        probability[i] = 0.0;
        mean[i] = 0.0;
        between[i] = 0.0;
    }

    probability[0] = histogram[0];

    for(int i = 1; i < 256; i++) {
        probability[i] = probability[i - 1] + histogram[i];
        mean[i] = mean[i - 1] + i * histogram[i];
    }

    threshold = 0;
    max_between = 0.0;

    for(int i = 0; i < 255; i++) {
        if(probability[i] != 0.0 && probability[i] != 1.0)
            between[i] = pow(mean[255] * probability[i] - mean[i], 2) / (probability[i] * (1.0 - probability[i]));
    else
        between[i] = 0.0;
        if(between[i] > max_between) {
            max_between = between[i];
            threshold = i;
        }
    }

    return threshold;
}



void images_ellipse(DMRRun *aRun, int ihd, int igr, int ipol, double *fit_info, TCanvas *c1)
{
  gStyle->SetOptFit(1111);
  //cout << ihd << " " << icl <<" " << ipol << endl;
  int nbinsX=0;
  int nbinsY=0;
  int cont=0;
  Double_t bkg_mean=0;
  Double_t bkg_sigma=0;
  Double_t sig_peak=0;
  int locmax,locmay,locmaz;
  int total_pixels;
  float histogram[256]={};
  float occurrence[256]={};
  float threshold_value=0;
  float pixX=0.02752;
  float pixY=-0.02752;
  double fit_result[3]={};

  //TF2 *bigaus = new TF2("fit_bigaus","[0]*x*x+2*[1]*x*y+2*[2]*x+2*[3]*y+[4]");
  //TF2 *bigaus = new TF2("fit_bigaus","bigaus+[6]");
  TF2 *bigaus = new TF2("fit_bigaus","bigaus");
  TF2 *biconst = new TF2("fit_biconst","[6]");
  TH1F *hspect = new TH1F("hspect","hs",255,0,255);
  TF1 *gaus = new TF1("fit_gaus","gaus");
  gaus->SetParameter(1,23);
  gaus->SetParameter(2,1);

  int dr=20;
  TH2F *h = new TH2F();
  TH2F *hh = new TH2F();
  TH2F *hhh = new TH2F();
  TGraph *gr = new TGraph();
  //cout << iv << " " << icl << " " << ipol << endl;
  DMRImageCl *im    = new DMRImageCl();
  im = aRun->GetGRIM(ihd,igr,ipol,dr);
  if(im){
    h = im->GetHist2();
    hhh = im->GetHist2();
    nbinsX = h->GetNbinsX();
    nbinsY = h->GetNbinsY();
    //total_pixels = (nbinsX-1)*(nbinsY-1); //-1 per la riga/colonna vuota
    
    for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	cont = h->GetBinContent(i+1,j+1);
	if(cont>0){
	  hspect->Fill(cont);
	  total_pixels++;
	  //occurrence[cont] = occurrence[cont] + 1;
	}
      }
    }

    
    
  if(hspect->GetEntries()!=0){
    double max_pos = hspect->GetBinCenter(hspect->GetMaximumBin());
    hspect->Fit(gaus,"Q","",max_pos-5,max_pos+5);
    //TF1 *gaus = hspect->GetFunction("gaus");
    bkg_mean = gaus->GetParameter(1);
    bkg_sigma =  gaus->GetParameter(2);
    sig_peak = h->GetMaximum();
  }
  else {
    bkg_mean=-1;
    sig_peak=-1;
  }
  
  }
  else{
    bkg_mean=-1;
    sig_peak=-1;
  }

  if(sig_peak!=-1){
    h->GetMaximumBin(locmax,locmay,locmaz);
    cout << locmax << " " << locmay << endl; 
  for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	cont = h->GetBinContent(i+1,j+1);
	if(cont>0){
	  //hspect->Fill(cont);
	  occurrence[cont] = occurrence[cont] + 1;
	}
      }
    }

  for(int i = 0; i <= 255; i++) {
    /* TAKES NUMBER OF OCCURRENCES OF A PARTICULAR PIXEL 
     * AND DIVIDES BY THE TOTAL NUMBER OF PIXELS YIELDING 
     * A RATIO */
    histogram[i] = (float) occurrence[i] / (float) total_pixels;
  }

  threshold_value = otsu_method(histogram, total_pixels);

   for(int i=0;i<nbinsX;i++){
      for(int j=0;j<nbinsY;j++){
	//if(h->GetBinContent(i+1,j+1)>threshold_value)h->SetBinContent(i+1,j+1,sig_peak);
        //if(h->GetBinContent(i+1,j+1)<threshold_value)h->SetBinContent(i+1,j+1,0);
	if(h->GetBinContent(i+1,j+1)<(bkg_mean+3*bkg_sigma))h->SetBinContent(i+1,j+1,0);
	//else h->SetBinContent(i+1,j+1,255);
      }
    }
   hh = (TH2F*)h->Clone();
   int index=0;
   for(int i=1;i<nbinsX-1;i++){
      for(int j=1;j<nbinsY-1;j++){
	//if(h->GetBinContent(i+1,j+1)>threshold_value)h->SetBinContent(i+1,j+1,sig_peak);
        if(h->GetBinContent(i,j)!=0 && h->GetBinContent(i,j+1)!=0 && h->GetBinContent(i,j+2)!=0 && h->GetBinContent(i+1,j)!=0 &&  h->GetBinContent(i+1,j+2)!=0 &&
	    h->GetBinContent(i+2,j)!=0 &&  h->GetBinContent(i+2,j+1)!=0 &&  h->GetBinContent(i+2,j+2)!=0)hh->SetBinContent(i+1,j+1,0);
      }
   }
   
   for(int i=1;i<nbinsX-1;i++){
     for(int j=1;j<nbinsY-1;j++){
       //if(h->GetBinContent(i,j)==0 && h->GetBinContent(i,j+1)==0 && h->GetBinContent(i,j+2)==0 && h->GetBinContent(i+1,j)==0 &&  h->GetBinContent(i+1,j+2)==0 &&
       //  h->GetBinContent(i+2,j)==0 &&  h->GetBinContent(i+2,j+1)==0 &&  h->GetBinContent(i+2,j+2)==0)hh->SetBinContent(i+1,j+1,0);	
       if(hh->GetBinContent(i+1,j+1)!=0){
	 hh->SetBinContent(i+1,j+1,255);
	 gr->SetPoint(index,hh->GetYaxis()->GetBinCenter(i+1),hh->GetYaxis()->GetBinCenter(j+1));
	 index++;
       }
     }
   }
   /*
   double xc = h->GetXaxis()->GetBinCenter(locmax); 
   double yc = h->GetYaxis()->GetBinCenter(locmay);
   cout << xc << " " << yc << endl;
   bigaus->SetParameter(1,locmax);
   bigaus->SetParameter(3,locmay);
   bigaus->SetParameter(2,6);
   bigaus->SetParameter(4,6);
   bigaus->SetParameter(6,18);
   //h->Fit(bigaus,"BQ");
   fit_info[2] = TMath::ATan((2*bigaus->GetParameter(5))/(bigaus->GetParameter(2)-bigaus->GetParameter(4)))/2.;
   if(bigaus->GetParameter(2)<bigaus->GetParameter(4)){
     fit_info[0] = bigaus->GetParameter(2);
     fit_info[1] = bigaus->GetParameter(4);
   }
   else {
     fit_info2[0] = bigaus->GetParameter(4);
     fit_info[1] = bigaus->GetParameter(2);
    }
   */
   

   //hh->Fit(bigaus,"BQ");
}

  fitEllipseTGraphDLSF(gr,fit_info);
  TFile *fit_file = TFile::Open("fitEll.root");
  TCanvas *c = (TCanvas *)fit_file->Get("c");
  /*TPad *pad = (TPad*)gPad;
  TObject *obj;
  TIter next(c->GetListOfPrimitives());
  while ((obj=next())) {
     gROOT->SetSelectedPad(pad);
     pad->GetListOfPrimitives()->Add(obj->Clone());
  }
  gPad->Modified();*/
  
  //fit_info[0]=fit_result[0];
  //fit_info[1]=fit_result[1];
  //fit_info[2]=fit_result[2];
  //cout << fit_info[0] << " " << fit_info[1] << " " << fit_info[2] << " " << threshold_value << endl;
  c1 = new TCanvas(Form("vid_%d_gr_%d",ihd,igr));
  c1->Divide(2,2);
  c1->cd(1);
  hhh->Draw("colz");
  c1->cd(2);
  h->Draw("colz");
  c1->cd(3);
  hh->Draw("colz");
  c1->cd(4);
  c->DrawClonePad();
  c1->Modified();
  c1->Update();
  c1->SaveAs(Form("ellipse_vid_%d_gr_%d.root",ihd,igr));
  delete bigaus;
  delete biconst;
  delete hspect;
  delete im;
  delete gaus;
  delete h;
  delete hh;
  //delete c1;
}



void myDatacard(char* datacard){
  ifstream myfile(datacard);
  vector<vector<string> > runDescription;
  if (myfile.is_open()){
    string line;
    // cout << line << endl;
    while (getline (myfile,line)){
      istringstream iss(line);
      string a;
      vector<string> tmpS;
      while(!iss.eof()){     
	iss>>a;
	tmpS.push_back(a);
	//cout << a << " ";
      }
      runDescription.push_back(tmpS);
      tmpS.clear();
    }
    myfile.close();
    
  }
  else cout << "Unable to open file: "<< datacard << endl; 

  cout << "SETTINGS from " << datacard << endl;
  for(int i = 0; i< runDescription.size();i++){
    int nstring  = runDescription[i].size();
    // cout << endl;
    for(int j = 0; j<nstring; j++){
      //cout << runDescription[i][j]<<endl;
      if(runDescription[i][j] == "pix_size_x") {
	x_pix_size = atof( runDescription[i][2].c_str());
	cout << "pixel size x [um] is " << x_pix_size << endl;
      }
      if(runDescription[i][j] == "pix_size_y") {
	y_pix_size = atof( runDescription[i][2].c_str());
	cout << "pixel size y [um] is "<< y_pix_size << endl;
      }
      if(runDescription[i][j] == "len_view_x") {
	len_view_x = atof( runDescription[i][2].c_str());
	cout << "view size x [um] is " << len_view_x << endl;
      }
      if(runDescription[i][j] == "len_view_y") {
	len_view_y = atof( runDescription[i][2].c_str());
	cout << "view size y [um] is " << len_view_y << endl;
      }
      if(runDescription[i][j] == "thr_ldust_br") {
	thr_ldust_br = atof( runDescription[i][2].c_str());
	cout << "threshold brigthness for large dust  (log10 scale) is " << thr_ldust_br << endl;
      }
      if(runDescription[i][j] == "thr_ldust_area") {
	thr_ldust_area = atof( runDescription[i][2].c_str());
	cout << "threshold number of pixel for large dust (log 10 scale) is " << thr_ldust_area << endl;
      }
      if(runDescription[i][j] == "fid_cut_par") {
	fid_cut_par = atof( runDescription[i][2].c_str());
	cout << "fiducial cut parameter for nearby large dust (unit of large dust area) is "<< fid_cut_par << endl;
      }
      if(runDescription[i][j] == "maxcl") {
	maxcl = atof( runDescription[i][2].c_str());
	cout << "max cluster number per view to start the  analysis (no dirty view control) is "<< maxcl << endl;
      }
      if(runDescription[i][j] == "cut_nofit") {
	cut_nofit = atof( runDescription[i][2].c_str());
	cout << "No fit cut for ellipticity equal to "<< cut_nofit << endl;
      }
      if(runDescription[i][j] == "cut_goodzone") {
	cut_goodzone = atof( runDescription[i][2].c_str());
	cout << "Boolean Good zone cut is "<< cut_goodzone << endl;
      }
      if(runDescription[i][j] == "cut_minor") {
	cut_minor = atof( runDescription[i][2].c_str());
	cout << "Lower cut on minor length is "<< cut_minor << " [um]" <<  endl;
      }
      if(runDescription[i][j] == "cut_ncl") {
	cut_ncl = atof( runDescription[i][2].c_str());
	cout << "Cut on number of merged cluster liked to a grain is "<< cut_ncl << endl;
      }
      if(runDescription[i][j] == "cut_mtrk") {
	cut_mtrk = atof( runDescription[i][2].c_str());
	cout << "Mtrk index is set to "<< cut_mtrk << endl;
      }
      if(runDescription[i][j] == "cut_npol") {
	cut_npol = atof( runDescription[i][2].c_str());
	cout << "Lower cut of number of polarization is "<< cut_npol << endl;
      } 
      if(runDescription[i][j] == "cut_isolated") {
	cut_isolated = atof( runDescription[i][2].c_str());
	cout << "Isolated index is set to "<< cut_isolated << endl;
      }
      if(runDescription[i][j] == "cut_bar_l") {
	cut_bar_l = atof( runDescription[i][2].c_str());
	cout << "Lower cut on barshift is "<< cut_bar_l << " [um]" << endl;
      }
      if(runDescription[i][j] == "cut_reg_l") {
	cut_reg_l = atof( runDescription[i][2].c_str());
	cout << "Lower cut on regularity index is "<< cut_reg_l << endl;
      }
      if(runDescription[i][j] == "cut_bar_u") {
	cut_bar_u = atof( runDescription[i][2].c_str());
	cout << "Upper cut on barshift is "<< cut_bar_u << " [um]" << endl;
      }
      if(runDescription[i][j] == "cut_reg_u") {
	cut_reg_u = atof( runDescription[i][2].c_str());
	cout << "Upper cut on regularity index is "<< cut_reg_u << endl;
      }
      if(runDescription[i][j] == "cut_view") {
	cut_view = atof( runDescription[i][2].c_str());
	cout << "Views removed "<< cut_view << endl;
      }

      if(runDescription[i][j] == "channel") {
	channel = atof( runDescription[i][2].c_str());
	cout << "channel "<< channel << endl;
      } 
    }
  }
}

void myData::Loop()
{
  
  //gROOT->ProcessLine("_file0->GetName()");
  //cout << _file0->GetName() << endl;
  /////////////  SETTINGS  ARGUMENT /////////////////////////////////////////////////        
  char datacard[100]="settings.mac";
  //cout << "Type the settings macro [settings.mac or settings_ini.mac or ls from shell]" << endl;
  //cout << "Enter some text here: ";
  // cin >> datacard;
  cout << "You are using: " << datacard << "\n\n";
  myDatacard(datacard);
  cout << "PAY ATTENTION 1" << endl;
  cout << "max number of polarization angles (excluding 180°) is "<< npol << endl;
  cout << "If the number of polarizations is changed please fix in 'myData_vx.C'" << endl;
  cout << "PAY ATTENTION 2" << endl;
  cout << "The fiducial cut parameter has to be checked in the histo 'r_dist_ldust'" << endl;
  cout << "If it is not correct please change 'fid_cut_par' in 'settings.mac'" << endl; 
  cout << "End of settings (enjoy the results)" << "\n\n";
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const char *file="dm_tracks_cl.dm.root";
  aRun = new DMRRun(file);
  aRun->SetFixEncoderFaults(1);
  view = aRun->GetView();

  //// LOG OUTPUT
  ofstream mybfcl("pred_bfcl.txt");    /// lista bfcl per ogni grano
  ofstream bfcl8("bfcl8copy.txt");     /// lista bfcl (solo 8pol) per funzioni immagini animate
  ofstream yandex("yandex_bfcl.txt");     /// lista bfcl (solo 8pol) per funzioni immagini animate
  ofstream cut8("bfcl8_with_cuts.txt");     /// lista bfcl (solo 8pol) per funzioni immagini animate con tagli

  /// HISTOGRAMS
  TH2F * hcutarea = new TH2F("cut_area","",len_view_x*100,-len_view_x/2.,len_view_x/2.,len_view_y*100,-len_view_y/2.,len_view_y/2.);  // mappa area_tagliata
  TH2F * hradius = new TH2F("radius","",totXbin,-totXbin/2.,totXbin/2.,totYbin,-totYbin/2.,totYbin/2.);  // mappa xy della proiezione x di phi_bar
  TH2F * hxmap = new TH2F("xmap","",totXbin,-totXbin/2.,totXbin/2.,totYbin,-totYbin/2.,totYbin/2.);  // mappa xy della proiezione x di phi_bar
  TH2F * hymap = new TH2F("ymap","",totYbin,-totYbin/2.,totYbin/2.,totYbin,-totYbin/2.,totYbin/2.);  // mappa xy della proiezione y di phi_bar
  TH2F * hrphimap = new TH2F("rphimap","",100,-pi/2.,pi/2.,100,0,0.1);  // mappa xy della di phi_bar ?
  TH2F * hphimap = new TH2F("phimap","",totYbin,-totYbin/2.,totYbin/2.,totYbin,-totYbin/2.,totYbin/2.);  // mappa xy della di phi_bar ?
  TH1F * hphicell = new TH1F("phimap_int","",100,-pi/2.,pi/2.);  // phi map integrata
  TH1F * hrdist = new TH1F("rdist","",256,0,20);  // distanza mutua tra grani dopo i cuts primari (ldust, nearby ldust, no_ell_fit, minor_cut, long_chains)
  TH1F * hrdist2 = new TH1F("rdist2","",128,0,1); // distanza tra i picchi dei grani npeaks
  TH1F * hrdist_ldust = new TH1F("rdist_ldust","",10000,0,500); // distanza di ogni grano no large dust dai grani di large dust (espressa in unità di raggio del large dust)
  TH1F * hmtrk = new TH1F("dist_mtrk","",10000,0,500); // distanza minima tra grani appartenenti a due microtracce differenti

  TH1F * hclx = new TH1F("clx","",800,-10,10);  // distanza mutua tra cluster di una collezione da otto [nm]
  TH1F * hcly = new TH1F("cly","",800,-10,10);  // distanza mutua tra cluster di una collezione da otto [nm]
  
  TCanvas *c1 = new TCanvas();
  TGraph * grNcl = new TGraph(); // ncl per view
  grNcl->SetName("ncl_per_view");

  //// DATA OUTPUT
  TFile * f_out = new TFile("debug6_test_grain.root","RECREATE");
  TTree * Tree_out = new TTree("tree1","data");
  

  Tree_out->Branch("eHeaderID",&eHeaderID,"eHeaderID/I");
  Tree_out->Branch("eViewID",&eViewID,"eViewID/I");
  Tree_out->Branch("eFlag",&eFlag,"eFlag/I");
  Tree_out->Branch("eGrainID",&eGrainID,"eGrainID/I");
  Tree_out->Branch("ePolID",&ePolID,"ePolID/I");
  Tree_out->Branch("eBfcID",&eBfcID,"eBfcID/I");
  Tree_out->Branch("eBfcPolID",&eBfcPolID,"eBfcPolID/I");
  Tree_out->Branch("eBfcGap",&eBfcGap,"eBfcGap/I");
  Tree_out->Branch("eBfcSigPeak",&eBfcSigPeak,"eBfcSigPeak/D");
  Tree_out->Branch("eBfcMeanBkg",&eBfcMeanBkg,"eBfcMeanBkg/D");
  Tree_out->Branch("eBfcWeigth",&eBfcWeigth,"eBfcWeigth/D");
  Tree_out->Branch("eSclSigPeak",&eSclSigPeak,"eSclSigPeak/D");
  Tree_out->Branch("eSclMeanBkg",&eSclMeanBkg,"eSclMeanBkg/D");
  Tree_out->Branch("eGrSigPeak",&eGrSigPeak,"eGrSigPeak/D");
  Tree_out->Branch("eGrMeanBkg",&eGrMeanBkg,"eGrMeanBkg/D");
  //Tree_out->Branch("eGrWeigth",&eGrWeigth,"eGrWeigth/D");
  Tree_out->Branch("eGrSclSigPeak",&eGrSclSigPeak,"eGrSclSigPeak/D");
  Tree_out->Branch("eGrSclMeanBkg",&eGrSclMeanBkg,"eGrSclMeanBkg/D");
  Tree_out->Branch("eNcl",&eNcl,"eNcl/I");
  Tree_out->Branch("eNgr",&eNgr,"eNgr/I");
  Tree_out->Branch("eNclFr",&eNclFr,"eNclFr/I");
  Tree_out->Branch("eGrainx",&eGrainx,"eGrainx/D");
  Tree_out->Branch("eGrainy",&eGrainy,"eGrainy/D");
  Tree_out->Branch("eGrainz",&eGrainz,"eGrainz/D");
  Tree_out->Branch("eClustx",&eClustx,"eClustx/D");
  Tree_out->Branch("eClusty",&eClusty,"eClusty/D");
  Tree_out->Branch("eClustz",&eClustz,"eClustz/D");
  Tree_out->Branch("eClustMin",&eClustMin,"eClustMin/D");
  Tree_out->Branch("eClustMaj",&eClustMaj,"eClustMaj/D");
  Tree_out->Branch("eClustEll",&eClustEll,"eClustEll/D");
  Tree_out->Branch("eClustPhi",&eClustPhi,"eClustPhi/D");
  Tree_out->Branch("eClustMaxPeak",&eClustMaxPeak,"eClustMaxPeak/D");
  Tree_out->Branch("eClustMeanBkg",&eClustMeanBkg,"eClustMeanBkg/D");
  Tree_out->Branch("ePuls",&ePuls,"ePuls/I");
  Tree_out->Branch("eGrainMin",&eGrainMin,"eGrainMin/D");
  Tree_out->Branch("eGrainMaj",&eGrainMaj,"eGrainMaj/D");
  Tree_out->Branch("eGrainEll",&eGrainEll,"eGrainEll/D");
  Tree_out->Branch("eGrainPhi",&eGrainPhi,"eGrainPhi/D");
  Tree_out->Branch("eGrainTheta",&eGrainTheta,"eGrainTheta/D");
  Tree_out->Branch("eEllPrjX",&eEllPrjX,"eEllPrjX/D");
  Tree_out->Branch("eEllPrjY",&eEllPrjY,"eEllPrjY/D");
  Tree_out->Branch("eGoodZone",&eGoodZone,"eGoodZone/B");
  Tree_out->Branch("eLargeDust",&eLargeDust,"eLargeDust/B");
  Tree_out->Branch("eClustTypeDust",&eClustTypeDust,"eClustTypeDust/B");
  Tree_out->Branch("eVolume",&eVolume,"eVolume/D");
  Tree_out->Branch("eBfcVolume",&eBfcVolume,"eBfcVolume/D");
  Tree_out->Branch("eFrBfVolume",&eFrBfVolume,"eFrBfVolume/D");
  Tree_out->Branch("eArea",&eArea,"eArea/D");
  Tree_out->Branch("eBfcArea",&eBfcArea,"eBfcArea/D");
  Tree_out->Branch("eFrBfArea",&eFrBfArea,"eFrBfArea/D");
  Tree_out->Branch("eGrainArea",&eGrainArea,"eGrainArea/D");
  Tree_out->Branch("eXView",&eXView,"eXView/D");
  Tree_out->Branch("eYView",&eYView,"eYView/D");
  Tree_out->Branch("eZ",&eZ,"eZ/D");
  Tree_out->Branch("eZlen",&eZlen,"eZlen/D");
  Tree_out->Branch("eNFrame",&eNFrame,"eNFrame/I");
  Tree_out->Branch("eSetFrame",&eSetFrame,"eSetFrame/I");
  Tree_out->Branch("eSetNCopy",&eSetNCopy,"eSetNCopy/I");
  Tree_out->Branch("eSetNStatic",&eSetNStatic,"eSetNStatic/I");
  Tree_out->Branch("eSetXRms",&eSetXRms,"eSetXRms/D");
  Tree_out->Branch("eSetYRms",&eSetYRms,"eSetYRms/D");
  Tree_out->Branch("eSetXBar",&eSetXBar,"eSetXBar/D");
  Tree_out->Branch("eSetYBar",&eSetYBar,"eSetYBar/D");
  Tree_out->Branch("eSetXMaxBar",&eSetXMaxBar,"eSetXMaxBar/D");
  Tree_out->Branch("eSetYMaxBar",&eSetYMaxBar,"eSetYMaxBar/D");
  Tree_out->Branch("eSetXMinBar",&eSetXMinBar,"eSetXMinBar/D");
  Tree_out->Branch("eSetYMinBar",&eSetYMinBar,"eSetYMinBar/D");
  Tree_out->Branch("eSetNpeaks",&eSetNpeaks,"eSetNpeaks/D");
  Tree_out->Branch("eSetNpeaksMax",&eSetNpeaksMax,"eSetNpeaksMax/D");
  Tree_out->Branch("eSetNpeaksNcopy",&eSetNpeaksNcopy,"eSetNpeaksNcopy/D");
  Tree_out->Branch("eSetNpeaksDist",&eSetNpeaksDist,"eSetNpeaksDist/D");
  Tree_out->Branch("eSetNpeaksPhi",&eSetNpeaksPhi,"eSetNpeaksPhi/D");
  Tree_out->Branch("eSetNpeaksDVol",&eSetNpeaksDVol,"eSetNpeaksDVol/D");
  Tree_out->Branch("eSetNpeaksDNpx",&eSetNpeaksDNpx,"eSetNpeaksDNpx/D");
  Tree_out->Branch("eSetNpeaksDBri",&eSetNpeaksDBri,"eSetNpeaksDBri/D");
  Tree_out->Branch("eSetNpeaksMaxPhiAmp",&eSetNpeaksMaxPhiAmp,"eSetNpeaksMaxPhiAmp/D");
  Tree_out->Branch("eSetGapZ",&eSetGapZ,"eSetGapZ/D");
  Tree_out->Branch("eSetGapZ2",&eSetGapZ2,"eSetGapZ2/D");
  Tree_out->Branch("eSetGap",&eSetGap,"eSetGap/D");
  Tree_out->Branch("eSetPath",&eSetPath,"eSetPath/D");
  Tree_out->Branch("eSetMeanPath",&eSetMeanPath,"eSetMeanPath/D");
  Tree_out->Branch("eSetRmsPath",&eSetRmsPath,"eSetRmsPath/D");
  Tree_out->Branch("eSetMaxPath",&eSetMaxPath,"eSetMaxPath/D");
  Tree_out->Branch("eSetMaxDist",&eSetMaxDist,"eSetMaxDist/D");
  Tree_out->Branch("eSetMaxPol1",&eSetMaxPol1,"eSetMaxPol1/D");
  Tree_out->Branch("eSetMaxPol2",&eSetMaxPol2,"eSetMaxPol2/D");
  Tree_out->Branch("eSetMaxBar",&eSetMaxBar,"eSetMaxBar/D");
  Tree_out->Branch("eSetMaxAmp",&eSetMaxAmp,"eSetMaxAmp/D");
  Tree_out->Branch("eSetRmsAmp",&eSetRmsAmp,"eSetRmsAmp/D");
  Tree_out->Branch("eSetMeanAmp",&eSetMeanAmp,"eSetMeanAmp/D");
  Tree_out->Branch("eSetPhi",&eSetPhi,"eSetPhi/D");
  Tree_out->Branch("eSetPhiRms",&eSetPhiRms,"eSetPhiRms/D");
  Tree_out->Branch("eSetPhiMean",&eSetPhiMean,"eSetPhiMean/D");
  Tree_out->Branch("eSetPhiBar",&eSetPhiBar,"eSetPhiBar/D");
  Tree_out->Branch("eSetPhiMaxAmp",&eSetPhiMaxAmp,"eSetPhiMaxAmp/D");
  Tree_out->Branch("eSetTheBar",&eSetTheBar,"eSetTheBar/D");
  Tree_out->Branch("eSetMeanBright",&eSetMeanBright,"eSetMeanBright/D");
  Tree_out->Branch("eSetBrAmp",&eSetBrAmp,"eSetBrAmp/D");
  Tree_out->Branch("eSetBrMaxPol",&eSetBrMaxPol,"eSetBrMaxPol/D");
  Tree_out->Branch("eSetBrMinPol",&eSetBrMinPol,"eSetBrMinPol/D");
  Tree_out->Branch("eSetBkgAmp",&eSetBkgAmp,"eSetBkgAmp/D");
  Tree_out->Branch("eSetPeakAmp",&eSetPeakAmp,"eSetPeakAmp/D");
  Tree_out->Branch("eIsolated",&eIsolated,"eIsolated/D");
  Tree_out->Branch("eSetNpxRms",&eSetNpxRms,"eSetNpxRms/D");
  Tree_out->Branch("eSetVolRms",&eSetVolRms,"eSetVolRms/D");
  Tree_out->Branch("eSetVolRatio",&eSetVolRatio,"eSetVolRatio/D");
  Tree_out->Branch("eMTrk",&eMTrk,"eMTrk/I");
  Tree_out->Branch("eMTID",&eMTID,"eMTID/I");
  Tree_out->Branch("eMTGr",&eMTGr,"eMTGr/I");
  Tree_out->Branch("eMTNfr",&eMTNfr,"eMTNfr/I");
  Tree_out->Branch("eMTPhi",&eMTPhi,"eMTPhi/D");
  Tree_out->Branch("eMTThe",&eMTThe,"eMThe/D"); 
  Tree_out->Branch("eMTLen",&eMTLen,"eMTLen/D");
  Tree_out->Branch("eMTVolDif",&eMTVolDif,"eMTVolDif/D");
  Tree_out->Branch("eChannel",&eChannel ,"eChannel/I");
  Tree_out->Branch("eDeltaPhi",&eDeltaPhi ,"eDeltaPhi/D");
  Tree_out->Branch("eGrFitMin",&eGrFitMin ,"eGrFitMin/D");
  Tree_out->Branch("eGrFitMaj",&eGrFitMaj ,"eGrFitMaj/D");
  Tree_out->Branch("eGrFitEll",&eGrFitEll ,"eGrFitEll/D");
  Tree_out->Branch("eGrFitPhi",&eGrFitPhi ,"eGrFitPhi/D");
  //Tree_out->Branch("eSetPhiRms",&eSetPhiRms,"eSetPhiRms/D");
  //Tree_out->Branch("eSetPhiMean",&eSetPhiMean,"eSetPhiMean/D");
  /////////////////////////////////////////////////////////

  /// READ INPUT FILE
  Long64_t nentries = fChain->GetEntries();
  cout << "nentries "<< nentries << endl;

  /////// TIPO DI SCANNING /////////  
  fChain->Draw("vid","","goff");
  Float_t vid_mean = TMath::Mean(fChain->GetSelectedRows(),fChain->GetV1());
  fChain->Draw("aid","","goff");
  Float_t aid_mean = TMath::Mean(fChain->GetSelectedRows(),fChain->GetV1());
  cout << aid_mean << " " << vid_mean << " " << viewID << endl;
  ////////////////////////// end
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // su tutte le View
    
    GetEntry(jentry);
    Long64_t ientry = LoadTree(jentry);
    
    //cout << aid << " " << vid << " " << viewID << endl;
    
    if(vid_mean>aid_mean)viewID=vid;
    else viewID=aid;

    hID=id;

    grNcl->SetPoint(viewID,viewID,ncl);
    hcutarea->Reset();
    
    if((ncl<maxcl || jentry==cut_view) /*&& flag==0*/){   /// cut per view troppo sporche (maxcl va settato)   evita i CRASH
      cout << viewID << " " << ncl << endl;
      area_scan += len_view_x*len_view_y;
      ///////// SEARCH ENCODING FAULTS
      //Double_t cl_x2[ncl];//=NULL;
      //Double_t cl_y2[ncl];//=NULL;
      Double_t *cl_x2 = new Double_t[ncl];
      Double_t *cl_y2 = new Double_t[ncl];

      //std::vector<double> cl_x2(ncl);
      //std::vector<double> cl_y2(ncl);
      //cl_x2.clear();
      //cl_y2.clear();

      for(int jn=0; jn<cl_;jn++){
	if(jn==0){
	  cl_x2[jn]=cl_x[jn];
	  cl_y2[jn]=cl_y[jn];
	}
	if(jn!=0){
	  if(fr_x[cl_ifr[jn]]!=fr_x[cl_ifr[0]]) cl_x2[jn]=cl_x[jn] - (fr_x[cl_ifr[jn]]-fr_x[cl_ifr[0]]);


	  else cl_x2[jn]=cl_x[jn];
	  
	  if(fr_y[cl_ifr[jn]]!=fr_y[cl_ifr[0]]) cl_y2[jn]=cl_y[jn] - (fr_y[cl_ifr[jn]]-fr_y[cl_ifr[0]]);
	  else cl_y2[jn]=cl_y[jn];
	}
      }
      ///////////////////////////

      //// INIZIALIZATIONS
      for(int in=0; in<gr_;in++){
	for(int jn=0;jn<npol;jn++){
	  xb_frbf[in][jn]=0;
	  yb_frbf[in][jn]=0;
	  vol_frbf[in][jn]=0;
	  npx_frbf[in][jn]=0;
	  same_frame[in][jn]=-1;
	  num_peaks[in][jn]=1;
	  ipol_gr[in][jn]=-1;
	  frbf_ent[in][jn]=0;
	}
	mtrk[gr_imt[in]]=false;
	mt_dif_vol[in]=-1;
	phi_mt[in]=-10;
	twopeak_dist[in]=-1;
	twopeak_dnpx[in]=-1;
	twopeak_dvol[in]=-1;
	twopeak_dbri[in]=-1;
	twopeak_phi[in]=10;
	bfc_gap[in]=-1;
	gr_180[in]=-1;
	num_tot_peaks[in]=0;

	gr_max_br_amp[in]=0;
	gr_max_peak_amp[in]=0;
	gr_mean_bkg_amp[in]=0;
	
      }
    
      /////////////////// PARTE 1: CARATTERIZZAZIONE DEI SINGOLI GRANI RICOSTRUITI ////////////////////////////////////////////////////////////////////
      
      for(int in=0; in<gr_;in++){
	//// inizialization

	//for(int jpol=0;jpol<npol;jpol++){
	//ipol_gr[in][jpol]=-1;
	//}

	//// INIZIALIZATIONS GRAINS
	tmp_gap=-1;
	tmp_gr=-1;    
	dist_cl_gr[in]=0;
	grain[in]=true;
	ctd[in]=false;
	ldust[in]=false;
	bfc_border[in]=false;
	index_clust=0;      
	n_ifr=1;
	//fiducial_cut=fid_cut_par*TMath::Sqrt((gr_npx[in]/gr_ncl[in])/(TMath::Pi()))*((x_pix_size-y_pix_size)/2.);
      
	//cout << fiducial_cut << endl;
	Int_t tcl_gr = gr_ncl[in]*npol;       
	Double_t *x_ldust = new Double_t[tcl_gr];
	Double_t *y_ldust = new Double_t[tcl_gr];
	Double_t *z_ldust = new Double_t[tcl_gr];
	Double_t *lx_ldust = new Double_t[tcl_gr];
	Double_t *ly_ldust = new Double_t[tcl_gr];
	Double_t *phi_ldust = new Double_t[tcl_gr];
	Double_t *area_ldust = new Double_t[tcl_gr];     
	TString chain;
	ostringstream frame;

	Double_t tmp_vol[npol]={-1,-1,-1,-1,-1,-1,-1,-1};
	Int_t bfc_ifr[npol]={-1,-1,-1,-1,-1,-1,-1,-1};
	Int_t bfc_zfr[npol]={-1,-1,-1,-1,-1,-1,-1,-1};
	Int_t tmp_same_frame[npol]={-1,-1,-1,-1,-1,-1,-1,-1};
	////////////////// end

            
	/// start cluster
	for(int jn=0; jn<cl_;jn++){ //loop sui clusters
	  //if(cl_igr[jn]==gr_id[in])cout  << gr_id[in] << " " << cl_igr[jn] <<  " " << cl_flags[jn] << " " << cl_ipol[jn] << endl;
	  ///////////////////////// start if cluster-grain
	  if(cl_igr[jn]==gr_id[in]  && cl_flags[jn]==0){  // clusters appartenenti all'iesimo grano e non mergiati, quindi collegati a specifica polarizzazione
	    ////// MAX BRIGHTNESS PER OGNI POLARIZZAZIONE  (CERCA IL BFC DEL GRANO PER OGNI POL) --> BFC definito come il più luminoso
	    if((double)cl_vol[jn]/cl_npx[jn]>tmp_vol[cl_ipol[jn]]){
		tmp_vol[cl_ipol[jn]]=(double)cl_vol[jn]/cl_npx[jn];
		ipol_gr[in][cl_ipol[jn]]=cl_id[jn];
		//cout << in << " " << " " << cl_ipol[jn] << " " << ipol_gr[in][cl_ipol[jn]] << " " << fr_iz[cl_ifr[ipol_gr[in][cl_ipol[jn]]]] <<  endl;
		//////////// FRAME DEL BEST FOCUS CLUSTER PER OGNI POL 
		bfc_ifr[cl_ipol[jn]]=cl_ifr[ipol_gr[in][cl_ipol[jn]]];              // frame
		bfc_zfr[cl_ipol[jn]]=fr_iz[cl_ifr[ipol_gr[in][cl_ipol[jn]]]];       // zeta
	      }
	    //if(in==159 && viewID==0 && cl_ipol[jn]==4)cout <<"bfc_ifr "<< bfc_ifr[cl_ipol[jn]] << " " <<(double)cl_vol[jn]/cl_npx[jn] << " "<< cl_id[jn]<< endl;
	    /////////////////////////////////////
	  
	  
	    ////////  MARGINI DEI GRANI   (LOOP SU TU TUTTI I CLUSTER COLLEGATI) --> i cluster di un grano formano un reticolo 
	    if(index_clust==0){
	      ld_area[in][0]=cl_x2[jn];
	      ld_area[in][1]=cl_y2[jn];
	      ld_area[in][2]=cl_x2[jn];
	      ld_area[in][3]=cl_y2[jn];
	    }
	    else{
	      if(cl_x2[jn]<ld_area[in][0])ld_area[in][0]=cl_x2[jn];
	      if(cl_y2[jn]<ld_area[in][1])ld_area[in][1]=cl_y2[jn];
	      if(cl_x2[jn]>ld_area[in][2])ld_area[in][2]=cl_x2[jn];
	      if(cl_y2[jn]>ld_area[in][3])ld_area[in][3]=cl_y2[jn];
	    }
	    index_clust++;	  
	  }  /// end if-cluster grain
	
	  ///// eGAP  (DISCREPANZA 0° E 180°)  /// DA VERIFICARE
	  if(cl_igr[jn]==gr_id[in]  && cl_flags[jn]==2){   // cluster del grano con polarizzazione 180° (ovvero ripetizione di 0°)
	    if((double)cl_vol[jn]/cl_npx[jn]>tmp_gap){   // prendo il più luminoso tra essi
	      tmp_gap=(double)cl_vol[jn]/cl_npx[jn];
	      gr_180[in]=cl_id[jn];  
	    }
	  }
	  //////////////////////// end if clust-grain 
	} // end cluster

      
	///////  CENTRO E RAGGIO DEL GRANO RICOSTRUITO
	grain_Ox[in]=(ld_area[in][2]+ld_area[in][0])/2.;
	grain_Oy[in]=(ld_area[in][3]+ld_area[in][1])/2.;
	grain_rx[in]=(ld_area[in][2]-ld_area[in][0])/2.;
	grain_ry[in]=(ld_area[in][3]-ld_area[in][1])/2.;
	grain_area[in]=(ld_area[in][2]-ld_area[in][0])*(ld_area[in][3]-ld_area[in][1]);
	fiducial_cut=fid_cut_par*TMath::Sqrt(grain_area[in]/TMath::Pi());//*((x_pix_size-y_pix_size)/2.);
	////////////////////////////////////////////////////

      
	////// MARGINI NEARBY LARGE DUST   (FATTO SUL BFC DEL GRANO)
	if(TMath::Log10((double)cl_vol[gr_ibfc[in]]/cl_npx[gr_ibfc[in]])>thr_ldust_br || TMath::Log10(cl_npx[gr_ibfc[in]])>thr_ldust_area){
	  ldust[in]=true;
	  cut_area[iLarge_dust][0]=grain_Ox[in]+x-fiducial_cut;
	  cut_area[iLarge_dust][1]=grain_Oy[in]+y-fiducial_cut;
	  cut_area[iLarge_dust][2]=grain_Ox[in]+x+fiducial_cut;
	  cut_area[iLarge_dust][3]=grain_Oy[in]+y+fiducial_cut;
	  min_area_x = TMath::Ceil((len_view_x/2. + cut_area[iLarge_dust][0]-x)*100) + 1;
	  min_area_y = TMath::Ceil((len_view_y/2. + cut_area[iLarge_dust][1]-y)*100) + 1;
	  max_area_x = TMath::Ceil((len_view_x/2. + cut_area[iLarge_dust][2]-x)*100) + 1;
	  max_area_y = TMath::Ceil((len_view_y/2. + cut_area[iLarge_dust][3]-y)*100) + 1;
	  cout << "pos "<< cut_area[iLarge_dust][0]-x << " " <<  cut_area[iLarge_dust][1]-y  << " " <<  cut_area[iLarge_dust][2]-x << " " <<  cut_area[iLarge_dust][3]-y << endl;
	  cout << "area_tagliata "<<min_area_x << " " << min_area_y << " " << max_area_x << " " << max_area_y << endl;
	  for(int i=min_area_x;i<(max_area_x+1);i++){
	    for(int j=min_area_y;j<(max_area_y+1);j++){
	      hcutarea->SetBinContent(i,j,1);
	    }
	  }
	  //cout << gr_x[in] << " " << cut_area[iLarge_dust][0] << " " << iLarge_dust << endl;
	  //area_cut += 2*fiducial_cut*fiducial_cut;
	  iLarge_dust++;
	}
	/////////// end MARGINI NEARBY LARGE DUST
      
	delete [] x_ldust;
	delete [] y_ldust;
	delete [] z_ldust;
	delete [] lx_ldust;
	delete [] ly_ldust;
	delete [] phi_ldust;
	delete [] area_ldust;

      
	//// eGoodZone (SOLO GRANI LONTANI DAI LARGE DUST)
	// su tutti i grani di una View
	if(!ldust[in]){
	  for(int jn=0;jn<iLarge_dust;jn++){
	    if(gr_x[in]+x>=cut_area[jn][0] && gr_x[in]+x<=cut_area[jn][2] && gr_y[in]+y>=cut_area[jn][1] && gr_y[in]+y<=cut_area[jn][3])
	      grain[in]=false;   // near large dust
	  }
	} else grain[in]=false;   // is a large dust
	//// end eGoodZone 


	/////////  NUMERO DI POLARIZZAZIONI IN UN GRANO 
	Int_t nCopy=npol;
	for(int jn=0;jn<npol;jn++){
	  if(ipol_gr[in][jn]==-1)nCopy--;
	}
	//cout << nCopy << endl;
	///// end
      
	///// GAP Z (VIBRATION OR ENCODER FAULTS AND BFCFR GAP) //////////////////////
	//Double_t zbfcl[nCopy];//={};	
	Double_t *zbfcl = new Double_t[nCopy];
	//std::vector<double> zbfcl(nCopy);
	Int_t iCopy=0;
	for(int jn=0;jn<npol;jn++){
	  if(ipol_gr[in][jn]!=-1){
	    zbfcl[iCopy]=fr_z[cl_ifr[ipol_gr[in][jn]]];
	    //cout << in << " " << " " << iCopy << " " << ipol_gr[in][jn] << " " << fr_z[cl_ifr[ipol_gr[in][jn]]] <<  " " << zbfcl[iCopy] << endl;
	    iCopy++;
	  }
	}
	if(nCopy!=0)gr_gap_z[in] = TMath::MaxElement(nCopy,zbfcl) - TMath::MinElement(nCopy,zbfcl);
	//if(nCopy!=0)gr_gap_z[in] = *max_element(zbfcl.begin(), zbfcl.end()) - *min_element(zbfcl.begin(), zbfcl.end());
	else gr_gap_z[in]=-100;
	//cout << gr_gap_z[in] << endl;
	////////////////////////////////////////////////////////////////
      
      
	///// CORRECTION OF BFC FRAME //////////////////////////////////  (correggo i frame con diverso zeta rispetto alla popolazione dominante)
	Int_t tmp_nfr[npol]={0,0,0,0,0,0,0,0};
	Int_t tmp_clfr[npol]={-1,-1,-1,-1,-1,-1,-1,-1};
	Int_t ifr=0;
      
	for(int jn=0;jn<npol;jn++){
	  //cout << in << " " << " " << jn << " " << ipol_gr[in][jn] << " " << fr_iz[cl_ifr[ipol_gr[in][jn]]] <<  endl;
	  bool match_fr=false;               // booleano
	  if(ipol_gr[in][jn]!=-1 && jn==0){  // se la polarizzazione è 0
	    tmp_clfr[ifr]=bfc_zfr[jn];       // il tmp frz è il bfczfr di 0
	    tmp_nfr[ifr]++;                  // è di partenza sempre 1 (numero di frame uguali)
	    //cout << in << " " << jn << " " << bfc_fr[jn] << " " << tmp_clfr[ifr] <<  " " << tmp_nfr[jn] << endl;
	    ifr++;                           // numero di frame tra i bfcfr (iframe diventa 1 (numero minimo di frame))
	  }
	  for(int kn=0;kn<ifr;kn++){         
	    if(ipol_gr[in][jn]!=-1 && jn>0 && bfc_zfr[jn]==tmp_clfr[kn]){  // solo pol diverse da 0 e con frame uguale al frame della pol jn-esima
	      tmp_nfr[kn]++;    // conta numero di frame uguali al frame della polarizzazione iesima, quando è uguale match_fr diventa true
	      match_fr=true;
	      //cout << in << " " << jn << " " << kn << " " << bfc_fr[jn] << " " << tmp_clfr[kn] << " " << tmp_nfr[kn] << endl;
	    }
	  }
	  if(ipol_gr[in][jn]!=-1 && jn>0 && !match_fr){   // solo se match_fr è falso, cioè se diverso dal frame 0, guardo i frame successivi
	    ifr++;                                        // incrementa
	    tmp_clfr[ifr]=bfc_zfr[jn];                    // nuovo frame nel tmp
	    tmp_nfr[ifr]++;                               // tmp frame è 1 (molteplicità minima)
	    //cout << in << " " << jn << " " << bfc_fr[jn] <<  " " << tmp_clfr[ifr] << " " << tmp_nfr[ifr] << endl;
	  }
	}

	Int_t fr_sort[npol]={-1,-1,-1,-1,-1,-1,-1,-1};   // fr_sort ordina in ordine decrescente per indice i molteplicità dei frame
	TMath::Sort(npol,tmp_nfr,fr_sort);
	//cout << "s1 " << tmp_nfr[fr_sort[0]] << " " << tmp_clfr[fr_sort[0]] << endl;
	//cout << "s2 " << tmp_nfr[fr_sort[1]] << " " << tmp_clfr[fr_sort[1]] << endl;
	//cout << "s3 " << tmp_nfr[fr_sort[2]] << " " << tmp_clfr[fr_sort[2]] << endl;
	Int_t nbfc_same_fr = TMath::MaxElement(npol,tmp_nfr);   // Indice di massima molteplicità tra i frame
	Int_t max_fr = tmp_clfr[fr_sort[0]];                 // cerca il frame associato all'indice di massima molteplicità

	//cout << nbfc_same_fr << " " << max_fr << endl;

	if((nCopy-nbfc_same_fr)/(nCopy/2.)<1 && nCopy>=5){  // start correction (algoritmo di correzione, solo per collezioni con almeno 5 elementi) 
	  /// start cluster
	  for(int jn=0; jn<cl_;jn++){	 
	    if(cl_igr[jn]==gr_id[in]  && cl_flags[jn]==0){   // cerco tra i cluster associati al grano
	      //cout <<in << " " << cl_ipol[jn] << " " <<  max_fr << " " << bfc_zfr[cl_ipol[jn]] << " " << fr_iz[cl_ifr[jn]] << endl;
	      if(fr_iz[cl_ifr[jn]]==max_fr && bfc_zfr[cl_ipol[jn]]!=max_fr){ // se il bfclfr non coincide con quello di massima molteplicità
		ipol_gr[in][cl_ipol[jn]]=cl_id[jn];                          // rimpiazza quest'ultimo al posto del precedente
		bfc_ifr[cl_ipol[jn]]=cl_ifr[ipol_gr[in][cl_ipol[jn]]];              // frame
		bfc_zfr[cl_ipol[jn]]=fr_iz[cl_ifr[ipol_gr[in][cl_ipol[jn]]]];
		bfc_gap[in]=1;                                               // segnalo che c'era un gap nei bfclfr
		//cout << max_fr << " " << bfc_zfr[cl_ipol[jn]] << " " << fr_iz[cl_ifr[jn]] << endl;
	      }
	    }	  
	    //////////////////////// end if clust-grain 
	  } // end cluster
      
	} ///// end correction
	if((nCopy-nbfc_same_fr)==0 && nCopy!=0)bfc_gap[in]=0;           // gap zero se la molteplicità massima è uguale al numero di copie
	//cout << "gfsdff "<< bfc_gap[in] << endl;
	////////////////////////////////////////////////////////////////   end correction bfcfr


	///// GAP Z (CROSS-CHECK BFC_GAP2) ---> controllo se la correzione ha avuto un effetto desiderato 
	//Double_t zbfcl2[nCopy];//={};
	Double_t *zbfcl2 = new Double_t[nCopy];
	//std::vector<double> zbfcl2(nCopy);
	iCopy=0;
	for(int jn=0;jn<npol;jn++){
	  if(ipol_gr[in][jn]!=-1){
	    zbfcl2[iCopy]=fr_z[cl_ifr[ipol_gr[in][jn]]];
	    //cout << in << " " << " " << iCopy << " " << ipol_gr[in][jn] << " " << fr_z[cl_ifr[ipol_gr[in][jn]]] <<  " " << zbfcl2[iCopy] << endl;
	    iCopy++;
	  }
	}
	//if(nCopy!=0)gr_gap_z2[in] = *max_element(zbfcl2.begin(), zbfcl2.end()) - *min_element(zbfcl2.begin(), zbfcl2.end());
	if(nCopy!=0)gr_gap_z2[in] = TMath::MaxElement(nCopy,zbfcl2) - TMath::MinElement(nCopy,zbfcl2);
	else gr_gap_z2[in]=-100;
	//cout <<"gap rec "<< gr_gap_z2[in] << endl;
	////////////////////////////////////////////////////////////////

	delete [] zbfcl;
	delete [] zbfcl2; 

	/// FRAME OF THE BEST FOCUS CLUSTER FOR EACH POLARIZATION  (PIU PRECISO SE NEL FR DEL BFC CI SONO PIU CLUSTER)
	for(int jn=0;jn<cl_;jn++){
	  if(cl_igr[jn]==gr_id[in]  && cl_flags[jn]==0 && cl_ifr[jn]==cl_ifr[ipol_gr[in][cl_ipol[jn]]] && ipol_gr[in][cl_ipol[jn]]!=-1){
	    xb_frbf[in][cl_ipol[jn]]+=cl_x2[jn];//*cl_vol[jn]/cl_npx[jn];
	    yb_frbf[in][cl_ipol[jn]]+=cl_y2[jn];//*cl_vol[jn]/cl_npx[jn];
	    vol_frbf[in][cl_ipol[jn]]+=cl_vol[jn];
	    npx_frbf[in][cl_ipol[jn]]+=cl_npx[jn];
	    frbf_ent[in][cl_ipol[jn]]++;
	    //cout << in << " " << cl_ipol[jn] << " " <<  xb_frbf[in][cl_ipol[jn]] << " " << yb_frbf[in][cl_ipol[jn]] << " " << cl_vol[jn] << " " << cl_npx[jn] << endl;
	  }
	  ////////////////////////////////////////
		
	  ////// start frame //////////////// (MOLTEPLICITA' BEST FOCUS CLUSTER FRAME - SAME_FRAME==2)
	  if(cl_igr[jn]==gr_id[in]  && cl_flags[jn]==0){
	    if(cl_ifr[jn]==bfc_ifr[cl_ipol[jn]] && jn!=ipol_gr[in][cl_ipol[jn]]){
	      same_frame[in][cl_ipol[jn]]=2;
	      num_peaks[in][cl_ipol[jn]]++;
	    }
	    
	    //if(viewID==0&&in==159)cout << in << " " << jn << " " << ipol_gr[in][cl_ipol[jn]] << " " <<  cl_ipol[jn]  << " " << bfc_ifr[cl_ipol[jn]] << " " << cl_ifr[jn] << " " << tmp_same_frame[cl_ipol[jn]] << " " << same_frame[in][cl_ipol[jn]] << endl;	  
	  }
	  //////////// end frame
	}
	/////////////////  end
           
	//////////// BARYCENTER BEST FOCUS CLUSTER FRAME  (BARICENTRO CALCOLATO NEL FRAME DEL BFC, PIU PRECISO DEL BAR DEL SOLO BFC)  
	for(int jn=0;jn<npol;jn++){
	  //cout << "2 "<< jn << " " << ipol_gr[in][jn] << " " << fr_iz[cl_ifr[ipol_gr[in][jn]]] << endl;
	  if(ipol_gr[in][jn]!=-1){
	    xb_frbf[in][jn]=xb_frbf[in][jn]/frbf_ent[in][jn];//]vol_frbf[in][jn]/npx_frbf[in][jn];
	    yb_frbf[in][jn]=yb_frbf[in][jn]/frbf_ent[in][jn];//]vol_frbf[in][jn]/npx_frbf[in][jn];
	    //cout << in << " " << jn << " " << ipol_gr[in][jn] << " " << xb_frbf[in][jn] <<  " " << vol_frbf[in][jn] << endl;
	  }
	  num_tot_peaks[in] = num_tot_peaks[in] + num_peaks[in][jn];
	  //if(same_frame[in][jn]==2)cout << "primo_debug "<<num_tot_peaks[in] << " " << num_peaks[in][jn] << endl;
	}
	/////////////////// end     
      }
      /////////////////// END 1

      
      /////////////////// PARTE 2: CARATTERIZZAZIONE DELLE MICROTRACCE  ////////////////////////////////////////////////////////////////////
    
      ///// DIFFERENZA IN VOLUME MAX E MIN TRA I GRANI DI UNA MICROTRACCIA E CALCOLO DI PHI
      for(int kn=0;kn<mt_;kn++){
	//tmp_mtrk_dist[kn]=1000;
	int imt=0;
	int tmp_mt_id=-1;
	Double_t *vol_mt = new Double_t[mt_ngr[kn]];
	Double_t *x_mt = new Double_t[mt_ngr[kn]];
	Double_t *y_mt = new Double_t[mt_ngr[kn]];
	for(int in=0;in<gr_;in++){
	  if(gr_imt[in]==kn){
	    vol_mt[imt]=gr_vol[in];
	    x_mt[imt]=gr_x[in];
	    y_mt[imt]=gr_y[in];
	    if(imt==0)tmp_mt_id=in;
	    //cout << jentry << " " << kn << " " << in << " " << vol_mt[imt] << " " << tmp_mt_id <<  " " << x_mt[imt] << " " << y_mt[imt] << endl;
	    imt++;
	  }
	}
	if(mt_ngr[kn]==2){   // coefficiente angolare se il numero di grani è 2
	  phi_mt[tmp_mt_id]=TMath::ATan((y_mt[1]-y_mt[0])/(x_mt[1]-x_mt[0]));
	  //cout << phi_mt[tmp_mt_id] << endl;
	}
	if(mt_ngr[kn]>2){   /// fit lineare se il numero di grani è maggiore di 2
	  TGraph *grmt = new TGraph(mt_ngr[kn],x_mt,y_mt);
	  grmt->Fit("pol1","Q");
	  TF1 *pol1 = grmt->GetFunction("pol1");
	  phi_mt[tmp_mt_id] = TMath::ATan(pol1->GetParameter(1));
	  pol1->ReleaseParameter(1);
	  //cout << phi_mt[tmp_mt_id] << endl;
	}
	mt_dif_vol[tmp_mt_id]=TMath::MaxElement(imt,vol_mt) - TMath::MinElement(imt,vol_mt); // differenza tra volume massimo e minimo dei grani in una microtraccia
	//cout << tmp_mt_id << " " << mt_dif_vol[tmp_mt_id] << endl;
	delete [] vol_mt;
	delete [] x_mt;
	delete [] y_mt;
	//delete grmt;
	//delete pol1;
      }
      //////////////////////// END 2
    
      
      /////////////////// PARTE 3: COLLEZIONI DI GRANI  ////////////////////////////////////////////////////////////////////////////////////////////////
      
      for(int in = 0;in<gr_; in++){

	if(gr_phi[in]<=(TMath::Pi()/2.))eGrainPhi=gr_phi[in];
	if(gr_phi[in]>(TMath::Pi()/2.) && gr_phi[in]<=(3*TMath::Pi()/2.))eGrainPhi=gr_phi[in] - TMath::Pi();
	if(gr_phi[in]>(3*TMath::Pi()/2.) && gr_phi[in]<=(2*TMath::Pi()))eGrainPhi=gr_phi[in] - 2*TMath::Pi();

	mybfcl  << "Id candidate: "<< gr_id[in] << endl;
	if(gr_imt[in]!=-1)mybfcl << "It's a microtrack " << endl;

	///////////////// ELEMENTI IN UNA COLLEZIONE
	gr_copy[in]=npol;     
	for(int jpol=0;jpol<npol;jpol++){
	  if(ipol_gr[in][jpol]==-1)gr_copy[in]--;
	}
	///////////////////////////////////////////

	//////////////// INIZIALIZATIONS
	index_pol=0;
	index_2pol=0;
	match=0;
	Int_t max_br_pol=0;
	Int_t min_br_pol=255;
	Double_t *gr_phi_pol = new Double_t[gr_copy[in]];
        Double_t *gr_min_pol = new Double_t[gr_copy[in]];
	Double_t *gr_maj_pol = new Double_t[gr_copy[in]];
	Double_t *gr_ell_pol = new Double_t[gr_copy[in]];
	Int_t ndelta= (gr_copy[in]*(gr_copy[in]-1))/2;
	//if(ndelta<2)ndelta=2;
	Int_t idelta=0;
	//Double_t *delta_phi_pol = new Double_t[gr_copy[in]];
        Double_t *delta_phi_pol = new Double_t[ndelta];
	Double_t *gr_x_pol = new Double_t[gr_copy[in]];
	Double_t *gr_y_pol = new Double_t[gr_copy[in]];
	Double_t *gr_x_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_y_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_npx_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_vol_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_bright_pol_bar = new Double_t[gr_copy[in]];
	Double_t *gr_z_pol = new Double_t[gr_copy[in]];
	Double_t *gr_angpol = new Double_t[gr_copy[in]];
	Double_t *gr_npx_pol = new Double_t[gr_copy[in]];
	Double_t *gr_pol_id = new Double_t[gr_copy[in]];
	Double_t *gr_amp = new Double_t[gr_copy[in]];
	Double_t *gr_step = new Double_t[npol];

	Double_t *bkg_mean = new Double_t[gr_copy[in]];
	Double_t *sig_peak = new Double_t[gr_copy[in]];
	Double_t *br_amp = new Double_t[gr_copy[in]];
	Double_t *gr_xy_dist_bar = new Double_t[gr_copy[in]];

        Double_t clx_rel_dist=0;
	Double_t cly_rel_dist=0;

	q_line[in]=0;
	m_line[in]=0;
	gr_max_amp[in]=0;
	gr_max_dist_bar[in]=0;
	gr_max_dist[in]=0;
	gr_path[in]=0;
	gr_fr[in]=0;
	gr_isolated[in]=-1;
	gr_same_pos[in]=0;
	gr_npeaks_npol[in]=0;
	gr_npeaks[in]=0;
      
	Double_t sort_gr[npol]={-1,-1,-1,-1,-1,-1,-1,-1};

	////////////////////////////////////////////////////

	
	/////////// CARATTERIZZAZIONE ELEMENTI DI UNA COLLEZIONE
	for(int jpol=0;jpol<npol;jpol++){
	  if(ipol_gr[in][jpol]!=-1){
	    gr_phi_pol[index_pol]=cl_phi[ipol_gr[in][jpol]];	    
	    gr_x_pol_bar[index_pol]=xb_frbf[in][jpol];
	    gr_y_pol_bar[index_pol]=yb_frbf[in][jpol];
	    gr_npx_pol_bar[index_pol]=npx_frbf[in][jpol];
	    gr_vol_pol_bar[index_pol]=vol_frbf[in][jpol];
	    gr_bright_pol_bar[index_pol]=vol_frbf[in][jpol]/npx_frbf[in][jpol];
	    gr_x_pol[index_pol]=cl_x2[ipol_gr[in][jpol]];
	    gr_y_pol[index_pol]=cl_y2[ipol_gr[in][jpol]];
	    gr_z_pol[index_pol]=cl_z[ipol_gr[in][jpol]];

	    gr_maj_pol[index_pol]=cl_lx[ipol_gr[in][jpol]];
	    gr_min_pol[index_pol]=cl_ly[ipol_gr[in][jpol]];
	    gr_ell_pol[index_pol]=cl_lx[ipol_gr[in][jpol]]/cl_ly[ipol_gr[in][jpol]];
	   
	    
	    gr_npx_pol[index_pol]=cl_npx[ipol_gr[in][jpol]];
	    gr_pol_id[index_pol]=cl_ipol[ipol_gr[in][jpol]];
	    gr_angpol[index_pol]=cl_pol[ipol_gr[in][jpol]];
	    if(cl_ipol[ipol_gr[in][jpol]]==0)gr_vol_ref[in]=vol_frbf[in][jpol];	  
	    if(same_frame[in][jpol]==2){
	      gr_isolated[in]=2;
	      gr_npeaks_npol[in]++;  //da implementare
	      gr_npeaks[in] = gr_npeaks[in] + num_peaks[in][jpol]-1; // -1 per non avere molteplicità 2
	    }
	    set_fr[in]=cl_ifr[gr_ibfc[in]];
	    gr_npx_sum[in] += npx_frbf[in][jpol];
	    gr_vol_sum[in] += vol_frbf[in][jpol];
	    tmp_ifr=fr_iz[cl_ifr[gr_ibfc[in]]];
	    //cout << index_pol << " " <<  xb_frbf[in][jpol] << " " << yb_frbf[in][jpol] << endl;
	    //if(gr_isolated[in]==-1)cout << index_pol << " " << gr_x_pol[index_pol] << " " << gr_y_pol[index_pol] << " " << gr_x_pol_bar[index_pol] << " " << gr_y_pol_bar[index_pol] <<  endl;
	    //if(gr_isolated[in]==2)cout << in << " " << index_pol << " " << ipol_gr[in][jpol] << " " << gr_x_pol_bar[index_pol] << " " << gr_y_pol_bar[index_pol] << " " << gr_x_pol[index_pol] << " " << gr_y_pol[index_pol] << endl;

	    mybfcl  << "Bfc candidates: "<< ipol_gr[in][jpol] << " " << gr_pol_id[index_pol] << " " << gr_x_pol[index_pol] << " " << gr_y_pol[index_pol] << " " <<  gr_z_pol[index_pol] << endl;
	    sort_gr[cl_ipol[ipol_gr[in][jpol]]]=ipol_gr[in][jpol];
	    
	    if(gr_copy[in]==npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && flag==0){
	      
	      //cout << "PRIMA" << endl;
	      //if(viewID>0 &&  viewID<10)
	      //images_info(aRun,id,ipol_gr[in][jpol],jpol,im_info);
	      //cout << "eccomiiiiiiiii " << jentry << " " << in << endl;
	      //if(jentry==1 && in==2 && index_pol==0)images_ellipse(aRun,1,7,0,im_info,c1);
	      //cout << jentry << " " << in << " " << im_info[0] << " " << im_info[1] << endl;
	      //cout << "DOPO" << endl;
	      bkg_mean[index_pol]=im_info[0];
	      sig_peak[index_pol]=im_info[1];
	      br_amp[index_pol]=im_info[1]-im_info[0];
	      if(sig_peak[index_pol]>max_br_pol){
		gr_max_br_pol[in]=gr_pol_id[index_pol];
		max_br_pol = sig_peak[index_pol];
	      }
	      if(sig_peak[index_pol]<min_br_pol){
		gr_min_br_pol[in]=gr_pol_id[index_pol];
		min_br_pol = sig_peak[index_pol];
	      }
	    }
	    index_pol++;
	  }
	}
	////////////////////////////////////////////////////////////////

	
	///// SET STATISTICS INFO /////////////////////
	if(index_pol>1){
	  gr_max_br_amp[in] = TMath::Mean(index_pol,br_amp);
	  //gr_phi_max[in] = TMath::MaxElement(index_pol,gr_phi_pol)-TMath::MinElement(index_pol,gr_phi_pol);
	  gr_x_rms[in] = TMath::RMS(index_pol,gr_x_pol,gr_bright_pol_bar);
	  gr_y_rms[in] = TMath::RMS(index_pol,gr_y_pol,gr_bright_pol_bar);
	  gr_z_rms[in] = TMath::RMS(index_pol,gr_z_pol,gr_bright_pol_bar);
	  gr_npx_rms_bar[in] = TMath::RMS(index_pol,gr_npx_pol_bar);
	  gr_vol_rms_bar[in] = TMath::RMS(index_pol,gr_vol_pol_bar);
	  gr_vol_mean_bar[in] = TMath::Mean(index_pol,gr_bright_pol_bar);       
	  gr_x_mean[in] = TMath::Mean(index_pol,gr_x_pol,gr_bright_pol_bar);
	  gr_y_mean[in] = TMath::Mean(index_pol,gr_y_pol,gr_bright_pol_bar);
	  //cout << gr_phi_rms[in] << endl;
	}
      
	///////// MAX DISTANCE E PHI PER BFCL E BARYCENTER IN BFCLFR
	for(int iset = 0;iset<index_pol;iset++){
	  /// clusters info
	  cl_x_pos[in][iset]=gr_x_pol[iset];
	  cl_y_pos[in][iset]=gr_y_pol[iset];
	  cl_z_pos[in][iset]=gr_z_pol[iset];
	  cl_min_bf[in][iset]=gr_min_pol[iset];
	  cl_maj_bf[in][iset]=gr_maj_pol[iset];
	  cl_ell_bf[in][iset]=gr_ell_pol[iset];
	  cl_phi_bf[in][iset]=gr_phi_pol[iset];
	  cl_bkg_mean[in][iset]=bkg_mean[iset];
	  cl_sig_peak[in][iset]=sig_peak[iset];

	  //// PRELIMINAR CORRECTION OF DISTORSION
	  //gr_x_pol_bar[iset] -= xcorr[iset]/1000.; //um
	  //gr_y_pol_bar[iset] -= ycorr[iset]/1000.; //um
	  gr_xy_dist_bar[iset] = TMath::Sqrt(TMath::Power(gr_x_pol[iset]-gr_x_mean[in],2) + TMath::Power(gr_y_pol[iset]-gr_y_mean[in],2));
	}

	gr_mean_dist_bar[in]  = TMath::Mean(npol,gr_xy_dist_bar);
	gr_min_maj[in]  = TMath::MinElement(npol,gr_maj_pol);
	/*
	for(int iset = 0;iset<index_pol;iset++){
	  //Double_t ang_bar = TMath::ATan((gr_y_pol_bar[iset]-gr_y_mean[in])/(gr_x_pol_bar[iset]-gr_x_mean[in]));
	  //Double_t gap_bar = (major_mean[4]-major_mean[iset]);
	  //Double_t gap_bar = TMath::Abs(gr_maj_pol[iset] - gr_min_maj[in]);
	  //Double_t xgap_bar = TMath::Abs(gap_bar*TMath::Cos(phi_elliptic[iset]));
	  //Double_t ygap_bar = TMath::Abs(gap_bar*TMath::Sin(phi_elliptic[iset]));
	  
	  //cout << iset << " " << gr_mean_dist_bar[in] << " " << gr_xy_dist_bar[iset] << " " << gap_bar << " " << xgap_bar << " " << ygap_bar << " " << ang_bar << endl;
	  //cout << iset << " " << gr_x_mean[in] << " " << gr_y_mean[in] << " " << gr_x_pol_bar[iset] << " " << gr_y_pol_bar[iset] << endl;
	  //if(gr_xy_dist_bar[iset]<gr_mean_dist_bar[in]){
	  //  if((gr_y_pol_bar[iset] - gr_y_mean[in])<0) gr_y_pol_bar[iset] -= ygap_bar;
	  //  if((gr_y_pol_bar[iset] - gr_y_mean[in])>0) gr_y_pol_bar[iset] += ygap_bar;
	  //  if((gr_x_pol_bar[iset] - gr_x_mean[in])<0) gr_x_pol_bar[iset] -= xgap_bar;
	  //  if((gr_x_pol_bar[iset] - gr_x_mean[in])>0) gr_x_pol_bar[iset] += xgap_bar;
	 // }
	  //if(gr_xy_dist_bar[iset]>gr_mean_dist_bar[in]){
	  //if((gr_y_pol_bar[iset] - gr_y_mean[in])<0) gr_y_pol_bar[iset] += ygap_bar;
	    //if((gr_y_pol_bar[iset] - gr_y_mean[in])>0) gr_y_pol_bar[iset] -= ygap_bar;
	    //if((gr_x_pol_bar[iset] - gr_x_mean[in])<0) gr_x_pol_bar[iset] += xgap_bar;
	    //if((gr_x_pol_bar[iset] - gr_x_mean[in])>0) gr_x_pol_bar[iset] -= xgap_bar;
	    //}
	  //gr_y_pol_bar[iset] -= xy_gap[iset];
	  if((gr_y_pol_bar[iset]-gr_y_mean[in])>0)gr_y_pol_bar[iset] -= xy_gap[iset];
	  if((gr_y_pol_bar[iset]-gr_y_mean[in])<0)gr_y_pol_bar[iset] += xy_gap[iset];
	  //cout << iset << " " << gap_bar << " " << xgap_bar << " " << ygap_bar << " " << gr_x_pol_bar[iset] << " " << gr_y_pol_bar[iset] << endl;
	  }*/
	
	for(int iset = 0;iset<index_pol;iset++){
	  /////
	  for(int jset = 0;jset<index_pol;jset++){
	    if(TMath::Sqrt(TMath::Power(gr_x_pol[iset]-gr_x_pol[jset],2)+TMath::Power(gr_y_pol[iset]-gr_y_pol[jset],2))>gr_max_dist[in]){  // test
	      gr_max_dist[in]=TMath::Sqrt(TMath::Power(gr_x_pol[iset]-gr_x_pol[jset],2)+TMath::Power(gr_y_pol[iset]-gr_y_pol[jset],2));
	      phi_set[in] = TMath::ATan((gr_y_pol[iset]-gr_y_pol[jset])/(gr_x_pol[iset]-gr_x_pol[jset]));
	      gr_x_maxbar[in]=gr_x_pol[iset];
	      gr_y_maxbar[in]=gr_y_pol[iset];
	      gr_x_minbar[in]=gr_x_pol[jset];
	      gr_y_minbar[in]=gr_y_pol[jset];
	      gr_maxpol1[in]=gr_angpol[iset];
	      gr_maxpol2[in]=gr_angpol[jset];
	    }
	    //// BFCLFR
	    if(TMath::Sqrt(TMath::Power(gr_x_pol_bar[iset]-gr_x_pol_bar[jset],2)+TMath::Power(gr_y_pol_bar[iset]-gr_y_pol_bar[jset],2))>gr_max_dist_bar[in]){
	      gr_max_dist_bar[in]=TMath::Sqrt(TMath::Power(gr_x_pol_bar[iset]-gr_x_pol_bar[jset],2)+TMath::Power(gr_y_pol_bar[iset]-gr_y_pol_bar[jset],2));
	      phi_set_bar[in] = TMath::ATan((gr_y_pol_bar[iset]-gr_y_pol_bar[jset])/(gr_x_pol_bar[iset]-gr_x_pol_bar[jset]));
	      the_set_bar[in] = TMath::ATan((TMath::Power((gr_y_pol_bar[iset]-gr_y_pol_bar[jset]),2)+TMath::Power((gr_x_pol_bar[iset]-gr_x_pol_bar[jset]),2))/(gr_z_pol[iset]-gr_z_pol[jset]));
	      m_line[in] = TMath::Tan(phi_set_bar[in]);                             // servono per l'ampiezza  (coefficiente angolare della retta)
	      q_line[in] = gr_y_pol_bar[jset] - m_line[in]*gr_x_pol_bar[jset];      // servono per l'ampiezza  (intercetta)
	      //cout << in << " " << iset << " " << cl_x_pos[in][iset] << " " << cl_y_pos[in][iset] << " " << cl_z_pos[in][iset] <<  endl;
	      //cout << in << " " << jset << " " << cl_x_pos[in][jset] << " " << cl_y_pos[in][jset] << " " << cl_z_pos[in][jset] <<  endl;
	      //cout << the_set_bar[in] << " " << gr_z_pol[iset]-gr_z_pol[jset] << endl;
	    }
	    //// MAX BRIGHTNESS AMPLITUDE
	    /*if(TMath::Abs((sig_peak[iset]-bkg_mean[iset])-(sig_peak[jset]-bkg_mean[jset]))>gr_max_br_amp[in]){
	      gr_max_br_amp[in] = TMath::Abs((sig_peak[iset]-bkg_mean[iset])-(sig_peak[jset]-bkg_mean[jset]));	      
	    }*/
	    //// MAX PEAK AMPLITUDE
	    if(TMath::Abs(sig_peak[iset] - sig_peak[jset])>gr_max_peak_amp[in]){
	      gr_max_peak_amp[in] = TMath::Abs(sig_peak[iset]-sig_peak[jset]);	      
	    }
	    //// MEAN BKG AMPLITUDE
	    if(TMath::Abs(bkg_mean[iset] - bkg_mean[jset])>gr_mean_bkg_amp[in]){
	      gr_mean_bkg_amp[in] = TMath::Abs(bkg_mean[iset]-bkg_mean[jset]);	      
	    }

	    /// RELATIVE DISTANCE
	    if(/*gr_copy[in]==8 &&*/ jset>iset){
	      clx_rel_dist = (gr_x_pol_bar[iset]-gr_x_pol_bar[jset])*1000.; //nm
	      cly_rel_dist = (gr_y_pol_bar[iset]-gr_y_pol_bar[jset])*1000.; //nm
	      if(abs(clx_rel_dist)==0 && abs(clx_rel_dist)==0) gr_same_pos[in]++; 
		//cout << iset << " " << jset << " " << gr_x_pol_bar[iset] << " " << gr_x_pol_bar[jset] << endl;
	      if(gr_copy[in]==8 && (clx_rel_dist!=0 && cly_rel_dist!=0)){
		hclx->Fill(clx_rel_dist);
		hcly->Fill(cly_rel_dist);
	      }

	      //// delta phi
	      //cout << iset << " " << jset << " " << index_2pol << " " << ndelta << endl;
	      delta_phi_pol[index_2pol] = TMath::Abs(gr_phi_pol[iset]-gr_phi_pol[jset]);
	      if(delta_phi_pol[index_2pol]>(TMath::Pi()/2.)) delta_phi_pol[index_2pol] = TMath::Pi()- TMath::Abs(delta_phi_pol[index_2pol]);
	      index_2pol++;
	    }
	    
	  }
	}

	
	///////////////////////////////////////////////////////
	
	////////////////////// MAX AMPLITUDE    (distanza massima di un bfcl dalla retta che collega i bfc che danno origine alla massima distanza)
	for(int iset = 0;iset<index_pol;iset++){
	   gr_amp[iset]=TMath::Abs((gr_y_pol_bar[iset]) - (m_line[in]*gr_x_pol_bar[iset] + q_line[in]))/(TMath::Sqrt(1+m_line[in]*m_line[in]));
	  if(TMath::Abs((gr_y_pol_bar[iset]) - (m_line[in]*gr_x_pol_bar[iset] + q_line[in]))/(TMath::Sqrt(1+m_line[in]*m_line[in]))>gr_max_amp[in]) gr_max_amp[in]=TMath::Abs((gr_y_pol_bar[iset]) - (m_line[in]*gr_x_pol_bar[iset] + q_line[in]))/(TMath::Sqrt(1+m_line[in]*m_line[in]));
	  //cout << in << " " << gr_maxpol1[in] << " " << gr_maxpol2[in] << " " << gr_pol_id[iset] << " " << q_line[in] << " " << m_line[in] << " " << gr_x_pol_bar[iset] << " " << gr_y_pol_bar[iset] << " " << gr_max_amp[in] << endl; 
	}
	///////////////////////////////
      
	///////// LUNGHEZZA DEL PERCORSO (LINEA CHIUSA) CALCOLATO A PARTIRE DA POL 0 (ha senso solo per numero di copie uguale a 8)
	for(int pid = 0;pid<npol; pid++){
	  gr_step[pid]=-1;
	  for(int iset = 0;iset<index_pol;iset++){
	    if(gr_pol_id[iset]==pid && match==0){   // solo per match=0 (ovvero il punto di origine che può essere anche diverso da pol=0°)
	      pre_step_x=gr_x_pol_bar[iset];
	      pre_step_y=gr_y_pol_bar[iset];
	      //cout << pid << " " << gr_pol_id[iset] << " " << gr_path[in] << endl;
	      first_step_x=gr_x_pol_bar[iset];
	      first_step_y=gr_y_pol_bar[iset];
	      match++;
	    }
	    if(gr_pol_id[iset]==pid && match>0 && pid>0){ /// punti successivi
	      gr_path[in] += TMath::Sqrt(TMath::Power(gr_x_pol_bar[iset]-pre_step_x,2)+TMath::Power(gr_y_pol_bar[iset]-pre_step_y,2));
	      if(gr_copy[in]==npol)gr_step[pid] = TMath::Sqrt(TMath::Power(gr_x_pol_bar[iset]-pre_step_x,2)+TMath::Power(gr_y_pol_bar[iset]-pre_step_y,2));
	      pre_step_x=gr_x_pol_bar[iset];
	      pre_step_y=gr_y_pol_bar[iset];
	      //cout << pid << " " << gr_pol_id[iset] << " " << gr_path[in] << " " << match <<  endl;
	      match++;
	    }
	  }
	  //cout << pid <<  " " << gr_path[in] << " " << match <<  endl;
	  if(match==npol){
	    gr_path[in]+=TMath::Sqrt(TMath::Power(first_step_x-pre_step_x,2)+TMath::Power(first_step_y-pre_step_y,2));  // chiudo la linea
	    if(gr_copy[in]==npol)gr_step[0] = TMath::Sqrt(TMath::Power(first_step_x-pre_step_x,2)+TMath::Power(first_step_y-pre_step_y,2));
	    //cout << pid <<  " " << gr_path[in] << " " << match <<  " " << gr_step[0] << endl;
	  }
	  //cout << pid <<  " " << gr_path[in] << " " << match <<  " " << gr_step[pid] << endl;
	}
	//////////////////////////////////////////// end 

	
	gr_phi_mean[in] = TMath::Mean(ndelta,delta_phi_pol);
	gr_phi_rms[in] = TMath::RMS(ndelta,delta_phi_pol);
	gr_phi_max[in] = TMath::MaxElement(ndelta,delta_phi_pol);
	gr_mean_path[in] = TMath::Mean(npol,gr_step);
	gr_rms_path[in] = TMath::RMS(npol,gr_step);
	gr_max_path[in] = TMath::MaxElement(npol,gr_step);
	gr_rms_amp[in] = TMath::RMS(index_pol,gr_amp);
	gr_mean_amp[in] = TMath::Mean(index_pol,gr_amp);

	delete [] gr_amp;
	delete [] gr_step;
	delete [] delta_phi_pol;
	delete [] gr_phi_pol;
	delete [] gr_x_pol;
	delete [] gr_y_pol;
	delete [] gr_z_pol;
	delete [] gr_min_pol;
	delete [] gr_maj_pol;
	delete [] gr_ell_pol;
	delete [] gr_npx_pol;
	delete [] gr_x_pol_bar;
	delete [] gr_y_pol_bar;
	delete [] gr_npx_pol_bar;
	delete [] gr_vol_pol_bar;
	delete [] gr_bright_pol_bar;
	delete [] gr_angpol;
	delete [] gr_pol_id;
	delete [] bkg_mean;
	delete [] sig_peak;
	delete [] br_amp;
	delete [] gr_xy_dist_bar;
	
	/////// OUTPUT ////////////////
	mybfcl << "View: "<< viewID << " - Grain ID: " << in << " - # Elements: " << index_pol << endl;
	mybfcl << "bfcl("<< hID << "," << viewID << "," << in << ",40,";
	if(gr_copy[in]==npol)bfcl8 << "bfcl("<< hID << "," << viewID << "," << in << ",40,";

	if(gr_copy[in]==npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated)yandex <<  hID << ","  << viewID << "," << in << ",";
	
	if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && gr_max_dist_bar[in]>cut_bar_l && gr_max_dist_bar[in]<cut_bar_u &&  gr_path[in]/gr_max_dist_bar[in]>cut_reg_l && gr_path[in]/gr_max_dist_bar[in]<cut_reg_u)cut8 << "bfcl("<< hID << "," << viewID << "," << in << ",40,";
	for(int k=0;k<npol;k++){
	  if(k<7){
	    mybfcl <<sort_gr[k] << ",";
	    if(gr_copy[in]==npol)bfcl8 <<sort_gr[k] << ",";
	    if(gr_copy[in]==npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated)yandex <<sort_gr[k] << ",";
	    if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && gr_max_dist_bar[in]>cut_bar_l && gr_max_dist_bar[in]<cut_bar_u &&  gr_path[in]/gr_max_dist_bar[in]>cut_reg_l && gr_path[in]/gr_max_dist_bar[in]<cut_reg_u)cut8 <<sort_gr[k] << ",";
	  }
	  if(k==7){
	    mybfcl <<sort_gr[k] << ")";
	    if(gr_copy[in]==npol)bfcl8 <<sort_gr[k] << ")";
	    if(gr_copy[in]==npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated)yandex <<sort_gr[k];
	    if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && gr_max_dist_bar[in]>cut_bar_l && gr_max_dist_bar[in]<cut_bar_u &&  gr_path[in]/gr_max_dist_bar[in]>cut_reg_l && gr_path[in]/gr_max_dist_bar[in]<cut_reg_u)cut8 <<sort_gr[k] << ")";
	  }
	}
	mybfcl << endl;
	if(gr_copy[in]==npol)bfcl8 << endl;
        if(gr_copy[in]==npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated)yandex << endl;
        if(gr_copy[in]==cut_npol && gr_lx[in]/gr_ly[in]!=cut_nofit && grain[in]==cut_goodzone && gr_ly[in]>=cut_minor && gr_ncl[in]<cut_ncl && gr_imt[in]==cut_mtrk && gr_isolated[in]==cut_isolated && gr_max_dist_bar[in]>cut_bar_l && gr_max_dist_bar[in]<cut_bar_u &&  gr_path[in]/gr_max_dist_bar[in]>cut_reg_l && gr_path[in]/gr_max_dist_bar[in]<cut_reg_u)cut8 << endl;
	///////////////////////////////////

	//cout << gr_max_br_amp[in] << " " << gr_max_peak_amp[in] << " " << gr_mean_bkg_amp[in] << endl;
      }
      ////////////////// END 3
     
      /////////////////// PARTE 4: CARATTERIZZAZIONE DEI SINGOLI GRANI RICOSTRUITI ////////////////////////////////////////////////////////////////////
      
      if(flag==0){
      for(int in=0;in<gr_;in++){

	//// INIZIALIZATIONS
	tmp_gr_dist=0;
	tmp_rdist=100;
	tmp_rdist2=100;
	tmp_rdist_ldust=100;
	gr_npeaks_max_phi[in]=-5;

	Int_t npk_index=0;
	Bool_t npeaks_ok=false;
	//if(gr_npeaks_npol[in]>2)npk_index=(gr_npeaks_npol[in]*(gr_npeaks_npol[in]-1))/2;
	if((num_tot_peaks[in]-gr_npeaks_npol[in])>1)npk_index=(num_tot_peaks[in]-gr_npeaks_npol[in])*(num_tot_peaks[in]-gr_npeaks_npol[in]-1)/2;
	Float_t * delta_npeaks_phi = new Float_t[npk_index];
	Float_t * gr_npeaks_phi = new Float_t[npk_index];	
	Int_t npeaks_index=0;
	Int_t npeaks_index2=0;
	if(gr_npeaks_npol[in]==0)delta_npeaks_phi[0]=-5;
	if(gr_npeaks_npol[in]==1)delta_npeaks_phi[0]=0;
	////////////////////////
      
	if(grain[in] && gr_lx[in]!=gr_ly[in] && gr_ncl[in]<=30 && gr_ly[in]>=0.100  && gr_isolated[in]==-1){ // SOLO GRANI ISOLATI CHE PASSANO I CUT
	  /// MINIMA DISTANZA 
	  for(int kn=0;kn<gr_;kn++){
	    if(grain[kn] && gr_lx[kn]!=gr_ly[kn] && gr_ncl[kn]<=30 && gr_ly[kn]>=0.100 &&  gr_isolated[kn]==-1 && kn!=in){ // SOLO GRANI ISOLATI CHE PASSANO I CUT
	      bfcl_rdist=TMath::Sqrt(TMath::Power(gr_x[kn]-gr_x[in],2)+TMath::Power(gr_y[kn]-gr_y[in],2)+TMath::Power(gr_z[kn]-gr_z[in],2));
	      if(bfcl_rdist<tmp_rdist){
		tmp_rdist=bfcl_rdist;
	      }	    
	    }	    	    
	  }

	  
	  for(int kbin=0;kbin<totXbin;kbin++){
	    for(int jbin=0;jbin<totYbin;jbin++){
	      if(gr_x[in]>(-totXbin/2. + kbin) && gr_x[in]<(-totXbin/2. + kbin+1) && gr_y[in]>(-totYbin/2. + jbin) && gr_y[in]<(-totYbin/2. + jbin+1)){
		new_phi_bar = TMath::ATan(1./TMath::Tan(phi_set_bar[in]));
		phi_projX[kbin][jbin]+=TMath::Sin(new_phi_bar)*gr_max_dist_bar[in];
		phi_projY[kbin][jbin]+=TMath::Sin(phi_set_bar[in])*gr_max_dist_bar[in];
		phi_projX_adi[kbin][jbin]+=TMath::Sin(new_phi_bar);
		phi_projY_adi[kbin][jbin]+=TMath::Sin(phi_set_bar[in]);
		ngr_cell[kbin][jbin]++;
	      }
	    }
	  }

	  hrdist->Fill(tmp_rdist);

	  
	}
	//////////////////////end

	if(grain[in] && gr_lx[in]!=gr_ly[in]  && gr_isolated[in]==2 && gr_ncl[in]<=30 && gr_ly[in]>=0.100 && gr_imt[in]==-1){ // GRANI CON DUE PICCHI (MTRK) CHE SUPERANO I CUT
	  /// DISTANZA TRA I DUE PICCHI
	  for(int jn=0;jn<cl_;jn++){
	    if(cl_igr[jn]==gr_id[in] &&  cl_ifr[jn]==cl_ifr[ipol_gr[in][cl_ipol[jn]]] && cl_flags[jn]==0 && jn!=ipol_gr[in][cl_ipol[jn]]){ // GRANI CON DUE PICCHI (MTRK) CHE SUPERANO I CUT 
	      bfcl_rdist2=TMath::Sqrt(TMath::Power(cl_x2[ipol_gr[in][cl_ipol[jn]]]-cl_x2[jn],2)+TMath::Power(cl_y2[ipol_gr[in][cl_ipol[jn]]]-cl_y2[jn],2));	      
	      if(bfcl_rdist2<tmp_rdist2){
		tmp_rdist2=bfcl_rdist2;
		phiang2 = TMath::ATan((cl_y2[ipol_gr[in][cl_ipol[jn]]]-cl_y2[jn])/(cl_x2[ipol_gr[in][cl_ipol[jn]]]-cl_x2[jn]));
		npeak_vol_dif = (cl_vol[ipol_gr[in][cl_ipol[jn]]] - cl_vol[jn]);   // differenza volume tra i due picchi
		npeak_npx_dif = (cl_npx[ipol_gr[in][cl_ipol[jn]]] - cl_npx[jn]);   // differenza area tra i due picchi
		npeak_bri_dif = ((double)cl_vol[ipol_gr[in][cl_ipol[jn]]]/cl_npx[ipol_gr[in][cl_ipol[jn]]] - (double)cl_vol[jn]/cl_npx[jn]);  // differenza luminosità tra i due picchi
		//cout << phiang2 << endl;
	      }
	      //cout << "debug in progress "<< npk_index << " " << gr_npeaks_npol[in] << " " << num_tot_peaks[in] << " " << npeaks_index << endl;
	      gr_npeaks_phi[npeaks_index] = TMath::ATan((cl_y2[ipol_gr[in][cl_ipol[jn]]]-cl_y2[jn])/(cl_x2[ipol_gr[in][cl_ipol[jn]]]-cl_x2[jn]));
	      npeaks_index++;
	    }
	  }

	  /////// Distanza angolare su npeaks elements
	  //npeaks_index=0;
	  for(int kn=0;kn<npeaks_index;kn++){
	    for(int ln=0;ln<npeaks_index;ln++){
	      if(ln > kn){
		delta_npeaks_phi[npeaks_index2] = TMath::Abs(gr_npeaks_phi[kn] - gr_npeaks_phi[ln]);
		if(delta_npeaks_phi[npeaks_index2]>(TMath::Pi()/2.)) delta_npeaks_phi[npeaks_index2] = TMath::Pi()- TMath::Abs(delta_npeaks_phi[npeaks_index2]);
		npeaks_index2++;
	      }
	    }
	  }
	  gr_npeaks_max_phi[in] = TMath::MaxElement(npeaks_index2,delta_npeaks_phi);
	  cout << jentry << " " << in << " " <<  gr_npeaks_max_phi[in] << endl;
	  delete [] gr_npeaks_phi;
	  delete [] delta_npeaks_phi;
	  /////////////////////////////////////////

	  
	  hrdist2->Fill(bfcl_rdist2);
	  twopeak_dist[in]=bfcl_rdist2;
	  twopeak_phi[in]=phiang2;
	  twopeak_dvol[in]=TMath::Abs(npeak_vol_dif);
	  twopeak_dnpx[in]=TMath::Abs(npeak_npx_dif);
	  twopeak_dbri[in]=TMath::Abs(npeak_bri_dif);
	  //cout << in << " "  << twopeak_dist[in] << endl;
	}
	////////////////// end
	
	//////  MIN DIST BETWEEN LARGE DUST AND NEARBY LARGE DUST
	if(ldust[in] && gr_lx[in]!=gr_ly[in]){  /// prendo solo i large dust
	  ///  DISTANZA DAL LARGE DUST  (espressa in unità di raggio medio del grano)
	  for(int jn=0;jn<gr_;jn++){
	    if(!ldust[jn] &&  gr_lx[jn]!=gr_ly[jn]) { ///// tutti i grani tranne i large dust
	      bfcl_rdist_ldust=(TMath::Sqrt(TMath::Power(gr_x[jn]-grain_Ox[in],2)+TMath::Power(gr_y[jn]-grain_Oy[in],2)))/(TMath::Sqrt(grain_area[in]/TMath::Pi()));//*((x_pix_size-y_pix_size)/2.));
	      //bfcl_rdist_ldust=(TMath::Sqrt(TMath::Power(gr_x[jn]-grain_Ox[in],2)+TMath::Power(gr_y[jn]-grain_Oy[in],2)))/(TMath::Sqrt((gr_npx[in]/gr_ncl[in])/(TMath::Pi()))*((x_pix_size-y_pix_size)/2.));
	      //cout << (TMath::Sqrt(TMath::Power(gr_x[jn]-grain_Ox[in],2)+TMath::Power(gr_y[jn]-grain_Oy[in],2))) << " " << (TMath::Sqrt((gr_npx[in]/gr_ncl[in])/(TMath::Pi()))*((x_pix_size-y_pix_size)/2.)) << endl;
	      hrdist_ldust->Fill(bfcl_rdist_ldust);
	      // cout << jentry << " " << ent << endl;
	    }
	  }
	}
	//////////////////// end
	//cout << in << " "  << twopeak_dist[in] << " " << twopeak_phi[in] << endl;
	
      } /// relative distance

      //////////////////////////

      vector<double> min_mtrk_dist;
      bool dist_not_taken = false;
      
      for(int in=0;in<mt_;in++){
	tmp_mtrk_dist=1000;
	dist_not_taken=false;
	
	if(mt_id[in]!=-1){
	  for(int kn=in;kn<mt_;kn++){
	    if(mt_id[kn]!=-1 && in!=kn ){ 
	      mtrk_rdist[0]=TMath::Sqrt(TMath::Power(mt_limX0[kn]-mt_limX0[in],2)+TMath::Power(mt_limY0[kn]-mt_limY0[in],2)+TMath::Power(mt_limZ0[kn]-mt_limZ0[in],2));
	      mtrk_rdist[1]=TMath::Sqrt(TMath::Power(mt_limX1[kn]-mt_limX0[in],2)+TMath::Power(mt_limY1[kn]-mt_limY0[in],2)+TMath::Power(mt_limZ1[kn]-mt_limZ0[in],2));
	      mtrk_rdist[2]=TMath::Sqrt(TMath::Power(mt_limX0[kn]-mt_limX1[in],2)+TMath::Power(mt_limY0[kn]-mt_limY1[in],2)+TMath::Power(mt_limZ0[kn]-mt_limZ1[in],2));
	      mtrk_rdist[3]=TMath::Sqrt(TMath::Power(mt_limX1[kn]-mt_limX1[in],2)+TMath::Power(mt_limY1[kn]-mt_limY1[in],2)+TMath::Power(mt_limZ1[kn]-mt_limZ1[in],2));
	      dist_mtrk = TMath::MinElement(4,mtrk_rdist);
	      //cout << jentry << " " << aid << " " << mt_id[in] << " " << mt_id[kn] << " " << dist_mtrk << endl;
	      hmtrk->Fill(dist_mtrk);
	      if(dist_mtrk<tmp_mtrk_dist){
		tmp_mtrk_dist=dist_mtrk;
		//tmp_mtrk_match=gr_imt[kn];
	      }	    
	    }	    	    
	  }
	  for(int jn=0;jn<min_mtrk_dist.size();jn++){
	    if(tmp_mtrk_dist==min_mtrk_dist.at(jn)) dist_not_taken=true;
	  }
	  if(!dist_not_taken && tmp_mtrk_dist!=1000){
	    min_mtrk_dist.push_back(tmp_mtrk_dist);
	    //hmtrk->Fill(tmp_mtrk_dist);
	    //cout << jentry << " " << aid << " " << mt_id[in] << " " << tmp_mtrk_dist << endl;
	  }

	  //for(int kn=0;kn<mt_;kn++){
	  //if(mt_id[kn]!=-1 && in!=kn && TMath::Abs(TMath::ATan(TMath::Tan(mt_phi[kn]))-TMath::ATan(TMath::Tan(mt_phi[in])))<0.79){
	  //  
	  //}
	  
	  
	}
      }
      
    } // flag ==0
      ///////////////////////////// END 4

      /////////////////// PARTE 5: DATASET (costruzione del Tree)  ////////////////////////////////////////////////////////////////////////////////////////////////
     
      for(int in = 0;in<gr_; in++){   /// indice sui grani
	  int clset[npol]={};
	  for(int jn=0; jn<npol;jn++){     /// indice sulle polarizzazioni
	    eChannel=channel;             // for smart comparison
	    eHeaderID=id;
	    eViewID=viewID;
	    eFlag=flag;
	    eGrainID=in;
	    eNgr=ngr;
	    ePolID=jn;
	    ePuls=gr_ncl[in];
	    eBfcPolID=ipol_gr[in][jn];
	    eBfcID=gr_ibfc[in];
	    eBfcGap=bfc_gap[in];
	    eNcl=ncl;
	    eNclFr=fr_ncl[cl_ifr[gr_ibfc[in]]];
	    eGrainx=gr_x[in];
	    eGrainy=gr_y[in];
	    eGrainz=gr_z[in];
	    eClustx=cl_x_pos[in][jn];
	    eClusty=cl_y_pos[in][jn];
	    eClustz=cl_z_pos[in][jn];
	    eClustMin=cl_min_bf[in][jn];
	    eClustMaj=cl_maj_bf[in][jn];
	    eClustEll=cl_ell_bf[in][jn];
	    eClustPhi=cl_phi_bf[in][jn];
	    eClustMeanBkg=cl_bkg_mean[in][jn];
	    eClustMaxPeak=cl_sig_peak[in][jn];
	    eSetBrAmp=gr_max_br_amp[in];
	    eSetBkgAmp=gr_mean_bkg_amp[in];
	    eSetPeakAmp=gr_max_peak_amp[in];
	    eSetBrMaxPol=gr_max_br_pol[in];
	    eSetBrMinPol=gr_min_br_pol[in];
	    eGrainMin=gr_ly[in];
	    eGrainMaj=gr_lx[in];
	    eGrainEll=gr_lx[in]/gr_ly[in];
	    if(gr_phi[in]<=(TMath::Pi()/2.))eGrainPhi=gr_phi[in];
	    if(gr_phi[in]>(TMath::Pi()/2.) && gr_phi[in]<=(3*TMath::Pi()/2.))eGrainPhi=gr_phi[in] - TMath::Pi();
	    if(gr_phi[in]>(3*TMath::Pi()/2.) && gr_phi[in]<=(2*TMath::Pi()))eGrainPhi=gr_phi[in] - 2*TMath::Pi();
	    eGrainTheta=gr_theta[in];
	    eVolume=gr_vol[in];
	    eBfcVolume=cl_vol[gr_ibfc[in]];
	    eFrBfVolume=vol_frbf[in][jn];
	    eArea=gr_npx[in];
	    eBfcArea=cl_npx[gr_ibfc[in]];
	    eFrBfArea=npx_frbf[in][jn];
	    eGrainArea=grain_area[in];
	    eXView=x;
	    eYView=y;
	    eZ=(zmin+zmax)*0.5;
	    eZlen=gr_lz[in];
	    eIsolated=gr_isolated[in];
	    eNFrame=gr_frLast[in]-gr_frFirst[in]+1;
	    eSetFrame=cl_ifr[gr_ibfc[in]];
	    eSetNCopy=gr_copy[in];
	    eSetNStatic=gr_same_pos[in];
	    eSetGapZ=gr_gap_z[in];
	    eSetGapZ2=gr_gap_z2[in];
	    //if(viewID==726 && (ipol_gr[in][jn]!=-1) && gr_180[in]!=-1)cout << eSetGap << " " << cl_y2[gr_180[in]] << " " << cl_y2[ipol_gr[in][jn]] << " " << cl_x2[gr_180[in]] << " " << cl_x2[ipol_gr[in][jn]] << endl;
	    if(ipol_gr[in][jn]!=-1)eSetGap=TMath::Sqrt(TMath::Power(cl_y2[gr_180[in]] - cl_y2[ipol_gr[in][jn]],2) + TMath::Power(cl_x2[gr_180[in]] - cl_x2[ipol_gr[in][jn]],2));
	    //cout << cl_y[gr_180[in]] << " " <<  cl_y[ipol_gr[in][jn]] << " " << cl_x[gr_180[in]] << " " <<  cl_x[ipol_gr[in][jn]] << endl;
	    eSetPath=gr_path[in];
	    eSetMeanPath=gr_mean_path[in];
	    eSetRmsPath=gr_rms_path[in];
	    eSetMaxPath=gr_max_path[in];
	    eSetMaxDist=gr_max_dist[in];
	    eSetNpeaksDist=twopeak_dist[in];
	    eSetNpeaksPhi=twopeak_phi[in];
	    eSetNpeaksDVol=twopeak_dvol[in];
	    eSetNpeaksDNpx=twopeak_dnpx[in];
	    eSetNpeaksDBri=twopeak_dbri[in];
	    eSetNpeaks=gr_npeaks[in]/gr_copy[in];
	    eSetNpeaksNcopy=gr_npeaks_npol[in];
	    eSetNpeaksMaxPhiAmp = gr_npeaks_max_phi[in];
	    eSetNpeaksMax = TMath::MaxElement(npol,num_peaks[in]);
	    eSetMaxPol1=gr_maxpol1[in];
	    eSetMaxPol2=gr_maxpol2[in];
	    eSetMaxBar=gr_max_dist_bar[in];
	    eSetMaxAmp=gr_max_amp[in];
	    eSetRmsAmp=gr_rms_amp[in];
	    eSetMeanAmp=gr_mean_amp[in];
	    eSetXRms=gr_x_rms[in];
	    eSetYRms=gr_y_rms[in];
	    eSetXBar=gr_x_mean[in];
	    eSetYBar=gr_y_mean[in];
	    eSetXMaxBar=gr_x_maxbar[in];
	    eSetYMaxBar=gr_y_maxbar[in];
	    eSetXMinBar=gr_x_minbar[in];
	    eSetYMinBar=gr_y_minbar[in];
	    eSetPhiRms=gr_phi_rms[in];
	    eSetPhiMean=gr_phi_mean[in];
	    eSetPhiMaxAmp=gr_phi_max[in];
	    eSetPhi=phi_set[in];
	    eSetPhiBar=phi_set_bar[in];
	    eSetTheBar=the_set_bar[in];
	    eSetMeanBright=gr_vol_mean_bar[in];
	    eSetNpxRms=gr_npx_rms_bar[in];
	    eSetVolRms=gr_vol_rms_bar[in];
	    if(gr_vol_ref[in]!=0)eSetVolRatio=vol_frbf[in][jn]/gr_vol_ref[in];
	    
	    eDeltaPhi = eSetPhiBar-eGrainPhi;
	    if(TMath::Abs(eDeltaPhi)>(TMath::Pi()/2.) && eDeltaPhi>0) eDeltaPhi = TMath::Pi()/2. - TMath::Abs(eDeltaPhi);
	    if(TMath::Abs(eDeltaPhi)>(TMath::Pi()/2.) && eDeltaPhi<0) eDeltaPhi = -(TMath::Pi()/2. - TMath::Abs(eDeltaPhi));
	  
	  
	    for(int kn=0;kn<gr_;kn++){
	      if(grain[in] && grain[kn] && kn!=in){
		rdist=TMath::Sqrt(TMath::Power(gr_x[in]-gr_x[kn],2)+TMath::Power(gr_y[in]-gr_y[kn],2)+TMath::Power(gr_z[in]-gr_z[kn],2));
		xydist=TMath::Sqrt(TMath::Power(gr_x[in]-gr_x[kn],2)+TMath::Power(gr_y[in]-gr_y[kn],2));
		hdist=TMath::Abs(gr_z[in]-gr_z[kn]);
		//hrdist->Fill(rdist);
		if(rdist<2){                // user setting
		ctd[in]=true;
		ctd[kn]=true;
		}
	      }
	    } 
	    eLargeDust=ldust[in];
	    eGoodZone=grain[in];
	    eShadow=shadow[in];
	    eClustTypeDust=ctd[in];

	    if(TMath::Abs(eGrainMaj*TMath::Cos(eGrainPhi))>TMath::Abs(eGrainMin*TMath::Sin(eGrainPhi)))eEllPrjX = TMath::Abs(eGrainMaj*TMath::Cos(eGrainPhi));
	    else eEllPrjX = TMath::Abs(eGrainMin*TMath::Sin(eGrainPhi));
	    if(TMath::Abs(eGrainMin*TMath::Cos(eGrainPhi))>TMath::Abs(eGrainMaj*TMath::Sin(eGrainPhi)))eEllPrjY = TMath::Abs(eGrainMin*TMath::Cos(eGrainPhi));
	    else eEllPrjY = TMath::Abs(eGrainMaj*TMath::Sin(eGrainPhi));

	    if(eGrainEll!=1 && eGoodZone && eGrainMin>0.1 && eGrainMin<0.2 && ePolID==0 && eFlag==0){
	      images_ellipse(aRun,eHeaderID,eGrainID,0,fit_info,c1);
	      eGrFitMin=fit_info[0];
	      eGrFitMaj=fit_info[1];
	      eGrFitEll=fit_info[1]/fit_info[0];
	      eGrFitPhi=fit_info[2];
	    }
	    else{
	      eGrFitMin=-1;
	      eGrFitMaj=-1;
	      eGrFitEll=-1;
	      eGrFitPhi=-1;
	    }
	    
	    //////// MICROTRACK //////////////////////
	    eMTrk=-1;
	    eMTGr=-1;
	    eMTNfr=-1;
	    eMTPhi=phi_mt[in];
	    eMTThe=0;
	    eMTLen=-1;
	    eMTChi2=-1;
	  eMTVolDif=mt_dif_vol[in];
	  eMTID=gr_imt[in];
      
	  if(gr_imt[in]!=-1){
	    if(!mtrk[gr_imt[in]]){
	      eMTrk=1;
	      eMTGr=mt_ngr[gr_imt[in]];
	      eMTNfr=mt_frLast[gr_imt[in]];	  	  
	      //if(mt_phi[gr_imt[in]]<=TMath::Pi())eMTPhi=mt_phi[gr_imt[in]];
	      //else eMTPhi=mt_phi[gr_imt[in]]-2*TMath::Pi();
	      eMTThe=mt_theta[gr_imt[in]];
	      eMTLen=mt_len[gr_imt[in]];
	      eMTChi2=mt_chi2[gr_imt[in]];
	      mtrk[gr_imt[in]]=true;	  
	    }
	    else eMTrk=11;
	  }

	  if(ipol_gr[in][jn]!=-1 && gr_copy[in]==npol && eGoodZone){
	    //scaling(aRun,viewID,ipol_gr[in][jn],jn,scale_result);
	    eBfcMeanBkg = scale_result[0];
	    eBfcSigPeak = scale_result[1];
	    eSclMeanBkg = scale_result[2];
	    eSclSigPeak = scale_result[3];	
	    if(jn==0){
	      scale_ref = eBfcSigPeak/eBfcMeanBkg;
	    }
	    eBfcWeigth = eBfcSigPeak/eBfcMeanBkg/scale_ref;
	    //cout << scale_result[0] << " " << scale_result[1] << endl;
	  }
	  else {
	    eBfcMeanBkg=-1;
	    eBfcSigPeak=-1;
	    eSclMeanBkg=-1;
	    eSclSigPeak=-1;
	    eBfcWeigth=-1;

	  }
      
	  if(jn==0 && eGoodZone && gr_copy[in]==npol){
	    TH2F* h_sum = new TH2F(Form("hSum_%d_%d",viewID,in),"hSum",100,0,50,100,0,50);
	    h_sum->Reset(0);
	    for(int kn=0;kn<npol;kn++){
	      clset[kn]=ipol_gr[in][kn];
	    }
	    //h_sum = merged_histo(aRun,view,viewID,in,clset);
	    //gr_spectrum(h_sum,gr_scale_result);
	    eGrMeanBkg = gr_scale_result[0];
	    eGrSigPeak = gr_scale_result[1];
	    delete h_sum;
	  }
	  Tree_out->Fill();      
	}
      }
      ////////////////// END 5

      delete [] cl_x2;
      delete [] cl_y2;
    
    } //end crash

    //// CUT AREA /////
    for(int i=0;i<hcutarea->GetNbinsX();i++){
      for(int j=0;j<hcutarea->GetNbinsY();j++){
	if(hcutarea->GetBinContent(i+1,j+1)!=0)area_cut++;
      }
    }
    cout << "area_tagliata "<<eViewID << " " << area_cut/10000. << endl;
  } // end su tutte le View
  
  for(int kbin=0;kbin<totXbin;kbin++){
    for(int jbin=0;jbin<totYbin;jbin++){
      phi_projY[kbin][jbin]/=ngr_cell[kbin][jbin];
      phi_projX[kbin][jbin]/=ngr_cell[kbin][jbin];
      phi_projY_adi[kbin][jbin]/=ngr_cell[kbin][jbin];
      phi_projX_adi[kbin][jbin]/=ngr_cell[kbin][jbin];
      radius=TMath::Sqrt(TMath::Power(phi_projY[kbin][jbin],2)+TMath::Power(phi_projX[kbin][jbin],2));
      if(ngr_cell[kbin][jbin]!=0){
	//cout << phi_projY[kbin][jbin] << " " << phi_projX[kbin][jbin] << endl;
	hxmap->SetBinContent(kbin+1,jbin+1,phi_projX_adi[kbin][jbin]);
	hymap->SetBinContent(kbin+1,jbin+1,phi_projY_adi[kbin][jbin]);
	hphimap->SetBinContent(kbin+1,jbin+1,TMath::ATan(phi_projY[kbin][jbin]/phi_projX[kbin][jbin]));
	hphicell->Fill(TMath::ATan(phi_projY[kbin][jbin]/phi_projX[kbin][jbin]));	
	hrphimap->Fill(TMath::ATan(phi_projY[kbin][jbin]/phi_projX[kbin][jbin]),radius);
	hradius->SetBinContent(kbin+1,jbin+1,radius);
      }
    }
  }

  ofstream log_run("log_run.txt");
  log_run << "area_scannata [um]^2: " << area_scan << endl;
  log_run << "area_tagliata [um]^2: " << area_cut/10000. << endl;
  log_run << "area_effettiva [um^2]: "<< area_scan - area_cut/1000. << endl;
  log_run.close();

  f_out->cd();
  
  grNcl->Write();
  hcutarea->Write();
  hrdist->Write();
  hrdist2->Write();
  hrdist_ldust->Write();
  hmtrk->Write();
  hxmap->Write();
  hymap->Write();
  hradius->Write();
  hphimap->Write();
  hrphimap->Write();
  hphicell->Write();
  hclx->Write();
  hcly->Write();
  Tree_out->Write();
  mybfcl.close();
  bfcl8.close();
  yandex.close();
  cut8.close();
  f_out->Close();
  
} // end Loop



int myrun(){
  vv.Loop();  
}




