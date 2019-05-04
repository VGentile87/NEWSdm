#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

#include "myData.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include <vector>
#include <iostream>

using namespace std;
const double pi=TMath::Pi();


//DMRView        *view;
//DMRImageCl     *image;

const char *file="dm_tracks_cl.dm.root";
DMRRun *aRun = new DMRRun(file);
//aRun->SetFixEncoderFaults(1);
DMRView *view = aRun->GetView();

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
Double_t fit_cluster[9]={};

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
Int_t ctd[kMaxgr]={};
Float_t cut_ctd[kMaxgr]={};
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
Float_t  gr_npeaks_dist_rms[kMaxgr]={};
Float_t  gr_npeaks_max_phi[kMaxgr]={};
Float_t  gr_npeaks_mean_phi[kMaxgr]={};
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
//Double_t phi_set_fit[kMaxgr]={};
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
Double_t mt_dif_br[kMaxgr]={};
Double_t mt_dif_npx[kMaxgr]={};
Double_t mt_dif_z[kMaxgr]={};
Int_t gr_mtrk[kMaxgr]={};
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



Double_t cl_fit_min[kMaxgr][npol]={};
Double_t cl_fit_maj[kMaxgr][npol]={};
Double_t cl_fit_ell[kMaxgr][npol]={};
Double_t cl_fit_phi[kMaxgr][npol]={};
Double_t cl_fit_x[kMaxgr][npol]={};
Double_t cl_fit_y[kMaxgr][npol]={};
Double_t gr_max_dist_fit[kMaxgr]={};
Double_t phi_set_fit[kMaxgr]={};
Double_t gr_fit_x_mean[kMaxgr]={};
Double_t gr_fit_y_mean[kMaxgr]={};
Double_t gr_fit_phi_mean[kMaxgr]={};
/////////////////////////

Int_t tcl_gr;
Double_t *x_ldust;
Double_t *y_ldust;
Double_t *z_ldust; 
Double_t *lx_ldust;
Double_t *ly_ldust; 
Double_t *phi_ldust;
Double_t *area_ldust;
TString chain;
ostringstream frame;
Double_t tmp_vol[npol];
Int_t bfc_ifr[npol];
Int_t bfc_zfr[npol];
Int_t tmp_same_frame[npol];


int imt;
int tmp_mt_id;
int tmp_br_mt;
Double_t *br_mt; 
Double_t *npx_mt;
Double_t *x_mt;
Double_t *y_mt;
Double_t *z_mt;


void first_initializer(int dim, int npol, int *gr_imt)
{
  for(int in=0; in<dim;in++){
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
    mt_dif_br[in]=-1;
    mt_dif_npx[in]=-1;
    mt_dif_z[in]=-1;
    gr_mtrk[in]=-1;
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
    gr_mean_bkg_amp[in]=-1;
    
    //cout << gr_mean_bkg_amp[in] << endl;
  }  
}


void second_initializer(int in)
{
  tmp_gap=-1;
  tmp_gr=-1;    
  dist_cl_gr[in]=0;
  grain[in]=true;
  //ctd[in]=0;
  ldust[in]=false;
  bfc_border[in]=false;
  index_clust=0;      
  n_ifr=1;
}


void gr_proc_definitions(int in, UInt_t *gr_ncl)
{
  tcl_gr = gr_ncl[in]*npol;       
  x_ldust = new Double_t[tcl_gr];
  y_ldust = new Double_t[tcl_gr];
  z_ldust = new Double_t[tcl_gr];
  lx_ldust = new Double_t[tcl_gr];
  ly_ldust = new Double_t[tcl_gr];
  phi_ldust = new Double_t[tcl_gr];
  area_ldust = new Double_t[tcl_gr];     
  
  for(int i=0;i<npol;i++){
    tmp_vol[i]=-1;
    bfc_zfr[i]=-1;
    tmp_same_frame[i]=-1;
  }
}

void gr_delete_obj()
{
  delete [] x_ldust;
  delete [] y_ldust;
  delete [] z_ldust;
  delete [] lx_ldust;
  delete [] ly_ldust;
  delete [] phi_ldust;
  delete [] area_ldust;
}


void mt_initialization(int mt_ngr)
{
  imt=0;
  tmp_mt_id=-1;
  tmp_br_mt=-1;
  br_mt = new Double_t[mt_ngr];
  npx_mt = new Double_t[mt_ngr];
  x_mt = new Double_t[mt_ngr];
  y_mt = new Double_t[mt_ngr];
  z_mt = new Double_t[mt_ngr];
}

void mt_delete_obj()
{
  delete [] br_mt;
  delete [] npx_mt;
  delete [] x_mt;
  delete [] y_mt;
  delete [] z_mt;
}
