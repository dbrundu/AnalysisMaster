/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/

/******************************
*  Utils 
******************************/

namespace toyMC {

inline double cos_decay_angle(hydra::Vector4R const& p, hydra::Vector4R const& q, hydra::Vector4R const& d){

  hydra::GReal_t pd = p*d;
  hydra::GReal_t pq = p*q;
  hydra::GReal_t qd = q*d;
  hydra::GReal_t mp2 = p.mass2();
  hydra::GReal_t mq2 = q.mass2();
  hydra::GReal_t md2 = d.mass2();

  return (pd * mq2 - pq * qd)
    / ::sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));
  }
        
        
/******************************
******************************/
        
inline double cos_decay_angle(TLorentzVector p, TLorentzVector q, TLorentzVector d) {

  double pd = p*d;
  double pq = p*q;
  double qd = q*d;
  double mp2 = p.Mag2();
  double mq2 = q.Mag2();
  double md2 = d.Mag2();

  return (pd * mq2 - pq * qd)
    / sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));

  }
        
/******************************
******************************/
        
inline double chi_plane_angle(hydra::Vector4R& d2, hydra::Vector4R& d3, hydra::Vector4R& h1) {

  hydra::Vector4R D = d2 + d3;
     
  hydra::Vector4R d1_perp = d2 - (D.dot(d2) / D.dot(D)) * D;
  hydra::Vector4R h1_perp = h1 - (D.dot(h1) / D.dot(D)) * D;

  // orthogonal to both D and d1_perp
  hydra::Vector4R d1_prime = D.cross(d1_perp);

  d1_perp  = d1_perp / d1_perp.d3mag();
  d1_prime = d1_prime / d1_prime.d3mag();

  hydra::GReal_t x, y;

  x = d1_perp.dot(h1_perp);
  y = d1_prime.dot(h1_perp);

  hydra::GReal_t chi = ::atan2(y, x);
      

  return chi;
}

/******************************
******************************/

inline double chi_plane_angle(TLorentzVector d2, TLorentzVector d3, TLorentzVector h1) {

  TVector3 d22 = d2.Vect();
  TVector3 d33 = d3.Vect();
  TVector3 h11 = h1.Vect();
  
  TVector3 D = d22 + d33;
     
  TVector3 d1_perp = d22 - (D.Dot(d22) / D.Dot(D)) * D;
  TVector3 h1_perp = h11 - (D.Dot(h11) / D.Dot(D)) * D;

  // orthogonal to both D and d1_perp
  TVector3 d1_prime = D.Cross(d1_perp);        // identify (normal to) the d2-d3 plane
  d1_perp = d1_perp * (1. / d1_perp.Mag());    // stay in the d2-d3 plane
  d1_prime = d1_prime * (1. / d1_prime.Mag());
  double x, y;
  x = d1_perp.Dot(h1_perp);
  y = d1_prime.Dot(h1_perp);
  double chi = atan2(y, x);
      
  return chi;
}

/******************************
******************************/

void BuildTree(TString dirname_string, TChain* treedata)
{
  const char *dirname    = dirname_string.Data();
  TSystemDirectory dir(dirname, dirname);

  TList *files = dir.GetListOfFiles();

  TIter next(files);

  TSystemFile *file;
  TString fname;
  const char *ext=".root";
  //cout << "  - Adding ntuples from "+dirname_string <<  endl;
  int k=0;
  while ((file=(TSystemFile*)next())) {
    fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
         treedata->Add(dirname_string+fname);
         k++;
      }
  }
  //return treedata;
}

/******************************
******************************/

void addDecays(hydra::Decays<4, hydra::device::sys_t > &Events_d, TChain *DecayTree){

  double D_M;
  double Dst_DTF_DiLepton_Mass_D0constr;
  double Dst_DTF_DiHadron_Mass_D0constr;
  double Dst_DTF_h0_PX_D0constr;
  double Dst_DTF_h0_PY_D0constr;
  double Dst_DTF_h0_PZ_D0constr;
  double Dst_DTF_h0_E_D0constr;
  double Dst_DTF_h1_PX_D0constr; 
  double Dst_DTF_h1_PY_D0constr;
  double Dst_DTF_h1_PZ_D0constr;
  double Dst_DTF_h1_E_D0constr;
  double Dst_DTF_l0_PX_D0constr; 
  double Dst_DTF_l0_PY_D0constr;
  double Dst_DTF_l0_PZ_D0constr;
  double Dst_DTF_l0_E_D0constr;
  double Dst_DTF_l1_PX_D0constr; 
  double Dst_DTF_l1_PY_D0constr;
  double Dst_DTF_l1_PZ_D0constr;
  double Dst_DTF_l1_E_D0constr;
  double weight_event;
  
  DecayTree->SetBranchAddress("D_M", &D_M);
  DecayTree->SetBranchAddress("Dst_DTF_DiLepton_Mass_D0constr" , &Dst_DTF_DiLepton_Mass_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_DiHadron_Mass_D0constr" , &Dst_DTF_DiHadron_Mass_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_h0_PX_D0constr"         , &Dst_DTF_h0_PX_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_h0_PY_D0constr"         , &Dst_DTF_h0_PY_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_h0_PZ_D0constr"         , &Dst_DTF_h0_PZ_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_h0_E_D0constr"          , &Dst_DTF_h0_E_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_h1_PX_D0constr"         , &Dst_DTF_h1_PX_D0constr); 
  DecayTree->SetBranchAddress("Dst_DTF_h1_PY_D0constr"         , &Dst_DTF_h1_PY_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_h1_PZ_D0constr"         , &Dst_DTF_h1_PZ_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_h1_E_D0constr"          , &Dst_DTF_h1_E_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_l0_PX_D0constr"         , &Dst_DTF_l0_PX_D0constr); 
  DecayTree->SetBranchAddress("Dst_DTF_l0_PY_D0constr"         , &Dst_DTF_l0_PY_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_l0_PZ_D0constr"         , &Dst_DTF_l0_PZ_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_l0_E_D0constr"          , &Dst_DTF_l0_E_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_l1_PX_D0constr"         , &Dst_DTF_l1_PX_D0constr); 
  DecayTree->SetBranchAddress("Dst_DTF_l1_PY_D0constr"         , &Dst_DTF_l1_PY_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_l1_PZ_D0constr"         , &Dst_DTF_l1_PZ_D0constr);
  DecayTree->SetBranchAddress("Dst_DTF_l1_E_D0constr"          , &Dst_DTF_l1_E_D0constr);
  DecayTree->SetBranchAddress("weight_event"                   , &weight_event);
  
  // =========>> Loop on events files
  for (Long64_t ievt=0; ievt<DecayTree->GetEntries(); ievt++) 
  {
  
    DecayTree->GetEntry(ievt);
    
    hydra::Vector4R pi_0(Dst_DTF_h0_E_D0constr, 
                         Dst_DTF_h0_PX_D0constr,
                         Dst_DTF_h0_PY_D0constr,
                         Dst_DTF_h0_PZ_D0constr);
                         
    hydra::Vector4R pi_1(Dst_DTF_h1_E_D0constr, 
                         Dst_DTF_h1_PX_D0constr,
                         Dst_DTF_h1_PY_D0constr,
                         Dst_DTF_h1_PZ_D0constr);
                         
    hydra::Vector4R mu_0(Dst_DTF_l0_E_D0constr, 
                         Dst_DTF_l0_PX_D0constr,
                         Dst_DTF_l0_PY_D0constr,
                         Dst_DTF_l0_PZ_D0constr);
                         
    hydra::Vector4R mu_1(Dst_DTF_l1_E_D0constr, 
                         Dst_DTF_l1_PX_D0constr,
                         Dst_DTF_l1_PY_D0constr,
                         Dst_DTF_l1_PZ_D0constr);
                         
    const hydra::Vector4R std_decay[4] = {pi_0, pi_1, mu_0, mu_1};
    
    Events_d.AddDecay(weight_event, std_decay);

  }

}

}
