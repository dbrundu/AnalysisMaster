/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/
 
 /*
 * generated_dataset.h
 *
 *  Created on: 19/11/2018
 *      Author: Davide Brundu
 */
 
#ifndef TOYMC_GENERATE_DATASET_H_
#define TOYMC_GENERATE_DATASET_H_

#include <external/std_libs.h>
#include <external/Hydra_libs.h>
#include <external/ROOT_libs.h>

#include <functions/BreitWignerLineShapes.h>
#include <functions/FourBodyAngularDist.h>
#include <functions/ConvertJohnsonParamters.h>


// toy MC utils
// contains functions in toyMC namespace
#include <toyMC/utils.h>


    static inline
    double decay_angle(hydra::Vector4R const& p, hydra::Vector4R const& q, hydra::Vector4R const& d) 
    {
        double pd = p*d;
        double pq = p*q;
        double qd = q*d;
        double mp2 = p.mass2();
        double mq2 = q.mass2();
        double md2 = d.mass2();

        return (pd * mq2 - pq * qd)
                / ::sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));
    }
    
    
    static inline
    double chi_angle(hydra::Vector4R const& d2, hydra::Vector4R const& d3, hydra::Vector4R const& h1) 
    {
        hydra::Vector4R D = d2 + d3;
     
        hydra::Vector4R d1_perp = d2 - (D.dot(d2) / D.dot(D)) * D; // d2 will be mu^+
        hydra::Vector4R h1_perp = h1 - (D.dot(h1) / D.dot(D)) * D;

        // orthogonal to both D and d1_perp
        hydra::Vector4R d1_prime = D.cross(d1_perp);

        d1_perp  = d1_perp / d1_perp.d3mag();
        d1_prime = d1_prime / d1_prime.d3mag();

        double x, y;

        x =  d1_perp.dot(h1_perp);
        y = d1_prime.dot(h1_perp);
      
        return ::atan2(y, x);
    }
    



namespace toyMC {


template<typename Backend, typename Container >
size_t generate_dataset(Backend const& system, 
                        const std::array<double, 3>& masses, 
                        Container& decays, 
                        const size_t nentries, 
                        const size_t bunch_size, 
                        bool Save=false) {

    const double D0_mass        = masses[0]; //1864.83
    const double h_mass         = masses[1]; //105.658
    const double mu_mass        = masses[2]; //pi=139.5706, K=493.677

    const double rho_mass   = 769.0;     
    const double phi_mass   = 1019.462;  
    const double rho_width  = 150.9;    
    const double phi_width  = 4.249;  


    hydra::Vector4R D0(D0_mass, 0.0, 0.0, 0.0);
    const double ph_masses[4]{h_mass, h_mass, mu_mass, mu_mass };

    // Create phase space object D0 -> pi pi mu mu
    hydra::PhaseSpace<4> phsp(ph_masses);
    

    auto mass_r1  = hydra::Parameter::Create().Name("MASS_RHO").Value(rho_mass).Error(0.0001).Limits(rho_mass*0.95,  rho_mass*1.05 );
    auto width_r1 = hydra::Parameter::Create().Name("WIDTH_RHO").Value(rho_width).Error(0.0001).Limits(rho_width*0.95, rho_width*1.05);
    
    auto mass_r2  = hydra::Parameter::Create().Name("MASS_PHI").Value(phi_mass).Error(0.0001).Limits(phi_mass*0.95,  phi_mass*1.05 );
    auto width_r2 = hydra::Parameter::Create().Name("WIDTH_PHI").Value(phi_width).Error(0.0001).Limits(phi_width*0.95, phi_width*1.05);
    
    auto coef_re1 = hydra::Parameter::Create().Name("A_RE1").Value(1).Error(0.001).Limits(0.95,1.05);
    auto coef_im1 = hydra::Parameter::Create().Name("A_IM1").Value(0).Error(0.001).Fixed();
    
    auto coef_re2 = hydra::Parameter::Create().Name("A_RE2").Value(5).Error(0.001).Limits(5*0.95,5*1.05);
    auto coef_im2 = hydra::Parameter::Create().Name("A_IM2").Value(0.0).Error(0.001).Fixed();

    hydra::BreitWignerLineShapes<hydra::PWave, hydra::PWave, hydra::SWave> 
        LineShape1(coef_re1, coef_im1, 
                   mass_r1, width_r1,
                   mass_r2, width_r2, 
                   D0_mass, h_mass, h_mass, mu_mass, mu_mass,
                   rho_mass, phi_mass, 
                   5.0);

    hydra::BreitWignerLineShapes<hydra::PWave,hydra::PWave,hydra::SWave> 
        LineShape2(coef_re2, coef_im2, 
                   mass_r1, width_r1,
                   mass_r1, width_r1, 
                   D0_mass, h_mass, h_mass, mu_mass, mu_mass,
                   rho_mass, rho_mass, 
                   5.0);

    double lambda_muplus  =  0.5;
    double lambda_muminus = -0.5;
    hydra::FourBodyAngularDistribution AngularDist(lambda_muplus, lambda_muminus, 1.0, 1.0);

    
    // Building complete amplitude
    auto Ampl = (LineShape1 + LineShape2)*AngularDist;
    
    // Lambda to perform |x|^2
    auto Norm = hydra::wrap_lambda( [] __hydra_dual__ ( unsigned int n, hydra::complex<double>* x){ return hydra::norm(x[0]); });
    
    auto resonance = hydra::compose(Norm, Ampl);



    std::srand((unsigned int)7531594562);

    // Allocate memory to hold the final states particles
    hydra::Decays<4, Backend > _data(bunch_size);

    auto start = std::chrono::high_resolution_clock::now();
    
    
    // Generate ph.space events in bunches and unweight until reach the correct size
    do {
        phsp.SetSeed(std::rand()); // std::chrono::system_clock::now().time_since_epoch().count() 

        //generate the final state particles
        phsp.Generate(D0, _data.begin(), _data.end());

        auto last = _data.Unweight(resonance, 1.0);

        decays.insert(decays.size()==0? decays.begin():decays.end(),
                _data.begin(), _data.begin()+last );

    } while(decays.size()<=nentries );
    
    decays.erase(decays.begin()+nentries, decays.end());

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed = end - start;

    //output
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "----- Toy MC Generation - Device --------"<< std::endl;
    std::cout << "| D0 -> h h mu mu"                        << std::endl;
    std::cout << "| Number of events :"<< nentries          << std::endl;
    std::cout << "| Time (ms)        :"<< elapsed.count()   << std::endl;
    std::cout << "-----------------------------------------"<< std::endl;





  if(Save){

    TH2D Dalitz_d2("Dalitz_d2", "Reweighted Sample [Device]; M(#mu #mu) [MeV/c^{2}]; M(#pi #pi) [MeV/c^{2}]", 
          100, 2*mu_mass, D0_mass - 2*h_mass,
          100, 2*h_mass, D0_mass - 2*mu_mass);
            
    TH1D cos_hel_pion_hist_d2("cos_hel_pion_hist_d2", "cos_hel_pion_hist_d2; cos(#theta_{#pi^{+}})", 50, -1, 1);
    TH1D cos_hel_muon_hist_d2("cos_hel_muon_hist_d2", "cos_hel_muon_hist_d2; cos(#theta_{#mu^{+}})", 50, -1, 1);
    TH1D chi_angle_hist_d2("chi_angle_hist_d2", "chi_angle_hist_d2; #phi", 50, -TMath::Pi(), TMath::Pi());
    TH1D dileptonmass_hist("dileptonmass_hist", "dileptonmass_hist; M(#mu #mu)", 100, 2*mu_mass, D0_mass - 2*h_mass);
    TH1D dihadronmass_hist("dihadronmass_hist", "dihadronmass_hist; M(#pi #pi)", 100, 2*h_mass, D0_mass - 2*mu_mass);
    
    cos_hel_pion_hist_d2.Sumw2();
    cos_hel_muon_hist_d2.Sumw2();
    chi_angle_hist_d2.Sumw2();
    dileptonmass_hist.Sumw2();
    dihadronmass_hist.Sumw2();
    
    int k_evt = 0;
    
    TRandom3 random;

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
    
    TFile toy_file("toy_file.root","RECREATE");
    TTree DecayTree("DecayTree","DecayTree");
    
    TBranch *D_M_br                            = DecayTree.Branch("D_M", &D_M);
    TBranch *Dst_DTF_DiLepton_Mass_D0constr_br = DecayTree.Branch("Dst_DTF_DiLepton_Mass_D0constr" , &Dst_DTF_DiLepton_Mass_D0constr);
    TBranch *Dst_DTF_DiHadron_Mass_D0constr_br = DecayTree.Branch("Dst_DTF_DiHadron_Mass_D0constr" , &Dst_DTF_DiHadron_Mass_D0constr);
    TBranch *Dst_DTF_h0_PX_D0constr_br         = DecayTree.Branch("Dst_DTF_h0_PX_D0constr"         , &Dst_DTF_h0_PX_D0constr);
    TBranch *Dst_DTF_h0_PY_D0constr_br         = DecayTree.Branch("Dst_DTF_h0_PY_D0constr"         , &Dst_DTF_h0_PY_D0constr);
    TBranch *Dst_DTF_h0_PZ_D0constr_br         = DecayTree.Branch("Dst_DTF_h0_PZ_D0constr"         , &Dst_DTF_h0_PZ_D0constr);
    TBranch *Dst_DTF_h0_E_D0constr_br          = DecayTree.Branch("Dst_DTF_h0_E_D0constr"          , &Dst_DTF_h0_E_D0constr);
    TBranch *Dst_DTF_h1_PX_D0constr_br         = DecayTree.Branch("Dst_DTF_h1_PX_D0constr"         , &Dst_DTF_h1_PX_D0constr);
    TBranch *Dst_DTF_h1_PY_D0constr_br         = DecayTree.Branch("Dst_DTF_h1_PY_D0constr"         , &Dst_DTF_h1_PY_D0constr);
    TBranch *Dst_DTF_h1_PZ_D0constr_br         = DecayTree.Branch("Dst_DTF_h1_PZ_D0constr"         , &Dst_DTF_h1_PZ_D0constr);
    TBranch *Dst_DTF_h1_E_D0constr_br          = DecayTree.Branch("Dst_DTF_h1_E_D0constr"          , &Dst_DTF_h1_E_D0constr);
    TBranch *Dst_DTF_l0_PX_D0constr_br         = DecayTree.Branch("Dst_DTF_l0_PX_D0constr"         , &Dst_DTF_l0_PX_D0constr);
    TBranch *Dst_DTF_l0_PY_D0constr_br         = DecayTree.Branch("Dst_DTF_l0_PY_D0constr"         , &Dst_DTF_l0_PY_D0constr);
    TBranch *Dst_DTF_l0_PZ_D0constr_br         = DecayTree.Branch("Dst_DTF_l0_PZ_D0constr"         , &Dst_DTF_l0_PZ_D0constr);
    TBranch *Dst_DTF_l0_E_D0constr_br          = DecayTree.Branch("Dst_DTF_l0_E_D0constr"          , &Dst_DTF_l0_E_D0constr);
    TBranch *Dst_DTF_l1_PX_D0constr_br         = DecayTree.Branch("Dst_DTF_l1_PX_D0constr"         , &Dst_DTF_l1_PX_D0constr);
    TBranch *Dst_DTF_l1_PY_D0constr_br         = DecayTree.Branch("Dst_DTF_l1_PY_D0constr"         , &Dst_DTF_l1_PY_D0constr);
    TBranch *Dst_DTF_l1_PZ_D0constr_br         = DecayTree.Branch("Dst_DTF_l1_PZ_D0constr"         , &Dst_DTF_l1_PZ_D0constr);
    TBranch *Dst_DTF_l1_E_D0constr_br          = DecayTree.Branch("Dst_DTF_l1_E_D0constr"          , &Dst_DTF_l1_E_D0constr);
    TBranch *weight_event_br                   = DecayTree.Branch("weight_event"                   , &weight_event);
    
    for( auto event : decays )
    {

        double weight        = hydra::get<0>(event); //weight = weight/weight_NR[k_evt];
        hydra::Vector4R pi1  = hydra::get<1>(event);
        hydra::Vector4R pi2  = hydra::get<2>(event);
        hydra::Vector4R mu1  = hydra::get<3>(event);
        hydra::Vector4R mu2  = hydra::get<4>(event);
        
        auto   util_pseudo_rndm = std::abs( pi1.get(1)-::floor(pi1.get(1)) );
        
        double M2_pipi     = (pi1+pi2).mass2();
        double M2_mumu     = (mu1+mu2).mass2();
        
        double M_pipi     = (pi1+pi2).mass();
        double M_mumu     = (mu1+mu2).mass();
        
        // Start to fill histograms
        Dalitz_d2.Fill(M_mumu, M_pipi);
        dileptonmass_hist.Fill(M_mumu);
        dihadronmass_hist.Fill(M_pipi);
        
        //if(util_pseudo_rndm>0.5){
            cos_hel_pion_hist_d2.Fill(cos_decay_angle(pi1+pi2+mu1+mu2 , pi1+pi2 , pi1)); 
            cos_hel_muon_hist_d2.Fill(cos_decay_angle(pi1+pi2+mu1+mu2 , mu1+mu2 , mu1)); 
            chi_angle_hist_d2.Fill(chi_angle(mu1,mu2,pi1));
        //}else{
        //    cos_hel_pion_hist_d2.Fill(cos_decay_angle(pi1+pi2+mu1+mu2 , pi1+pi2 , pi2)); 
        //    cos_hel_muon_hist_d2.Fill(cos_decay_angle(pi1+pi2+mu1+mu2 , mu1+mu2 , mu2)); 
        //    chi_angle_hist_d2.Fill(chi_angle(mu1,mu2,pi2));
        //}
        
        
        // Start to fill tree
        
        weight_event = weight;
        D_M          = random.Gaus((pi1+pi2+mu1+mu2).mass(), 10);
        
        Dst_DTF_DiLepton_Mass_D0constr = M_mumu;
        Dst_DTF_DiHadron_Mass_D0constr = M_pipi;
        
        Dst_DTF_h0_PX_D0constr = pi1.get(1);
        Dst_DTF_h0_PY_D0constr = pi1.get(2);
        Dst_DTF_h0_PZ_D0constr = pi1.get(3);
        Dst_DTF_h0_E_D0constr  = pi1.get(0);
        
        Dst_DTF_h1_PX_D0constr = pi2.get(1);
        Dst_DTF_h1_PY_D0constr = pi2.get(2);
        Dst_DTF_h1_PZ_D0constr = pi2.get(3);
        Dst_DTF_h1_E_D0constr  = pi2.get(0);
        
        Dst_DTF_l0_PX_D0constr = mu1.get(1); 
        Dst_DTF_l0_PY_D0constr = mu1.get(2);
        Dst_DTF_l0_PZ_D0constr = mu1.get(3);
        Dst_DTF_l0_E_D0constr  = mu1.get(0);
        
        Dst_DTF_l1_PX_D0constr = mu2.get(1); 
        Dst_DTF_l1_PY_D0constr = mu2.get(2);
        Dst_DTF_l1_PZ_D0constr = mu2.get(3);
        Dst_DTF_l1_E_D0constr  = mu2.get(0);
        
        DecayTree.Fill();
        ++k_evt;
    } //end for loop


    gStyle->SetOptStat(1);
    toy_file.cd();

    TCanvas Dalitz_canvas_d2("Dalitz_canvas_d2", "Dalitz_canvas_d2", 1200, 800);
    Dalitz_d2.Smooth();
    Dalitz_d2.Draw("COLZ");
    Dalitz_d2.Write();
    Dalitz_canvas_d2.SaveAs("GEN_Dalitz.pdf");
    
    TCanvas cos_hel_pion_canvas_d2("cos_hel_pion_canvas_d2", "cos_hel_pion_canvas_d2", 1200, 800);
    cos_hel_pion_hist_d2.Draw();
    cos_hel_pion_hist_d2.Write();
    cos_hel_pion_canvas_d2.SaveAs("GEN_costheta_pi.pdf");
    
    TCanvas cos_hel_muon_canvas_d2("cos_hel_muon_canvas_d2", "cos_hel_muon_canvas_d2", 1200, 800);
    cos_hel_muon_hist_d2.Draw();
    cos_hel_muon_hist_d2.Write();
    cos_hel_muon_canvas_d2.SaveAs("GEN_costheta_mu.pdf");
    
    TCanvas phi_canvas_d2("phi_canvas_d2", "phi_canvas_d2", 1200, 800);
    chi_angle_hist_d2.Draw();
    chi_angle_hist_d2.Write();
    phi_canvas_d2.SaveAs("GEN_chiangle2.pdf");
    
    TCanvas dileptonmass_canv("dileptonmass_canv", "dileptonmass_canv", 1200, 800);
    dileptonmass_hist.Draw();
    dileptonmass_hist.Write();
    dileptonmass_canv.SaveAs("GEN_dilepton.pdf");
    
    TCanvas dihadronmass_canv("dihadronmass_canv", "dihadronmass_canv", 1200, 800);
    dihadronmass_hist.Draw();
    dihadronmass_hist.Write();
    dihadronmass_canv.SaveAs("GEN_dihadron.pdf");

    TDirectory *cd_PiPiMuMu = toy_file.mkdir("DstD2PiPiMuMu");
    cd_PiPiMuMu->cd();
    DecayTree.Write();
    toy_file.Close();
  }// end save scope
  
  return decays.size();
}




} // end namespace

#endif /* TOYMC_GENERATE_DATASET_H_ */
