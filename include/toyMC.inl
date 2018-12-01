/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/
/*
 * main.inl
 *
 *  Created on: 
 *      Author: 
 */

#ifndef TOYMC_INL_
#define TOYMC_INL_

#include <external/std_libs.h>
#include <external/Hydra_libs.h>
#include <external/ROOT_libs.h>
#include <tclap/CmdLine.h>
#include <toyMC/generate_dataset.h>



void inline parse_cmdline(int argc, char *argv[], size_t& nentries){
    try {

        TCLAP::CmdLine cmd("Command line arguments for ", '=');

        TCLAP::ValueArg<size_t> EArg("n", "number-of-events","Number of events", true, 10e6, "size_t");
        cmd.add(EArg);

        // Parse the argv array.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        nentries = EArg.getValue();

    }
    catch (TCLAP::ArgException &e)  {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}



int main(int argc, char** argv)
{

    size_t nentries = 0;
    
    parse_cmdline(argc, argv, nentries);
    
    const double D0_mass    = 1864.83;
    const double mu_mass    = 105.658;  
    const double pi_mass    = 139.57061;
    const double K_mass     = 493.677;
    
    
    // Allocate memory to hold the final states particles
    hydra::Decays<4, hydra::host::sys_t> Events_h;
    
    toyMC::generate_dataset(hydra::device::sys, 
                            {D0_mass, pi_mass, mu_mass}, 
                            Events_h,
                            nentries, 
                            2*nentries,
                            true);

    
    std::cout << std::endl << "<======= Generated Dataset =======>"<< std::endl;
    for( size_t i=0; i<5; i++ ) std::cout << Events_h[i] << std::endl;
    
    
    // Smearing of D_M data for fitting
    hydra::host::vector<double>  D_M_data_h(Events_h.size());
    hydra::Random<> Generator( std::chrono::system_clock::now().time_since_epoch().count() );
    Generator.Gauss(D0_mass, 10.0, D_M_data_h.begin(), D_M_data_h.end());
    
}





#endif /* TOYMC_INL_ */

