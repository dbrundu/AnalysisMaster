/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/

/*
 * BreitWignerLineShapes.h
 *
 *  Created on: 13/11/2018
 *      Author: Davide Brundu
 */

#ifndef BREITWIGNERLINESHAPES_H_
#define BREITWIGNERLINESHAPES_H_

#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/detail/utility/CheckValue.h>
#include <hydra/Parameter.h>
#include <hydra/Tuple.h>
#include <hydra/Complex.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/BlattWeisskopfFunctions.h>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <cmath>

// Insert it in the namespace "hydra"
namespace hydra {


/**
 * \ingroup common_functions
 *
 * \class BreitWignerLineShapes
 *
 * Breit-Wigner line shape for 4 body double-resonant decays \f$ A -> r1 r2 , r1-> a b , r2 -> c d \f$ ,
 * where r1 and r2 are the intermediate resonances and \f$ a, b\f$ and \f$c ,d \f$ are the final states.
 * The Lineshapes are relativistic BreitWigners.
 *
 * \f[

 *  R(m_{a,b}|m_1,\Lambda_1 , m_{c,d}|m_2,\Lambda_2) =  BW(m_{a,b}|m_1,\Lambda_1) \times BW(m_{c,d}|m_2,\Lambda_2) 
 *                                                      B'_{L_r}(d, q_0, q)(\frac{q}{q_r})^{L_r}|_{r1} \times 
 *                                                      B'_{L_r}(d, q_0, q)(\frac{q}{q_r})^{L_r}|_{r2} \times
 *                                                      B'_{L_A}(d, p_0, p)(\frac{p}{m_A})^{L_A}|_{mother}
 * 
 * \f]
 *
 * where Breit-Wigner amplitude is given by:
 *
 *\f[ BW(m_{ab}|m_0,\Lambda_0)= \frac{1}{m_0^2 - m_{ab}^2 - im_0\Lambda(m_{ab})} \f]
 *
 *and
 *\f[  \Lambda(m_{ab}) = \Lambda_0(\frac{q}{q_0})^{2L_{r}+1}\frac{m_0}{m}B'_{L_r}(d, q_0, q)\f]
 *
 *@tparam Resonance1Wave hydra::Wave resonance decay vertex wave
 *@tparam MotherWave hydra::Wave mother particle decay vertex wave
 *
 * @param mass_R1,2 resonances mass. -> hydra::Parameter to be fitted
 * @param width_R1,2 resonance width. -> hydra::Parameter to be fitted
 * @param mother_mass resonances mother mass. -> hydra::Parameter to be fitted
 * @param daugther1_mass mass of daugther 1 of resonance 1
 * @param daugther2_mass mass of daugther 2 of resonance 1
 * @param daugther3_mass mass of daugther 1 of resonance 2
 * @param daugther4_mass mass of daugther 2 of resonance 2
 * @param res1_mass_nominal : nominal mass of resonance 1 (to be used as bachelor quantity)
 * @param res2_mass_nominal : nominal mass of resonance 2 (to be used as bachelor quantity)
 * @param radi decay vertex radio.


 */
template<Wave Resonance1Wave, Wave Resonance2Wave, Wave MotherWave=SWave, unsigned int ArgIndex=0>
class BreitWignerLineShapes : public BaseFunctor<BreitWignerLineShapes< Resonance1Wave,Resonance2Wave,MotherWave,ArgIndex>, hydra::complex<double>, 6>
{
        using BaseFunctor<BreitWignerLineShapes<Resonance1Wave,Resonance2Wave,MotherWave,ArgIndex>, hydra::complex<double>, 6>::_par;

private:
        double fDaughter1Mass;
        double fDaughter2Mass;
        double fDaughter3Mass;
        double fDaughter4Mass;
        double fRes1NominalMass;
        double fRes2NominalMass;
        double fMotherMass;
        double fRadi;


public:

        BreitWignerLineShapes()=delete;


        BreitWignerLineShapes(hydra::Parameter const& c_re,    hydra::Parameter const& c_im,
                          hydra::Parameter const& mass_R1, hydra::Parameter const& width_R1,
                          hydra::Parameter const& mass_R2, hydra::Parameter const& width_R2,
                          double mother_mass,              
                          double daugther1_mass,     double daugther2_mass,
                          double daugther3_mass,     double daugther4_mass,
                          double res1_mass_nominal,  double res2_mass_nominal, //used as bachelor masses in BW
                          double radi):
                BaseFunctor<BreitWignerLineShapes<Resonance1Wave,Resonance2Wave,MotherWave,ArgIndex>, hydra::complex<double>, 6>{c_re, c_im, mass_R1, width_R1, mass_R2, width_R2},
                fDaughter1Mass(daugther1_mass),
                fDaughter2Mass(daugther2_mass),
                fDaughter3Mass(daugther3_mass),
                fDaughter4Mass(daugther4_mass),
                fRes1NominalMass(res1_mass_nominal),
                fRes2NominalMass(res2_mass_nominal),
                fMotherMass(mother_mass),
                fRadi(radi)
        {}

        __hydra_host__  __hydra_device__
        BreitWignerLineShapes(BreitWignerLineShapes<Resonance1Wave,Resonance2Wave,MotherWave,ArgIndex>  const& other):
                BaseFunctor<BreitWignerLineShapes<Resonance1Wave,Resonance2Wave,MotherWave,ArgIndex>, hydra::complex<double>, 6>(other),
                fDaughter1Mass(other.GetDaughter1Mass()),
                fDaughter2Mass(other.GetDaughter2Mass()),
                fDaughter3Mass(other.GetDaughter3Mass()),
                fDaughter4Mass(other.GetDaughter4Mass()),
                fRes1NominalMass(other.GetRes1NominalMass()),
                fRes2NominalMass(other.GetRes2NominalMass()),
                fMotherMass(other.GetMotherMass()),
                fRadi(other.GetRadi())
                {}

        __hydra_host__  __hydra_device__
        BreitWignerLineShapes<Resonance1Wave,Resonance2Wave,MotherWave,ArgIndex>&
        operator=(BreitWignerLineShapes<Resonance1Wave,Resonance2Wave,MotherWave,ArgIndex>  const& other)
        {
                if(this==&other) return  *this;

                BaseFunctor<BreitWignerLineShapes<Resonance1Wave,Resonance2Wave,MotherWave,ArgIndex>,
                        hydra::complex<double>, 6>::operator=(other);

                fDaughter1Mass   = other.GetDaughter1Mass();
                fDaughter2Mass   = other.GetDaughter2Mass();
                fDaughter3Mass   = other.GetDaughter3Mass();
                fDaughter4Mass   = other.GetDaughter4Mass();
                fRes1NominalMass = other.GetRes1NominalMass();
                fRes2NominalMass = other.GetRes2NominalMass();
                fMotherMass= other.GetMotherMass();
                fRadi= other.GetRadi();

                 return  *this;
        }

        __hydra_host__  __hydra_device__ inline
        double GetDaughter1Mass() const {
                return fDaughter1Mass;
        }

        __hydra_host__  __hydra_device__ inline
        void SetDaughter1Mass(double daughter1Mass) {
                fDaughter1Mass = daughter1Mass;
        }

        __hydra_host__  __hydra_device__ inline
        double GetDaughter2Mass() const {
                return fDaughter2Mass;
        }

        __hydra_host__  __hydra_device__ inline
        void SetDaughter2Mass(double daughter2Mass) {
                fDaughter2Mass = daughter2Mass;
        }
        
        __hydra_host__  __hydra_device__ inline
        double GetDaughter3Mass() const {
                return fDaughter3Mass;
        }

        __hydra_host__  __hydra_device__ inline
        void SetDaughter3Mass(double daughter3Mass) {
                fDaughter3Mass = daughter3Mass;
        }
        
        
        __hydra_host__  __hydra_device__ inline
        double GetDaughter4Mass() const {
                return fDaughter4Mass;
        }

        __hydra_host__  __hydra_device__ inline
        void SetDaughter4Mass(double daughter4Mass) {
                fDaughter4Mass = daughter4Mass;
        }

        __hydra_host__  __hydra_device__ inline
        double GetRes1NominalMass() const {
                return fRes1NominalMass;
        }

        __hydra_host__  __hydra_device__ inline
        void SetRes1NominalMass(double Res1NominalMass) {
                fRes1NominalMass = Res1NominalMass;
        }
        
        
        __hydra_host__  __hydra_device__ inline
        double GetRes2NominalMass() const {
                return fRes2NominalMass;
        }

        __hydra_host__  __hydra_device__ inline
        void SetRes2NominalMass(double Res2NominalMass) {
                fRes2NominalMass = Res2NominalMass;
        }

        __hydra_host__  __hydra_device__ inline
        double GetMotherMass() const {
                return fMotherMass;
        }

        __hydra_host__  __hydra_device__ inline
        void SetMotherMass(double motherMass) {
                fMotherMass = motherMass;
        }

        __hydra_host__  __hydra_device__ inline
        double GetRadi() const {
                return fRadi;
        }

        __hydra_host__  __hydra_device__ inline
        void SetRadi(double radi) {
                fRadi = radi;
        }


        __hydra_host__ __hydra_device__ inline
        hydra::complex<double> Evaluate(unsigned int n, hydra::Vector4R* p) const {
        
        hydra::Vector4R p1 = p[0];
        hydra::Vector4R p2 = p[1];
        hydra::Vector4R p3 = p[2];
        hydra::Vector4R p4 = p[3];

        const double inv_mass_12 = (p1+p2).mass() ;
        const double inv_mass_34 = (p3+p4).mass() ;

        const double coeff_re   = _par[0];
        const double coeff_im   = _par[1];
        const double res1_mass  = _par[2];
        const double res1_width = _par[3];
        const double res2_mass  = _par[4];
        const double res2_width = _par[5];

        return LineShapes(inv_mass_12, inv_mass_34, coeff_re, coeff_im, res1_mass, res1_width, res2_mass, res2_width);

        }
        
        template<typename T>
        __hydra_host__ __hydra_device__ inline
        hydra::complex<double> Evaluate(T p) const {
        
        hydra::Vector4R p1 = get<0>(p);
        hydra::Vector4R p2 = get<1>(p);
        hydra::Vector4R p3 = get<2>(p);
        hydra::Vector4R p4 = get<3>(p);

        const double inv_mass_12 = (p1+p2).mass() ;
        const double inv_mass_34 = (p3+p4).mass() ;

        const double coeff_re   = _par[0];
        const double coeff_im   = _par[1];
        const double res1_mass  = _par[2];
        const double res1_width = _par[3];
        const double res2_mass  = _par[4];
        const double res2_width = _par[5];

        return LineShapes(inv_mass_12, inv_mass_34, coeff_re, coeff_im, res1_mass, res1_width, res2_mass, res2_width);

        }



private:

         template<Wave ResWave>
           __hydra_host__ __hydra_device__  inline
         double  Width( const double m, const double resonance_mass, const double resonance_width,
                           const double  p0, const double  p) const {

                 const double  B = BarrierFactor<ResWave>( fRadi, p0,  p);

                   return resonance_width*\
                                  pow<double, 2*ResWave+1>(p/p0)*\
                                  (resonance_mass/m)*\
                                  B*B;

           }

         __hydra_host__ __hydra_device__   inline
         hydra::complex<double> LineShapes(const double m12, 
                                           const double m34, 
                                           const double coeff_re,
                                           const double coeff_im,
                                           const double res1_mass_par,
                                           const double res1_width_par,
                                           const double res2_mass_par, 
                                           const double res2_width_par ) const {

                 const double p0 = pmf( fMotherMass,   res1_mass_par,  fRes2NominalMass); // mother -> Res1 Res2 (Res1 as resonance, Res2 as bachelor !)
                 const double q0 = pmf( res1_mass_par, fDaughter1Mass, fDaughter2Mass);   // Res1   -> dau1 dau2
                 const double r0 = pmf( res2_mass_par, fDaughter3Mass, fDaughter4Mass);   // Res2   -> dau3 dau4

                 const double p  = pmf( fMotherMass, m12,            fRes2NominalMass);  // mother -> Res1 Res2 (Res1 as resonance, Res2 as bachelor !)
                 const double q  = pmf( m12,         fDaughter1Mass, fDaughter2Mass);    // Res1   -> dau1 dau2
                 const double r  = pmf( m34,         fDaughter3Mass, fDaughter4Mass);    // Res2   -> dau3 dau4

                 const double  B_R1     = BarrierFactor<Resonance1Wave>( fRadi, q0,  q);
                 const double  B_R2     = BarrierFactor<Resonance2Wave>( fRadi, r0,  r);
                 const double  B_mother = BarrierFactor<MotherWave>( fRadi, p0,  p);

                 const double width_1 = Width<Resonance1Wave>( m12, res1_mass_par, res1_width_par, q0, q);
                 const double width_2 = Width<Resonance2Wave>( m34, res2_mass_par, res2_width_par, r0, r);

                 hydra::complex<double> coeff(coeff_re, coeff_im);

                 hydra::complex<double> numerator(  B_R1 * pow<double,Resonance1Wave>(q/q0) *\
                                                    B_R2 * pow<double,Resonance2Wave>(r/r0) *\
                                                    B_mother * pow<double,MotherWave>(p/p0), 0);

                 hydra::complex<double> denominator_res1(m12*m12 - res1_mass_par*res1_mass_par,  -res1_mass_par*width_1);
                 hydra::complex<double> denominator_res2(m34*m34 - res2_mass_par*res2_mass_par,  -res2_mass_par*width_2);

                 return hydra::complex<double>( (coeff*numerator) / (denominator_res1*denominator_res2) ) ;

         }


};



}  // namespace hydra


#endif /* BREITWIGNERLINESHAPES_H_ */
