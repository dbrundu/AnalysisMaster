/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/

/*
 * TheoreticalFourBodyAngularDist.h
 *
 *  Created on: 16/11/2018
 *      Author: Davide Brundu
 */

#ifndef THEORETICALFOURBODYANGULARDISTRIBUTION_H_
#define THEORETICALFOURBODYANGULARDISTRIBUTION_H_

#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/detail/utility/CheckValue.h>
#include <hydra/Parameter.h>
#include <hydra/Tuple.h>
#include <hydra/Complex.h>
#include <hydra/functions/Math.h>
#include <hydra/functions/Utils.h>
//#include <hydra/functions/CosHelicityAngle.h>
//#include <hydra/functions/PlanesDeltaAngle.h>
//#include <hydra/functions/BlattWeisskopfFunctions.h>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <cmath>


//Put it inside hydra namespace for convenience
namespace hydra {


class TheoreticalFourBodyAngularDist: public hydra::BaseFunctor<TheoreticalFourBodyAngularDist, double, 12>
{

using BaseFunctor<TheoreticalFourBodyAngularDist, double, 12>::_par;
typedef enum {I_1=0, I_2, I_3, I_4, I_5, I_6, I_7, I_8, I_9} I_num;

public:

    TheoreticalFourBodyAngularDist() = delete;

    TheoreticalFourBodyAngularDist(hydra::Parameter const& HL_0_re,    hydra::Parameter const& HL_0_im,
                                   hydra::Parameter const& HL_par_re,  hydra::Parameter const& HL_par_im,
                                   hydra::Parameter const& HL_perp_re, hydra::Parameter const& HL_perp_im,
                                   hydra::Parameter const& HR_0_re,    hydra::Parameter const& HR_0_im,
                                   hydra::Parameter const& HR_par_re,  hydra::Parameter const& HR_par_im,
                                   hydra::Parameter const& HR_perp_re, hydra::Parameter const& HR_perp_im):
              BaseFunctor<TheoreticalFourBodyAngularDist, hydra::complex<double>, 12>{HL_0_re, HL_0_im, HL_par_re, HL_par_im, 
                                                                                      HL_perp_re, HL_perp_im, HR_0_re, HR_0_im,
                                                                                      HR_par_re, HR_par_im, HR_perp_re, HR_perp_im}

    __hydra_dual__
    TheoreticalFourBodyAngularDist(const TheoreticalFourBodyAngularDist& other):
              BaseFunctor<TheoreticalFourBodyAngularDist, double, 12>(other)
    {}

    __hydra_dual__  
    TheoreticalFourBodyAngularDist&
    operator=( TheoreticalFourBodyAngularDist const& other)
    {
        if(this==&other) return *this;
        BaseFunctor<TheoreticalFourBodyAngularDist, double, 12>::operator=(other);
        return *this;
    }
    



    __hydra_dual__  inline
    double Evaluate(unsigned int n, hydra::Vector4R* p)  const 
    {

        hydra::Vector4R p1 = p[0];
        hydra::Vector4R p2 = p[1];
        hydra::Vector4R p3 = p[2];
        hydra::Vector4R p4 = p[3];

        const double  theta_pi = ::acos(decay_angle( (p1+p2+p3+p4), (p1+p2), p1 ));
        const double  theta_l  = ::acos(decay_angle( (p1+p2+p3+p4), (p3+p4), p3 ));
        const double  chi      = chi_angle(p3, p4, p1);
        
        const hydra::complex<double>    HL_0(_par[0] , _par[1]);
        const hydra::complex<double>  HL_par(_par[2] , _par[3]);
        const hydra::complex<double> HL_perp(_par[4] , _par[5]);
        const hydra::complex<double>    HR_0(_par[6] , _par[7]);
        const hydra::complex<double>  HR_par(_par[8] , _par[9]);
        const hydra::complex<double> HR_perp(_par[10], _par[11]);
        
        double result = 0;
        for(I_num i=I_1; i<=I_9; i++)
        {
           result += C_fun(theta_l, chi, i) *\
                     I_fun(HL_0, HL_par, HL_perp, HR_0, HR_par, HR_perp, theta_pi, (I_num)i);
        }

        return result;
    }
    
    
    template<typename T>
    __hydra_dual__  inline
    double Evaluate(T& p)  const 
    {

        hydra::Vector4R p1 = get<0>(p);
        hydra::Vector4R p2 = get<1>(p);
        hydra::Vector4R p3 = get<2>(p);
        hydra::Vector4R p4 = get<3>(p);
        
        const double  theta_pi = ::acos(decay_angle( (p1+p2+p3+p4), (p1+p2), p1 ));
        const double  theta_l  = ::acos(decay_angle( (p1+p2+p3+p4), (p3+p4), p3 ));
        const double  chi      = chi_angle(p3, p4, p1);
        
        const hydra::complex<double>    HL_0(_par[0] , _par[1]);
        const hydra::complex<double>  HL_par(_par[2] , _par[3]);
        const hydra::complex<double> HL_perp(_par[4] , _par[5]);
        const hydra::complex<double>    HR_0(_par[6] , _par[7]);
        const hydra::complex<double>  HR_par(_par[8] , _par[9]);
        const hydra::complex<double> HR_perp(_par[10], _par[11]);
        
        double result = 0;
        for(int i=I_1; i<=I_9; i++)
        {
           result += C_fun(theta_l, chi, i) *\
                     I_fun(HL_0, HL_par, HL_perp, HR_0, HR_par, HR_perp, theta_pi, (I_num)i);
        }

        return result;
        
    }




private:

    __hydra_dual__ inline
    double decay_angle(Vector4R const& p, Vector4R const& q, Vector4R const& d) const 
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
    
    
    __hydra_dual__ inline
    double chi_angle(Vector4R const& d2, Vector4R const& d3, Vector4R const& h1) const 
    {
        hydra::Vector4R D = d2 + d3;
     
        hydra::Vector4R d1_perp = d2 - (D.dot(d2) / D.dot(D)) * D; // d2 will be mu^+
        hydra::Vector4R h1_perp = h1 - (D.dot(h1) / D.dot(D)) * D;

        // orthogonal to both D and d1_perp
        hydra::Vector4R d1_prime = D.cross(d1_perp);

        d1_perp  = d1_perp / d1_perp.d3mag();
        d1_prime = d1_prime / d1_prime.d3mag();in 

        double x, y;

        x =  d1_perp.dot(h1_perp);
        y = d1_prime.dot(h1_perp);
      
        return ::atan2(y, x);
    }
    

    
    __hydra_dual__ inline
    double I_fun(const hydra::complex<double>& HL_0,
                                 const hydra::complex<double>& HL_par,
                                 const hydra::complex<double>& HL_perp,
                                 const hydra::complex<double>& HR_0, 
                                 const hydra::complex<double>& HR_par,
                                 const hydra::complex<double>& HR_perp,
                                 const double& theta_pi, I_num N){
        double Res = 0;
        switch(N)
        {
            case I_1: 
                Res = (1./16.)*(hydra::norm(HL_0) + hydra::norm(HR_0) + 1.5*pow(sin(theta_pi),2) *\
                      (hydra::norm(HL_perp) + hydra::norm(HL_par) + hydra::norm(HR_perp) + hydra::norm(HR_par) ) ); break;
            case I_2: 
                Res = (-1./16.)*(hydra::norm(HL_0) + hydra::norm(HR_0) - 0.5*pow(sin(theta_pi),2) *\
                      (hydra::norm(HL_perp) + hydra::norm(HL_par) + hydra::norm(HR_perp) + hydra::norm(HR_par) ) ); break;
            case I_3: 
                Res = (1./16.)*pow(sin(theta_pi),2)*(hydra::norm(HL_perp) - hydra::norm(HL_par) + hydra::norm(HR_perp) - hydra::norm(HR_par)); break; 

            case I_4: 
                Res = (-1./8.)*sin(theta_pi)*( (HL_0*hydra::conj(HL_par)).real() + (HR_0*hydra::conj(HR_par)).real() ); break; 
                
            case I_5: 
                Res = (-1./4.)*sin(theta_pi)*( (HL_0*hydra::conj(HL_perp)).real() - (HR_0*hydra::conj(HR_perp)).real() ); break; 
                
            case I_6: 
                Res = (1./4.)*pow(sin(theta_pi),2)*( (HL_par*hydra::conj(HL_perp)).real() - (HR_par*hydra::conj(HR_perp)).real() ); break; 
                
            case I_7: 
                Res = (-1./4.)*sin(theta_pi)*( (HL_0*hydra::conj(HL_par)).imag() - (HR_0*hydra::conj(HR_par)).imag() ); break; 
                
            case I_8: 
                Res = (-1./8.)*sin(theta_pi)*( (HL_0*hydra::conj(HL_perp)).imag() + (HR_0*hydra::conj(HR_perp)).imag() ); break; 

            case I_9: 
                Res = (1./8.)*pow(sin(theta_pi),2)*( (HL_perp*hydra::conj(HL_par)).imag() + (HR_perp*hydra::conj(HR_par)).imag() ); break; 
        }
        
        return Res;
    };
    
    
    
    __hydra_dual__ inline
    double C_fun(const double& theta_l, const double& chi, I_num N){
    
        double Res = 0;
        
        switch(N)
        {
            case I_1: 
                Res = 1; break;
            case I_2: 
                Res = cos(2*theta_l); break;
            case I_3: 
                Res = pow(sin(theta_l),2) * cos(2*chi); break; 
            case I_4: 
                Res = sin(2*theta_l)*cos(chi); break; 
            case I_5: 
                Res = sin(theta_l)*cos(chi); break; 
            case I_6: 
                Res = cos(theta_l); break; 
            case I_7: 
                Res = sin(theta_l)*sin(chi); break; 
            case I_8: 
                Res = sin(2*theta_l)*sin(chi); break; 
            case I_9: 
                Res = pow(sin(theta_l),2)*sin(2*chi); break; 
        }
        
        return Res;

    
    }
    




};

} // namespace hydra

#endif /* THEORETICALFOURBODYANGULARDISTRIBUTION_H_ */
