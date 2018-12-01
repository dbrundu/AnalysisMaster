/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/

/*
 * FourBodyAngularDist.h
 *
 *  Created on: 16/11/2018
 *      Author: Davide Brundu
 */

#ifndef FOURBODYANGULARDISTRIBUTION_H_
#define FOURBODYANGULARDISTRIBUTION_H_

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


class FourBodyAngularDistribution: public hydra::BaseFunctor<FourBodyAngularDistribution, hydra::complex<double>, 0>
{


public:

    FourBodyAngularDistribution() = delete;

    FourBodyAngularDistribution(double Lambda_muplus, double Lambda_muminus,
                        double j_r1, double j_r2):
              BaseFunctor<FourBodyAngularDistribution, hydra::complex<double>, 0>{},
              fLambda_muplus(Lambda_muplus),
              fLambda_muminus(Lambda_muminus),
              fJ_r1(j_r1), 
              fJ_r2(j_r2)
    {}


    __hydra_dual__
    FourBodyAngularDistribution(FourBodyAngularDistribution const& other):
              BaseFunctor<FourBodyAngularDistribution, hydra::complex<double>, 0>(other),
              fLambda_muplus(other.GetLambda_muplus()),
              fLambda_muminus(other.GetLambda_muminus()),
              fJ_r1(other.GetJr1()),
              fJ_r2(other.GetJr2())
    {}

    __hydra_dual__  
    FourBodyAngularDistribution&
    operator=( FourBodyAngularDistribution const& other)
    {
        if(this==&other) return *this;
        BaseFunctor<FourBodyAngularDistribution, hydra::complex<double>, 0>::operator=(other);
        fLambda_muplus  = other.GetLambda_muplus();
        fLambda_muminus = other.GetLambda_muminus();
        fJ_r1 = other.GetJr1();
        fJ_r2 = other.GetJr2();

        return *this;
    }
    
    __hydra_dual__ inline
    double GetJr1() const 
    {
        return fJ_r1;
    }
    
    __hydra_dual__ inline
    double GetJr2() const 
    {
        return fJ_r2;
    }
    
    __hydra_dual__ inline
    double GetLambda_muplus() const 
    {
        return fLambda_muplus;
    }
    
    __hydra_dual__ inline
    double GetLambda_muminus() const 
    {
        return fLambda_muminus;
    }
    
    __hydra_dual__ inline
    void SetLambda_muplus(double Lambda_muplus)
    {
        fLambda_muplus = Lambda_muplus;
    }
    
    __hydra_dual__ inline
    void SetLambda_muminus(double Lambda_muminus)
    {
        fLambda_muminus = Lambda_muminus;
    }



    __hydra_dual__  inline
    hydra::complex<double> Evaluate(unsigned int n, hydra::Vector4R* p)  const 
    {

        hydra::Vector4R p1 = p[0];
        hydra::Vector4R p2 = p[1];
        hydra::Vector4R p3 = p[2];
        hydra::Vector4R p4 = p[3];

        double  theta_pi = ::acos(decay_angle( (p1+p2+p3+p4), (p1+p2), p1 ));
        double  theta_mu = ::acos(decay_angle( (p1+p2+p3+p4), (p3+p4), p3 ));
        double       chi = chi_angle(p3, p4, p1);
        double   Dlambda = fLambda_muplus - fLambda_muminus;
        double lambda_r1 = fJ_r1;

        hydra::complex<double> r(0.0 , 0.0);

        for(unsigned int i=0; i<(2*fJ_r1+1); i++)
        {
           r += hydra::wigner_d_matrix(fJ_r1, lambda_r1, 0.0,     theta_pi) *\
                hydra::wigner_d_matrix(fJ_r2, lambda_r1, Dlambda, theta_mu) *\
                hydra::complex<double>(::cos(-1*lambda_r1*chi) , ::sin(-1*lambda_r1*chi));
           lambda_r1 = lambda_r1-1.0;
        }

        return r;
    }
    
    
    template<typename T>
    __hydra_dual__  inline
    hydra::complex<double> Evaluate(T& p)  const 
    {

        hydra::Vector4R p1 = get<0>(p);
        hydra::Vector4R p2 = get<1>(p);
        hydra::Vector4R p3 = get<2>(p);
        hydra::Vector4R p4 = get<3>(p);
        
        double  theta_pi = ::acos(decay_angle( (p1+p2+p3+p4), (p1+p2), p1 ));
        double  theta_mu = ::acos(decay_angle( (p1+p2+p3+p4), (p3+p4), p3 ));
        double       chi = chi_angle(p3, p4, p1);
        double   Dlambda = fLambda_muplus - fLambda_muminus;
        double lambda_r1 = fJ_r1;
        
        hydra::complex<double> r(0.0, 0.0);
        for(unsigned int i=0; i<(2*fJ_r1+1); i++)
        {
           r += hydra::wigner_d_matrix(fJ_r1, lambda_r1, 0.0,     theta_pi) *\
                hydra::wigner_d_matrix(fJ_r2, lambda_r1, Dlambda, theta_mu) *\
                hydra::complex<double>(::cos(-1*lambda_r1*chi) , ::sin(-1*lambda_r1*chi));
           lambda_r1 = lambda_r1-1.0;
        }

        return r;
        
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
        d1_prime = d1_prime / d1_prime.d3mag();

        double x, y;

        x =  d1_perp.dot(h1_perp);
        y = d1_prime.dot(h1_perp);
      
        return ::atan2(y, x);
    }
    
    double fLambda_muplus;
    double fLambda_muminus;
    double fJ_r1;
    double fJ_r2;



};

} // namespace hydra

#endif /* FOURBODYANGULARDISTRIBUTION_H_ */
