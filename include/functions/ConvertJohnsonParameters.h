/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/

#ifndef CONVERTJOHNSONPARAMETERS_H_
#define CONVERTJOHNSONPARAMETERS_H_

#include <iostream>


class ConvertJohnsonParameters 

{

public:

    ConvertJohnsonParameters() = delete;

    ConvertJohnsonParameters(double mean, double width, double nu, double tau):
              _mean(mean),
              _width(width),
              _nu(nu),
              _tau(tau)
    {
      _UpdateNewPar();
    }

    ConvertJohnsonParameters(std::array<double, 4>&& old_par):
              _mean(old_par[0]),
              _width(old_par[1]),
              _nu(old_par[2]),
              _tau(old_par[3])
    {
      _UpdateNewPar();
    }

    ConvertJohnsonParameters& operator=( const ConvertJohnsonParameters& other)
    {
        if(this==&other) return *this;
        _mean  = other.GetMean();
        _width = other.GetWidth();
        _nu    = other.GetNu();
        _tau   = other.GetTau();
        _UpdateNewPar();
        return *this;
    }

    inline double GetMean()  const { return _mean; }
    
    inline double GetWidth() const { return _width; }
    
    inline double GetNu()    const { return _nu; }
    
    inline double GetTau()   const { return _tau; }
    
    inline double GetGamma() const { return _gamma; }
    
    inline double GetDelta() const { return _delta; }
    
    inline double GetLamda() const { return _lambda; }
    
    inline double GetXi()   const { return _xi; }
    
    
    inline double SetMean(double& mean)  const { 
      _mean = mean;
      _UpdateNewPar(); 
    }
    
    inline double SetWidth(const double& width) const { 
       _width = width; 
       _UpdateNewPar(); 
    }
    
    inline double SetNu(const double& nu) const { 
       _nu = nu; 
       _UpdateNewPar(); 
    }
    
    inline double SetTau(const double& tau) const { 
       _tau = tau; 
       _UpdateNewPar(); 
    }
    


private: 
    
    void inline _UpdateNewPar(){
        _A     = 0.5 *(exp( _tau * _tau)-1);
        _B     = (exp( _tau * _tau) * cosh(2 * (-1*_nu *_tau)) + 1);
        _gamma  = -1*_nu;
        _delta  = 1./_tau;
        _lambda = _width/sqrt(_A*_B); 
        _xi     = _mean + _lambda*exp(0.5*pow(_delta,-2))*sinh(_gamma/_delta);
    }
    
    double _mean;
    double _width;
    double _nu;
    double _tau;
    double _cc;
    double _A;
    double _B;
    double _gamma;
    double _delta;
    double _lambda;
    double _xi;



};

#endif //CONVERTJOHNSONPARAMETERS_H_
