#ifndef X2B_H
#define X2B_H

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>

#include "kit.h" 
#include "poly-2b.h"
#include "x2b.h"

////////////////////////////////////////////////////////////////////////////////
namespace x2b {

void cart_to_vars(const double* xyz, double* v, double& s)
{
    const double* Rg  = xyz;

    const double* Cb  = xyz + 3;
    const double* Ob1 = xyz + 6;
    const double* Ob2 = xyz + 9;

    using kit::distance;

    v[0]  = var_intra(d0_intra_CO, m_k_CO_intra, distance(  Cb, Ob1));
    v[1]  = var_intra(d0_intra_CO, m_k_CO_intra, distance(  Cb, Ob2));
    v[2]  = var_intra(d0_intra_OO, m_k_OO_intra, distance( Ob1, Ob2));

    v[3]  = var_inter(d0_inter_RC, m_k_RC_inter, distance( Rg,  Cb));
    v[4]  = var_inter(d0_inter_RO, m_k_RO_inter, distance( Rg, Ob1));
    v[5]  = var_inter(d0_inter_RO, m_k_RO_inter, distance( Rg, Ob2));

    s = 1.0;
} 

double value(const double xyz[18])
{
    double v[num_vars], s;
    cart_to_vars(xyz, v, s);

    double E_poly(0);

    {   
        double mono[num_linear_params+1];
        eval(v, mono);
		
        for (size_t n = 0; n < num_linear_params; ++n) 
            E_poly += poly[n]*mono[n];
    }

    return s*E_poly;
}

}
#endif // X2B_H

////////////////////////////////////////////////////////////////////////////////
