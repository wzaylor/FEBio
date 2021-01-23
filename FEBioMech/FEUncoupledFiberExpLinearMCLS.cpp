/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEUncoupledFiberExpLinearMCLS.h"
#include <stdlib.h>
#include <limits>
#ifdef HAVE_GSL
#include "gsl/gsl_sf_expint.h"
#endif

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEUncoupledFiberExpLinearMCLS, FEElasticFiberMaterialUC);
	ADD_PARAMETER(m_c3  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c3");
	ADD_PARAMETER(m_c4  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c4");
	ADD_PARAMETER(m_c5  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c5");
	ADD_PARAMETER(m_lam1, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
	ADD_PARAMETER(m_fiber, "fiber");
	ADD_PARAMETER(m_epsf, "epsilon_scale");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEUncoupledFiberExpLinearMCLS::FEUncoupledFiberExpLinearMCLS(FEModel* pfem) : FEElasticFiberMaterialUC(pfem)
{
	m_c3 = m_c4 = m_c5 = 0;
	m_lam1 = 1;

	m_epsf = 0.0;
}

//-----------------------------------------------------------------------------
//! Fiber material stress
mat3ds FEUncoupledFiberExpLinearMCLS::DevFiberStress(FEMaterialPoint &mp, const vec3d& at)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d f_Ai = pt.m_F; // **MCLS** Note that this is the inverse deformation gradient. 
	double j = pt.m_J; // **MCLS** Note that this is the determinant of the inverse deformation gradient. 
	double J = 1.0 / j;
	double J_13 = pow(J, -1.0 / 3.0);
	// double twoJi = 2.0*Ji;

	// **MCLS** Calculate lambda
	mat3ds b_ij = pt.LeftCauchyGreenMCLS();
	mat3ds c_ij = b_ij.inverse();
	mat3ds aXa = dyad(at); // a_ij = at_i at_j
	double i4 = c_ij.dotdot(aXa); // i4 = (at_i f_Ai G_AB f_Bj at_j), where f_Ai G_AB f_Bj is the Finger deformation tensor
	double lambda = 1./sqrt(i4); // I4 = 1/i4 = lambda^2
	double lambda_tilde = J_13*lambda;

	// invariant I4
	// double I4 = lamd*lamd;

	// strain energy derivative
	double dfdl = 0;
	if (lambda_tilde >= 1)
	{
		double lambda_i = 1.0 / lambda_tilde;
		// double Wl;
		if (lambda < m_lam1)
		{
			dfdl = lambda_i*m_c3*(exp(m_c4*(lambda_tilde - 1)) - 1);
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1)) - 1) - m_c5*m_lam1;
			dfdl = m_c5 + 0.5*c6*lambda_i;
		}
	}
	else
	{
		dfdl = 0;
	}

	// Calculate the Cauchy stress
	mat3ds identity(1,1,1,0,0,0);

	mat3ds cauchyStress = lambda_tilde*lambda_tilde*((1./3.)*identity - aXa); // dI4/dc_ij, note that this is I4_tilde

	cauchyStress *= dfdl; // df/dlambda
	cauchyStress *= 0.5*(1./lambda_tilde); // dlambda/dI4, note that this is I4_tilde
	cauchyStress *= -2.*j; // This is part of the Cauchy stress equation.

	return cauchyStress;
}

//-----------------------------------------------------------------------------
//! Fiber material tangent
tens4ds FEUncoupledFiberExpLinearMCLS::DevFiberTangent(FEMaterialPoint &mp, const vec3d& at)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d f_Ai = pt.m_F; // **MCLS** Note that this is the inverse deformation gradient. 
	mat3d F_iA = f_Ai.inverse(); // **MCLS** The deformation gradient.
	mat3ds V; // **MCLS** The left stretch tensor in the polar decomposition of F_iA = V^ij g_jk R^k_A 
	mat3d R; // **MCLS** The rotation tensor in the polar decomposition of F_iA = V^ij g_jk R^k_A 
	F_iA.left_polar(V, R); // Perform the polar decomposition.
	vec3d Vat = V*at; // The left stretch tensor multiplied by at, v^i = V^ij g_jk at^k
	double j = pt.m_J; // **MCLS** Note that this is the determinant of the inverse deformation gradient. 
	double J = 1.0 / j;
	double J_13 = pow(J, -1.0 / 3.0);
	// double twoJi = 2.0*Ji;

	// **MCLS** Calculate lambda
	mat3ds b_ij = pt.LeftCauchyGreenMCLS();
	mat3ds c_ij = b_ij.inverse();
	mat3ds aXa = dyad(at); // a_ij = at_i at_j
	double i4 = c_ij.dotdot(aXa); // i4 = (at_i f_Ai G_AB f_Bj at_j), where f_Ai G_AB f_Bj is the Finger deformation tensor
	double lambda = 1./sqrt(i4); // I4 = 1/i4 = lambda^2
	double lambda_tilde = J_13*lambda;
	double lambda_tilde43 = pow(lambda_tilde, -4./3.0);
	double I4tilde = lambda_tilde*lambda_tilde;


	// invariant I4
	// double I4 = lamd*lamd;

	// strain energy derivative
	double dfdl = 0; // (df_2/dlambda_tilde) where f_2 is defined in Weiss et al. 1996
	double dfdll = 0; // (df_2/dlambda_tilde) where f_2 is defined in Weiss et al. 1996
	if (lambda_tilde >= 1)
	{
		double lambda_i = 1.0 / lambda_tilde;
		// double Wl;
		if (lambda < m_lam1)
		{
			dfdl = m_c3*lambda_i*(exp(m_c4*(lambda_tilde - 1.) - 1.));
			dfdll = -lambda_i*lambda_i*m_c3*(exp(m_c4*(lambda_tilde - 1.) - 1.));
			dfdll += lambda_i*m_c3*m_c4*exp(m_c4*(lambda_tilde - 1.));

		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1)) - 1) - m_c5*m_lam1;
			dfdl = m_c5 + c6*lambda_i;
			dfdll = -c6*lambda_i*lambda_i;
			
		}
	}
	else
	{
		dfdl = 0;
		dfdll = 0;
	}

	// --- calculate tangent ---

	// The fourth order tensors used to calculate the material stiffness matrix
	tens4ds bXb = dyad1s(b_ij); // The tensor product a^ijkl = b^ij*b^kl
	tens4ds b_ikjl = dyad4s(b_ij); // 0.5*(b^ik*b^jl + b^il*b^jk)
	mat3ds aVVa_ij = dyad(Vat); // The tensor product x^ij = (V^ik g_kl at^l)*(V^jm g_mn at^n)
	mat3ds dI4dcij = I4tilde*((1./3.)*b_ij - aVVa_ij); // dI4/dc_ij = I4*[(1/3)*b_ij - (V^ik g_kl at^l)*(V^jm g_mn at^n)]
	tens4ds dI4dcij_X_dI4dcab = dyad1s(dI4dcij); // fourth order tensor: dI4/dc_ij dI4/dc_ab


	// Initialize the 6x6 matrices.
	tens4ds dI4dcij_term; // The term that is multiplied by the fourth order tensor (dI4/dc_ij)(dI4/dc_ab)
	tens4ds dI4dcijcab_term; // The term that is multiplied by the foruth order tensor (d^2I4/dc_ij dc_ab)
	tens4ds stiffnessMatrix;

	dI4dcij_term = 0.25*(1/(lambda_tilde*lambda_tilde))*dfdll*dI4dcij_X_dI4dcab;
	dI4dcij_term += -0.25*lambda_tilde43*dfdl*dI4dcij_X_dI4dcab;

	tens4ds d2I4_dcijdcab = (1./I4tilde)*dI4dcij_X_dI4dcab - (1./3.)*I4tilde*b_ikjl; // The second partial derivative if I4_tilde w.r.t. c_ij. Note that this is multiplied because dI4dcij_X_dI4dcab has an I4^2 (tilde) term which should just be I4 (tilde)

	dI4dcijcab_term = 0.5*dfdl*(1./lambda_tilde)*d2I4_dcijdcab;

	stiffnessMatrix = dI4dcij_term + dI4dcijcab_term;
	return stiffnessMatrix;
}

//-----------------------------------------------------------------------------
//! Fiber material strain energy density
double FEUncoupledFiberExpLinearMCLS::DevFiberStrainEnergyDensity(FEMaterialPoint &mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// strain energy density
	double sed = 0.0;
#ifdef HAVE_GSL
	if (lamd >= 1)
	{
		if (lamd < m_lam1)
		{
			sed = m_c3*(exp(-m_c4)*
				(gsl_sf_expint_Ei(m_c4*lamd) - gsl_sf_expint_Ei(m_c4))
				- log(lamd));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1 - 1)) - 1) - m_c5*m_lam1;
			sed = m_c5*(lamd - 1) + c6*log(lamd);
		}
	}
#endif
	// --- active contraction contribution to sed is zero ---

	return sed;
}
