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
#include "FETransIsoVerondaWestmannMCLS.h"

// define the material parameters
BEGIN_FECORE_CLASS(FETransIsoVerondaWestmannMCLS, FEUncoupledMaterialMCLS)
	ADD_PARAMETER(      m_c1  , "c1");
	ADD_PARAMETER(      m_c2  , "c2");
	ADD_PARAMETER(m_fib.m_c3  , "c3");
	ADD_PARAMETER(m_fib.m_c4  , "c4");
	ADD_PARAMETER(m_fib.m_c5  , "c5");
	ADD_PARAMETER(m_fib.m_lam1, "lam_max");
	ADD_PARAMETER(m_fib.m_fiber, "fiber");

	ADD_PROPERTY(m_ac, "active_contraction", FEProperty::Optional);
END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// FETransIsoVerondaWestmann
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
FETransIsoVerondaWestmannMCLS::FETransIsoVerondaWestmannMCLS(FEModel* pfem) : FEUncoupledMaterialMCLS(pfem), m_fib(pfem)
{
	m_ac = 0;
	m_fib.SetParent(this);
}

//-----------------------------------------------------------------------------
mat3ds FETransIsoVerondaWestmannMCLS::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double j = pt.m_J; // **MCLS** This is the determinant of the inverse deformation gradient (i.e. j = 1/J)
	double J = 1.0/j; // **MCLS** The determinant of the deformation gradient (J = 1/j).
	double j_23 = pow(j, 2.0/3.0); // **MCLS** The inverse jacobian raised to the 2/3 power
	double j_43 = pow(j, 4.0/3.0); // **MCLS** The inverse jacobian raised to the 4/3 power

	// calculate left Cauchy-Green tensor 
	// **MCLS** Calculate the left Cauchy-Green tensor (originally this calculated the deviatoric left Cauchy-Green tensor)
	mat3ds b_ij = pt.LeftCauchyGreenMCLS();
	mat3ds identity(1,1,1,0,0,0);

	// calculate square of b_ij
	mat3ds b_ij_squared = b_ij.sqr();

	// The first invariant of b_ij_tilde
	double I1_tilde = j_23*b_ij.tr();

	// Invariants of b_ij are equal to the invariants of C_AB (where C_AB is the right Cauchy-Green deformation tensor)
	double I1 = b_ij.tr();
	double I2 = 0.5*(I1*I1 - b_ij_squared.tr());

	// **MCLS** The partial derivative of \tilde{I1} with respect to g_ij (the metric on the deformed configuration)
	mat3ds partial_I1_cij = j_23*(b_ij - (1./3.)*I1*identity);

	// **MCLS** The partial derivative of \tilde{I2} with respect to g_ij (the metric on the deformed configuration)
	mat3ds partial_I2_cij = j_43*(I1*b_ij - b_ij_squared - (2./3.)*I2*identity);

	// **MCLS** The partial derivative of F with respect to I1_tilde and I2_tilde
	double dfdI1 = m_c1*m_c2*exp(m_c2*(I1_tilde - 3));
	double dfdI2 = -0.5*m_c1*m_c2;

	// Calculate deviatoricCauchyStress = -2*j*[(dW/dI1)(dI1/dg_ij) + (dW/dI2)(dI2/dg_ij)]
	// Note that I1 and I2 in the above comment are \tilde{I1} and \tilde{I2}.
	mat3ds deviatoricCauchyStress = 2.*j*(dfdI1*partial_I1_cij + dfdI2*partial_I2_cij);

	// calculate the passive fiber stress
	mat3ds fiberDeviatoricCauchyStress = m_fib.DevStress(mp);

	// **MCLS** This material doesn't account for active fiber stress
	// calculate the active fiber stress (if provided)
	// if (m_ac) fs += m_ac->FiberStress(m_fib.FiberVector(mp), pt);

	return deviatoricCauchyStress + fiberDeviatoricCauchyStress;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4dmm FETransIsoVerondaWestmannMCLS::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// **MCLS** Define the determinant of the inverse deformation gradient (f_Ai).
	// pt.m_J is calculated using pt.m_F, where pt.m_F is the inverse deformation gradient (f_Ai)
	// pt.m_F is the inverse deformation gradient because because we are taking the original known node positions to be in the deformed configuration, and taking the displacements to be defining the node positions in the reference configuration.
	// Therefore, pt.m_F would more accurately be called pt.m_f, however this would require updating much more of the code.
	double det_f_Ai_23 = pow(pt.m_J, 2./3.);
	double det_f_Ai_43 = pow(pt.m_J, 4./3.);

	// calculate left Cauchy-Green tensor: B = F*Ft
	mat3ds b_ij = pt.LeftCauchyGreenMCLS();
	mat3ds b_ij_squared = b_ij.sqr(); // b^ik g_kl b^lj
	mat3ds b_ij_cubed; // (b^3)^ij = b^ik g_kl b^lm g_mn b^nj
	// Define b_ij_cubed
	b_ij_cubed.xx() = b_ij.xx()*b_ij_squared.xx() + b_ij.xy()*b_ij_squared.xy() + b_ij.xz()*b_ij_squared.xz();
	b_ij_cubed.yy() = b_ij.xy()*b_ij_squared.xy() + b_ij.yy()*b_ij_squared.yy() + b_ij.yz()*b_ij_squared.yz();
	b_ij_cubed.zz() = b_ij.xz()*b_ij_squared.xz() + b_ij.yz()*b_ij_squared.yz() + b_ij.zz()*b_ij_squared.zz();
	b_ij_cubed.xy() = b_ij.xx()*b_ij_squared.xy() + b_ij.xy()*b_ij_squared.yy() + b_ij.xz()*b_ij_squared.yz();
	b_ij_cubed.yz() = b_ij.xy()*b_ij_squared.xz() + b_ij.yy()*b_ij_squared.yz() + b_ij.yz()*b_ij_squared.zz();
	b_ij_cubed.xz() = b_ij.xx()*b_ij_squared.xz() + b_ij.xy()*b_ij_squared.yz() + b_ij.xz()*b_ij_squared.zz();

	// Tensor invariants.
	double I1 = b_ij.tr();
	double I2 = 0.5*(I1*I1 - b_ij_squared.tr());
	// The first invariant of b_ij_tilde
	double I1_tilde = det_f_Ai_23*b_ij.tr();

	// Other terms
	mat3ds identity(1,1,1,0,0,0);

	// Initialize the stiffness matrix (i.e. the tangent matrix)
	tens4dmm stiffnessMatrix;

	// Second order tensors that are used to make the fourth order tensors
	mat3ds dI1dcij = (1./3.)*b_ij*I1 - b_ij_squared; // (1/3)b^ij*I1 - b^ik g_kl b^lj NOTE: The J^(-1/3) is not included here.
	mat3ds dI1dgij = b_ij*I1 - (1./3.)*I1*identity; // b^ij -(1/3)I1 g^ij NOTE: The J^(-1/3) is not included here.

	// The fourth order tensors used to calculate the material stiffness matrix
	tens4dmm dI1dgdc = dyad1mm(dI1dgij, dI1dcij); // The tensor product a^ijkl = (b^ij -(1/3)I1 g^ij)((1/3)b^kl*I1 - b^km g_mn b^nl) NOTE: The J^(-2/3) is not included here.
	tens4dmm b_ijkl = dyad1mm(b_ij, b_ij); // The tensor product a^ijkl = b^ij*b^kl
	tens4dmm gb_ijkl = dyad1mm(identity, b_ij); // The tensor product a^ijkl = g^ij*b^kl
	tens4dmm gbb_ijkl = dyad1mm(identity, b_ij_squared); // The tensor product a^ijkl = g^ij b^ka g_ab b^bl
	tens4ds b_ijkl_sym = dyad4s(b_ij); // 0.5*(b^ik b^jl + b^jk b^il)
	// ---
	tens4dmm bbb_ijkl = dyad1mm(b_ij_squared, b_ij); // a^ijkl = b^kl(b^ia g_ab b^bj)
	tens4ds bb_ikjl = dyad4s(b_ij_squared, b_ij); // a^ijkl = 0.5(b^ik b^ml + b^mk b^il)g_mn b^nj + 0.5(b^jk b^ml + b^mk b^jl)g_mn b^ni
	tens4dmm gbbb_ijkl = dyad1mm(identity, b_ij_cubed); // The tensor product a^ijkl = g^ij b^ka g_ab b^bm g_mn b^ml

	// Initialize the terms that are multiplied by c1 and c2
	tens4dmm dI1dgdc_term; // a_ijkl = d^2(I1_tilde)/dg_ij dc_kl
	tens4dmm dI2dgdc_term;
	dI1dgdc_term.zero();
	dI2dgdc_term.zero();

	// Populate the c1 term
	dI1dgdc_term += (1./3.)*b_ijkl - (1./9.)*I1*gb_ijkl;
	dI1dgdc_term += (1./3.)*gbb_ijkl - b_ijkl_sym;

	// Populate the c2 term
	dI2dgdc_term += (2./3.)*I1*b_ijkl;
	dI2dgdc_term += -(2./3.)*bbb_ijkl;
	dI2dgdc_term += -(4./9.)*I2*gb_ijkl;

	dI2dgdc_term += -bbb_ijkl.transpose(); // The transpose turns bb_ijkl into a^ijkl = b^ij(b^ja g_ab b^bl)
	dI2dgdc_term += -I1*b_ijkl_sym;
	dI2dgdc_term += bb_ikjl;
	dI2dgdc_term += (2./3.)*I1*gbb_ijkl;
	dI2dgdc_term += -(2./3.)*gbbb_ijkl;

	stiffnessMatrix = m_c1*m_c2*m_c2*exp(m_c2*(I1_tilde - 3))*det_f_Ai_43*dI1dgdc;
	stiffnessMatrix += m_c1*m_c2*exp(m_c2*(I1_tilde - 3))*det_f_Ai_23*dI1dgdc_term;
	stiffnessMatrix += -0.5*m_c1*m_c2*det_f_Ai_43*dI2dgdc_term;

	// add the passive fiber stiffness
	stiffnessMatrix += m_fib.DevTangent(mp);

	// **MCLS** This material doesn't account for active fiber stress
	// add the active fiber stiffness
	// if (m_ac) c += m_ac->FiberStiffness(m_fib.FiberVector(mp), mp);

	return stiffnessMatrix;
}

//-----------------------------------------------------------------------------
double FETransIsoVerondaWestmannMCLS::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B.sqr();
    
	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());
    
	// calculate sed
    double sed = m_c1*(exp(m_c2*(I1-3))-1) - m_c1*m_c2*(I2-3)/2;
    
	// add the fiber strain energy density
	sed += m_fib.DevStrainEnergyDensity(mp);
    
	return sed;
}
