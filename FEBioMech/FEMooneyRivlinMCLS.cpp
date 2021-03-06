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
#include "FEMooneyRivlinMCLS.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEMooneyRivlinMCLS, FEUncoupledMaterial)
	ADD_PARAMETER(m_c1, FE_RANGE_GREATER(0.0), "c1");
	ADD_PARAMETER(m_c2, "c2");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Calculate the deviatoric stress
mat3ds FEMooneyRivlinMCLS::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get material parameters
	double c1 = m_c1(mp);
	double c2 = m_c2(mp);

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

	// Invariants of b_ij are equal to the invariants of C_AB (where C_AB is the right Cauchy-Green deformation tensor)
	double I1 = b_ij.tr();

	// **MCLS** The partial derivative of \tilde{I1} with respect to c_ij (the Finger deformation tensor, c_ij = f_Ai G_AB f_Bj)
	mat3ds partial_I1_cij = j_23*(-b_ij + (1./3.)*I1*identity);

	// **MCLS** The partial derivative of \tilde{I2} with respect to c_ij (the Finger deformation tensor, c_ij = f_Ai G_AB f_Bj)
	mat3ds partial_I2_cij = j_43*b_ij_squared;
	partial_I2_cij += -j_43*(1./3.)*b_ij_squared.tr()*identity;

	partial_I2_cij += -j_43*I1*b_ij;
	partial_I2_cij += j_43*(1./3.)*I1*I1*identity;

	// Calculate deviatoricCauchyStress = -2*j*[(dW/dI1)(dI1/dc_ij) + (dW/dI2)(dI2/dc_ij)]
	// Note that I1 and I2 in the above comment are \tilde{I1} and \tilde{I2}.
	// This definition is similar to Sansour 1993, where W=rho_0*psi, and rho_0 is the density in the reference configuration, and j = rho_t/rho_0
	mat3ds deviatoricCauchyStress = -2.*j*(c1*partial_I1_cij + c2*partial_I2_cij);

	return deviatoricCauchyStress;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric tangent
tens4ds FEMooneyRivlinMCLS::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get material parameters
	double c1 = m_c1(mp);
	double c2 = m_c2(mp);

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

	// Initialize the stiffness matrix (i.e. the tangent matrix)
	tens4ds stiffnessMatrix;

	// The fourth order tensors used to calculate the material stiffness matrix
	tens4ds bXb = dyad1s(b_ij); // The tensor product a^ijkl = b^ij*b^kl
	tens4ds b_ikjl = dyad4s(b_ij); // 0.5*(b^ik*b^jl + b^il*b^jk)
	tens4ds b2b_ijlk = dyad4s(b_ij_squared, b_ij); // 0.5{[b^ma b^ib + b^mb b^ia] b^jn g_mn + b^mi[b^ja b^nb + b^jb b^na]}
	tens4ds b2Xb_ijkl = dyad1s(b_ij_squared, b_ij); // a^ijab = b^ab (b^im g_mn b^nj) + b^ij (b^am g_mn b^nb)
	tens4ds b2Xb2_ijkl = dyad1s(b_ij_squared); // The tensor product a^ijkl = b^im g_mn b^nj b^kr g_rs b^sl
	tens4ds b3Xb_ijkl = dyad1s(b_ij_cubed, b_ij); // The sum of tensor products a^ijkl = (b^3)^ij b^kl + b^ij (b^3)^kl
	tens4ds b3b_ijkl = dyad4s(b_ij_cubed, b_ij); // a^ijab = 0.5{[b^ra b^ib + b^rb b^ia]b^js g_sm b^mn g_rn + [b^ja b^sb + b^jb b^sa]b^ri g_sm b^mn g_rn}
	tens4ds b2_ijkl = dyad4s(b_ij_squared); // a^ijkl = 0.5 b^ri b^js g_sm [b^mk b^nl + b^ml b^nk]

	tens4ds c1Term; // The terms that are multiplied by c1 in the material stiffness matrix
	tens4ds c2Term; // The terms that are multiplied by c1 in the material stiffness matrix

	c1Term = -(1./3.)*b2Xb_ijkl;
	c1Term += (1./9.)*b_ij.tr()*bXb; // Note that bXb already has the 0.5 incorporated.
	c1Term += b2b_ijlk;
	c1Term += -(1./3.)*b_ij.tr()*b_ikjl;

	c2Term = (2./3.)*b3Xb_ijkl;
	c2Term += -(2./9.)*b_ij_squared.tr()*bXb; // Note that bXb already has the 0.5 incorporated.
	c2Term += -(2./3.)*b_ij.tr()*b2Xb_ijkl;
	c2Term += (2./9.)*b_ij.tr()*b_ij.tr()*bXb; // Note that bXb already has the 0.5 incorporated.

	c2Term += -b3b_ijkl;
	c2Term += -b2_ijkl;

	c2Term += (1./3.)*b_ij_squared.tr()*b_ikjl;

	c2Term += b_ij.tr()*b2b_ijlk;
	c2Term += b2Xb2_ijkl; // Note that b2Xb2_ijkl already has the 0.5 incorporated.is multiplied by 0.5 when it is defined above.

	c2Term += -(1./3.)*b_ij.tr()*b_ij.tr()*b_ikjl;

	// Use the terms above to calculate the stiffnessMatrix
	stiffnessMatrix = c1*det_f_Ai_23*c1Term;
	stiffnessMatrix += c2*det_f_Ai_43*c2Term;

	return stiffnessMatrix;
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEMooneyRivlinMCLS::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// get material parameters
	double c1 = m_c1(mp);
	double c2 = m_c2(mp);	
	
	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B.sqr();
    
	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());
    
	//
	// W = C1*(I1 - 3) + C2*(I2 - 3)
	//
    double sed = c1*(I1-3) + c2*(I2-3);
    
    return sed;
}
