/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FENeoHookeanMCLS.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FENeoHookeanMCLS, FEElasticMaterial)
	ADD_PARAMETER2(m_E, FE_PARAM_DOUBLE, FE_RANGE_GREATER   (      0.0), "E");
	ADD_PARAMETER2(m_v, FE_PARAM_DOUBLE, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
mat3ds FENeoHookeanMCLS::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double j = pt.m_J; // **MCLS** This is the determinant of the inverse deformation gradient (i.e. j = 1/J)
	double J = 1.0/j; // **MCLS** The determinant of the deformation gradient (J = 1/j).
	double lnJ = log(J);

	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreenMCLS();

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// Identity
	mat3dd I(1);

	// calculate stress
	mat3ds s = (b - I)*(mu*j) + I*(lam*lnJ*j);

	return s;
}

//-----------------------------------------------------------------------------
// NOTE:: This differs from the original code because the input 'D' is an input that is defined by reference.
// This function does not return a tens4ds
void FENeoHookeanMCLS::Tangent(FEMaterialPoint& mp, double D[6][6])
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	double detF = pt.m_J;

	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreenMCLS();
	double identity[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// double lam1 = lam / detF;
	// double mu1  = (mu - lam*log(detF)) / detF;
	
	int indices[6][2] = {{0,0}, {1,1}, {2,2}, {0,1}, {1,2}, {0,2}}; // The indices we need to setup the Voigt form of the Material stiffness matrix.
	int i, j, k, l; // The indices that the values in 'indices' will be assigned to.

	for (int ij_index = 0; ij_index<6; ++ij_index)
    {
		i = indices[ij_index][0];
		j = indices[ij_index][1];
		for (int kl_index = 0; kl_index<6; ++kl_index)
		{
			k = indices[kl_index][0];
			l = indices[kl_index][1];
			D[ij_index][kl_index] = -mu*0.5*(b(i,k)*b(l,j) + b(i,l)*b(k,j)) - 0.5*lam*b(k,l)*identity[i][j];
		}
	}
//            materialStiffnessMatrix[ij_index, kl_index] = -muLame*0.5*(c_ij_inv[i,k]*c_ij_inv[l,j] + c_ij_inv[i,l]*c_ij_inv[k,j]) - 0.5*lambdaLame*c_ij_inv[k,l]*identity[i,j]

	return;
}

//-----------------------------------------------------------------------------
// This is here to the project will build appropriately
tens4ds FENeoHookeanMCLS::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	double detF = pt.m_J;

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	double lam1 = lam / detF;
	double mu1  = (mu - lam*log(detF)) / detF;
	
	double D[6][6] = {0};
	D[0][0] = lam1+2.*mu1; D[0][1] = lam1       ; D[0][2] = lam1       ;
	D[1][0] = lam1       ; D[1][1] = lam1+2.*mu1; D[1][2] = lam1       ;
	D[2][0] = lam1       ; D[2][1] = lam1       ; D[2][2] = lam1+2.*mu1;
	D[3][3] = mu1;
	D[4][4] = mu1;
	D[5][5] = mu1;

	return tens4ds(D);
}

//-----------------------------------------------------------------------------
double FENeoHookeanMCLS::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	double lnJ = log(J);
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = b.tr();
	
	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);
	
	double sed = mu*((I1-3)/2.0 - lnJ)+lam*lnJ*lnJ/2.0;
	
	return sed;
}
