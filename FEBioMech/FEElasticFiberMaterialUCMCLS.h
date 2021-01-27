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



#pragma once
#include "FEUncoupledMaterialMCLS.h"

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterialUCMCLS : public FEUncoupledMaterialMCLS
{
public:
    FEElasticFiberMaterialUCMCLS(FEModel* pfem);

	// Get the fiber direction (in global coordinates) at a material point
	vec3d FiberVector(FEMaterialPoint& mp);

	// calculate stress in fiber direction a0
	virtual mat3ds DevFiberStress(FEMaterialPoint& mp, const vec3d& a0) = 0;

	// Spatial tangent
	virtual tens4dmm DevFiberTangent(FEMaterialPoint& mp, const vec3d& a0) = 0;

	//! Strain energy density
	virtual double DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0) = 0;

public:
	// These are made private since fiber materials should implement the functions above instead. 
	// The functions can still be reached when a fiber material is used in an elastic mixture. 
	// In those cases the fiber vector is taken from the first column of Q. 
	mat3ds DevStress(FEMaterialPoint& mp) final { return DevFiberStress(mp, FiberVector(mp)); }
	tens4dmm DevTangent(FEMaterialPoint& mp) final { return DevFiberTangent(mp, FiberVector(mp)); }
	double DevStrainEnergyDensity(FEMaterialPoint& mp) final { return DevFiberStrainEnergyDensity(mp, FiberVector(mp)); }

public:
	FEParamVec3		m_fiber;	//!< fiber orientation

	double	m_epsf;

	DECLARE_FECORE_CLASS();
};
