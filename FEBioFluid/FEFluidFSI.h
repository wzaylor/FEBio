//
//  FEFluidFSI.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/13/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEFluidFSI_hpp
#define FEFluidFSI_hpp

#include "FEBioMech/FEElasticMaterial.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! Biphasic material point class.
//
class FEFSIMaterialPoint : public FEMaterialPoint
{
public:
    //! constructor
    FEFSIMaterialPoint(FEMaterialPoint* pt);
    
    //! create a shallow copy
    FEMaterialPoint* Copy();
    
    //! data serialization
    void Serialize(DumpStream& ar);
    
    //! Data initialization
    void Init();
    
public:
    // FSI material data
    vec3d       m_w;      //!< fluid flux relative to solid
    vec3d       m_aw;     //!< material time derivative of m_wt
    double      m_Jdot;   //!< time derivative of solid volume ratio
};

//-----------------------------------------------------------------------------
//! Base class for FluidFSI materials.

class FEFluidFSI : public FEMaterial
{
public:
    FEFluidFSI(FEModel* pfem);
    
    // returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData() override;
    
    // Get the elastic component (overridden from FEMaterial)
    FEElasticMaterial* GetElasticMaterial() override { return m_pSolid->GetElasticMaterial(); }
    
    //! performs initialization
    bool Init() override;
    
public:
    FEFluid* Fluid() { return m_pFluid; }
    FEElasticMaterial* Solid() { return m_pSolid; }
    
private: // material properties
    FEPropertyT<FEElasticMaterial>			m_pSolid;	//!< pointer to elastic solid material
    FEPropertyT<FEFluid>                    m_pFluid;	//!< pointer to fluid material
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidFSI_hpp */