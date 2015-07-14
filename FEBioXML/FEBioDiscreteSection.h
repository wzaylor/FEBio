#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
class FEBioDiscreteSection : public FEBioFileSection
{
public:
	FEBioDiscreteSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseSpringSection(XMLTag& tag);
};
