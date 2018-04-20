#pragma once
#include "DumpStream.h"

//-----------------------------------------------------------------------------
class FETimeInfo
{
public:
	FETimeInfo();

	FETimeInfo(double time, double tinc);

	void Serialize(DumpStream& ar);

public:
	double	currentTime;		//!< current time value
	double	timeIncrement;		//!< current time step (difference between this time and previous one)
	int		currentIteration;	//!< iteration number

	// HHT time integration parameters
	double	alpha;
	double	beta;
	double	gamma;
	double  alphaf;
	double  alpham;
};