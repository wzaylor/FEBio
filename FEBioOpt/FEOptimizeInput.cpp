#include "stdafx.h"
#include "FEOptimizeData.h"
#include "FELMOptimizeMethod.h"
#include "FEPowellOptimizeMethod.h"
#include "FEScanOptimizeMethod.h"
#include "FEConstrainedLMOptimizeMethod.h"
#include "FEOptimizeInput.h"
#include <FECore/log.h>

//=============================================================================
InvalidVariableName::InvalidVariableName(const char* sz)
{
	strcpy(szname, sz);
}

//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool FEOptimizeInput::ReadParameter(XMLTag& tag, FEParameterList& pl)
{
	// see if we can find this parameter
	FEParam* pp = pl.Find(tag.Name());
	if (pp == 0) return false;

	if (pp->dim() == 1)
	{
		switch (pp->type())
		{
		case FE_PARAM_DOUBLE: tag.value(pp->value<double>()); break;
		case FE_PARAM_INT: tag.value(pp->value<int   >()); break;
		case FE_PARAM_BOOL: tag.value(pp->value<bool  >()); break;
		case FE_PARAM_STRING: tag.value(pp->cvalue()); break;
		default:
			assert(false);
			return false;
		}
	}
	else
	{
		switch (pp->type())
		{
		case FE_PARAM_INT: tag.value(pp->pvalue<int   >(), pp->dim()); break;
		case FE_PARAM_DOUBLE: tag.value(pp->pvalue<double>(), pp->dim()); break;
		default:
			assert(false);
			return false;
		}
	}

	int nattr = tag.m_natt;
	for (int i = 0; i<nattr; ++i)
	{
		const char* szat = tag.m_att[i].m_szatt;
		if (strcmp(szat, "lc") == 0)
		{
			int lc = atoi(tag.m_att[i].m_szatv) - 1;
			if (lc < 0) throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
			switch (pp->type())
			{
			case FE_PARAM_BOOL:
			case FE_PARAM_INT: pp->SetLoadCurve(lc); break;
			case FE_PARAM_DOUBLE: pp->SetLoadCurve(lc, pp->value<double>()); break;
			case FE_PARAM_VEC3D: pp->SetLoadCurve(lc, pp->value<vec3d >()); break;
			default:
				assert(false);
			}
		}
		else
		{
			felog.printf("WARNING: attribute \"%s\" of parameter \"%s\" ignored (line %d)\n", szat, tag.Name(), tag.m_ncurrent_line - 1);
		}
	}

	return true;
}

//=============================================================================
// FEOptimizeInput
//=============================================================================

//-----------------------------------------------------------------------------
//! Read the data from the xml input file
//!
bool FEOptimizeInput::Input(const char* szfile, FEOptimizeData* pOpt)
{
	// try to open the file
	XMLReader xml;
	if (xml.Open(szfile) == false) return false;

	// find the root element
	XMLTag tag;
	if (xml.FindTag("febio_optimize", tag) == false) return false;

	// check for the version attribute
	const char* szversion = tag.AttributeValue("version", true);
	if ((szversion == 0) || (strcmp(szversion, "1.0") != 0))
	{
		fprintf(stderr, "\nFATAL ERROR: Invalid version number for febio_optimize!\n\n");
		return false;
	}

	FEOptimizeData& opt = *pOpt;

	// parse the file
	try
	{
		++tag;
		bool bret = true;
		do
		{
			if      (tag == "Task"       ) bret = ParseTask(tag, opt);
			else if (tag == "Options"    ) bret = ParseOptions(tag, opt);
			else if (tag == "Objective"  ) bret = ParseObjective(tag, opt);
			else if (tag == "Parameters" ) bret = ParseParameters(tag, opt);
			else if (tag == "Constraints") bret = ParseConstraints(tag, opt);
			else if (tag == "LoadData"   ) bret = ParseLoadData(tag, opt);
			else throw XMLReader::InvalidTag(tag);

			if (bret == false) return false;

			// go to the next tag
			++tag;
		} while (!tag.isend());
	}
	catch (InvalidVariableName e)
	{
		fprintf(stderr, "FATAL ERROR: the variable %s is not recognized.\n\n", e.szname);
		return false;
	}
	catch (NothingToOptimize)
	{
		fprintf(stderr, "FATAL ERROR: there is nothing to optimize.\n\n");
		return false;
	}
	catch (XMLReader::Error& e)
	{
		fprintf(stderr, "FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
		return false;
	}
	catch (...)
	{
		fprintf(stderr, "FATAL ERROR: an exception occured in the optimize routine.\n\n");
		return false;
	}
	xml.Close();

	return true;
}

//-----------------------------------------------------------------------------
//! Read the optimizer section of the input file
bool FEOptimizeInput::ParseOptions(XMLTag& tag, FEOptimizeData& opt)
{
	FEOptimizeMethod* popt = 0;
	const char* szt = tag.AttributeValue("type", true);
	if (szt == 0) popt = new FELMOptimizeMethod;
	else
	{
		if (strcmp(szt, "levmar") == 0) popt = new FELMOptimizeMethod;
		else if (strcmp(szt, "powell") == 0) popt = new FEPowellOptimizeMethod;
		else if (strcmp(szt, "scan") == 0) popt = new FEScanOptimizeMethod;
#ifdef HAVE_LEVMAR
		else if (strcmp(szt, "constrained levmar") == 0) popt = new FEConstrainedLMOptimizeMethod;
#endif
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
	}

	// get the parameter list
	FEParameterList& pl = popt->GetParameterList();

	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if (ReadParameter(tag, pl) == false)
			{
				if (tag == "log_level")
				{
					char szval[256];
					tag.value(szval);
					if (strcmp(szval, "LOG_DEFAULT") == 0) {} // don't change the plot level
					else if (strcmp(szval, "LOG_NEVER") == 0) popt->m_loglevel = Logfile::LOG_NEVER;
					else if (strcmp(szval, "LOG_FILE_ONLY") == 0) popt->m_loglevel = Logfile::LOG_FILE;
					else if (strcmp(szval, "LOG_SCREEN_ONLY") == 0) popt->m_loglevel = Logfile::LOG_SCREEN;
					else if (strcmp(szval, "LOG_FILE_AND_SCREEN") == 0) popt->m_loglevel = Logfile::LOG_FILE_AND_SCREEN;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (tag == "print_level")
				{
					char szval[256];
					tag.value(szval);
					if (strcmp(szval, "PRINT_ITERATIONS") == 0) popt->m_print_level = PRINT_ITERATIONS;
					else if (strcmp(szval, "PRINT_VERBOSE") == 0) popt->m_print_level = PRINT_VERBOSE;
				}
				else throw XMLReader::InvalidTag(tag);
			}
			++tag;
		} while (!tag.isend());
	}

	opt.SetSolver(popt);

	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimizeInput::ParseTask(XMLTag& tag, FEOptimizeData& opt)
{
	opt.m_pTask = fecore_new<FECoreTask>(FETASK_ID, tag.szvalue(), &opt.GetFEM());
	if (opt.m_pTask == 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
//! Read the objectives section of the input file
bool FEOptimizeInput::ParseObjective(XMLTag &tag, FEOptimizeData& opt)
{
	FEModel& fem = opt.GetFEM();

	const char* sztype = tag.AttributeValue("type");
	if (sztype == 0) return false;

	if (strcmp(sztype, "data-fit") == 0)
	{
		FEDataFitObjective* obj = new FEDataFitObjective(&fem);
		opt.SetObjective(obj);

		++tag;
		do
		{
			if (tag == "fnc")
			{
				tag.value(obj->m_szname);

				// get the loadcurve for this objective function
				tag.AttributeValue("lc", obj->m_nlc);
				obj->m_nlc--;
			}
			else throw XMLReader::InvalidTag(tag);

			++tag;
		} while (!tag.isend());
	}
	else
	{
		FEMinimizeObjective* obj = new FEMinimizeObjective(&fem);
		opt.SetObjective(obj);

		++tag;
		do
		{
			if (tag == "var")
			{
				const char* szname = tag.AttributeValue("name");
				if (szname == 0) return false;

				double d[2] = {0};
				tag.value(d, 2);

				if (obj->AddFunction(szname, d[0]) == false) throw XMLReader::InvalidAttributeValue(tag, "name", szname);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Read the variables section of the input file
bool FEOptimizeInput::ParseParameters(XMLTag& tag, FEOptimizeData& opt)
{
	FEModel& fem = opt.GetFEM();

	// read the parameters
	++tag;
	do
	{
		if (tag == "param")
		{
			FEModelParameter* var = new FEModelParameter(&fem);

			// get the variable name
			const char* sz = tag.AttributeValue("name");
			if (sz == 0) throw InvalidVariableName("[Unknown]");

			if (var->SetParameter(sz) == 0) throw InvalidVariableName(sz);

			// set initial values and bounds
			double d[4] = { 0, 0, 0, 1 };
			tag.value(d, 4);
			var->SetValue(d[0]);
			var->MinValue() = d[1];
			var->MaxValue() = d[2];
			var->ScaleFactor() = d[3];

			// add the variable
			opt.AddInputParameter(var);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	} while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
//! Parse the Constraints
bool FEOptimizeInput::ParseConstraints(XMLTag &tag, FEOptimizeData &opt)
{
	int NP = opt.InputParameters();
	if ((NP > OPT_MAX_VAR) || (NP < 2)) throw XMLReader::InvalidTag(tag);

	double v[OPT_MAX_VAR + 1];
	++tag;
	do
	{
		if (tag == "constraint")
		{
			int m = tag.value(v, OPT_MAX_VAR + 1);
			if (m != NP + 1) throw XMLReader::InvalidValue(tag);

			OPT_LIN_CONSTRAINT con;
			for (int i = 0; i<NP; ++i) con.a[i] = v[i];
			con.b = v[NP];

			opt.AddLinearConstraint(con);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
//! Read the load data section of the input file
bool FEOptimizeInput::ParseLoadData(XMLTag &tag, FEOptimizeData& opt)
{
	FEObjectiveFunction& obj = opt.GetObjective();

	++tag;
	do
	{
		if (tag == "loadcurve")
		{
			FELoadCurve::INTFUNC ntype = FELoadCurve::LINEAR;
			FELoadCurve::EXTMODE nextm = FELoadCurve::CONSTANT;

			// get the (optional) type
			const char* szt = tag.AttributeValue("type", true);
			if (szt)
			{
				if (strcmp(szt, "step") == 0) ntype = FELoadCurve::STEP;
				else if (strcmp(szt, "linear") == 0) ntype = FELoadCurve::LINEAR;
				else if (strcmp(szt, "smooth") == 0) ntype = FELoadCurve::SMOOTH;
				else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
			}

			// get the optional extend mode
			const char* szm = tag.AttributeValue("extend", true);
			if (szm)
			{
				if (strcmp(szm, "constant") == 0) nextm = FELoadCurve::CONSTANT;
				else if (strcmp(szm, "extrapolate") == 0) nextm = FELoadCurve::EXTRAPOLATE;
				else if (strcmp(szm, "repeat") == 0) nextm = FELoadCurve::REPEAT;
				else if (strcmp(szm, "repeat offset") == 0) nextm = FELoadCurve::REPEAT_OFFSET;
				else throw XMLReader::InvalidAttributeValue(tag, "extend", szt);
			}

			// see if the user wants to read the data from a text file
			const char* szf = tag.AttributeValue("import", true);
			if (szf)
			{
				// make sure this tag is a leaf
				if ((tag.isempty() == false) || (tag.isleaf() == false)) throw XMLReader::InvalidValue(tag);

				// read the data form a text file
				FILE* fp = fopen(szf, "rt");
				if (fp == 0) throw XMLReader::InvalidAttributeValue(tag, "import", szf);
				vector< pair<double, double> > data;
				char szline[256] = { 0 };
				do
				{
					fgets(szline, 255, fp);
					double t, v;
					int n = sscanf(szline, "%lg%lg", &t, &v);
					if (n == 2)
					{
						data.push_back(pair<double, double>(t, v));
					}
					else break;
				} while ((feof(fp) == 0) && (ferror(fp) == 0));

				fclose(fp);

				// create the load curve
				const int nlp = (const int)data.size();
				FELoadCurve* plc = new FELoadCurve;
				plc->Create(nlp);
				plc->SetInterpolation(ntype);
				plc->SetExtendMode(nextm);

				// TODO: This is a hack. Find a better way
				if (dynamic_cast<FEDataFitObjective*>(&obj))
				{
					(dynamic_cast<FEDataFitObjective*>(&obj))->AddLoadCurve(plc);
				}
				

				// set the load points
				for (int i = 0; i<nlp; ++i)
				{
					plc->LoadPoint(i).time = data[i].first;
					plc->LoadPoint(i).value = data[i].second;
				}
			}
			else
			{
				// count how many points we have
				XMLTag t(tag); ++t;
				int nlp = 0;
				while (!t.isend()) { ++nlp; ++t; }

				// create the loadcurve
				FELoadCurve* plc = new FELoadCurve;
				plc->Create(nlp);
				plc->SetInterpolation(ntype);
				plc->SetExtendMode(nextm);

				// TODO: This is a hack. Find a better way
				if (dynamic_cast<FEDataFitObjective*>(&obj))
				{
					(dynamic_cast<FEDataFitObjective*>(&obj))->AddLoadCurve(plc);
				}

				// read the points
				double d[2];
				++tag;
				for (int i = 0; i<nlp; ++i)
				{
					tag.value(d, 2);
					plc->LoadPoint(i).time = d[0];
					plc->LoadPoint(i).value = d[1];

					++tag;
				}
			}
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	} while (!tag.isend());

	return true;
}