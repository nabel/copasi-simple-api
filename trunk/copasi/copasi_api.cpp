//std
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <limits> //get max and min for double

//copasi
#define COPASI_MAIN 1
#include "copasi_api.h"
#include "copasi/copasi.h"
#include "copasi/report/CCopasiRootContainer.h"
#include "copasi/CopasiDataModel/CCopasiDataModel.h"
#include "copasi/model/CModel.h"
#include "copasi/model/CCompartment.h"
#include "copasi/model/CMetab.h"
#include "copasi/model/CReaction.h"
#include "copasi/model/CChemEq.h"
#include "copasi/model/CModelValue.h"
#include "copasi/function/CFunctionDB.h"
#include "copasi/function/CFunction.h"
#include "copasi/function/CEvaluationTree.h"
#include "copasi/report/CReport.h"
#include "copasi/report/CReportDefinition.h"
#include "copasi/report/CReportDefinitionVector.h"
#include "copasi/trajectory/CTrajectoryTask.h"
#include "copasi/trajectory/CTrajectoryMethod.h"
#include "copasi/trajectory/CTrajectoryProblem.h"
#include "copasi/scan/CScanTask.h"
#include "copasi/scan/CScanMethod.h"
#include "copasi/scan/CScanProblem.h"
#include "copasi/trajectory/CTimeSeries.h"
#include "copasi/steadystate/CSteadyStateTask.h"
#include "copasi/steadystate/CSteadyStateProblem.h"
#include "copasi/steadystate/CMCATask.h"
#include "copasi/steadystate/CMCAMethod.h"
#include "copasi/elementaryFluxModes/CFluxMode.h"
#include "copasi/elementaryFluxModes/CEFMTask.h"
#include "copasi/elementaryFluxModes/CEFMProblem.h"
#include "copasi/commandline/COptions.h"
#include "copasi/report/CCopasiContainer.h"
#include "copasi/parameterFitting/CFitTask.h"
#include "copasi/parameterFitting/CFitMethod.h"
#include "copasi/parameterFitting/CFitProblem.h"
#include "copasi/parameterFitting/CFitItem.h"
#include "copasi/parameterFitting/CExperimentSet.h"
#include "copasi/parameterFitting/CExperiment.h"
#include "copasi/parameterFitting/CExperimentObjectMap.h"
#include "copasi/report/CKeyFactory.h"

//genetic algorithm (used for optimization)
#include "GASStateGA.h"
#include "GA1DArrayGenome.h"

//libstruct
#include "libstructural.h"
#include "matrix.h"

//regex
#include "boost/regex.hpp"

//parse math (used for optimization)
#include "muParserDef.h"
#include "muParser.h"
#include "muParserInt.h"
extern "C"
{
	#include "mtrand.h"
}

//Antimony lib
extern "C"
{
	#define LIB_EXPORTS 1
	#include "src/antimony_api.h"
}

using namespace LIB_STRUCTURAL;
using namespace LIB_LA;
using namespace std;

//using macro instead of this variable (as in original copasi code) causes some issues in visual studio
unsigned C_INT32 C_INVALID_INDEX = std::numeric_limits< unsigned C_INT32 >::max();
//#ifdef _WIN32
	static double NaN = std::numeric_limits<double>::quiet_NaN();
//#else
	//static double NaN = 0.0/0.0;
//#endif

static double _EPSILON = 1E-3;

double cSetEpsilon(double eps)
{
	if (eps > 0)
		_EPSILON = eps;
	return _EPSILON;
}

//this "wrapper" struct is used to store pointer to
//either a compartment, species, reaction, or parameter
//it also stores the copasi's unique name and the normal human readable name
struct CopasiPtr
{
	string name;
	string key;
	CMetab * species;
	CCompartment * compartment;
	CReaction * reaction;
	CModelValue * param;
	string assignmentRule;
	bool unused;
};


/***************************
  All of these global variables
  are related hashing names
  and return types
***************************/
typedef map< string, CopasiPtr > CCMap;
static list< CCMap* > hashTablesToCleanup;
static list< copasi_model > copasi_modelsToCleanup;
static map< void*,  map<string,bool> > returnTypeFilters;

/*****************************
  These functions are used for
   string replacements
******************************/
static int rename(string& target, const string& oldname,const string& newname0);
static int replaceSubstring(string& target, const string& oldname,const string& newname0);
static boost::regex sbmlPowFunction("pow\\s*\\(\\s*([^,]+)\\s*,\\s*([^,]+)\\s*\\)", boost::regex::perl);

/*****************************
  These functions are used for
   hashing
******************************/
template <typename T1, typename T2>
bool contains(map<T1,T2> * hash, const T1 & s)
{
	return hash && (hash->find(s) != hash->end());
}

bool contains(const string& str, const string & s)
{
	return str.find(s) != string::npos;
}

template <typename T1, typename T2>
T2 & getHashValue(map<T1,T2> * hash, const T1 & s)
{
	return (*hash)[s];
}

template <typename T1, typename T2>
void hashInsert(map<T1,T2> * hash, const T1 & s, T2 v)
{
	(*hash)[s] = v;
}

int indexOf( list<string>& lst, const string & s)
{
	int k = 0;
	for (list<string>::iterator i=lst.begin(); i != lst.end(); i++, ++k)
		if ((*i) == s)
			return k;
	return -1;
}

int indexOf( vector<string>& lst, const string & s)
{
	int k = 0;
	for (vector<string>::iterator i=lst.begin(); i != lst.end(); i++, ++k)
		if ((*i) == s)
			return k;
	return -1;
}

double string_to_double( const std::string& s )
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

list<string> splitString(const string& seq, const string& delim);

/*****************************
       Copasi helper function
******************************/

void copasi_init()
{
	CCopasiRootContainer::init(0, NULL);
}

void copasi_end()
{
	for (list<CCMap*>::iterator i = hashTablesToCleanup.begin(); i != hashTablesToCleanup.end(); i++)
		delete (*i);

	list< copasi_model > models = copasi_modelsToCleanup;
	copasi_modelsToCleanup.clear();

	for (list<copasi_model>::iterator i = models.begin(); i != models.end(); i++)
		cRemoveModel(*i);

	CCopasiRootContainer::destroy();
}

//enable assignment rule substitutions, e.g. if A = sin(time), then B = 2*A  becomes B = 2*sin(time)
std::string cSetAssignmentRuleHelper(copasi_model , const char * );

static int DO_CLEANUP_ASSIGNMENT_RULES = 1;

void cEnableAssignmentRuleReordering()
{
	DO_CLEANUP_ASSIGNMENT_RULES = 1;
}

void cDisableAssignmentRuleReordering()
{
	DO_CLEANUP_ASSIGNMENT_RULES = 0;
}

//populate the hash table
void populate_hash(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel || !hash) return;

	CCopasiVectorNS < CCompartment > & compartments = pModel->getCompartments();
	CCopasiVector< CMetab > & species = pModel->getMetabolites();
	CCopasiVectorNS < CReaction > & reacs = pModel->getReactions();
	CCopasiVectorN< CModelValue > & params = pModel->getModelValues();

	for (int i=0; i < compartments.size(); ++i)
		if (compartments[i])
		{
			CopasiPtr copasiPtr = {
				(std::string)compartments[i]->getCN(),
				compartments[i]->getKey(),
				0,
				compartments[i],
				0,
				0,
				"",
				false};

			hashInsert(hash,   compartments[i]->getObjectName(),		copasiPtr );
			hashInsert(hash,   compartments[i]->getSBMLId(),		copasiPtr );
		}

	for (int i=0; i < species.size(); ++i)
		if (species[i])
		{
			CopasiPtr copasiPtr = {
				(std::string)species[i]->getCN(),
				species[i]->getKey(),
				species[i],
				0,
				0,
				0,
				"",
				false};

			hashInsert(hash,  species[i]->getObjectName(),	  copasiPtr );
			hashInsert(hash,  species[i]->getSBMLId(),	  copasiPtr );

			if (species[i]->getCompartment())
				hashInsert(hash,  species[i]->getCompartment()->getObjectName() + string("_") + species[i]->getObjectName(), copasiPtr);
				hashInsert(hash,  species[i]->getCompartment()->getCN() + string("_") + species[i]->getSBMLId(), copasiPtr);

		}

	for (int i=0; i < reacs.size(); ++i)
		if (reacs[i])
		{
			CopasiPtr copasiPtr = {
				(std::string)reacs[i]->getCN(),
				reacs[i]->getKey(),
				0,
				0,
				reacs[i],
				0,
				"",
				false};

			hashInsert(hash,   reacs[i]->getObjectName(),		copasiPtr );
			hashInsert(hash,   reacs[i]->getSBMLId(),		copasiPtr );
		}

	for (int i=0; i < params.size(); ++i)
		if (params[i])
		{
			CopasiPtr copasiPtr = {
				(std::string)params[i]->getCN(),
				params[i]->getKey(),
				0,
				0,
				0,
				params[i],
				"",
				false};

			hashInsert(hash,   params[i]->getObjectName(),	 copasiPtr );
			hashInsert(hash,   params[i]->getSBMLId(),	 copasiPtr );
		}
}

int copasi_cleanup_assignments(copasi_model model)
{
	CCMap * hash = (CCMap*)(model.qHash);
	if (!hash) return 0;

	CMetab* pSpecies = 0;
	vector<string> names, assignments;
	names.reserve(hash->size());
	assignments.reserve(hash->size());

	for (CCMap::iterator i = hash->begin(); i != hash->end(); i++)
		if ( (*i).second.species && !(*i).second.assignmentRule.empty())
		{
			names.push_back( (*i).first );
			assignments.push_back( (*i).second.assignmentRule );
		}
		else
		{
			names.push_back( string() );
			assignments.push_back( string() );
		}

	int retval = 1;
	bool replace_needed = (bool)(DO_CLEANUP_ASSIGNMENT_RULES);

	while (replace_needed)
	{
		replace_needed = false;
		for (CCMap::iterator i = hash->begin(); i != hash->end(); i++)
			if ((*i).second.species && !(*i).second.assignmentRule.empty())
			{
				for (int j=0; j < names.size(); ++j)
					if (!names[j].empty() &&
						names[j] != (*i).first &&
						contains((*i).second.assignmentRule, names[j]))
					{
						if (rename((*i).second.assignmentRule, names[j], assignments[j]))
							replace_needed = true;
					}
			}
	}

	for (CCMap::iterator i = hash->begin(); i != hash->end(); i++)
		if ((*i).second.species)
			if ((*i).second.assignmentRule.empty())
			{
				if ((*i).second.unused)
				{
					(*i).second.species->setStatus(CModelEntity::FIXED); //unused species
				}
			}
			else
			{
                std::string s = cSetAssignmentRuleHelper(model, (*i).second.assignmentRule.c_str());
                if (s.length() < 1)
                {
                    (*i).second.species->setStatus(CModelEntity::REACTIONS);
                }
                else
                {
                    (*i).second.species->setStatus(CModelEntity::ASSIGNMENT);
                    retval = retval & (*i).second.species->setExpression(s);
                }
			}
	return retval;
}

/**********************************************************************************
 The following function is a "hack" for avoiding repeated memory allocations when calling
  c_createMatrix. For example, suppose a function is calling getSpeciesConcentrations
  inside a loop, then it can use setStorageMatrix to prevent re-allocation of new memory
  during every iteration (for speedup). The calling function must be extremely careful.
  The unsetStorageMatrix must be called at the end of the loop.
***********************************************************************************/
static c_matrix * _StorageMatrix = NULL;
void setStorageMatrix(c_matrix * m)
{
	_StorageMatrix = m;
}

void unsetStorageMatrix()
{
	_StorageMatrix = NULL;
}

c_matrix efficiently_createMatrix(int r, int c)
{
	c_matrix m;
	if (_StorageMatrix && _StorageMatrix->rows == r && _StorageMatrix->cols == c)
		return (*_StorageMatrix);

	m = c_createMatrix(r,c);

	//if (!_StorageMatrix)
		//_StorageMatrix = &m;

	return m;
}

/* The Main API functions */

void cRemoveModel(copasi_model model)
{
	//remove from list
	for (list<copasi_model>::iterator i=copasi_modelsToCleanup.begin(); i != copasi_modelsToCleanup.end(); i++)
		if ((*i).CopasiDataModelPtr == model.CopasiDataModelPtr)
		{
			copasi_model m = { (void*)NULL, (void*)NULL, (void*)NULL, (char*)NULL, (char*)NULL};
			(*i) = m;
		}

	//delete model
	if (model.errorMessage)
		free(model.errorMessage);
	if (model.warningMessage)
		free(model.warningMessage);
	if (model.CopasiDataModelPtr)
		CCopasiRootContainer::removeDatamodel((CCopasiDataModel*)model.CopasiDataModelPtr);
}

void clearcopasi_model(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel || !pDataModel || !hash) return;

	CopasiPtr p;

	for (CCMap::iterator i = hash->begin(); i != hash->end(); i++)
	{
		p = (*i).second;
		if (p.species)
			pModel->remove(p.species);
		else
		if (p.param)
			pModel->remove(p.param);
		else
		if (p.compartment)
			pModel->remove(p.compartment);
		else
		if (p.reaction)
			pModel->remove(p.reaction);
	}

	hash->clear();
}

// ------------------------------------------------------------------
// Create model group
// ------------------------------------------------------------------


copasi_model cCreateModel(const char * name)
{
	copasi_init();

	CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();
	CModel* pModel = pDataModel->getModel();
	CCMap * qHash = new CCMap();
	copasi_model m = { (void*)(pModel) , (void*)(pDataModel), (void*)(qHash), (char*)(NULL), (char*)(NULL)};

	hashTablesToCleanup.push_back(qHash);
	copasi_modelsToCleanup.push_back(m);

	pModel->setSBMLId( string(name) );
	pModel->setObjectName( string(name) );
	//pModel->setTimeUnit(CModel::s);
	//pModel->setVolumeUnit(CModel::microl);
	//pModel->setQuantityUnit(CModel::nMol);
	pModel->setTimeUnit(CModel::dimensionlessTime);
	pModel->setVolumeUnit(CModel::dimensionlessVolume);
	pModel->setQuantityUnit(CModel::dimensionlessQuantity);

	cCreateVariable(m, "time", "time");

	return m;
}

void cCreateSpecies(copasi_compartment compartment, const char* name, double iv)
{
	CModel* pModel = (CModel*)(compartment.CopasiModelPtr);
	CCompartment* pCompartment = (CCompartment*)(compartment.CopasiCompartmentPtr);
	CCMap * hash = (CCMap*)(compartment.qHash);
	CMetab* pSpecies;

	if (!pModel || !hash || !pCompartment) return;
	if (contains(hash, string(name)))
	{
		pSpecies = getHashValue(hash,string(name)).species;
		if (pSpecies)
		{
			pSpecies->setConcentration(iv);
			//pSpecies->setValue(iv);
			//pSpecies->setInitialValue(iv);
			pSpecies->setInitialConcentration(iv);
		}
		return;
	}

	pSpecies = pModel->createMetabolite(name, pCompartment->getObjectName(), iv, CMetab::REACTIONS);
	pSpecies->setConcentration(iv);
	//pSpecies->setValue(iv);
	//pSpecies->setInitialValue(iv);
	pSpecies->setInitialConcentration(iv);

	CopasiPtr copasiPtr = {
			pSpecies->getCN(),
			pSpecies->getKey(),
			pSpecies,
			0,
			0,
			0,
			"",
			true};

	hashInsert(hash,
				(pCompartment->getObjectName() + string("_") + string(name)),
				copasiPtr
				);

	hashInsert(hash, string(name), copasiPtr); //for speedy lookup
}

copasi_compartment cCreateCompartment(copasi_model model, const char* name, double volume)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCompartment* pCompartment = pModel->createCompartment(name, volume);
	CCMap * hash = (CCMap*)(model.qHash);
	copasi_compartment c = { 0, 0 , 0};

	if (!pModel || !hash) return c;

	c.CopasiModelPtr = (void*)(pModel);
	c.qHash = model.qHash;

	if (contains(hash,string(name)))
	{
		if (getHashValue(hash, string(name)).compartment)
			c.CopasiCompartmentPtr = (void*)(getHashValue(hash, string(name)).compartment);

		return c;
	}
	else
	{
		c.CopasiCompartmentPtr = (void*)(pCompartment);
	}

	CopasiPtr copasiPtr = {
			pCompartment->getCN(),
			pCompartment->getKey(),
			0,
			pCompartment,
			0,
			0,
			"",
			false};

	hashInsert(hash, string(name),copasiPtr); //for speedy lookup

	return c;
}

int cSetValue(copasi_model model, const char * name, double value)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);

	if (!hash) return 0;

	if (!contains(hash,s))
	{
		cSetGlobalParameter(model,name,value);
		return 0;
	}

	CopasiPtr & p = getHashValue(hash,s);

	if (p.compartment)
	{
		p.compartment->setInitialValue(value);
		p.compartment->setValue(value);
		return 1;
	}
	else
	if (p.species)
	{
		p.species->setConcentration(value);
		p.species->setValue(value);
		p.species->setInitialValue(value);
		p.species->setInitialConcentration(value);
		return 1;
	}
	else
	if (p.param)
	{
		p.param->setInitialValue(value);
		p.param->setValue(value);
		return 1;
	}

	cSetGlobalParameter(model,name,value);
	return 0;
}

double cGetValue(copasi_model model, const char * name)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	double value = 0.0;

	if (!hash || !contains(hash,s))
	{
		return value;
	}

	CopasiPtr & p = getHashValue(hash,s);

	if (p.compartment)
	{
		value = p.compartment->getValue();
	}
	else
	if (p.species)
	{
		value = p.species->getConcentration();
	}
	else
	if (p.param)
	{
		value = p.param->getInitialValue();
	}

	return value;
}

void cSetVolume(copasi_model model, const char * name, double vol)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CCompartment* pVol = NULL;

	if (!hash) return;

	if (contains(hash,s) &&
		(pVol = getHashValue(hash,s).compartment))
	{
		pVol->setInitialValue(vol);
		pVol->setValue(vol);
	}
}

void cSetSpeciesConcentration(copasi_model model, const char * name, double conc)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CMetab* pSpecies = NULL;

	if (!hash) return;

	if (contains(hash,s) &&
		(pSpecies = getHashValue(hash,s).species))
	{
		pSpecies->setConcentration(conc);
	}
}

void cSetInitialConcentration(copasi_model model, const char * name, double conc)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CMetab* pSpecies = NULL;

	if (!hash) return;

	if (contains(hash,s) &&
		(pSpecies = getHashValue(hash,s).species))
	{
		pSpecies->setInitialConcentration(conc);
	}
}

void cSetSpeciesAmount(copasi_model model, const char * name, double amnt)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CMetab* pSpecies = NULL;

	if (!hash) return;

	if (contains(hash,s) &&
		(pSpecies = getHashValue(hash,s).species))
	{
		pSpecies->setValue(amnt);
		pSpecies->setInitialValue(amnt);
	}
}

int cSetGlobalParameter(copasi_model model, const char * name, double value)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CModelValue * pValue = NULL;

	if (!hash || !pModel) return 0;

	if (contains(hash,s) &&
		(pValue = getHashValue(hash,s).param))
	{
		pValue->setInitialValue(value);
		pValue->setValue(value);
		return 1;
	}

	//parameter not found, so create it
	if (!pValue)
	{
		pValue = pModel->createModelValue(string(name),value);
		pValue->setInitialValue(value);

		CopasiPtr copasiPtr = {
				pValue->getCN(),
				pValue->getKey(),
				0,
				0,
				0,
				pValue,
				"",
				true};

		hashInsert(hash, s, copasiPtr); //for speedy lookup
	}

	return 0;
}

void cSetSpeciesType(copasi_model model, const char * name, int isBoundary)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CMetab* pSpecies = NULL;

	if (!hash) return;

	if (contains(hash,s) &&
		(pSpecies = getHashValue(hash,s).species))
	{
		if (isBoundary)
			pSpecies->setStatus(CModelEntity::FIXED);
		else
			pSpecies->setStatus(CModelEntity::REACTIONS);
	}
}

void cCreateSpecies(CModel * pModel, CCMap * hash, string s)
{
	if (!pModel || !hash) return;

	CCopasiVectorNS < CCompartment > & compartments = pModel->getCompartments();
		if (compartments.size() > 0 && compartments[0] != NULL)
		{
			CCompartment* pCompartment = compartments[0];
			if (pCompartment)
			{
				string s(pCompartment->getObjectName());
				if	(contains(hash,s) &&
					getHashValue(hash,s).compartment)
					{
						copasi_compartment c = { (void*)getHashValue(hash,s).compartment, pModel, hash };
						cCreateSpecies(c,s.c_str(),0.0);
					}
			}
		}
}

int cSetAssignmentRule(copasi_model model, const char * name, const char * formula)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	int i;
	bool retval=true;

	if (!pModel || !hash) return 0;

	if (!contains(hash,s))
	{
		cCreateSpecies(pModel, hash, s);
	}

	if (contains(hash,s) && getHashValue(hash,s).species)
	{
		CopasiPtr & p = (*hash)[s];
		p.assignmentRule = string(formula);
		p.assignmentRule = boost::regex_replace(p.assignmentRule, sbmlPowFunction, string("((\\1)^(\\2))"));
		return 1;
	}
	return 0;
}

std::string cSetAssignmentRuleHelper(copasi_model model, const char * formula)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	int i;
	bool retval=true;

	if (!pModel || !hash) return string();

	if (formula)
	{
		CFunction pFunction;
		string qFormula(formula);
		if (pFunction.setInfix(string(formula)))
		{
			CFunctionParameters& variables = pFunction.getVariables();
			CFunctionParameter* pParam;

			for (i=0; i < variables.size(); ++i)
			{
				pParam = variables[i];

				string s0(pParam->getObjectName());

				if (s0 == string("time") ||
					  s0 == string("Time") ||
			     	  s0 == string("TIME"))
				{
					string s1("<");
						s1 += pModel->getValueReference()->getCN();
						s1 += string(">");
					rename(qFormula,s0,s1);
				}
				else
				{
					if (!contains(hash,s0))
						cSetGlobalParameter(model,pParam->getObjectName().c_str(),1.0);
					if (contains(hash,s0))
					{
					 	string s1("<");
							s1 += getHashValue(hash,s0).name;
							s1 += string(">");
						rename(qFormula,s0,s1);
					}
				}
			}
		}
		return qFormula;
	}

	return string();
}

int cCreateVariable(copasi_model model, const char * name, const char * formula)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!hash || !pModel) return 0;

	CModelValue* pModelValue;
	string qname(name);

	if (contains(hash,qname))
	{
			CopasiPtr & ptr = getHashValue(hash,qname);
			if (ptr.species)
				return cSetAssignmentRule(model, name, formula);
			if (ptr.param)
				pModelValue = ptr.param;
			else
				return 0;
	}
	else
	{
		pModelValue = pModel->createModelValue(string(name), 0.0);
	}
	pModelValue->setStatus(CModelValue::ASSIGNMENT);
	int i;
	bool retval = true;

	CFunction pFunction;
	string qFormula(formula);

	if (pFunction.setInfix(string(formula)))
	{
		CFunctionParameters& variables = pFunction.getVariables();
		CFunctionParameter* pParam;

		for (i=0; i < variables.size(); ++i)
		{
			pParam = variables[i];

			string s0(pParam->getObjectName());
			if (s0 == string("time") ||
				  s0 == string("Time") ||
		     	  s0 == string("TIME"))
			{
				string s1("<");
					s1 += pModel->getValueReference()->getCN();
					s1 += string(">");
				rename(qFormula,s0,s1);
			}
			else
			{
				if (!contains(hash,s0))
					cSetGlobalParameter(model,pParam->getObjectName().c_str(),1.0);

				if (contains(hash,s0))
				{
				 	string s1("<");
						s1 += getHashValue(hash,s0).name;
						s1 += string(">");
					rename(qFormula,s0,s1);
				}
			}
		}
	}

	string sFormula( qFormula );

	retval = retval & pModelValue->setInitialExpression(sFormula);
	retval = retval & pModelValue->setExpression(sFormula);

	CopasiPtr copasiPtr = {
			pModelValue->getCN(),
			pModelValue->getKey(),
			0,
			0,
			0,
			pModelValue,
			"",
			true};

	hashInsert(hash, qname, copasiPtr); //for speedy lookup

	return (int)retval;
}

int cCreateEvent(copasi_model model, const char * name, const char * trigger, const char * variable, const char * formula)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!hash || !pModel)
	{
		return 0;
	}

	int i;
	bool retval = true;

	if (!contains(hash,string(variable)))
	{
		cSetGlobalParameter(model,variable,1.0);
	}

	if (!contains(hash,string(variable))) return 0;

	CopasiPtr & ptr = getHashValue(hash,string(variable));

	if (!ptr.species && !ptr.param) return 0;

	CEvent * pEvent = pModel->createEvent(string(name));

	CFunction pFunction;
	string qFormula(trigger);
	replaceSubstring(qFormula,">"," gt ");
	replaceSubstring(qFormula,"<"," lt ");
	replaceSubstring(qFormula,">="," ge ");
	replaceSubstring(qFormula,"<="," le ");
	replaceSubstring(qFormula,"="," eq ");

	if (pFunction.setInfix(qFormula))  //parse trigger
	{
		CFunctionParameters& variables = pFunction.getVariables();
		CFunctionParameter* pParam;

		for (i=0; i < variables.size(); ++i)
		{
			pParam = variables[i];

			string s0(pParam->getObjectName());
			if (s0 == string("time") ||
				  s0 == string("Time") ||
		     	  s0 == string("TIME"))
			{
				string s1("<");
					s1 += pModel->getValueReference()->getCN();
					s1 += string(">");
				rename(qFormula,s0,s1);
			}
			else
			{
				if (!contains(hash,s0))
					cSetGlobalParameter(model,s0.c_str(),1.0);
				if (contains(hash,s0))
				{
				 	string s1("<");
						s1 += getHashValue(hash,s0).name;
						s1 += string(">");
					rename(qFormula,s0,s1);
				}
			}
		}
	}
	else
	{
		retval = false;
	}

	CExpression * expression = new CExpression(name,pModel);
	retval = retval & expression->setInfix(qFormula);
	pEvent->setTriggerExpressionPtr(expression);   //set trigger

	qFormula = string(formula);

	if (pFunction.setInfix(string(formula)))   //parse response expression
	{
		CFunctionParameters& variables = pFunction.getVariables();
		CFunctionParameter* pParam;

		for (i=0; i < variables.size(); ++i)
		{
			pParam = variables[i];

			string s0(pParam->getObjectName().c_str());
			if (s0 == string("time") ||
				  s0 == string("Time") ||
			 	  s0 == string("TIME"))
			{
				string s1("<");
					s1 += pModel->getValueReference()->getCN();
					s1 += string(">");
				rename(qFormula,s0,s1);
			}
			else
			{
				if (!contains(hash,s0))
					cSetGlobalParameter(model,s0.c_str(),1.0);

				if (contains(hash,s0))
				{
				 	string s1("<");
						s1 += getHashValue(hash,s0).name;
						s1 += string(">");
					rename(qFormula,s0,s1);
				}
			}
		}
	}
	else
	{
		return 0;
	}

	CCopasiVectorN< CEventAssignment > & assignments = pEvent->getAssignments();
	CEventAssignment * assgn = new CEventAssignment(string(name) + string("_assgn"),pModel);
	if (ptr.species)
		retval = retval & assgn->setTargetKey(ptr.species->getKey());   //set target
	else
		retval = retval & assgn->setTargetKey(ptr.param->getKey());

	expression = new CExpression(name,pModel);
	retval = retval & expression->setInfix(qFormula);
	assgn->setExpressionPtr(  expression );   //set expression
	assignments.add(assgn);

	return (int)retval;
}

copasi_reaction cCreateReaction(copasi_model model, const char* name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel || !hash)
	{
		copasi_reaction r = { 0, 0, 0 };
		return r;
	}

	CReaction* pReaction = pModel->createReaction(name);

	copasi_reaction r = { (void*)(pReaction), (void*)(pModel), (void*)hash };

	CopasiPtr copasiPtr = {
			pReaction->getCN(),
			pReaction->getKey(),
			0,
			0,
			pReaction,
			0,
			"",
			false};

	string qname(name);
	hashInsert(hash, qname, copasiPtr); //for speedy lookup

	return r;
}

void cAddReactant(copasi_reaction reaction, const char * species, double stoichiometry)
{
	CReaction* pReaction = (CReaction*)(reaction.CopasiReactionPtr);
	CCMap * hash = (CCMap*)(reaction.qHash);

	if (!pReaction || !hash)
	{
		return;
	}

	CMetab* pSpecies = NULL;

	string s(species);
	if (!contains(hash, s))
	{
		cCreateSpecies((CModel*)reaction.CopasiModelPtr, (CCMap*)reaction.qHash, s);
	}

	if (contains(hash,s))
	{
		CopasiPtr & p = getHashValue(hash,s);
		if (pSpecies = p.species)
		{
			CChemEq* pChemEq = &pReaction->getChemEq();
			pChemEq->addMetabolite(pSpecies->getKey(), stoichiometry, CChemEq::SUBSTRATE);
			p.unused = false;
		}
	}

	if (pSpecies && pSpecies->getCompartment())
	{
		s = pSpecies->getCompartment()->getObjectName() + string("_") + s;
		if (contains(hash,s))
		{
			CopasiPtr & p = getHashValue(hash,s);
			if (pSpecies = p.species)
				p.unused = false;
		}
	}
}

void cAddProduct(copasi_reaction reaction, const char * species, double stoichiometry)
{
	CReaction* pReaction = (CReaction*)(reaction.CopasiReactionPtr);
	CCMap * hash = (CCMap*)(reaction.qHash);
	CMetab* pSpecies = NULL;

	if (!pReaction || !hash || !reaction.CopasiModelPtr) return;

	string s(species);
	if (!contains(hash, s))
	{
		cCreateSpecies((CModel*)reaction.CopasiModelPtr, (CCMap*)reaction.qHash, s);
	}

	if (contains(hash,s))
	{
		CopasiPtr & p = getHashValue(hash,s);
		if (pSpecies = p.species)
		{
			CChemEq* pChemEq = &pReaction->getChemEq();
			pChemEq->addMetabolite(pSpecies->getKey(), stoichiometry, CChemEq::PRODUCT);
			p.unused = false;
		}
	}

	if (pSpecies && pSpecies->getCompartment())
	{
		s = pSpecies->getCompartment()->getObjectName() + string("_") + s;
		if (contains(hash,s))
		{
			CopasiPtr & p = getHashValue(hash,s);
			if (pSpecies = p.species)
				p.unused = false;
		}
	}
}

int cSetReactionRate(copasi_reaction reaction, const char * formula)
{
	int i,j,k;
	CReaction* pReaction = (CReaction*)(reaction.CopasiReactionPtr);
	CCMap * hash = (CCMap*)(reaction.qHash);
	CModel* pModel = (CModel*)(reaction.CopasiModelPtr);
	CFunctionDB* pFunDB = CCopasiRootContainer::getFunctionList();

	if (!pReaction || !pModel || !hash) return 0;

	if (pFunDB)
	{
		string rateLawName(pReaction->getObjectName() + string("_rate_law")); //existing rate law

		CFunction * pFunction = dynamic_cast<CFunction*>(pFunDB->findFunction(rateLawName));
		if (pFunction)
			return (int)(pReaction->setFunction(pFunction)) - 1;

		CKinFunction* pKinFunction = new CKinFunction(rateLawName);
		pFunDB->add(pKinFunction, true);
		pFunction = pKinFunction;//dynamic_cast<CFunction*>(pFunDB->findFunction(rateLawName));

		if (!pFunction)
			return 0;

		pFunction->setReversible(TriFalse);

		int retval = 0;

		string formula2(formula);
		//formula2.replace(sbmlPowFunction, QString("((\\1)^(\\2))"));
		formula2 = boost::regex_replace(formula2, sbmlPowFunction, string("((\\1)^(\\2))"));

		if (pFunction->setInfix(string(formula2.c_str())))
		{
			retval = (int)(pReaction->setFunction(pFunction));
			CFunctionParameters& variables = pFunction->getVariables();
			CFunctionParameter* pParam;

			for (i=0; i < variables.size(); ++i)
			{
				pParam = variables[i];

				string s(pParam->getObjectName().c_str());

				if (s == string("Time") || s == string("TIME") )
					s = string("time");

				if (!contains(hash,s))
				{
					copasi_model model = { (void*)(pModel) , (void*)(0), (void*)(hash), (char*)(NULL), (char*)(NULL)};
					cSetGlobalParameter(model,pParam->getObjectName().c_str(),1.0);
				}

				if (contains(hash,s))
				{
					CopasiPtr & p = getHashValue(hash,s);
					if (p.compartment)
					{
						pParam->setUsage(CFunctionParameter::VOLUME);
						pReaction->setParameterMapping(pParam->getObjectName(), p.compartment->getKey());
					}
					else
					if (p.species)
					{
						pParam->setUsage(CFunctionParameter::MODIFIER);
						const CCopasiVector < CChemEqElement > & substrates = pReaction->getChemEq().getSubstrates();
						for (k =0; k < substrates.size(); ++k)
							if (substrates[k]->getMetabolite() == p.species)
							{
								pParam->setUsage(CFunctionParameter::SUBSTRATE);
								break;
							}
						pReaction->setParameterMapping(pParam->getObjectName(), p.species->getKey());
					}
					else
					if (p.param)
					{
						pParam->setUsage(CFunctionParameter::PARAMETER);
						pReaction->setParameterMapping(pParam->getObjectName(), p.param->getKey());
					}
				}
			}

			pFunction->compile();

			return retval;
		}
	}

	return 0;
}

void cCompileModel(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);

	if (!pModel) return;

	copasi_cleanup_assignments(model);

	CCopasiVectorNS < CCompartment > & compartments = pModel->getCompartments();
	CCopasiVector< CMetab > & species = pModel->getMetabolites();
	CCopasiVectorN< CModelValue > & params = pModel->getModelValues();
	const CCopasiObject* pObject = NULL;
	set<const CCopasiObject*> changedObjects;

	for (int i=0; i < compartments.size(); ++i)
		if (compartments[i])
		{
			pObject = compartments[i]->getObject(CCopasiObjectName("Reference=InitialVolume"));
			if (pObject)
				changedObjects.insert(pObject);
		}

	for (int i=0; i < species.size(); ++i)
		if (species[i])
		{
			pObject = species[i]->getObject(CCopasiObjectName("Reference=InitialConcentration"));
			if (pObject)
				changedObjects.insert(pObject);
		}

	for (int i=0; i < params.size(); ++i)
		if (params[i])
		{
			pObject = params[i]->getObject(CCopasiObjectName("Reference=Value"));
			if (pObject)
				changedObjects.insert(pObject);

			pObject = params[i]->getObject(CCopasiObjectName("Reference=InitialValue"));
			if (pObject)
				changedObjects.insert(pObject);
		}

	// compile needs to be done before updating all initial values for
	// the model with the refresh sequence
	pModel->compileIfNecessary(NULL);

	// now that we are done building the model, we have to make sure all
	// initial values are updated according to their dependencies
	vector<Refresh*> refreshes = pModel->buildInitialRefreshSequence(changedObjects);

	vector<Refresh*>::iterator it2 = refreshes.begin(), endit2 = refreshes.end();

	while (it2 != endit2)
	{
		// call each refresh
		(**it2)();
		++it2;
	}
}

c_matrix simulate(copasi_model model, double startTime, double endTime, int numSteps, CCopasiMethod::SubType method)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the trajectory task object
	CTrajectoryTask* pTask = dynamic_cast<CTrajectoryTask*>(TaskList["Time-Course"]);
	// if there isn’t one
	if (pTask == NULL)
	{
		// create a new one
		pTask = new CTrajectoryTask();
		// remove any existing trajectory task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Time-Course");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CState state = pModel->getInitialState();

	CCopasiMessage::clearDeque();

	if (startTime >= endTime)
		endTime += startTime;

	if (pTask && pTask->setMethodType(method))
	{
		//set the start and end time, number of steps, and save output in memory
		CTrajectoryProblem* pProblem=(CTrajectoryProblem*)pTask->getProblem();
		pProblem->setModel(pModel);
		pTask->setScheduled(true);
		pProblem->setStepNumber(numSteps);
		pProblem->setDuration(endTime-startTime);
		pDataModel->getModel()->setInitialTime(startTime);
		pProblem->setTimeSeriesRequested(true);
		pTask->setUpdateModel(true);
		try
		{
			pTask->initialize(CCopasiTask::ONLY_TIME_SERIES, pDataModel, NULL);
			pTask->process(true);
			pTask->restore();
		}
		catch(...)
		{
			cerr << "Error. Running the simulation failed." << endl;
			// check if there are additional error messages
			if (CCopasiMessage::size() > 0)
			{
				// print the messages in chronological order
				cerr << CCopasiMessage::getAllMessageText(true);
			}
			pTask = NULL;
		}
	}

	if (pTask)
	{
		const CTimeSeries & timeSeries = pTask->getTimeSeries();
		int rows = timeSeries.getRecordedSteps(),
			  cols = (1+pModel->getNumMetabs());//timeSeries.getNumVariables();
		int i,j,k;

		c_matrix output = efficiently_createMatrix(rows, cols);
		list<string> colnames;

		for (j=1; j < cols; ++j)
			colnames.push_back( timeSeries.getTitle(j) );

		colnames.sort();
		colnames.push_front(timeSeries.getTitle(0).c_str());

		j = 0;
		for (list<string>::iterator it=colnames.begin(); j < cols && it != colnames.end(); ++j, it++)
			c_setColumnName( output, j, (*it).c_str() );

		for (j=0; j < cols; ++j)
		{
			k = indexOf(colnames,timeSeries.getTitle(j));
			for (i=0; i < rows; ++i)
				c_setMatrixValue( output, i, k, timeSeries.getConcentrationData(i,j) );
		}

		pModel->setInitialState(state);
		return output;
	}
	return c_createMatrix(0,0);
}

c_matrix cSimulateDeterministic(copasi_model model, double startTime, double endTime, int numSteps)
{
	return simulate(model,startTime,endTime,numSteps,CCopasiMethod::deterministic);
}

// STUB: NEEDS TO BE IMPLEMENTED
double cOneStep(copasi_model model, double timeStep)
{
	return 0.0;
}

c_matrix cSimulateTauLeap(copasi_model model, double startTime, double endTime, int numSteps)
{
	return simulate(model,startTime,endTime,numSteps,CCopasiMethod::tauLeap);
}

c_matrix cSimulateStochastic(copasi_model model, double startTime, double endTime, int numSteps)
{
	return simulate(model,startTime,endTime,numSteps,CCopasiMethod::stochastic);
}

c_matrix cSimulateHybrid(copasi_model model, double startTime, double endTime, int numSteps)
{
	return simulate(model,startTime,endTime,numSteps,CCopasiMethod::hybridLSODA);
}

static int _SBML_LEVEL = 2;
static int _SBML_VERSION = 2;
void cSetSBMLLevelAndVersion(int level, int version)
{
	_SBML_LEVEL = level;
	_SBML_VERSION = version;
}

void cWriteSBMLFile(copasi_model model, const char * filename)
{
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (pDataModel)
		pDataModel->exportSBML(filename, true, _SBML_LEVEL, _SBML_VERSION);
}

void cWriteAntimonyFile(copasi_model model, const char * filename)
{
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (pDataModel)
	{
		pDataModel->exportSBML(filename, true, _SBML_LEVEL, _SBML_VERSION);
		loadSBMLFile(filename);
		writeAntimonyFile(filename,NULL);
	}
}

copasi_model cReadAntimonyString(const char * model)
{
	loadString(model); //load Antimony
	const char * s = getSBMLString("__main");  //Antimony -> SBML (at worst, an empty model)
	copasi_model m = cReadSBMLString(s);
	freeAll(); //free Antimony
	return m;
}

copasi_model cReadAntimonyFile(const char * filename)
{
	loadFile(filename); //load Antimony
	const char * s = getSBMLString("__main");  //Antimony -> SBML (at worst, an empty model)
	copasi_model m = cReadSBMLString(s);
	freeAll(); //free Antimony
	return m;
}

//this helper class can load either sbml file or string, since both are identical except for one line
static copasi_model cReadSBML_helper(const char * sbml, bool isFile)
{
	copasi_init();

	CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();
	CModel* pModel = 0;
	CCMap * qHash = 0;
	char * error = NULL;
	char * warning = NULL;
	string s;
	CCopasiMessage::Type type;

	try
	{
		if (isFile)  //here, the code is different depending on file or string input
			pDataModel->importSBML(sbml); //SBML file -> COPASI
		else
			pDataModel->importSBMLFromString(sbml); //SBML string -> COPASI

		s = CCopasiMessage::getAllMessageText();
		type = CCopasiMessage::getHighestSeverity();
		pModel = pDataModel->getModel();
		qHash = new CCMap();
	}
	catch(...)
	{
		s = CCopasiMessage::getAllMessageText();
		type = CCopasiMessage::ERROR;
	}

	int len = s.length();
	if (len > 1)
	{
		char * msg = (char*)malloc((1+len) * sizeof(char));
		if (msg)
		{
			for (int i=0; i < len; ++i) msg[i] = s[i];
			msg[len-1] = 0;
		}

		//error or warning?
		if (type == CCopasiMessage::EXCEPTION || type == CCopasiMessage::ERROR)
			error = msg;
		else
			warning = msg;
	}

	copasi_model m = { (void*)(pModel) , (void*)(pDataModel), (void*)(qHash), (char*)(error), (char*)warning};
	if (pModel && qHash)
	{
		hashTablesToCleanup.push_back( qHash );
		copasi_modelsToCleanup.push_back(m);
		populate_hash(m);
	}
	return m;
}

copasi_model cReadSBMLFile(const char * filename)
{
	return cReadSBML_helper(filename, true);
}

copasi_model cReadSBMLString(const char * sbml)
{
	return cReadSBML_helper(sbml, false);
}

void cResetState(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	if (pModel)
		pModel->setState( pModel->getInitialState() );
}

c_matrix cGetJacobian(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	double epsilon = _EPSILON;
	CMatrix< C_FLOAT64 > jacobian(species.size(), species.size()); 
	pModel->calculateJacobian(jacobian, epsilon, epsilon);

	c_matrix J = efficiently_createMatrix(species.size(),species.size());
	for (int i=0; i < species.size(); ++i)
		for (int j=0; j < species.size(); ++j)
			c_setMatrixValue(J, i, j, jacobian[i][j]);
	
	return J;

	//cCompileModel(model);
	/*
	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the steady state task object
	CSteadyStateTask* pTask = dynamic_cast<CSteadyStateTask*>(TaskList["Steady-State"]);
	// if there isn’t one
	if (pTask == NULL)
	{
		// create a new one
		pTask = new CSteadyStateTask();
		// remove any existing steady state task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Steady-State");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		// now we run the actual trajectory
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when computing steady state." << endl;
		return c_createMatrix(0,0);
	}

	const CArrayAnnotation* pAJ = pTask->getJacobianAnnotated();
	//const CEigen & cGetEigenvalues() const;

	if (pAJ && pAJ->dimensionality() == 2)
	{
		vector<unsigned int> index(2);
		const vector<string>& annotations = pAJ->getAnnotationsString(1);

		int n = annotations.size();
		c_matrix J = efficiently_createMatrix(n,n);

		for (int i=0; i < J.rows; ++i)
		{
			c_setRowName(J, i, annotations[i].c_str());
			c_setColumnName(J, i, annotations[i].c_str());
		}

		for (int i=0; i < n; ++i)
		{
			index[0] = i;
			for (int j=0; j < n; ++j)
			{
				index[1] = j;
				c_setMatrixValue(J, i, j, (*pAJ->array())[index]);
			}
		}

		return J;
	}*/

	return c_createMatrix(0,0);
}

c_matrix cGetSteadyStateUsingSimulation(copasi_model model, int maxiter)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);

	//cCompileModel(model);

    int iter = 0;
    double err = 2.0, eps = _EPSILON, time = 10.0;

   	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	CTrajectoryTask* pTask = dynamic_cast<CTrajectoryTask*>(TaskList["Time-Course"]);
	// if there isn’t one
	if (pTask == NULL)
	{
		pTask = new CTrajectoryTask();
		TaskList.remove("Time-Course");
		TaskList.add(pTask, true);
	}

	CState state = pModel->getInitialState();

	if (pTask)
		pTask->setUpdateModel(true);

    while (iter < maxiter && err > eps)
    {
        ++iter;
        time *= 2.0;

	   CCopasiMessage::clearDeque();

	    if (pTask && pTask->setMethodType(CCopasiMethod::deterministic))
	    {
		    CTrajectoryProblem* pProblem=(CTrajectoryProblem*)pTask->getProblem();
		    pProblem->setModel(pModel);
		    pTask->setScheduled(true);
		    pProblem->setStepNumber(int(time * 2.0));
		    pProblem->setDuration(time);
		    pDataModel->getModel()->setInitialTime(0.0);
		    pProblem->setTimeSeriesRequested(true);
		    try
		    {
			    pTask->initialize(CCopasiTask::ONLY_TIME_SERIES, pDataModel, NULL);
			    pTask->process(true);
                pTask->restore();
		    }
		    catch(...)
		    {
			    cerr << CCopasiMessage::getAllMessageText(true);
			    pTask = NULL;
		    }
	    }

        if (pTask)
        {
            const CTimeSeries & timeSeries = pTask->getTimeSeries();
            int cols = (pModel->getNumMetabs());
            int j = timeSeries.getRecordedSteps() - 1;
            double diff;
            err = 0.0;
            if (j < 1)
                err = eps * 2.0;
            else
            {
                for (int i=1; i <= cols; ++i)
                {
                    diff = timeSeries.getConcentrationData(j,i) - timeSeries.getConcentrationData(j-1,i);
                    err += diff * diff;
                }
                err /= cols;
            }

            if (err < eps)
            {
				list<string> colnames;
				for (int i=0; i < cols; ++i)
                	colnames.push_back( timeSeries.getTitle(i+1) );
				colnames.sort();
                c_matrix output = efficiently_createMatrix(cols, 1);
				int k;
				list<string>::iterator it=colnames.begin();
                for (int i=0; i < cols && it != colnames.end(); ++i, it++)
                {
					k = indexOf(colnames, timeSeries.getTitle(i+1));
            		c_setRowName( output, i, (*it).c_str() );
                    c_setMatrixValue( output, k, 0, timeSeries.getConcentrationData(j,i+1) );
                }
				pModel->setInitialState(state);
                return output;
            }
        }
	}

    c_matrix m = c_createMatrix(0,0);
	return m;
}

c_matrix cGetSteadyState(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);

	//cCompileModel(model);

	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	CTrajectoryTask* pTrajTask = dynamic_cast<CTrajectoryTask*>(TaskList["Time-Course"]);
	// if there isn’t one
	if (pTrajTask == NULL)
	{
		// create a new one
		pTrajTask = new CTrajectoryTask();
		// remove any existing trajectory task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Time-Course");
		// add the new time course task to the task list
		TaskList.add(pTrajTask, true);
	}

	CCopasiMessage::clearDeque();
	CState state = pModel->getInitialState();

	if (pTrajTask && pTrajTask->setMethodType(CCopasiMethod::deterministic))
	{
		//set the start and end time, number of steps, and save output in memory
		CTrajectoryProblem* pProblem=(CTrajectoryProblem*)pTrajTask->getProblem();
		pProblem->setModel(pModel);
		pTrajTask->setScheduled(true);
		pProblem->setStepNumber(10);
		pProblem->setDuration(100.0);
		pDataModel->getModel()->setInitialTime(0.0);
		pProblem->setTimeSeriesRequested(true);
		pTrajTask->setUpdateModel(true);
		try
		{
			pTrajTask->initialize(CCopasiTask::ONLY_TIME_SERIES, pDataModel, NULL);
			pTrajTask->process(true);
			pTrajTask->restore();
		}
		catch(...)
		{
			cerr << CCopasiMessage::getAllMessageText(true);
			pTrajTask = NULL;
		}
	}

	CSteadyStateTask* pTask = dynamic_cast<CSteadyStateTask*>(TaskList["Steady-State"]);

	if (pTask == NULL)
	{
		pTask = new CSteadyStateTask();
		TaskList.remove("Steady-State");
		TaskList.add(pTask, true);
	}

	try
	{
		pTask->setUpdateModel(true);
		pTask->initialize(CCopasiTask::OUTPUT, pDataModel, NULL);
		pTask->process(true);
		pTask->restore();
	}
	catch (...)
	{
		cerr << "Error when computing steady state." << endl;
		return c_createMatrix(0,0);
	}

	CCopasiMessage::clearDeque();

	c_matrix m = cGetFloatingSpeciesConcentrations(model);
	pModel->setInitialState(state);

	return m;
}

c_matrix cGetEigenvalues(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the steady state task object
	CSteadyStateTask* pTask = dynamic_cast<CSteadyStateTask*>(TaskList["Steady-State"]);
	// if there isn’t one
	if (!pTask)
	{
		// create a new one
		pTask = new CSteadyStateTask();
		// remove any existing steady state task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Steady-State");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		// now we run the actual trajectory
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when computing steady state." << endl;
		return c_createMatrix(0,0);
	}

	const CEigen & eigen = pTask->getEigenValues();
	const CVector< C_FLOAT64 > & im = eigen.getI(),
													& re = eigen.getR();

	c_matrix E = efficiently_createMatrix(im.size(),2);

	c_setColumnName(E, 0, "real\0");
	c_setColumnName(E, 1, "imaginary\0");
	for (int i=0; i < im.size() && i < re.size(); ++i)
	{
		c_setMatrixValue(E, i,0,re[i]);
		c_setMatrixValue(E, i,1,im[i]);
	}

	return E;
}

c_matrix cGetUnscaledElasticities(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the MCA task object
	CMCATask* pTask = dynamic_cast<CMCATask*>(TaskList["Metabolic Control Analysis"]);
	// if there isn’t one
	if (!pTask)
	{
		// create a new one
		pTask = new CMCATask();
		// remove any existing steady state task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Metabolic Control Analysis");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when performing MCA" << endl;
		return c_createMatrix(0,0);
	}

	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());

	if (!mcaMethod) return c_createMatrix(0,0);

	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getUnscaledElasticities();
	const CArrayAnnotation * annot = mcaMethod->getUnscaledElasticitiesAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	c_matrix M = efficiently_createMatrix(rows, cols);

	for (int i=0; i < rows; ++i)
		c_setRowName(M, i, rownames[i].c_str());

	for (int i=0; i < cols; ++i)
		c_setColumnName(M, i, colnames[i].c_str());

	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			c_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

c_matrix cGetUnscaledConcentrationControlCoeffs(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the MCA task object
	CMCATask* pTask = dynamic_cast<CMCATask*>(TaskList["Metabolic Control Analysis"]);
	// if there isn’t one
	if (!pTask)
	{
		// create a new one
		pTask = new CMCATask();
		// remove any existing steady state task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Metabolic Control Analysis");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when performing MCA" << endl;
		return c_createMatrix(0,0);
	}

	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());

	if (!mcaMethod) return c_createMatrix(0,0);

	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getUnscaledConcentrationCC();
	const CArrayAnnotation * annot = mcaMethod->getUnscaledConcentrationCCAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	c_matrix M = efficiently_createMatrix(rows, cols);

	for (int i=0; i < rows; ++i)
		c_setRowName(M, i, rownames[i].c_str());

	for (int i=0; i < cols; ++i)
		c_setColumnName(M, i, colnames[i].c_str());

	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			c_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

c_matrix cGetUnscaledFluxControlCoeffs(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the MCA task object
	CMCATask* pTask = dynamic_cast<CMCATask*>(TaskList["Metabolic Control Analysis"]);
	// if there isn’t one
	if (!pTask)
	{
		// create a new one
		pTask = new CMCATask();
		// remove any existing steady state task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Metabolic Control Analysis");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when performing MCA" << endl;
		return c_createMatrix(0,0);
	}

	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());

	if (!mcaMethod) return c_createMatrix(0,0);

	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getUnscaledFluxCC();
	const CArrayAnnotation * annot = mcaMethod->getUnscaledFluxCCAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	c_matrix M = efficiently_createMatrix(rows, cols);

	for (int i=0; i < rows; ++i)
		c_setRowName(M, i, rownames[i].c_str());

	for (int i=0; i < cols; ++i)
		c_setColumnName(M, i, colnames[i].c_str());

	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			c_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

c_matrix cGetScaledElasticities(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the MCA task object
	CMCATask* pTask = dynamic_cast<CMCATask*>(TaskList["Metabolic Control Analysis"]);
	// if there isn’t one
	if (!pTask)
	{
		// create a new one
		pTask = new CMCATask();
		// remove any existing steady state task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Metabolic Control Analysis");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when performing MCA" << endl;
		return c_createMatrix(0,0);
	}

	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());

	if (!mcaMethod) return c_createMatrix(0,0);

	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getScaledElasticities();
	const CArrayAnnotation * annot = mcaMethod->getScaledElasticitiesAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	c_matrix M = efficiently_createMatrix(rows, cols);

	for (int i=0; i < rows; ++i)
		c_setRowName(M, i, rownames[i].c_str());

	for (int i=0; i < cols; ++i)
		c_setColumnName(M, i, colnames[i].c_str());

	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			c_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

c_matrix cGetScaledConcentrationConcentrationCoeffs(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the MCA task object
	CMCATask* pTask = dynamic_cast<CMCATask*>(TaskList["Metabolic Control Analysis"]);
	// if there isn’t one
	if (!pTask)
	{
		// create a new one
		pTask = new CMCATask();
		// remove any existing steady state task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Metabolic Control Analysis");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when performing MCA" << endl;
		return c_createMatrix(0,0);
	}

	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());

	if (!mcaMethod) return c_createMatrix(0,0);

	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getScaledConcentrationCC();
	const CArrayAnnotation * annot = mcaMethod->getScaledConcentrationCCAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	c_matrix M = efficiently_createMatrix(rows, cols);

	for (int i=0; i < rows; ++i)
		c_setRowName(M, i, rownames[i].c_str());

	for (int i=0; i < cols; ++i)
		c_setColumnName(M, i, colnames[i].c_str());

	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			c_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

c_matrix cGetScaledFluxControlCoeffs(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the MCA task object
	CMCATask* pTask = dynamic_cast<CMCATask*>(TaskList["Metabolic Control Analysis"]);
	// if there isn’t one
	if (!pTask)
	{
		// create a new one
		pTask = new CMCATask();
		// remove any existing steady state task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Metabolic Control Analysis");
		// add the new time course task to the task list
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when performing MCA" << endl;
		return c_createMatrix(0,0);
	}

	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());

	if (!mcaMethod) return c_createMatrix(0,0);

	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getScaledFluxCC();
	const CArrayAnnotation * annot = mcaMethod->getScaledFluxCCAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	size_t rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	c_matrix M = efficiently_createMatrix(rows, cols);

	for (int i=0; i < rows; ++i)
		c_setRowName(M, i, rownames[i].c_str());

	for (int i=0; i < cols; ++i)
		c_setColumnName(M, i, colnames[i].c_str());

	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			c_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

/** Sloten from TinkerCell  **/

static int rename(string& target, const string& oldname,const string& newname0)
{
	if (oldname == newname0 || target == newname0) return 0;
	string newname(newname0);

	boost::regex_replace(newname, boost::regex("[^A-Za-z0-9_]",boost::regex::perl), string("_@@@_"));
	//newname.replace(QRegExp("[^A-Za-z0-9_]"),QString("_@@@_"));

	boost::regex regexp1(string("^") + oldname + string("$"),boost::regex::perl),  //just old name
		regexp2(string("^") + oldname + string("([^A-Za-z0-9_\\]]).*"),boost::regex::perl),  //oldname+(!letter/num)
		regexp3(string(".*([^A-Za-z0-9_\\.=\\[])") + oldname + string("$"),boost::regex::perl), //(!letter/num)+oldname
		regexp4(string(".*([^A-Za-z0-9_\\.=\\[])") + oldname + string("([^A-Za-z0-9_\\]]).*"),boost::regex::perl); //(!letter/num)+oldname+(!letter/num)

	boost::regex regexp1r(string("^") + oldname + string("$"),boost::regex::perl),  //just old name
		regexp2r(string("^") + oldname + string("([^A-Za-z0-9_\\]])"),boost::regex::perl),  //oldname+(!letter/num)
		regexp3r(string("([^A-Za-z0-9_\\.=\\[])") + oldname + string("$"),boost::regex::perl), //(!letter/num)+oldname
		regexp4r(string("([^A-Za-z0-9_\\.=\\[])") + oldname + string("([^A-Za-z0-9_\\]])"),boost::regex::perl); //(!letter/num)+oldname+(!letter/num)

	int retval = 0;

	if (boost::regex_match(target, regexp1))
	{
		retval = 1;
		target = newname;
	}

	while (boost::regex_match(target, regexp2))
	{
		retval = 1;
		target = boost::regex_replace(target, regexp2r,newname+string("\\1"));
	}

	while (boost::regex_match(target, regexp3))
	{
		retval = 1;
		target = boost::regex_replace(target, regexp3r, string("\\1")+newname);
	}

	while (boost::regex_match(target, regexp4))
	{
		retval = 1;
		target = boost::regex_replace(target, regexp4r, string("\\1")+newname+string("\\2"));
	}

	target = boost::regex_replace(target, boost::regex(newname),newname0);
	return retval;
}

c_matrix cGetReducedStoichiometryMatrix(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	CCopasiVector< CMetab > & species = pModel->getMetabolitesX();
	CCopasiVectorNS < CReaction > & reacs = pModel->getReactions();
	CMatrix < C_FLOAT64 > stoi = pModel->getRedStoi();

	c_matrix N = c_createMatrix( stoi.numRows(), stoi.numCols() );

	for  (int i=0; i < N.rows && i < species.size(); ++i)
		if (species[i])
			c_setRowName(N, i, species[i]->getObjectName().c_str());

	for  (int i=0; i < N.cols && i < reacs.size(); ++i)
		if (reacs[i])
			c_setColumnName(N, i, reacs[i]->getObjectName().c_str());

	for  (int i=0; i < N.rows; ++i)
		for  (int j=0; j < N.cols; ++j)
			c_setMatrixValue(N, i, j, stoi(i,j));

	return N;
}

c_matrix cGetFullStoichiometryMatrix(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	CCopasiVector< CMetab > & species = pModel->getMetabolites();
	CCopasiVectorNS < CReaction > & reacs = pModel->getReactions();
	CMatrix < C_FLOAT64 > stoi = pModel->getStoi();

	c_matrix N = efficiently_createMatrix( stoi.numRows(), stoi.numCols() );

	for  (int i=0; i < N.rows && i < species.size(); ++i)
		if (species[i])
			c_setRowName(N, i, species[i]->getObjectName().c_str());

	for  (int i=0; i < N.cols && i < reacs.size(); ++i)
		if (reacs[i])
			c_setColumnName(N, i, reacs[i]->getObjectName().c_str());

	for  (int i=0; i < N.rows; ++i)
		for  (int j=0; j < N.cols; ++j)
			c_setMatrixValue(N, i, j, stoi(i,j));

	return N;
}

c_matrix cGetElementaryFluxModes(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return c_createMatrix(0,0);
	//cCompileModel(model);

	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();

	CEFMTask* pTask = 0;

	if (TaskList["Elementary Flux Modes"])
		pTask = dynamic_cast<CEFMTask*>(TaskList["Elementary Flux Modes"]);

	if (!pTask)
	{
		pTask = new CEFMTask();
		TaskList.remove("Elementary Flux Modes");
		TaskList.add(pTask, true);
	}

	CCopasiMessage::clearDeque();

	try
	{
		pTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when computing EFM" << endl;
		return c_createMatrix(0,0);
	}

	const vector< CFluxMode > & fluxModes = pTask->getFluxModes();
	CEFMProblem* pProblem = dynamic_cast<CEFMProblem*>(pTask->getProblem());

	if (!pProblem)
		return c_createMatrix(0,0);

	vector< const CReaction * > & reactions = pProblem->getReorderedReactions();
	c_matrix M = efficiently_createMatrix( reactions.size() , fluxModes.size() );
	for (int i=0; i < reactions.size(); ++i)
		c_setRowName(M, i, reactions[i]->getObjectName().c_str());

	for (int i=0; i < fluxModes.size(); ++i)
	{
		CFluxMode::const_iterator itMode = fluxModes[i].begin();
		CFluxMode::const_iterator endMode = fluxModes[i].end();
		for (; itMode != endMode; ++itMode)
			c_setMatrixValue( M, itMode->first, i, itMode->second);
	}
	return M;
}

void cFitModelToData(copasi_model model, const char * filename, c_matrix params, const char * method)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel || !pDataModel || !hash) return;

	string delim;

	ifstream file(filename);
	list<string> words;
	int numlines=0;

	if (file.is_open())
	{
		string line;
		getline(file, line);

		if (line[0] == '#')
			line = line.substr(1,line.length());

		if (contains(line,"\t"))
			delim = "\t";
		else
		if (contains(line,","))
			delim = ",";
		else
		if (contains(line,";"))
			delim = ";";
		else
			delim = " ";

		words = splitString(line, delim);

		while (file.good())
		{
			getline(file, line);
			++numlines;
		}

		file.close();
	}

std::cout << "delim = " << delim << "\n";

	//find the species from the header of the data file
	vector<CMetab*> targetSpecies(words.size(), (CMetab*)0);
	CopasiPtr copasiPtr;
	list<string>::iterator itr = words.begin();
	for (int i=0; i < words.size() && itr != words.end(); ++i, ++itr)
	{
		string s(*itr);
		remove_if(s.begin(), s.end(), ::isspace);
		if (contains(hash,s))
		{
			copasiPtr = getHashValue(hash,s);
			if (copasiPtr.species)
			{
				targetSpecies[i] = copasiPtr.species;
				cout << i << "  =  " << s << endl;
			}
		}
		else
			cout << "'" << s << "' not found \n";
	}

	//get the target parameters
	vector< CModelValue* > targetParams(params.rows, (CModelValue*)0);
	for (int i=0; i < params.rows; ++i)
	{
		string rowname(c_getRowName(params, i));
		if (contains(hash,rowname))
		{
			copasiPtr = getHashValue(hash,rowname);
			if (copasiPtr.param && copasiPtr.param->getStatus() != CModelValue::ASSIGNMENT)
			{
				targetParams[i] = copasiPtr.param;
				//cout << "good  " << i << "  =  " << rowname << endl;
			}
		}
	}

	// get the task object
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();


	// get the optim task object
	CFitTask* pFitTask = dynamic_cast<CFitTask*>(TaskList["Parameter Estimation"]);

	// if there isn’t one
	if (pFitTask == NULL)
	{
		// create a new one
		pFitTask = new CFitTask(CCopasiTask::parameterFitting,pModel);
		// remove any existing task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Parameter Estimation");
		// add the new task to the task list
		TaskList.add(pFitTask, true);
	}

	//set method
	string sMethod( method ), s;

	for (int i=0; s.size() > 1; ++i)
	{
        s = string(CCopasiMethod::SubTypeName[i]);
	    if (sMethod.compare( s ) == 0)
		{
			pFitTask->setMethodType((CCopasiMethod::SubType)i);
			sMethod = string("");
			break;
		}
	}

	// the method in a fit task is an instance of COptMethod or a subclass
	COptMethod* pFitMethod = dynamic_cast<COptMethod*>(pFitTask->getMethod());
	// the object must be an instance of COptMethod or a subclass
	CFitProblem* pFitProblem = dynamic_cast<CFitProblem*>(pFitTask->getProblem());
	pFitProblem->setModel(pModel);
	CExperimentSet* pExperimentSet = dynamic_cast<CExperimentSet*>(pFitProblem->getParameter("Experiment Set"));
	// first experiment (we only have one here)
	CExperiment* pExperiment = new CExperiment(pDataModel);
	// tell COPASI where to find the data
	// reading data from string is not possible with the current C++ API
	pExperiment->setFileName(filename);
	pExperiment->setSeparator(delim); //tab-delimited
	pExperiment->setFirstRow(1);
	pExperiment->setLastRow(numlines-1);
	pExperiment->setHeaderRow(0);
	pExperiment->setExperimentType(CCopasiTask::timeCourse);
	pExperiment->setNumColumns(targetSpecies.size());
	CExperimentObjectMap* pObjectMap = &pExperiment->getObjectMap();

	//assign index for time
	pObjectMap->setNumCols(targetSpecies.size());
	pObjectMap->setRole(0, CExperiment::time);
	const CCopasiObject* pTimeReference = pModel->getObject(CCopasiObjectName("Reference=Time"));
	pObjectMap->setObjectCN(0, pTimeReference->getCN());

	// now we tell COPASI which column contain the concentrations of metabolites and belong to dependent variables
	int k;
	CMetab * pMetab;
	cout << "num = " << targetSpecies.size() << endl;
	for (int i=1; i < targetSpecies.size(); ++i)
	{
		k = i;
		pMetab = targetSpecies[i];
		if (pMetab)
		{
			pObjectMap->setRole( k , CExperiment::dependent );
			const CCopasiObject* pParticleReference = pMetab->getObject(CCopasiObjectName("Reference=Concentration"));
			pObjectMap->setObjectCN(k, pParticleReference->getCN());
		}
	}

	pExperimentSet->addExperiment(*pExperiment);
	// addExperiment makes a copy, so we need to get the added experiment
	delete pExperiment;
	pExperiment = pExperimentSet->getExperiment(0);
	// assign optimization parameters
	// get the list where we have to add the fit items
	CCopasiParameterGroup* pOptimizationItemGroup = dynamic_cast<CCopasiParameterGroup*>(pFitProblem->getParameter("OptimizationItemList"));

	// define CFitItem for each param
	for (int i=0; i < targetParams.size(); ++i)
		if (targetParams[i])
		{
			const CCopasiObject * pParameterReference = targetParams[i]->getObject(CCopasiObjectName("Reference=Value"));
			CFitItem* pFitItem = new CFitItem(pDataModel);
			pFitItem->setObjectCN(pParameterReference->getCN());
			pFitItem->setStartValue(c_getMatrixValue(params,i,0));

			stringstream ss1, ss2;
			ss1 << "0.0";//c_getMatrixValue(params,i,1);
			ss1 << "2.0"; //c_getMatrixValue(params,i,2);
			std::cout << "set lower bound " << pFitItem->setLowerBound(CCopasiObjectName("0.0"));
			std::cout << "set upper bound " << pFitItem->setUpperBound(CCopasiObjectName("1.0"));
			pOptimizationItemGroup->addParameter(pFitItem);
		}

	try
	{
		// initialize the fit task
		// we want complete output (HEADER, BODY and FOOTER)
		bool result = pFitTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		if (result == true)
		{
			// running the task for this example will probably take some time
			result = pFitTask->process(true);
		}
	}
	catch (...)
	{
		// failed
		return;
	}
	pFitTask->restore();

	//assign optimized values back into the model
	for (int i=0; i < params.rows && i < pFitProblem->getSolutionVariables().size(); ++i)
	{
		double x = pFitProblem->getSolutionVariables()[i];
		cSetValue(model, c_getRowName(params,i), x);
		c_setMatrixValue(params, i, 0, x);
		cout << c_getRowName(params,i) << " = " << x << endl;
	}
}

/*****************************************************************
   GENETIC ALGORITHM BASED OPTIMIZATION -- HELPER FUNCTIONS
******************************************************************/

typedef GA1DArrayGenome<float> RealGenome;

struct GAData
{
	mu::Parser * parser;
	double * parserValues;
	copasi_model * model;
	c_matrix * data;
	c_matrix * params;
};

static void InitializeGenome(GAGenome & x)
{
	RealGenome & g = (RealGenome &)x;
	GAData * data = (GAData*)(g.geneticAlgorithm()->userData());
	c_matrix * params = data->params;
	for (int i=0; i < g.size() && i < params->rows; ++i)
		g.gene(i,0) = c_getMatrixValue(*params,i,1) + mtrand() * (c_getMatrixValue(*params,i,2) - c_getMatrixValue(*params,i,1));
}

static float EuclideanDistance(const GAGenome & c1, const GAGenome & c2)
{
  const RealGenome & a = (RealGenome &)c1;
  const RealGenome & b = (RealGenome &)c2;

  float x=0.0;
  for(int i=0; i < b.length() && i < a.length(); ++i)
	  x += (a.gene(i) - b.gene(i))*(a.gene(i) - b.gene(i));

  return (float)(x);
}

static float MatrixDistance(c_matrix * data1, c_matrix * data2)
{
  int n = 0;
  float x=0.0, total=0.0;
  for(int i=1; i < data1->cols && i < data2->cols; ++i)
  	for(int j=0; j < data1->rows && j < data2->rows; ++j)
  	{
  	   x = c_getMatrixValue(*data1,i,j) - c_getMatrixValue(*data2,i,j);
	   total += (x*x);
	   ++n;
	}

  return total/(float)(n);
}

static float ObjectiveForFittingTimeSeries(GAGenome & x)
{
	RealGenome & g = (RealGenome &)x;

	if (!g.geneticAlgorithm())
		return numeric_limits<float>::max();

	GAData * pData = (GAData*)(g.geneticAlgorithm()->userData());
	copasi_model * model = pData->model;
	c_matrix * data = pData->data;
	c_matrix * params = pData->params;
	double retval;

	for (int i=0; i < params->rows && i < g.length(); ++i)
		cSetValue( *model, c_getRowName(*params,i), g.gene(i) );

	double start = c_getMatrixValue(*data, 0, 0),
			   end = c_getMatrixValue(*data, data->rows-1, 0);
	c_matrix output = cSimulateDeterministic(*model, start, end, data->rows);

	retval = MatrixDistance(data, &output);

	c_deleteMatrix(output);

	for (int i=0; i < params->rows; ++i)
	{
		if (g.gene(i) < c_getMatrixValue(*params, i, 1) ||
			 g.gene(i) > c_getMatrixValue(*params, i, 2))
		 {
		 	retval = numeric_limits<float>::max();
		 	break;
		 }
	}
	return retval;
}

static float ObjectiveForFittingSteadyStateData(GAGenome & x)
{
	RealGenome & g = (RealGenome &)x;

	GAData * pData = (GAData*)(g.geneticAlgorithm()->userData());
	copasi_model * model = pData->model;
	c_matrix * data = pData->data;
	c_matrix * params = pData->params;

	for (int i=0; i < params->rows && i < g.length(); ++i)
		cSetValue( *model, c_getRowName(*params,i), g.gene(i) );

	double retval;

	c_matrix output = efficiently_createMatrix(data->rows, data->cols);
	const char * name = c_getRowName(*data,0);

	for (int i=0; i < data->rows; ++i)
	{
		cSetValue(*model, name, c_getMatrixValue(*data, i, 0));
		c_setMatrixValue(output, i, 0, c_getMatrixValue(*data, i, 0));
		c_matrix ss = cGetSteadyState(*model);
		for (int j=0; j < ss.rows; ++j)
			for (int k=0; k < data->cols; ++k)
				if (string(c_getRowName(ss, j)) == string(c_getColumnName(*data,k)))
				{
					c_setMatrixValue(output, i, k+1, c_getMatrixValue(ss, j, 0));
					break;
				}
		c_deleteMatrix(ss);
	}

	retval = MatrixDistance(data, &output);
	c_deleteMatrix(output);

	for (int i=0; i < params->rows; ++i)
	{
		if (g.gene(i) < c_getMatrixValue(*params, i, 1) ||
			 g.gene(i) > c_getMatrixValue(*params, i, 2))
		 {
		 	retval = numeric_limits<float>::max();
		 	break;
		 }
	}

	return retval;
}

static float ObjectiveForMaximizingFormula(GAGenome & x)
{
	RealGenome & g = (RealGenome &)x;

	GAData * pData = (GAData*)(g.geneticAlgorithm()->userData());
	mu::Parser * parser = pData->parser;
	copasi_model * model = pData->model;
	c_matrix * data = pData->data;
	c_matrix * params = pData->params;
	double retval;

	for (int i=0; i < params->rows && i < g.length(); ++i)
		cSetValue( *model, c_getRowName(*params,i), g.gene(i) );

	if (parser)
	{
		c_matrix ss = cGetSteadyState(*model);
		for (int i=0; i < ss.rows; ++i)
			pData->parserValues[i] = c_getMatrixValue(ss, i, 0);

		retval = parser->Eval();
		c_deleteMatrix(ss);
	}
	else
	{
		retval = numeric_limits<float>::max();
	}

	for (int i=0; i < params->rows; ++i)
	{
		if (g.gene(i) < c_getMatrixValue(*params, i, 1) ||
			 g.gene(i) > c_getMatrixValue(*params, i, 2))
		 {
		 	retval = numeric_limits<float>::max();
		 	break;
		 }
	}

	return retval;
}

/**************************************************
   GENETIC ALGORITHM BASED OPTIMIZATION
***************************************************/
static int _OPTIMIZATION_MAX_ITER = 100;
static int _OPTIMIZATION_POPULATION_SIZE = 1000;
static double _OPTIMIZATION_MUTATION_RATE = 0.2;
static double _OPTIMIZATION_CROSSOVER_RATE = 0.8;

c_matrix cOptimize(copasi_model model, const char * objective, c_matrix params)
{
	ifstream file(objective);

	mu::Parser parser;
	GAData pData;
	c_matrix data;

	pData.parser = 0;
	pData.model = &model;
	pData.data = 0;
	pData.params = &params;

	list<string> words;
	int numlines=0;
	if (file.is_open())
	{
		string delim("\t");

		string line;
		getline(file, line);
		if (line[0] == '#')
			line.erase(0,1);

		if (contains(line, ","))
			delim = string(",");

		words = splitString(line , delim);

		if (!words.empty())
		{
			while (file.good())
			{
				getline(file,line);
				++numlines;
			}
		}

		file.close();

		if (words.empty() || numlines < 1)
		{
			return c_createMatrix(0,0);
		}
		else
		{
			ifstream file(objective);
			if (file.is_open())
			{
				data = c_createMatrix(numlines, words.size());
				pData.data = &data;
				getline(file, line);
				int i=0;
				double d;

				i = 0;
				for (list<string>::iterator it=words.begin(); it != words.end(); it++, ++i)
				{
					c_setColumnName(data, i, (*it).c_str());
				}

				i = 0;
				while (file.good() && i < numlines)
				{
					getline(file, line);
					words = splitString( line, delim);
					list<string>::iterator it = words.begin();
					for (int j=0; it != words.end() && j < data.cols; ++j, it++)
					{
						d = string_to_double((*it));
						c_setMatrixValue(data, i, j, d); //set data
					}
					++i;
				}
				file.close();
			}
		}
	}
	else
	{
		c_matrix ss = cGetSteadyState(model);
		double * array = new double[ss.rows];
		for (int i=0; i < ss.rows; ++i)
		{
			double * dp = &(array[i]);
			parser.DefineVar(c_getRowName(ss,i), dp);   //add all the model variables
		}
		pData.parserValues = array;
		parser.SetExpr(objective);

		try
		{
			parser.Eval();  //just checking if the eq. is valid
		}
		catch(...)
		{
			delete pData.parserValues;
			return c_createMatrix(0,0);
		}
	}

	GAPopulation pop;

	if (pData.parser)
	{
		RealGenome genome( params.rows , &ObjectiveForMaximizingFormula );
		genome.initializer(&InitializeGenome);
		GASteadyStateGA ga(genome);
		ga.userData(&pData);
		pop = ga.population();
		GASharing dist(EuclideanDistance);
		ga.scaling(dist);
		ga.pReplacement(1.0);
		ga.maximize();
		ga.populationSize(_OPTIMIZATION_POPULATION_SIZE);
		ga.nGenerations(_OPTIMIZATION_MAX_ITER);
		ga.pMutation(_OPTIMIZATION_MUTATION_RATE);
		ga.pCrossover(_OPTIMIZATION_CROSSOVER_RATE);
		ga.initialize();
		int k=0;
		while (ga.done() != gaTrue)
		{
			ga.step();
			//cout << "gen " << ++k << "\n";
		}
		//ga.evolve();
		pop = ga.population();
		pop.order(GAPopulation::HIGH_IS_BEST);
		pop.sort(gaTrue);
		delete pData.parserValues;
	}
	else
	if (pData.data)
	{
		if (string(c_getColumnName(data,0)) == string("time") ||
  			 string(c_getColumnName(data,0)) == string("Time") ||
			 string(c_getColumnName(data,0)) == string("TIME"))
		{
			RealGenome genome( params.rows , &ObjectiveForFittingTimeSeries );
			genome.initializer(&InitializeGenome);
			GASteadyStateGA ga(genome);
			ga.userData(&pData);
			pop = ga.population();
			GASharing dist(EuclideanDistance);
			ga.scaling(dist);
			ga.pReplacement(1.0);
			ga.maximize();
			ga.populationSize(_OPTIMIZATION_POPULATION_SIZE);
			ga.nGenerations(_OPTIMIZATION_MAX_ITER);
			ga.pMutation(_OPTIMIZATION_MUTATION_RATE);
			ga.pCrossover(_OPTIMIZATION_CROSSOVER_RATE);
			ga.minimize();
			ga.initialize();
			int k=0;
			while (ga.done() != gaTrue)
			{
				ga.step();
				//cout << "\n\ngen " << ++k << "\n";
			}
			pop = ga.population();
			pop.order(GAPopulation::LOW_IS_BEST);
			pop.sort(gaTrue);
		}
		else
		{
			RealGenome genome( params.rows , &ObjectiveForFittingSteadyStateData );
			genome.initializer(&InitializeGenome);
			GASteadyStateGA ga(genome);
			ga.userData(&pData);
			pop = ga.population();
			GASharing dist(EuclideanDistance);
			ga.scaling(dist);
			ga.pReplacement(1.0);
			ga.maximize();
			ga.populationSize(_OPTIMIZATION_POPULATION_SIZE);
			ga.nGenerations(_OPTIMIZATION_MAX_ITER);
			ga.pMutation(_OPTIMIZATION_MUTATION_RATE);
			ga.pCrossover(_OPTIMIZATION_CROSSOVER_RATE);
			ga.minimize();
			ga.initialize();
			int k=0;
			while (ga.done() != gaTrue)
			{
				ga.step();
				//cout << "gen " << ++k << "\n";
			}
			pop = ga.population();
			pop.order(GAPopulation::LOW_IS_BEST);
			pop.sort(gaTrue);
		}
		c_deleteMatrix(data);
	}
	else
	{
		return c_createMatrix(0,0);
	}

	c_matrix result = efficiently_createMatrix( pop.size(), params.rows );

	for (int i=0; i < pop.size(); ++i)
	{
		RealGenome & g = (RealGenome &)(pop.individual(i));
		for (int j=0; j < params.rows; ++j)
			c_setMatrixValue(result, i, j, g.gene(j));
	}

	return result;
}

void cSetOptimizerSize(int n)
{
	_OPTIMIZATION_POPULATION_SIZE = n;
}

void cSetOptimizerIterations(int n)
{
	_OPTIMIZATION_MAX_ITER = n;
}

void cSetOptimizerMutationRate(double q)
{
	_OPTIMIZATION_MUTATION_RATE = q;
}

void cSetOptimizerCrossoverRate(double q)
{
	_OPTIMIZATION_CROSSOVER_RATE = q;
}

/* LIBSTRUCTURAL */

c_matrix convertFromDoubleMatrix(DoubleMatrix& matrix, vector< string > &rowNames, vector< string > &colNames)
{
	c_matrix m = efficiently_createMatrix(matrix.numRows(), matrix.numCols());

	for (int i=0; i < m.rows && i < rowNames.size(); ++i)
		c_setRowName(m, i, rowNames[i].c_str());

	for (int i=0; i < m.cols && i < colNames.size(); ++i)
		c_setColumnName(m, i, colNames[i].c_str());

	for (int i=0; i < m.rows; ++i)
		for (int j=0; j < m.cols; ++j)
			c_setMatrixValue(m, i, j, matrix(i,j));

	return m;
}

void convertToDoubleMatrix(c_matrix m, DoubleMatrix & matrix, vector< string > &rowNames, vector< string > &colNames)
{
	matrix.resize(m.rows, m.cols);

	rowNames.resize(m.rows);
	colNames.resize(m.cols);

	for (int i=0; i < m.rows && i < rowNames.size(); ++i)
		rowNames[i] = string(c_getRowName(m, i));

	for (int i=0; i < m.cols && i < colNames.size(); ++i)
		colNames[i] = string(c_getColumnName(m, i));

	for (int i=0; i < m.rows; ++i)
		for (int j=0; j < m.cols; ++j)
			matrix(i,j) = c_getMatrixValue(m, i, j);
}

c_matrix cGetGammaMatrix(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return c_createMatrix(0,0);

	//get stoichiometry
	c_matrix c_N = cGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;
	DoubleMatrix N;

	convertToDoubleMatrix( c_N , N, rowNames, colNames );

	//use libstructural
	LibStructural * instance = LibStructural::getInstance();
	instance->loadStoichiometryMatrix (N);

	CCopasiVector< CMetab > & metabolites = pModel->getMetabolites();
	vector<double> iv(metabolites.size(), 0);  //species concentrations
	for (int i=0; i < metabolites.size(); ++i)
		if (metabolites[i] != NULL)
			iv[i] = metabolites[i]->getInitialConcentration();

	instance->loadSpecies (rowNames, iv);
	instance->loadReactionNames (colNames);
	instance->analyzeWithQR();

	DoubleMatrix * matrix = instance->getGammaMatrix();

	rowNames.clear();
	colNames.clear();
	instance->getGammaMatrixLabels(rowNames, colNames);

	//convert
	c_matrix m = convertFromDoubleMatrix(*matrix, rowNames, colNames);

	//cleanup
	//delete matrix;
	c_deleteMatrix(c_N);

	return m;
}

c_matrix cGetKMatrix(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return c_createMatrix(0,0);

	//get stoichiometry
	c_matrix c_N = cGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;
	DoubleMatrix N;

	convertToDoubleMatrix( c_N , N, rowNames, colNames );

	//use libstructural
	LibStructural * instance = LibStructural::getInstance();
	instance->loadStoichiometryMatrix (N);

	CCopasiVector< CMetab > & metabolites = pModel->getMetabolites();
	vector<double> iv(metabolites.size(), 0);  //species concentrations
	for (int i=0; i < metabolites.size(); ++i)
		if (metabolites[i] != NULL)
			iv[i] = metabolites[i]->getInitialConcentration();

	instance->loadSpecies (rowNames, iv);
	instance->loadReactionNames (colNames);
	instance->analyzeWithQR();

	DoubleMatrix * matrix = instance->getKMatrix();

	rowNames.clear();
	colNames.clear();
	instance->getKMatrixLabels(rowNames, colNames);

	//convert
	c_matrix m = convertFromDoubleMatrix(*matrix, rowNames, colNames);

	//cleanup
	//delete matrix;
	c_deleteMatrix(c_N);

	return m;
}

c_matrix cGetLinkMatrix(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return c_createMatrix(0,0);

	//get stoichiometry
	c_matrix c_N = cGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;
	DoubleMatrix N;

	convertToDoubleMatrix( c_N , N, rowNames, colNames );

	//use libstructural
	LibStructural * instance = LibStructural::getInstance();
	instance->loadStoichiometryMatrix (N);

	CCopasiVector< CMetab > & metabolites = pModel->getMetabolites();
	vector<double> iv(metabolites.size(), 0);  //species concentrations
	for (int i=0; i < metabolites.size(); ++i)
		if (metabolites[i] != NULL)
			iv[i] = metabolites[i]->getInitialConcentration();

	instance->loadSpecies (rowNames, iv);
	instance->loadReactionNames (colNames);
	instance->analyzeWithQR();

	DoubleMatrix * matrix = instance->getLinkMatrix();

	rowNames.clear();
	colNames.clear();
	instance->getLinkMatrixLabels(rowNames, colNames);

	//convert
	c_matrix m = convertFromDoubleMatrix(*matrix, rowNames, colNames);

	//cleanup
	//delete matrix;
	c_deleteMatrix(c_N);

	return m;
}

c_matrix cGetK0Matrix(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return c_createMatrix(0,0);

	//get stoichiometry
	c_matrix c_N = cGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;
	DoubleMatrix N;

	convertToDoubleMatrix( c_N , N, rowNames, colNames );

	//use libstructural
	LibStructural * instance = LibStructural::getInstance();
	instance->loadStoichiometryMatrix (N);

	CCopasiVector< CMetab > & metabolites = pModel->getMetabolites();
	vector<double> iv(metabolites.size(), 0);  //species concentrations
	for (int i=0; i < metabolites.size(); ++i)
		if (metabolites[i] != NULL)
			iv[i] = metabolites[i]->getInitialConcentration();

	instance->loadSpecies (rowNames, iv);
	instance->loadReactionNames (colNames);
	instance->analyzeWithQR();

	DoubleMatrix * matrix = instance->getK0Matrix();

	rowNames.clear();
	colNames.clear();
	instance->getK0MatrixLabels(rowNames, colNames);

	//convert
	c_matrix m = convertFromDoubleMatrix(*matrix, rowNames, colNames);

	//cleanup
	//delete matrix;
	c_deleteMatrix(c_N);

	return m;
}

c_matrix cGetL0Matrix(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return c_createMatrix(0,0);

	//get stoichiometry
	c_matrix c_N = cGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;
	DoubleMatrix N;

	convertToDoubleMatrix( c_N , N, rowNames, colNames );

	//use libstructural
	LibStructural * instance = LibStructural::getInstance();
	instance->loadStoichiometryMatrix (N);

	CCopasiVector< CMetab > & metabolites = pModel->getMetabolites();
	vector<double> iv(metabolites.size(), 0);  //species concentrations
	for (int i=0; i < metabolites.size(); ++i)
		if (metabolites[i] != NULL)
			iv[i] = metabolites[i]->getInitialConcentration();

	instance->loadSpecies (rowNames, iv);
	instance->loadReactionNames (colNames);
	instance->analyzeWithQR();

	DoubleMatrix * matrix = instance->getL0Matrix();

	rowNames.clear();
	colNames.clear();
	instance->getL0MatrixLabels(rowNames, colNames);

	//convert
	c_matrix m = convertFromDoubleMatrix(*matrix, rowNames, colNames);

	//cleanup
	//delete matrix;
	c_deleteMatrix(c_N);

	return m;
}

list<string> splitString(const string& seq, const string& _1cdelim)
{
	bool keeptoken = false, _removews = true;
	list<string> L;

	string delims = _1cdelim;
	string STR;
	if(delims.empty()) delims = "\n\r";
	if(_removews) delims += " ";

	string::size_type pos=0, LEN = seq.size();
	while(pos < LEN ){
	STR=""; // Init/clear the STR token buffer
	// remove any delimiters including optional (white)spaces
	while( (delims.find(seq[pos]) != string::npos) && (pos < LEN) ) ++pos;
	// leave if @eos
	if(pos==LEN) return L;
	// Save token data
	while( (delims.find(seq[pos]) == string::npos) && (pos < LEN) ) STR += seq[pos++];
	// put valid STR buffer into the supplied list
	//cout << "[" << STR << "]";
	if( ! STR.empty() ) L.push_back(STR);
	}
	return L;
}


int replaceSubstring(std::string& s,const std::string& from, const std::string& to)
{
	int cnt = -1;

	if(from != to && !from.empty())
	{
		string::size_type pos1(0);
		string::size_type pos2(0);
		const string::size_type from_len(from.size());
		const string::size_type to_len(to.size());
		cnt = 0;

		while((pos1 = s.find(from, pos2)) != std::string::npos)
		{
			s.replace(pos1, from_len, to);
			pos2 = pos1 + to_len;
			++cnt;
		}
	}

	return cnt;
}

c_matrix cGetReactionRates(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);

	if (!pModel) return c_createMatrix(0,0);

	const CCopasiVectorNS < CReaction > & reactions = pModel->getReactions();

	c_matrix res  = efficiently_createMatrix(1, reactions.size());

	for (int i=0; i < reactions.size(); ++i)
		if (reactions[i])
		{
			c_setColumnName(res, i, reactions[i]->getObjectName().c_str());
			c_setMatrixValue(res, 0, i, reactions[i]->calculateFlux());
		}

	return res;
}

double cGetReactionRate(copasi_model model, const char * name)
{
	return cGetFlux(model, name);
}

c_matrix cGetReactionRatesEx(copasi_model model, c_matrix conc)
{
	int i,j;
	c_matrix savedMatrix = c_createMatrix(conc.rows, conc.cols);
	for (i=0; i < conc.rows; ++i)
		for (j=0; j < conc.cols; ++j)
			c_setMatrixValue( savedMatrix, i, j, c_getMatrixValue(conc, i, j));
	cSetValues(model,conc);

	double * temp = conc.values;
	conc.values = savedMatrix.values;

	cSetValues(model,conc);

	conc.values = temp;
	c_deleteMatrix(savedMatrix);

	c_matrix rates = cGetReactionRates(model);

	cResetState(model);

	return rates;
}

c_matrix cGetFloatingSpeciesConcentrations(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return c_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	int n = 0;
	for (int i=0; i < species.size(); ++i)
		if (species[i] &&
			(species[i]->getStatus() == CModelEntity::ODE || species[i]->getStatus() == CModelEntity::REACTIONS))
			++n;

	c_matrix res  = efficiently_createMatrix(n,1);

	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] &&
			(species[i]->getStatus() == CModelEntity::ODE || species[i]->getStatus() == CModelEntity::REACTIONS))
		{
			c_setMatrixValue(res, j, 0, species[i]->getConcentration());
			c_setRowName(res, j, species[i]->getObjectName().c_str());
			++j;
		}

	return res;
}

c_matrix cGetBoundarySpecies(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return c_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	int n = 0;
	for (int i=0; i < species.size(); ++i)
		if (species[i] &&
			(species[i]->getStatus() == CModelEntity::FIXED))
			++n;

	c_matrix res  = efficiently_createMatrix(n,1);

	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] &&
			(species[i]->getStatus() == CModelEntity::FIXED))
		{
			c_setMatrixValue(res, j, 0, species[i]->getConcentration());
			c_setRowName(res, j, species[i]->getObjectName().c_str());
			++j;
		}

	return res;
}

int cGetNumberOfSpecies(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return 0;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	int n = species.size();
	return n;
}

int cGetNumberOfFloatingSpecies(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return 0;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	int n = 0;
	for (int i=0; i < species.size(); ++i)
		if (species[i] &&
			(species[i]->getStatus() == CModelEntity::ODE || species[i]->getStatus() == CModelEntity::REACTIONS))
			++n;
	return n;
}

int cGetNumberOfBoundarySpecies(copasi_model model)
{
		CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return 0;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	int n = 0;
	for (int i=0; i < species.size(); ++i)
		if (species[i] &&
			(species[i]->getStatus() == CModelEntity::FIXED))
			++n;
	return n;
}
c_matrix cGetFloatingSpeciesIntitialConcentrations (copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return c_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	int n = 0;
	for (int i=0; i < species.size(); ++i)
		if (species[i] &&
			(species[i]->getStatus() == CModelEntity::ODE || species[i]->getStatus() == CModelEntity::REACTIONS))
			++n;

	c_matrix res  = efficiently_createMatrix(n,1);

	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] &&
			(species[i]->getStatus() == CModelEntity::ODE || species[i]->getStatus() == CModelEntity::REACTIONS))
		{
			c_setMatrixValue(res, j, 0, species[i]->getInitialConcentration());
			c_setRowName(res, j, species[i]->getObjectName().c_str());
			++j;
		}

	return res;
}

void cSetFloatingSpeciesIntitialConcentrations (copasi_model model, c_matrix sp)
{
	if (sp.rows > sp.cols)  //row vector or column vector (lets allow both)
	{
		for (int i=0; i < sp.rows; ++i)
			cSetInitialConcentration(model, c_getRowName(sp,i), c_getMatrixValue(sp, i, 0));
	}
	else
	{
		for (int i=0; i < sp.cols; ++i)
			cSetInitialConcentration(model, c_getColumnName(sp,i), c_getMatrixValue(sp, 0, i));
	}
}

void cSetBoundarySpeciesConcentrations (copasi_model model, c_matrix d)
{
	cSetValues(model, d);
}

c_matrix cGetAllSpecies(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return c_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	c_matrix res  = efficiently_createMatrix(1, species.size());

	for (int i=0; i < species.size(); ++i)
		if (species[i])
		{
			c_setColumnName(res, i, species[i]->getObjectName().c_str());
			c_setMatrixValue(res, 0, i, species[i]->getConcentration());
		}

	return res;
}

c_matrix cGetCompartments(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return c_createMatrix(0,0);

	const CCopasiVectorNS< CCompartment > & compartments = pModel->getCompartments();

	c_matrix res  = efficiently_createMatrix(1, compartments.size());

	for (int i=0; i < compartments.size(); ++i)
		if (compartments[i])
		{
			c_setColumnName(res, i, compartments[i]->getObjectName().c_str());
			c_setMatrixValue(res, 0, i, compartments[i]->getValue());
		}
	return res;
}

c_matrix cGetAmounts(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel) return c_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	c_matrix res  = efficiently_createMatrix(1, species.size());

	for (int i=0; i < species.size(); ++i)
		if (species[i] && species[i]->getCompartment())
		{
			c_setColumnName(res, i, species[i]->getObjectName().c_str());
			c_setMatrixValue(res, 0, i,
					CMetab::convertToNumber( species[i]->getConcentration(), *species[i]->getCompartment(), pModel ));
		}

	return res;
}

double cGetConcentration(copasi_model model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);

	if (!pModel || !contains(hash, s)) return -1.0;
	CopasiPtr & p = getHashValue(hash, s);

	if (!p.species || !p.species->getCompartment()) return -1.0;

	return p.species->getConcentration();
}

double cGetAmount(copasi_model model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);

	if (!pModel || !contains(hash, s)) return -1.0;
	CopasiPtr & p = getHashValue(hash, s);

	if (!p.species || !p.species->getCompartment()) return -1.0;

	return CMetab::convertToNumber( p.species->getConcentration(), *p.species->getCompartment(), pModel );
}

c_matrix cGetRatesOfChangeEx(copasi_model model, c_matrix m)
{
	cSetValues(model, m);
	c_matrix m2 = cGetRatesOfChange(model);
	cResetState(model);
	return m2;
}

c_matrix cGetRatesOfChange(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);

	if (!pModel)
		return c_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	int n = pModel->getState().getNumVariable();

	if (n > species.size()) n = species.size();
	c_matrix res  = efficiently_createMatrix(1, n);

	for (int i=0; i < species.size(); ++i)
		if (species[i] && species[i]->getCompartment())
		{
			c_setColumnName(res, i, species[i]->getObjectName().c_str());
			species[i]->refreshRate();
			c_setMatrixValue(res, 0, i, species[i]->getRate());
		}
	return res;
}

double cGetRateOfChange(copasi_model model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);

	if (!pModel || !contains(hash, s)) return 0.0;
	CopasiPtr & p = getHashValue(hash, s);

	if (!p.species || !p.species->getCompartment()) return 0.0;

	p.species->refreshRate();
	return p.species->getRate();
}

double cGetFlux(copasi_model model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s (name);

	if (!pModel || !contains(hash, s)) return NaN;

	CopasiPtr & p = getHashValue(hash, s);

	if (!p.reaction) return NaN;

	return p.reaction->calculateFlux();
}

double cGetParticleFlux(copasi_model model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s (name);

	if (!pModel || !contains(hash, s)) return NaN;

	CopasiPtr & p = getHashValue(hash, s);

	if (!p.reaction) return NaN;

	return p.reaction->calculateParticleFlux();
}


// ------------------------------------------------------------------
// Parameter Group
// ------------------------------------------------------------------

int sSetGlobalParameter(copasi_model model, const char * name, double value)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CModelValue * pValue = NULL;

	if (!hash || !pModel) return 0;

	if (contains(hash,s) &&
		(pValue = getHashValue(hash,s).param))
	{
		pValue->setInitialValue(value);
		pValue->setValue(value);
		return 1;
	}

	//parameter not found, so create it
	if (!pValue)
	{
		pValue = pModel->createModelValue(string(name),value);
		pValue->setInitialValue(value);

		CopasiPtr copasiPtr = {
				pValue->getCN(),
				pValue->getKey(),
				0,
				0,
				0,
				pValue,
				"",
				true};

		hashInsert(hash, s, copasiPtr); //for speedy lookup
	}

	return 0;
}

int cGetNumberOfGlobalParameters (copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	if (!pModel) return 0;

	CCopasiVectorN< CModelValue > & params = pModel->getModelValues();
	return params.size();
}

c_matrix cGetGlobalParameters (copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel || !hash) return c_createMatrix(0,0);

	CCopasiVectorN< CModelValue > & params = pModel->getModelValues();

	c_matrix m = efficiently_createMatrix(params.size(),1);

	for (int i=0; i < params.size(); ++i)
		if (params[i])
		{
			c_setRowName(m, i, params[i]->getObjectName().c_str());
			c_setMatrixValue(m, i, 0, params[i]->getInitialValue());
		}
	return m;
}

void cSetValues (copasi_model model, c_matrix gp)
{
	if (gp.rows > gp.cols)  //row vector or column vector (lets allow both)
	{
		for (int i=0; i < gp.rows; ++i)
			cSetValue(model, c_getRowName(gp,i), c_getMatrixValue(gp, i, 0));
	}
	else
	{
		for (int i=0; i < gp.cols; ++i)
			cSetValue(model, c_getColumnName(gp,i), c_getMatrixValue(gp, 0, i));
	}
}

void cSetGlobalParameterValues (copasi_model model, c_matrix gp)
{
	cSetValues(model, gp);
}

void cSetCompartmentVolumes (copasi_model model, c_matrix v)
{
	cSetValues(model, v);
}

void cSetSpeciesConcentrations (copasi_model model, c_matrix v)
{
	cSetValues(model, v);
}

void cSetSpeciesAmounts (copasi_model model, c_matrix v)
{
	cSetValues(model, v);
}

int cGetNumberOfCompartments (copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);

	if (!pModel) return 0;

	const CCopasiVectorNS< CCompartment > & compartments = pModel->getCompartments();
	return compartments.size();
}

int cGetNumberOfReactions (copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);

	if (!pModel) return 0;

	CCopasiVectorNS < CReaction > & reacs = pModel->getReactions();
	return reacs.size();
}

static int GET_CC = 0;
static int GET_ELAS = 1;
static int GET_RATES = 2;
static int GET_DERIV = 3;

//this function is a generic function used to get different MCA related calculations from simulated data
c_matrix cGetXFromTimeCourse(copasi_model model, c_matrix results, int type)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel || !hash) return results;

	vector<CMetab*> metabs;
	metabs.resize(results.cols-1, (CMetab*)NULL);
	for (int i=1; i < results.cols; ++i)
	{
		string s(c_getColumnName(results, i));
		if (contains(hash, s))
		{
			CopasiPtr & p = getHashValue(hash, s);
			metabs[i-1] = p.species;
		}
	}

	c_matrix m;
	if (type == GET_ELAS)
		m = cGetScaledElasticities(model);
	else
	if (type == GET_CC)
		m = cGetScaledConcentrationConcentrationCoeffs(model);
	else
	if (type == GET_RATES)
		m = cGetReactionRates(model);
	else
		m = cGetRatesOfChange(model);

	setStorageMatrix(&m); //for speed-up.. might be risky on multithreaded apps

	int i,j,k;
	//get column names
	list<string> colnames;
	for (i=0; i < m.rows; ++i)
	{
		char c0[] = {0};
		char * c = (char*)c_getRowName(m, i);
		if (c == 0)
			c = c0;
		string rowname( c );

		for (j=0; j < m.cols; ++j)
		{
			c = (char*)c_getColumnName(m, j);
			if (c == 0)
				c = c0;
			string colname( c );

			if (type == GET_ELAS)
				colnames.push_back( string("ec_") + rowname + string("_") + colname); //naming convention: es_x_y
			else
			if (type == GET_CC)
				colnames.push_back( string("cc_") + rowname + string("_") + colname); //naming convention: cc_x_y
			else
			if (m.rows > 1)
				colnames.push_back( rowname + string("_") + colname ); //
			else
			if (type == GET_DERIV)
				colnames.push_back( colname + string("'") ); //same
			else
				colnames.push_back( colname); //same
		}
	}

	c_matrix output = c_createMatrix(results.rows, 1 + colnames.size()); //the matrix with calculations for each row in the time series

	j = 0;
	c_setColumnName( output, 0, c_getColumnName(results,0) );  //first column name remains same, e.g. "time"

	for (list<string>::iterator it=colnames.begin(); j < (1+output.cols) && it != colnames.end(); ++j, it++) //set column names
		c_setColumnName( output, 1+j, (*it).c_str() );

	//the main loop -- make calculations for each row
	for (i=0; i < results.rows; ++i)
	{
		c_setMatrixValue( output, i, 0,
				c_getMatrixValue(results,i,0) );  //copy the first column, presumably the independent variable (e.g. time)

		for (j=1; j < results.cols; ++j) //update metab concentrations
			if (metabs[j-1])
			{
				double v = c_getMatrixValue(results, i, j);
				metabs[j-1]->setConcentration(v);
				metabs[j-1]->setValue(v);
			}

		if (type == GET_ELAS)
			m = cGetScaledElasticities(model);
		else
		if (type == GET_CC)
			m = cGetScaledConcentrationConcentrationCoeffs(model);
		else
		if (type == GET_RATES)
			m = cGetReactionRates(model);
		else
			m = cGetRatesOfChange(model);

		//c_printOutMatrix(m);

		int l = 0;
		for (j=0; j <  m.rows; ++j)
			for (k=0; k <  m.cols; ++k)
			{
				c_setMatrixValue( output, i, 1+l, c_getMatrixValue(m, j, k) );
				++l;
			}
	}

	cResetState(model);
	unsetStorageMatrix();
	return output;
}

c_matrix cGetCCFromTimeCourse(copasi_model model, c_matrix results)
{
	return cGetXFromTimeCourse(model, results, GET_CC);
}

c_matrix cGetElasticitiesFromTimeCourse(copasi_model model, c_matrix results)
{
	return cGetXFromTimeCourse(model, results, GET_ELAS);
}

c_matrix cGetReactionRatesFromTimeCourse(copasi_model model, c_matrix results)
{
	return cGetXFromTimeCourse(model, results, GET_RATES);
}

c_matrix cGetDerivativesFromTimeCourse(copasi_model model, c_matrix results)
{
	return cGetXFromTimeCourse(model, results, GET_DERIV);
}

c_matrix cFilterTimeCourseResults(copasi_model model, c_matrix results, c_strings names)
{
	c_matrix rates = cGetReactionRatesFromTimeCourse(model, results);
	c_matrix derivs = cGetDerivativesFromTimeCourse(model, results);
	c_matrix ccs = cGetCCFromTimeCourse(model, results);
	c_matrix elas = cGetCCFromTimeCourse(model, results);

	int numcols = 0;
	//these vectors will be used to identify which columns to include in the final output
	vector<bool> resultCols, ratesCols, derivCols, ccCols, elasCols;
	vector<string> snames;
	snames.resize(names.length);
	for (int i=0; i < names.length; ++i)
		snames[i] = string( c_getString(names, i) );  //convert C strings to std::string for comparison purpose

	if (names.length == 0) //include everything
	{
		resultCols.resize(results.cols,true);
		ratesCols.resize(rates.cols,true);
		derivCols.resize(derivs.cols,true);
		ccCols.resize(ccs.cols,true);
		elasCols.resize(elas.cols,true);
		numcols = results.cols + rates.cols + derivs.cols + ccs.cols + elas.cols;
	}
	else
	{
		numcols = 1; //independent variable (e.g time)
		ratesCols.resize(rates.cols,false);
		for (int i=0; i < rates.cols; ++i)
			for (int j=0; j < names.length; ++j) //go through all the desired names
				if (snames[j].compare( c_getColumnName(rates, i) ) == 0) //names match?
				{
					ratesCols[i] = true;
					++numcols;
					break;
				}

		derivCols.resize(results.cols,false);
		for (int i=0; i < derivs.cols; ++i)
			for (int j=0; j < names.length; ++j) //go through all the desired names
				if (snames[j].compare( c_getColumnName(derivs, i) ) == 0) //names match?
				{
					derivCols[i] = true;
					++numcols;
					break;
				}

		ccCols.resize(results.cols,false);
		for (int i=0; i < ccs.cols; ++i)
			for (int j=0; j < names.length; ++j) //go through all the desired names
				if (snames[j].compare( c_getColumnName(ccs, i) ) == 0) //names match?
				{
					ccCols[i] = true;
					++numcols;
					break;
				}

		elasCols.resize(results.cols,false);
		for (int i=0; i < elas.cols; ++i)
			for (int j=0; j < names.length; ++j) //go through all the desired names
				if (snames[j].compare( c_getColumnName(elas, i) ) == 0) //names match?
				{
					elasCols[i] = true;
					++numcols;
					break;
				}

		resultCols.resize(results.cols,false);
		if (results.cols > 0)
			resultCols[0] = true;

		for (int i=1; i < results.cols; ++i)
			for (int j=0; j < names.length; ++j) //go through all the desired names
				if (snames[j].compare( c_getColumnName(results, i) ) == 0) //names match?
				{
					resultCols[i] = true;
					++numcols;
					break;
				}
	}

	c_matrix output = c_createMatrix(results.rows, numcols);
	numcols = 0;
	for (int j=0; j < results.cols; ++j)
		if (resultCols[j])
		{
			c_setColumnName(output, numcols, c_getColumnName(results, j));
			for  (int i=0; i < output.rows; ++i)
				c_setMatrixValue(output, i, numcols, c_getMatrixValue(results,i,j));
			++numcols;
		}

	for (int j=0; j < rates.cols; ++j)
		if (ratesCols[j])
		{
			c_setColumnName(output, numcols, c_getColumnName(rates, j));
			for  (int i=0; i < output.rows; ++i)
				c_setMatrixValue(output, i, numcols, c_getMatrixValue(rates,i,j));
			++numcols;
		}

	for (int j=0; j < derivs.cols; ++j)
		if (derivCols[j])
		{
			c_setColumnName(output, numcols, c_getColumnName(derivs, j));
			for  (int i=0; i < output.rows; ++i)
				c_setMatrixValue(output, i, numcols, c_getMatrixValue(derivs,i,j));
			++numcols;
		}

	for (int j=0; j < ccs.cols; ++j)
		if (ccCols[j])
		{
			c_setColumnName(output, numcols, c_getColumnName(ccs, j));
			for  (int i=0; i < output.rows; ++i)
				c_setMatrixValue(output, i, numcols, c_getMatrixValue(ccs,i,j));
			++numcols;
		}

	for (int j=0; j < elas.cols; ++j)
		if (elasCols[j])
		{
			c_setColumnName(output, numcols, c_getColumnName(elas, j));
			for  (int i=0; i < output.rows; ++i)
				c_setMatrixValue(output, i, numcols, c_getMatrixValue(elas,i,j));
			++numcols;
		}

	return output;
}
/*
c_matrix cSimulationParameterScan(copasi_model model, const char * param, double start, double end, int numSteps, double startTime, double endTime, int nstep)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	CReportDefinitionVector* pReports = pDataModel->getReportDefinitionList();
	// create a new report definition object
	CReportDefinition* pReport = pReports->createReportDefinition("Report", "Output for timecourse");
	// set the task type for the report definition to timecourse
	pReport->setTaskType(CCopasiTask::timeCourse);
	// we don’t want a table
	pReport->setIsTable(false);
	//pReport->setIsTable(true);
	// the entries in the output should be seperated by a ", "
	//pReport->setSeparator(CCopasiReportSeparator("\t"));
	// we need a handle to the header and the body
	// the header will display the ids of the metabolites and "time" for
	// the first column
	// the body will contain the actual timecourse data
	std::vector<CRegisteredObjectName>* pHeader = pReport->getHeaderAddr();
	std::vector<CRegisteredObjectName>* pBody = pReport->getBodyAddr();
	pBody->push_back(CCopasiObjectName(pDataModel->getModel()->getCN() + ",Reference=Time"));
	pBody->push_back(CRegisteredObjectName(pReport->getSeparator().getCN()));
	pHeader->push_back(CCopasiStaticString("time").getCN());
	pHeader->push_back(pReport->getSeparator().getCN());
	unsigned int i, iMax = pModel->getMetabolites().size();
	for (i = 0;i < iMax;++i)
	{
		CMetab* pMetab = pModel->getMetabolites()[i];
		// we don’t want output for FIXED metabolites right now
		if (pMetab->getStatus() != CModelEntity::FIXED)
		{
			// we want the concentration oin the output
			// alternatively, we could use "Reference=Amount" to get the
			// particle number
			pBody->push_back(pMetab->getObject(CCopasiObjectName("Reference=Concentration"))->getCN());
			// after each entry, we need a seperator
			pBody->push_back(pReport->getSeparator().getCN());
			// add the corresponding id to the header
			pHeader->push_back(CCopasiStaticString(pMetab->getSBMLId()).getCN());
			// and a seperator
			pHeader->push_back(pReport->getSeparator().getCN());
		}
	}
	if (iMax > 0)
	{
		// delete the last separator
		// since we don’t need one after the last element on each line
		if ((*pBody->rbegin()) == pReport->getSeparator().getCN())
		{
			pBody->erase(--pBody->end());
		}
		if ((*pHeader->rbegin()) == pReport->getSeparator().getCN())
		{
			pHeader->erase(--pHeader->end());
		}
	}
	// get the task list
	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the trajectory task object
	CTrajectoryTask* pTrajectoryTask = dynamic_cast<CTrajectoryTask*>(TaskList["Time-Course"]);
	// if there isn’t one
	if (pTrajectoryTask == NULL)
	{
		// create a new one
		pTrajectoryTask = new CTrajectoryTask();
		// remove any existing trajectory task just to be sure since in
		// theory only the cast might have failed above
		TaskList.remove("Time-Course");
		// add the new time course task to the task list
		TaskList.add(pTrajectoryTask, true);
	}
	// run a stochastic time course
	pTrajectoryTask->setMethodType(CCopasiMethod::stochastic);
	// pass a pointer of the model to the problem
	pTrajectoryTask->getProblem()->setModel(pDataModel->getModel());
	// we don’t want the trajectory task to run by itself, but we want to
	// run it from a scan, so we deactivate the standalone trajectory task
	pTrajectoryTask->setScheduled(false);
	// get the problem for the task to set some parameters
	CTrajectoryProblem* pProblem = dynamic_cast<CTrajectoryProblem*>(pTrajectoryTask->getProblem());
	// simulate 100 steps
	pProblem->setStepNumber(nstep);
	// start at time 0
	pDataModel->getModel()->setInitialTime(0.0);
	// simulate a duration of 10 time units
	pProblem->setDuration(10);
	// tell the problem to actually generate time series data
	pProblem->setTimeSeriesRequested(true);
	// now we set up the scan
	CScanTask* pScanTask = dynamic_cast<CScanTask*>(TaskList["Scan"]);
	if (pScanTask == NULL)
	{
		// create a new scan task
		pScanTask = new CScanTask();
		// just to be on the save side, delete any existing scan task
		TaskList.remove("Scan");
		// add the new scan task
		TaskList.add(pScanTask, true);
	}
	// get the problem
	CScanProblem* pScanProblem = dynamic_cast<CScanProblem*>(pScanTask->getProblem());
	assert(pScanProblem != NULL);
	// set the model for the problem
	pScanProblem->setModel(pDataModel->getModel());
	// actiavate the task so that is is run
	// if the model is saved and passed to CopasiSE
	pScanTask->setScheduled(true);
	// set the report for the task
	pScanTask->getReport().setReportDefinition(pReport);
	// set the output file for the report
	pScanTask->getReport().setTarget("example4.txt");
	// don’t append to an existing file, but overwrite
	pScanTask->getReport().setAppend(false);
	// tell the scan that we want to make a scan over a trajectory task
	pScanProblem->setSubtask(CCopasiTask::timeCourse);
	// we just want to run the timecourse task a number of times, so we
	// create a repeat item with 100 repeats
	pScanProblem->createScanItem(CScanProblem::SCAN_REPEAT, 100);
	// we want the output from the trajectory task
	pScanProblem->setOutputInSubtask(true);
	// we don’t want to set the initial conditions of the model to the end
	// state of the last run
	pScanProblem->setAdjustInitialConditions(false);
	try
	{
		// initialize the trajectory task
		// we want complete output (HEADER, BODY and FOOTER)
		pScanTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
		// now we run the actual trajectory
		pScanTask->process(true);
	}
	catch (...)
	{
		std::cerr << "Error. Running the scan failed." << std::endl;
		// check if there are additional error messages
		if (CCopasiMessage::size() > 0)
		{
		// print the messages in chronological order
		std::cerr << CCopasiMessage::getAllMessageText(true);
		}
		return;
	}
	// restore the state of the trajectory
	pScanTask->restore();
}
*/

double runif(double min, double max)
{
	double d = mtrand();
	d *= (max - min);
	d += min;
	return d;
}

static int _COPASI_OPTIM_METHOD=10;
COPASIAPIEXPORT int cSetOptimizationMethod(const char * method)
{
    string s = string(CCopasiMethod::SubTypeName[0]),
              sMethod(method);

    for (int i=0; s.size() > 1; ++i)
	{
        s = string(CCopasiMethod::SubTypeName[i]);
	    if (sMethod.compare( s ) == 0)
		{
			_COPASI_OPTIM_METHOD = i;
			return 1;
		}
	}

	return 0;
}

double cMinimizeOrMaximize(copasi_model model, const char * formula, int minimize)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel || !pDataModel || !hash) return 0;

	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
	// get the optimization task
	COptTask* pOptTask = dynamic_cast<COptTask*>(TaskList["Optimization"]);

	//set method
	string s;
    pOptTask->setMethodType((CCopasiMethod::SubType)_COPASI_OPTIM_METHOD);

	// next we need to set subtask type on the problem
	COptProblem* pOptProblem = dynamic_cast<COptProblem*>(pOptTask->getProblem());
	pOptProblem->setModel(pModel);
	pOptProblem->setSubtaskType(CCopasiTask::steadyState);
	CSteadyStateTask* pTask = dynamic_cast<CSteadyStateTask*>(TaskList["Steady-State"]);
	CSteadyStateProblem* pProblem = dynamic_cast<CSteadyStateProblem*>(pTask->getProblem());
    //pProblem->setDuration(1);
    s = cSetAssignmentRuleHelper(model, formula);
    // now we set the objective function in the problem
    pOptProblem->setObjectiveFunction(s);
    // now we create the optimization items
    // i.e. the model elements that have to be changed during the optimization
    // in order to get to the optimal solution

    CCopasiVectorN< CModelValue > & params = pModel->getModelValues();

    for (int i=0; i < params.size(); ++i)
    {
        CModelValue * p = params[i];
        double x = p->getInitialValue();
        COptItem* pOptItem = &pOptProblem->addOptItem(CCopasiObjectName(p->getObject(CCopasiObjectName("Reference=InitialValue"))->getCN()));
        pOptItem->setStartValue(x);

        stringstream ss1,ss2;
        ss1 << (x*0.1);
        ss2 << (x*10.0);

        pOptItem->setLowerBound(CCopasiObjectName(ss1.str()));
        pOptItem->setUpperBound(CCopasiObjectName(ss2.str()));
    }

    COptMethod* pOptMethod = dynamic_cast<COptMethod*>(pOptTask->getMethod());
    CCopasiParameter* pParameter = pOptMethod->getParameter("Iteration Limit");
    if (pParameter != NULL)
        pParameter->setValue((long int)2000);
   pParameter = pOptMethod->getParameter("Tolerance");
   if (pParameter != NULL)
        pParameter->setValue(1.0e-5);

    try
    {
        // initialize the trajectory task
        pOptTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
        // now we run the actual trajectory
        pOptTask->process(true);
    }
    catch (...)
    {
        //std::cerr << "Error. Running the optimization failed." << std::endl;
        // check if there are additional error messages
        pOptTask->restore();
        return 0.0;
    }
    double bestValue = pOptProblem->getSolutionValue();
    const CVector< C_FLOAT64 > & solution = pOptProblem->getSolutionVariables();

    for (int i=0; i < params.size() && i < solution.size(); ++i)
    {
        CModelValue * p = params[i];
        p->setValue(solution[i]);
        p->setInitialValue(solution[i]);
        cout << solution[i] << "\n";
    }

    return bestValue;
}

double cMinimize(copasi_model model, const char * formula)
{
	return cMinimizeOrMaximize(model, formula, 0);
}

double cMaximize(copasi_model model, const char * formula)
{
	return cMinimizeOrMaximize(model, formula, 1);
}

int main1()
{
// initialize the backend library
CCopasiRootContainer::init(0, NULL);
assert(CCopasiRootContainer::getRoot() != NULL);
// create a new datamodel
CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();
assert(CCopasiRootContainer::getDatamodelList()->size() == 1);
// first we load a simple model
try
{
// load the model
pDataModel->importSBMLFromString(MODEL_STRING);
// we clear the message stack to get rid of all error messages
// this is important because otherwise the initialization of the
// fit task will fail
CCopasiMessage::clearDeque();
}
catch (...)
{
std::cerr << "Error while importing the model." << std::endl;
exit(1);
}
// now we need to run some time course simulation to get data to fit
// against
// get the trajectory task object
CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();
// get the optimization task
CTrajectoryTask* pTrajectoryTask = dynamic_cast<CTrajectoryTask*>(TaskList["Time-Course"]);
assert(pTrajectoryTask != NULL);
// if there isn’t one
if (pTrajectoryTask == NULL)
{
// create a new one
pTrajectoryTask = new CTrajectoryTask();
// add the new time course task to the task list
pDataModel->getTaskList()->add(pTrajectoryTask);
}
// run a deterministic time course
pTrajectoryTask->setMethodType(CCopasiMethod::deterministic);
// pass a pointer of the model to the problem
pTrajectoryTask->getProblem()->setModel(pDataModel->getModel());
// activate the task so that it will be run when the model is saved
// and passed to CopasiSE
pTrajectoryTask->setScheduled(true);
// get the problem for the task to set some parameters
CTrajectoryProblem* pProblem = dynamic_cast<CTrajectoryProblem*>(pTrajectoryTask-> getProblem());
// simulate 4000 steps
pProblem->setStepNumber(4000);
// start at time 0
pDataModel->getModel()->setInitialTime(0.0);
// simulate a duration of 400 time units
pProblem->setDuration(400);
// tell the problem to actually generate time series data
pProblem->setTimeSeriesRequested(true);
bool result = true;
try
{
// now we run the actual trajectory
pTrajectoryTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
result = pTrajectoryTask->process(true);
}
catch (...)
{
std::cerr <<
"Error. Running the time course simulation failed." << std::endl;
// check if there are additional error messages
if (CCopasiMessage::size() > 0)
{
// print the messages in chronological order
std::cerr << CCopasiMessage::getAllMessageText(true) << std::endl;
}
// clean up the library
CCopasiRootContainer::destroy();
exit(1);
}
if (result == false)
{
std::cerr << "An error occured while running the time course simulation." << std:: endl;
// check if there are additional error messages
if (CCopasiMessage::size() > 0)
{
// print the messages in chronological order
std::cerr << CCopasiMessage::getAllMessageText(true) << std::endl;
}
// clean up the library
CCopasiRootContainer::destroy();
exit(1);
}
// restore the trajectory task state
pTrajectoryTask->restore();
// we write the data to a file and add some noise to it
// This is necessary since COPASI can only read experimental data from
// file.
const CTimeSeries* pTimeSeries = &pTrajectoryTask->getTimeSeries();
// we simulated 4000 steps, including the initial state, this should be
// 4001 step in the timeseries
assert(pTimeSeries->getRecordedSteps() == 4001);
int i, iMax = pTimeSeries->getNumVariables();
// there should be four variables, the three metabolites and time
assert(iMax == 5);
int lastIndex = pTimeSeries->getRecordedSteps() - 1;
// open the file
// we need to remember in which order the variables are written to file
// since we need to specify this later in the parameter fitting task
std::set<int> indexSet;
std::vector<CMetab*> metabVector;
// write the header
// the first variable in a time series is a always time, for the rest
// of the variables, we use the SBML id in the header
double random = 0.0;
try
{
std::ofstream os("fakedata_example6.txt", std::ios_base::out | std::ios_base::trunc);
os << ("# time ");
for (i = 1; i < iMax; ++i)
{
std::string key = pTimeSeries->getKey(i);
CCopasiObject* pObject = CCopasiRootContainer::getKeyFactory()->get(key);
assert(pObject != NULL);
// only write header data or metabolites
if (dynamic_cast<const CMetab*>(pObject))
{
os << ", ";
os << pTimeSeries->getSBMLId(i, pDataModel);
CMetab* pM = dynamic_cast<CMetab*>(pObject);
indexSet.insert(i);
metabVector.push_back(pM);
}
}
os << std::endl;
double data = 0.0;
for (i = 0; i < lastIndex; ++i)
{
int j;
std::ostringstream s;
for (j = 0; j < iMax; ++j)
{
// we only want to write the data for metabolites
// the compartment does not interest us here
if (j == 0 || indexSet.find(j) != indexSet.end())
{
// write the data with some noise (+-5% max)
random = ((double)rand()) / (double)RAND_MAX;
data = pTimeSeries->getConcentrationData(i, j);
// don’t add noise to the time
if (j != 0)
{
data += data * (random * 0.1 - 0.05);
}
s << data;
s << ", ";
}
}
// remove the last two characters again
os << s.str().substr(0, s.str().length() - 2);
os << std::endl;
}
os.close();
}
catch (const std::exception& e)
{
std::cerr << "Error. Could not write time course data to file." << std::endl;
std::cout << e.what() << std::endl;
exit(1);
}
CModel* pModel = pDataModel->getModel();
// now we change the parameter values to see if the parameter fitting
// can really find the original values
random = (double)rand() / (double)RAND_MAX * 10.0;
CCopasiVectorN< CModelValue > & params = pModel->getModelValues();
// the parameter of a irreversible mass action is called k1
params[0]->setInitialValue(random);
// we know that it is an irreversible mass action, so there is one
// parameter
params[1]->setInitialValue(random);
CFitTask* pFitTask = dynamic_cast<CFitTask*>(TaskList["Parameter Estimation"]);
assert(pFitTask != NULL);
// the method in a fit task is an instance of COptMethod or a subclass of
// it.
COptMethod* pFitMethod = dynamic_cast<COptMethod*>(pFitTask->getMethod());
assert(pFitMethod != NULL);
// the object must be an instance of COptMethod or a subclass thereof
// (CFitMethod)
CFitProblem* pFitProblem = dynamic_cast<CFitProblem*>(pFitTask->getProblem());
assert(pFitProblem != NULL);
CExperimentSet* pExperimentSet = dynamic_cast<CExperimentSet*>(pFitProblem->getParameter ("Experiment Set"));
assert(pExperimentSet != NULL);
// first experiment (we only have one here)
CExperiment* pExperiment = new CExperiment(pDataModel);
assert(pExperiment != NULL);
// tell COPASI where to find the data
// reading data from string is not possible with the current C++ API
pExperiment->setFileName("fakedata_example6.txt");
// we have to tell COPASI that the data for the experiment is a komma
// separated list (the default is TAB separated)
pExperiment->setSeparator(",");
// the data start in row 1 and goes to row 4001
pExperiment->setFirstRow(1);
assert(pExperiment->getFirstRow() == 1);
pExperiment->setLastRow(4001);
assert(pExperiment->getLastRow() == 4001);
pExperiment->setHeaderRow(1);
assert(pExperiment->getHeaderRow() == 1);
pExperiment->setExperimentType(CCopasiTask::timeCourse);
assert(pExperiment->getExperimentType() == CCopasiTask::timeCourse);
pExperiment->setNumColumns(4);
assert(pExperiment->getNumColumns() == 4);
CExperimentObjectMap* pObjectMap = &pExperiment->getObjectMap();
assert(pObjectMap != NULL);
result = pObjectMap->setNumCols(4);
assert(result == true);
result = pObjectMap->setRole(0, CExperiment::time);
assert(result == true);
assert(pObjectMap->getRole(0) == CExperiment::time);
assert(pModel != NULL);
const CCopasiObject* pTimeReference = pModel->getObject(CCopasiObjectName("Reference=Time"));
assert(pTimeReference != NULL);
pObjectMap->setObjectCN(0, pTimeReference->getCN());
// now we tell COPASI which column contain the concentrations of
// metabolites and belong to dependent variables
pObjectMap->setRole(1, CExperiment::dependent);
CMetab* pMetab = metabVector[0];
assert(pMetab != NULL);
cout << "here 1\n";
const CCopasiObject* pParticleReference = pMetab->getObject(CCopasiObjectName("Reference=Concentration"));
assert(pParticleReference != NULL);
pObjectMap->setObjectCN(1, pParticleReference->getCN());
pObjectMap->setRole(2, CExperiment::dependent);
pMetab = metabVector[1];
assert(pMetab != NULL);
cout << "here 2\n";
pParticleReference = pMetab->getObject(CCopasiObjectName("Reference=Concentration"));
assert(pParticleReference != NULL);
pObjectMap->setObjectCN(2, pParticleReference->getCN());
pObjectMap->setRole(3, CExperiment::dependent);
pMetab = metabVector[2];
assert(pMetab != NULL);

pParticleReference = pMetab->getObject(CCopasiObjectName("Reference=Concentration"));
assert(pParticleReference != NULL);
pObjectMap->setObjectCN(3, pParticleReference->getCN());
pExperimentSet->addExperiment(*pExperiment);
assert(pExperimentSet->getExperimentCount() == 1);
// addExperiment makes a copy, so we need to get the added experiment
// again
delete pExperiment;
pExperiment = pExperimentSet->getExperiment(0);
assert(pExperiment != NULL);
// now we have to define the two fit items for the two local parameters
CModelValue* pParameter = params[0];
assert(pParameter != NULL);
// get the list where we have to add the fit items
CCopasiParameterGroup* pOptimizationItemGroup = dynamic_cast<CCopasiParameterGroup*>( pFitProblem->getParameter("OptimizationItemList"));
assert(pOptimizationItemGroup != NULL);
// define a CFitItem
const CCopasiObject* pParameterReference = pParameter->getObject(CCopasiObjectName("Reference=Value"));
assert(pParameterReference != NULL);
CFitItem* pFitItem1 = new CFitItem(pDataModel);
pFitItem1->setObjectCN(pParameterReference->getCN());
assert(pFitItem1 != NULL);
pFitItem1->setStartValue(4.0);
pFitItem1->setLowerBound(CCopasiObjectName("0.00001"));
pFitItem1->setUpperBound(CCopasiObjectName("10"));
// add the fit item
pOptimizationItemGroup->addParameter(pFitItem1);
pParameter = params[1];
assert(pParameter != NULL);
// define a CFitItem
pParameterReference = pParameter->getObject(CCopasiObjectName("Reference=Value"));
assert(pParameterReference != NULL);
CFitItem* pFitItem2 = new CFitItem(pDataModel);
pFitItem2->setObjectCN(pParameterReference->getCN());
assert(pFitItem2 != NULL);
pFitItem2->setStartValue(4.0);
pFitItem2->setLowerBound(CCopasiObjectName("0.00001"));
pFitItem2->setUpperBound(CCopasiObjectName("10"));
// add the fit item
pOptimizationItemGroup->addParameter(pFitItem2);
result = true;
try
{
// initialize the fit task
// we want complete output (HEADER, BODY and FOOTER)
result = pFitTask->initialize(CCopasiTask::OUTPUT_COMPLETE, pDataModel, NULL);
if (result == true)
{
// running the task for this example will probably take some time
std::cout << "This can take some time..." << std::endl;
result = pFitTask->process(true);
}
}
catch (...)
{
std::cerr << "Error. Parameter fitting failed." << std::endl;
// clean up the library
CCopasiRootContainer::destroy();
exit(1);
}
pFitTask->restore();
assert(result == true);
// assert that there are two optimization items
assert(pFitProblem->getOptItemList().size() == 2);
// the order should be the order in whih we added the items above
COptItem* pOptItem1 = pFitProblem->getOptItemList()[0];
COptItem* pOptItem2 = pFitProblem->getOptItemList()[1];
// the actual results are stored in the fit problem
assert(pFitProblem->getSolutionVariables().size() == 2);
std::cout << "value for " << pOptItem1->getObject()->getCN() << ": " << pFitProblem-> getSolutionVariables()[0] << std::endl;
std::cout << "value for " << pOptItem2->getObject()->getCN() << ": " << pFitProblem-> getSolutionVariables()[1] << std::endl;
// depending on the noise, the fit can be quite bad, so we are a litle
// relaxed here (we should be within 3% of the original values)
assert((fabs(pFitProblem->getSolutionVariables()[0] - 0.03) / 0.03) < 3e-2);
assert((fabs(pFitProblem->getSolutionVariables()[1] - 0.004) / 0.004) < 3e-2);
// clean up the library
CCopasiRootContainer::destroy();
}
const char* MODEL_STRING = "\
<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<!-- Created by COPASI version 4.5.30 (Debug) on 2009-03-30 08:01 with libSBML version\
3.3.2. -->\n\
<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">\n\
 <model metaid=\"COPASI1\" id=\"Model_1\" name=\"New Model\">\n\
<listOfUnitDefinitions>\n\
<unitDefinition id=\"volume\">\n\
<listOfUnits>\n\
<unit kind=\"litre\" scale=\"-3\"/>\n\
</listOfUnits>\n\
</unitDefinition>\n\
<unitDefinition id=\"substance\">\n\
<listOfUnits>\n\
<unit kind=\"mole\" scale=\"-3\"/>\n\
</listOfUnits>\n\
</unitDefinition>\n\
</listOfUnitDefinitions>\n\
<listOfParameters>\n\
<parameter id=\"k1\" name=\"k1\" value=\"0.004\"/>\n\
<parameter id=\"k1\" name=\"k1\" value=\"0.03\"/>\n\
</listOfParameters>\n\
<listOfCompartments>\n\
<compartment id=\"compartment_1\" name=\"compartment\" size=\"1\"/>\n\
</listOfCompartments>\n\
<listOfSpecies>\n\
<species id=\"species_1\" name=\"A\" compartment=\"compartment_1\" initialConcentration=\"5\"/>\n\
<species id=\"species_2\" name=\"B\" compartment=\"compartment_1\" initialConcentration=\"0\"/>\n\
<species id=\"species_3\" name=\"C\" compartment=\"compartment_1\" initialConcentration=\"0\"/>\n\
</listOfSpecies>\n\
<listOfReactions>\n\
<reaction id=\"reaction_1\" name=\"reaction\" reversible=\"false\">\n\
<listOfReactants>\n\
<speciesReference species=\"species_1\"/>\n\
</listOfReactants>\n\
<listOfProducts>\n\
<speciesReference species=\"species_2\"/>\n\
</listOfProducts>\n\
<kineticLaw>\n\
<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n\
<apply>\n\
<times/>\n\
<ci> compartment_1 </ci>\n\
<ci> k1 </ci>\n\
<ci> species_1 </ci>\n\
</apply>\n\
</math>\n\
</kineticLaw>\n\
</reaction>\n\
<reaction id=\"reaction_2\" name=\"reaction_1\" reversible=\"false\">\n\
<listOfReactants>\n\
<speciesReference species=\"species_2\"/>\n\
</listOfReactants>\n\
<listOfProducts>\n\
<speciesReference species=\"species_3\"/>\n\
</listOfProducts>\n\
<kineticLaw>\n\
<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n\
<apply>\n\
<times/>\n\
<ci> compartment_1 </ci>\n\
<ci> k1 </ci>\n\
<ci> species_2 </ci>\n\
</apply>\n\
</math>\n\
</kineticLaw>\n \
</reaction>\n \
</listOfReactions> \
</model>\n \
</sbml>";

