//std
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <map>
#include <limits> //get max and min for double

//copasi
#include "copasiSBWApi.h"
#include "copasi_api.h"
#define COPASI_MAIN 1
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
	//#define LIB_EXPORTS 1
	#include "src/antimony_api.h"
}

using namespace LIB_STRUCTURAL;
using namespace LIB_LA;
using namespace std;

// Do we need this??
//unsigned C_INT32 C_INVALID_INDEX = std::numeric_limits< unsigned C_INT32 >::max();

struct sb_CopasiPtr 
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

typedef map< string, sb_CopasiPtr > CCMap;
static int rename(string& target, const string& oldname,const string& newname0);
static int replaceSubstring(string& target, const string& oldname,const string& newname0);
static list< CCMap* > hashTablesToCleanup;
static list< copasiModel > copasiModelsToCleanup;

static boost::regex sb_stupidPowFunction("pow\\s*\\(\\s*([^,]+)\\s*,\\s*([^,]+)\\s*\\)", boost::regex::perl);

bool sb_contains(CCMap * hash, const string & s)
{
	return hash && (hash->find(s) != hash->end());
}

bool sb_contains(const string& str, const string & s)
{
	return str.find(s) != string::npos;
}

sb_CopasiPtr & sb_getHashValue(CCMap * hash, const string & s)
{
	return (*hash)[s];
}

void sb_hashInsert(CCMap * hash, const string & s, sb_CopasiPtr v)
{
	(*hash)[s] = v;
}

int sb_indexOf( list<string>& lst, const string & s)
{
	int k = 0;
	for (list<string>::iterator i=lst.begin(); i != lst.end(); i++, ++k)
		if ((*i) == s)
			return k;
	return -1;
}

int sb_indexOf( vector<string>& lst, const string & s)
{
	int k = 0;
	for (vector<string>::iterator i=lst.begin(); i != lst.end(); i++, ++k)
		if ((*i) == s)
			return k;
	return -1;
}

double sb_string_to_double( const std::string& s )
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
} 

list<string> sb_splitString(const string& seq, const string& _1cdelim);

void sCopasiInit()
{
	CCopasiRootContainer::init(0, NULL);
}

void sCopasiEnd()
{
	for (list<CCMap*>::iterator i = hashTablesToCleanup.begin(); i != hashTablesToCleanup.end(); i++)
		delete (*i);

	list< copasiModel > models = copasiModelsToCleanup;
	copasiModelsToCleanup.clear();

	for (list<copasiModel>::iterator i = models.begin(); i != models.end(); i++)
		sRemoveModel(*i);

	CCopasiRootContainer::destroy();
}

int sSetAssignmentRuleHelper(copasiModel , CMetab * , const char * );

static int DO_CLEANUP_ASSIGNMENT_RULES = 1;

void sEnableAssignmentRuleReordering()
{
	DO_CLEANUP_ASSIGNMENT_RULES = 1;
}

void sDisableAssignmentRuleReordering()
{
	DO_CLEANUP_ASSIGNMENT_RULES = 0;
}

int sCopasiCleanupAssignments(copasiModel model)
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
						sb_contains((*i).second.assignmentRule, names[j]))
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
					(*i).second.species->setStatus(CModelEntity::FIXED); //unused species
			}
			else
			{
				retval = retval * sSetAssignmentRuleHelper(model, (*i).second.species, (*i).second.assignmentRule.c_str());
			}
	return retval;
}

void sRemoveModel(copasiModel model)
{
	//remove from list
	for (list<copasiModel>::iterator i=copasiModelsToCleanup.begin(); i != copasiModelsToCleanup.end(); i++)
		if ((*i).CopasiDataModelPtr == model.CopasiDataModelPtr)
		{
			copasiModel m = { (void*)NULL, (void*)NULL, (void*)NULL, (char*)NULL };
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

void clearCopasiModel(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel || !pDataModel || !hash) return;
	
	sb_CopasiPtr p;

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


copasiModel sCreateModel(const char * name)
{
	sCopasiInit();

	CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();
	CModel* pModel = pDataModel->getModel();
	CCMap * qHash = new CCMap();
	copasiModel m = { (void*)(pModel) , (void*)(pDataModel), (void*)(qHash), (char*)(NULL), (char*)(NULL) };
	
	hashTablesToCleanup.push_back(qHash);
	copasiModelsToCleanup.push_back(m);

	pModel->setSBMLId( string(name) );
	pModel->setObjectName( string(name) );
	//pModel->setTimeUnit(CModel::s);
	//pModel->setVolumeUnit(CModel::microl);
	//pModel->setQuantityUnit(CModel::nMol);
	pModel->setTimeUnit(CModel::dimensionlessTime);
	pModel->setVolumeUnit(CModel::dimensionlessVolume);
	pModel->setQuantityUnit(CModel::dimensionlessQuantity);
	
	sCreateVariable(m, "time", "time");
	
	return m;
}

void sCreateSpecies(copasiCompartment compartment, const char* name, double iv)
{
	CModel* pModel = (CModel*)(compartment.CopasiModelPtr);
	CCompartment* pCompartment = (CCompartment*)(compartment.CopasiCompartmentPtr);
	CCMap * hash = (CCMap*)(compartment.qHash);
	CMetab* pSpecies;
	
	if (!pModel || !hash || !pCompartment) return;	
	if (sb_contains(hash, string(name)))
	{
		pSpecies = sb_getHashValue(hash,string(name)).species;
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

	sb_CopasiPtr copasiPtr = { 
			pSpecies->getCN(),
			pSpecies->getKey(),
			pSpecies,
			0,
			0,
			0,
			"",
			true};

	sb_hashInsert(hash, 
				(pCompartment->getObjectName() + string("_") + string(name)),
				copasiPtr
				);

	sb_hashInsert(hash, string(name), copasiPtr); //for speedy lookup
}

copasiCompartment sCreateCompartment(copasiModel model, const char* name, double volume)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCompartment* pCompartment = pModel->createCompartment(name, volume);
	CCMap * hash = (CCMap*)(model.qHash);
	copasiCompartment c = { 0, 0 , 0};
	
	if (!pModel || !hash) return c;
	c.CopasiModelPtr = (void*)(pModel);
	c.qHash = model.qHash;
	
	if (sb_contains(hash,string(name)))
	{
		if (sb_getHashValue(hash, string(name)).compartment)
			c.CopasiCompartmentPtr = (void*)(sb_getHashValue(hash, string(name)).compartment);
		
		return c;
	}
	else
	{
		c.CopasiCompartmentPtr = (void*)(pCompartment);
	}	
	
	sb_CopasiPtr copasiPtr = { 
			pCompartment->getCN(),
			pCompartment->getKey(),
			0,
			pCompartment,
			0,
			0,
			"",
			false};

	sb_hashInsert(hash, string(name),copasiPtr); //for speedy lookup
	
	return c;
}


int sSetGlobalParameter(copasiModel model, const char * name, double value)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CModelValue * pValue = NULL;
	
	if (!hash || !pModel) return 0;
	
	if (sb_contains(hash,s) && 
		(pValue = sb_getHashValue(hash,s).param))
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
	
		sb_CopasiPtr copasiPtr = {
				pValue->getCN(),
				pValue->getKey(),
				0,
				0,
				0,
				pValue,
				"",
				true};

		sb_hashInsert(hash, s, copasiPtr); //for speedy lookup
	}
	
	return 0;
}


int sSetValue(copasiModel model, const char * name, double value)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	
	if (!hash) return 0;
	
	if (!sb_contains(hash,s))
	{
		sSetGlobalParameter(model,name,value);
		return 0;
	}

	sb_CopasiPtr p = sb_getHashValue(hash,s);
	
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
	
	sSetGlobalParameter(model,name,value);
	return 0;
}

void sSetVolume(copasiModel model, const char * name, double vol)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CCompartment* pVol = NULL;
	
	if (!hash) return;
	
	if (sb_contains(hash,s) && 
		(pVol = sb_getHashValue(hash,s).compartment))
	{
		pVol->setInitialValue(vol);
		pVol->setValue(vol);
	}
}

void sSetConcentration(copasiModel model, const char * name, double conc)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CMetab* pSpecies = NULL;
	
	if (!hash) return;
	
	if (sb_contains(hash,s) && 
		(pSpecies = sb_getHashValue(hash,s).species))
	{
		pSpecies->setConcentration(conc);
	}
}

void sSetInitialConcentration(copasiModel model, const char * name, double conc)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CMetab* pSpecies = NULL;
	
	if (!hash) return;
	
	if (sb_contains(hash,s) && 
		(pSpecies = sb_getHashValue(hash,s).species))
	{
		pSpecies->setInitialConcentration(conc);
	}
}

void sSetAmount(copasiModel model, const char * name, double amnt)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CMetab* pSpecies = NULL;
	
	if (!hash) return;
	
	if (sb_contains(hash,s) && 
		(pSpecies = sb_getHashValue(hash,s).species))
	{
		pSpecies->setValue(amnt);
		pSpecies->setInitialValue(amnt);
	}
}

void sSetSpeciesType(copasiModel model, const char * name, int isBoundary)
{
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	CMetab* pSpecies = NULL;
	
	if (!hash) return;
	
	if (sb_contains(hash,s) && 
		(pSpecies = sb_getHashValue(hash,s).species))
	{
		double iv = pSpecies->getInitialConcentration();

		if (isBoundary)
			pSpecies->setStatus(CModelEntity::FIXED);
		else
			pSpecies->setStatus(CModelEntity::REACTIONS);
		
		pSpecies->setConcentration(iv);
		pSpecies->setValue(iv);
		pSpecies->setInitialValue(iv);
		pSpecies->setInitialConcentration(iv);
	}
}

void sCreateSpecies(CModel * pModel, CCMap * hash, string s)
{
	if (!pModel || !hash) return;

	CCopasiVectorNS < CCompartment > & compartments = pModel->getCompartments();
		if (compartments.size() > 0 && compartments[0] != NULL)
		{
			CCompartment* pCompartment = compartments[0];
			if (pCompartment)
			{
				string s(pCompartment->getObjectName());
				if	(sb_contains(hash,s) &&
					sb_getHashValue(hash,s).compartment)
					{
						copasiCompartment c = { (void*)sb_getHashValue(hash,s).compartment, pModel, hash };
						sCreateSpecies(c,s.c_str(),0.0);
					}
			}
		}
}


int sSetAssignmentRule(copasiModel model, const char * name, const char * formula)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	int i;
	bool retval=true;
	
	if (!pModel || !hash) return 0;
	
	if (!sb_contains(hash,s))
	{
		sCreateSpecies(pModel, hash, s);
	}

	if (sb_contains(hash,s) && sb_getHashValue(hash,s).species)
	{
		sb_CopasiPtr & p = (*hash)[s];
		p.assignmentRule = string(formula);
		boost::regex_replace(p.assignmentRule, sb_stupidPowFunction, string("((\\1)^(\\2))"));
		//cout << p.assignmentRule.c_str() << "\n";
		return 1;
	}
	return 0;
}

int sSetAssignmentRuleHelper (copasiModel model, CMetab* pSpecies, const char * formula)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	int i;
	bool retval=true;
	
	if (!pModel || !hash || !pSpecies) return 0;
	
	if (formula)
	{
		pSpecies->setStatus(CModelEntity::ASSIGNMENT);
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
					if (!sb_contains(hash,s0))
						sSetGlobalParameter(model,pParam->getObjectName().c_str(),1.0);
					if (sb_contains(hash,s0))
					{
					 	string s1("<");
							s1 += sb_getHashValue(hash,s0).name;
							s1 += string(">");
						//std::cout << "sub " << s1 << "  for " << s0 << "\n";
						rename(qFormula,s0,s1);
					}
				}
			}
		}

		//std::cout << qFormula << "\n";
		retval = retval & pSpecies->setExpression(qFormula);
	}
	else
		pSpecies->setStatus(CModelEntity::REACTIONS);
	
	return (int)retval;
}

int sCreateVariable(copasiModel model, const char * name, const char * formula)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!hash || !pModel) return 0;

	CModelValue* pModelValue;
	string qname(name);

	if (sb_contains(hash,qname))
	{
			sb_CopasiPtr ptr = sb_getHashValue(hash,qname);
			if (ptr.species)
				return sSetAssignmentRule(model, name, formula);
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
				if (!sb_contains(hash,s0))
					sSetGlobalParameter(model,pParam->getObjectName().c_str(),1.0);				
				
				if (sb_contains(hash,s0))
				{
				 	string s1("<");
						s1 += sb_getHashValue(hash,s0).name;
						s1 += string(">");
					rename(qFormula,s0,s1);
				}
			}
		}
	}

	string sFormula( qFormula );
	
	retval = retval & pModelValue->setInitialExpression(sFormula);
	retval = retval & pModelValue->setExpression(sFormula);
	
	sb_CopasiPtr copasiPtr = {
			pModelValue->getCN(),
			pModelValue->getKey(),
			0,
			0,
			0,
			pModelValue,
			"",
			true};

	sb_hashInsert(hash, qname, copasiPtr); //for speedy lookup
	
	return (int)retval;
}

int sCreateEvent(copasiModel model, const char * name, const char * trigger, const char * variable, const char * formula)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!hash || !pModel)
	{
		return 0;
	}
	
	int i;
	bool retval = true;

	if (!sb_contains(hash,string(variable)))
	{
		sSetGlobalParameter(model,variable,1.0);
	}
	
	if (!sb_contains(hash,string(variable))) return 0;

	sb_CopasiPtr ptr = sb_getHashValue(hash,string(variable));
	
	if (!ptr.species && !ptr.param) return 0;

	CEvent * pEvent = pModel->createEvent(string(name));

	CFunction pFunction;
	string qFormula(trigger);
	replaceSubstring(qFormula,">","gt");
	replaceSubstring(qFormula,"<","lt");
	replaceSubstring(qFormula,">=","ge");
	replaceSubstring(qFormula,"<=","le");
	replaceSubstring(qFormula,"=","eq");

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
				if (!sb_contains(hash,s0))
					sSetGlobalParameter(model,s0.c_str(),1.0);
				if (sb_contains(hash,s0))
				{
				 	string s1("<");
						s1 += sb_getHashValue(hash,s0).name;
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
				if (!sb_contains(hash,s0))
					sSetGlobalParameter(model,s0.c_str(),1.0);

				if (sb_contains(hash,s0))
				{
				 	string s1("<");
						s1 += sb_getHashValue(hash,s0).name;
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

copasiReaction sCreateReaction(copasiModel model, const char* name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel || !hash)
	{
		copasiReaction r = { 0, 0, 0 };
		return r;
	}
	
	CReaction* pReaction = pModel->createReaction(name);
	
	copasiReaction r = { (void*)(pReaction), (void*)(pModel), (void*)hash };
	
	sb_CopasiPtr copasiPtr = { 
			pReaction->getCN(),
			pReaction->getKey(),
			0,
			0,
			pReaction,
			0,
			"",
			false};

	string qname(name);
	sb_hashInsert(hash, qname, copasiPtr); //for speedy lookup

	return r;
}

void sAddReactant(copasiReaction reaction, const char * species, double stoichiometry)
{
	CReaction* pReaction = (CReaction*)(reaction.CopasiReactionPtr);
	CCMap * hash = (CCMap*)(reaction.qHash);
	
	if (!pReaction || !hash)
	{
		return;
	}

	CMetab* pSpecies = NULL;

	string s(species);
	if (!sb_contains(hash, s))
	{
		sCreateSpecies((CModel*)reaction.CopasiModelPtr, (CCMap*)reaction.qHash, s);
	}
	
	sb_CopasiPtr p = sb_getHashValue(hash,s);
	if (sb_contains(hash,s) && (pSpecies = p.species))
	{
		CChemEq* pChemEq = &pReaction->getChemEq();
		pChemEq->addMetabolite(pSpecies->getKey(), stoichiometry, CChemEq::SUBSTRATE);
		p.unused = false;
	}
}

void sAddProduct(copasiReaction reaction, const char * species, double stoichiometry)
{
	CReaction* pReaction = (CReaction*)(reaction.CopasiReactionPtr);
	CCMap * hash = (CCMap*)(reaction.qHash);
	CMetab* pSpecies = NULL;
	
	if (!pReaction || !hash) return;
	
	string s(species);
	if (!sb_contains(hash, s))
	{
		sCreateSpecies((CModel*)reaction.CopasiModelPtr, (CCMap*)reaction.qHash, s);
	}

	sb_CopasiPtr p = sb_getHashValue(hash,s);
	if (sb_contains(hash,s) && (pSpecies = p.species))
	{
		CChemEq* pChemEq = &pReaction->getChemEq();
		pChemEq->addMetabolite(pSpecies->getKey(), stoichiometry, CChemEq::PRODUCT);
		p.unused = false;
	}
}

int sSetReactionRate(copasiReaction reaction, const char * formula)
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
		//formula2.replace(sb_stupidPowFunction, QString("((\\1)^(\\2))"));
		boost::regex_replace(formula2, sb_stupidPowFunction, string("((\\1)^(\\2))"));
		
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
				
				if (!sb_contains(hash,s))
				{
					copasiModel model = { (void*)(pModel) , (void*)(0), (void*)(hash), (char*)(NULL), (char*)(NULL) };
					sSetGlobalParameter(model,pParam->getObjectName().c_str(),1.0);				
				}
				
				if (sb_contains(hash,s))
				{
					sb_CopasiPtr p = sb_getHashValue(hash,s);
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

void sCompileModel(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel) return;
	sCopasiCleanupAssignments(model);
	
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

sb_matrix sSimulate (copasiModel model, double startTime, double endTime, int numSteps, CCopasiMethod::SubType method)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);

	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	//sCompileModel(model);

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
	
		sb_matrix output = sb_createMatrix(rows, cols);
		list<string> colnames;

		for (j=1; j < cols; ++j)
			colnames.push_back( timeSeries.getTitle(j) );

		colnames.sort();
		colnames.push_front(timeSeries.getTitle(0).c_str());

		j = 0;
		for (list<string>::iterator it=colnames.begin(); j < cols && it != colnames.end(); ++j, it++)
			sb_setColumnName( output, j, (*it).c_str() );
	
		for (j=0; j < cols; ++j)
		{
			k = sb_indexOf(colnames,timeSeries.getTitle(j));
			for (i=0; i < rows; ++i)
				   sb_setMatrixValue( output, i, k, timeSeries.getConcentrationData(i,j) );
		}
		return output;
	}
	return sb_createMatrix(0,0);
}

sb_matrix sSimulateDeterministic(copasiModel model, double startTime, double endTime, int numSteps)
{
	return sSimulate(model,startTime,endTime,numSteps,CCopasiMethod::deterministic);
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT double cOneStep(copasiModel model, double timeStep)
{
	return 0.0;
}

sb_matrix sSimulateTauLeap(copasiModel model, double startTime, double endTime, int numSteps)
{
	return sSimulate(model,startTime,endTime,numSteps,CCopasiMethod::tauLeap);
}

sb_matrix sSimulateStochastic(copasiModel model, double startTime, double endTime, int numSteps)
{
	return sSimulate(model,startTime,endTime,numSteps,CCopasiMethod::stochastic);
}

sb_matrix sSimulateHybrid(copasiModel model, double startTime, double endTime, int numSteps)
{
	return sSimulate(model,startTime,endTime,numSteps,CCopasiMethod::hybridLSODA);
}

void sWriteSBMLFile(copasiModel model, const char * filename)
{
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	printf ("Stage A1\n");
	if (pDataModel) {
	   printf ("Stage B\n");
		 pDataModel->exportSBML(filename, true, 2, 3);
		 printf ("Stage C\n");
 }
}

void sWriteAntimonyFile(copasiModel model, const char * filename)
{
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
  printf ("Stage 1\n"); 

	if (pDataModel)
	{
    printf ("Stage 2\n"); 
		pDataModel->exportSBML(filename, true, 2, 3);
		printf ("Stage 2\n"); 
		loadSBMLFile(filename);
    printf ("Stage 3\n"); 
		writeAntimonyFile(filename,NULL);
	}
}

copasiModel sReadAntimonyString(const char * model)
{
	loadString(model); //load Antimony
	const char * s = getSBMLString("__main");  //Antimony -> SBML (at worst, an empty model)
	copasiModel m = sReadSBMLString(s);
	freeAll(); //free Antimony
	return m;
}

copasiModel sReadAntimonyFile(const char * filename)
{
	loadFile(filename); //load Antimony
	const char * s = getSBMLString("__main");  //Antimony -> SBML (at worst, an empty model)
	copasiModel m = sReadSBMLString(s);
	freeAll(); //free Antimony
	return m;
}

copasiModel sReadSBMLFile(const char * filename)
{
	sCopasiInit();
	
	CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();
	CModel* pModel = 0;
	CCMap * qHash = 0;	
	char * error = NULL;
	char * warning = NULL;
	string s;
	CCopasiMessage::Type type;

	try 
	{
		pDataModel->importSBML(filename); //SBML -> COPASI
		s = CCopasiMessage::getAllMessageText();
		type = CCopasiMessage::getHighestSeverity();
		pModel = pDataModel->getModel();
		qHash = new CCMap();	
	}
	catch(...)
	{
		s = CCopasiMessage::getAllMessageText();
		type = CCopasiMessage::EXCEPTION;
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

	copasiModel m = { (void*)(pModel) , (void*)(pDataModel), (void*)(qHash), (char*)(error), (char*)warning };
	if (pModel && qHash)
	{
		hashTablesToCleanup.push_back( qHash );
		copasiModelsToCleanup.push_back(m);
	}
	return m;
}

copasiModel sReadSBMLString(const char * sbml)
{
	sCopasiInit();
	
	CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();
	CModel* pModel = 0;
	CCMap * qHash = 0;	
	char * error = NULL;
	string s;
	try 
	{
		pDataModel->importSBMLFromString(sbml); //SBML -> COPASI	
		s = CCopasiMessage::getAllMessageText();
		pModel = pDataModel->getModel();
		qHash = new CCMap();	
	}
	catch(...)
	{
		s = CCopasiMessage::getAllMessageText();
	}

	int len = s.length();
	if (len > 1)
	{
		error = (char*)malloc((1+len) * sizeof(char));
		if (error)
		{
			for (int i=0; i < len; ++i) error[i] = s[i];
			error[len-1] = 0;
		}
	}
	copasiModel m = { (void*)(pModel) , (void*)(pDataModel), (void*)(qHash), (char*)(error) };
	if (pModel && qHash)
	{
		hashTablesToCleanup.push_back(qHash);
		copasiModelsToCleanup.push_back(m);
	}
	return m;
}

sb_matrix sGetJacobian(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}
	
	const CArrayAnnotation* pAJ = pTask->getJacobianAnnotated();
	//const CEigen & sGetEigenvalues() const;
	
	if (pAJ && pAJ->dimensionality() == 2)
	{
		vector<unsigned int> index(2);
		const vector<string>& annotations = pAJ->getAnnotationsString(1);
		
		int n = annotations.size();
		sb_matrix J = sb_createMatrix(n,n);
		
		for (int i=0; i < J.rows; ++i)
		{
			sb_setRowName(J, i, annotations[i].c_str());
			sb_setColumnName(J, i, annotations[i].c_str());
		}
		
		for (int i=0; i < n; ++i)
		{
			index[0] = i;
			for (int j=0; j < n; ++j)
			{
				index[1] = j;
				sb_setMatrixValue(J, i, j, (*pAJ->array())[index]);
			}
		}
		
		return J;
	}

	return sb_createMatrix(0,0);
}

sb_matrix sGetSteadyStateUsingSimulation(copasiModel model, int maxiter)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	
	sCompileModel(model);

    int iter = 0;
    double err = 2.0, eps = 0.01, time = 10.0;

   	CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();

    while (iter < maxiter && err > eps)
    {
        ++iter;
        time *= 2.0;

	    CTrajectoryTask* pTask = dynamic_cast<CTrajectoryTask*>(TaskList["Time-Course"]);
	    // if there isn’t one
	    if (pTask == NULL)
	    {
		    pTask = new CTrajectoryTask();
		    TaskList.remove("Time-Course");
		    TaskList.add(pTask, true);
	    }
	
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
                //pTask->restore();
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
                sb_matrix output = sb_createMatrix(cols, 1);
				int k;
				list<string>::iterator it=colnames.begin();
                for (int i=0; i < cols && it != colnames.end(); ++i, it++)
                {
					k = sb_indexOf(colnames, timeSeries.getTitle(i+1));
            		sb_setRowName( output, i, (*it).c_str() );
                    sb_setMatrixValue( output, k, 0, timeSeries.getConcentrationData(j,i+1) );
                }
                return output;
            }
        }
	}

    sb_matrix m = sb_createMatrix(0,0);
	return m;
}

sb_matrix sGetSteadyState(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	
	sCompileModel(model);

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
	
	if (pTrajTask && pTrajTask->setMethodType(CCopasiMethod::deterministic))
	{
		//set the start and end time, number of steps, and save output in memory
		CTrajectoryProblem* pProblem=(CTrajectoryProblem*)pTrajTask->getProblem();
		pProblem->setModel(pModel);
		pTrajTask->setScheduled(true);
		pProblem->setStepNumber(10);
		pProblem->setDuration(10.0);
		pDataModel->getModel()->setInitialTime(0.0);
		pProblem->setTimeSeriesRequested(true);
		try
		{
			pTrajTask->initialize(CCopasiTask::ONLY_TIME_SERIES, pDataModel, NULL);
			pTrajTask->process(true);
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
		pTask->initialize(CCopasiTask::OUTPUT, pDataModel, NULL);
		pTask->process(true);
	}
	catch (...)
	{
		cerr << "Error when computing steady state." << endl;
		return sb_createMatrix(0,0);
	}

	pTrajTask = dynamic_cast<CTrajectoryTask*>(TaskList["Time-Course"]);
	// if there isn’t one
	if (pTrajTask == NULL)
	{
		pTrajTask = new CTrajectoryTask();
		TaskList.remove("Time-Course");
		TaskList.add(pTrajTask, true);
	}
	
	CCopasiMessage::clearDeque();
	
	if (pTrajTask && pTrajTask->setMethodType(CCopasiMethod::deterministic))
	{
		//set the start and end time, number of steps, and save output in memory
		CTrajectoryProblem* pProblem=(CTrajectoryProblem*)pTrajTask->getProblem();
		pProblem->setModel(pModel);
		pTrajTask->setScheduled(true);
		pProblem->setStepNumber(10);
		pProblem->setDuration(10.0);
		pDataModel->getModel()->setInitialTime(0.0);
		pProblem->setTimeSeriesRequested(true);
		try
		{
			pTrajTask->initialize(CCopasiTask::ONLY_TIME_SERIES, pDataModel, NULL);
			pTrajTask->process(true);
			//pTrajTask->restore();
		}
		catch(...)
		{
			cerr << CCopasiMessage::getAllMessageText(true);
			pTrajTask = NULL;
		}
	}
	
	if (pTrajTask)
	{
		const CTimeSeries & timeSeries = pTrajTask->getTimeSeries();
		int rows = (pModel->getNumMetabs());
		int i,j,k;

		sb_matrix output = sb_createMatrix(rows, 1);

		list<string> rownames;
		for (i=0; i < rows; ++i)
        	rownames.push_back( timeSeries.getTitle(i+1) );
		rownames.sort();
		j = timeSeries.getRecordedSteps() - 1;	

		list<string>::iterator it=rownames.begin();
        for (i=0; i < rows && it != rownames.end(); ++i, it++)
        {
			k = sb_indexOf(rownames ,  timeSeries.getTitle(i+1) );
    		sb_setRowName( output, i, (*it).c_str() );
            sb_setMatrixValue( output, k, 0, timeSeries.getConcentrationData(j,i+1) );
        }

		return output;
	}

	return sGetSteadyStateUsingSimulation(model,10);
}

sb_matrix sGetEigenvalues(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}

	const CEigen & eigen = pTask->getEigenValues();
	const CVector< C_FLOAT64 > & im = eigen.getI(), 
													& re = eigen.getR();

	sb_matrix E = sb_createMatrix(im.size(),2);
	
	sb_setColumnName(E, 0, "real\0");
	sb_setColumnName(E, 1, "imaginary\0");
	for (int i=0; i < im.size() && i < re.size(); ++i)
	{
		sb_setMatrixValue(E, i,0,re[i]);
		sb_setMatrixValue(E, i,1,im[i]);
	}
	
	return E;
}

sb_matrix sGetUnscaledElasticities(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}
	
	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());
	
	if (!mcaMethod) return sb_createMatrix(0,0);
	
	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getUnscaledElasticities();
	const CArrayAnnotation * annot = mcaMethod->getUnscaledElasticitiesAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	sb_matrix M = sb_createMatrix(rows, cols);
	
	for (int i=0; i < rows; ++i)
		sb_setRowName(M, i, rownames[i].c_str());
	
	for (int i=0; i < cols; ++i)
		sb_setColumnName(M, i, colnames[i].c_str());
	
	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			sb_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

sb_matrix sGetUnscaledConcentrationControlCoeffs(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}
	
	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());
	
	if (!mcaMethod) return sb_createMatrix(0,0);
	
	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getUnscaledConcentrationCC();
	const CArrayAnnotation * annot = mcaMethod->getUnscaledConcentrationCCAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	sb_matrix M = sb_createMatrix(rows, cols);
	
	for (int i=0; i < rows; ++i)
		sb_setRowName(M, i, rownames[i].c_str());
	
	for (int i=0; i < cols; ++i)
		sb_setColumnName(M, i, colnames[i].c_str());
	
	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			sb_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

sb_matrix sGetUnscaledFluxControlCoeffs(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}
	
	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());
	
	if (!mcaMethod) return sb_createMatrix(0,0);
	
	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getUnscaledFluxCC();
	const CArrayAnnotation * annot = mcaMethod->getUnscaledFluxCCAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	sb_matrix M = sb_createMatrix(rows, cols);
	
	for (int i=0; i < rows; ++i)
		sb_setRowName(M, i, rownames[i].c_str());
	
	for (int i=0; i < cols; ++i)
		sb_setColumnName(M, i, colnames[i].c_str());
	
	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			sb_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

sb_matrix sGetScaledElasticities(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}
	
	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());
	
	if (!mcaMethod) return sb_createMatrix(0,0);
	
	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getScaledElasticities();
	const CArrayAnnotation * annot = mcaMethod->getScaledElasticitiesAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	sb_matrix M = sb_createMatrix(rows, cols);
	
	for (int i=0; i < rows; ++i)
		sb_setRowName(M, i, rownames[i].c_str());
	
	for (int i=0; i < cols; ++i)
		sb_setColumnName(M, i, colnames[i].c_str());
	
	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			sb_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

sb_matrix sGetScaledConcentrationConcentrationCoeffs(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}
	
	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());
	
	if (!mcaMethod) return sb_createMatrix(0,0);
	
	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getScaledConcentrationCC();
	const CArrayAnnotation * annot = mcaMethod->getScaledConcentrationCCAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	int rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();

	sb_matrix M = sb_createMatrix(rows, cols);
	
	for (int i=0; i < rows; ++i)
		sb_setRowName(M, i, rownames[i].c_str());
	
	for (int i=0; i < cols; ++i)
		sb_setColumnName(M, i, colnames[i].c_str());
	
	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			sb_setMatrixValue(M, i, j, cmatrix(i,j));
		}
	return M;
}

sb_matrix sGetScaledFluxControlCoeffs(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}
	
	CMCAMethod * mcaMethod = dynamic_cast<CMCAMethod*>(pTask->getMethod());
	
	if (!mcaMethod) return sb_createMatrix(0,0);
	
	const CMatrix<C_FLOAT64> & cmatrix = mcaMethod->getScaledFluxCC();
	const CArrayAnnotation * annot = mcaMethod->getScaledFluxCCAnn();
	const vector<string>& rownames = annot->getAnnotationsString(1),
												   & colnames = annot->getAnnotationsString(0);

	size_t rows = cmatrix.numRows(), cols = cmatrix.numCols();
	if (rows > rownames.size()) rows = rownames.size();
	if (cols > colnames.size()) cols = colnames.size();
	
	sb_matrix M = sb_createMatrix(rows, cols);

	for (int i=0; i < rows; ++i)
		sb_setRowName(M, i, rownames[i].c_str());
	
	for (int i=0; i < cols; ++i)	
		sb_setColumnName(M, i, colnames[i].c_str());
	
	for (int i=0; i < rows; ++i)
		for (int j=0; j < cols; ++j)
		{
			sb_setMatrixValue(M, i, j, cmatrix(i,j));
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

sb_matrix sGetReducedStoichiometryMatrix(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
	CCopasiVector< CMetab > & species = pModel->getMetabolitesX();
	CCopasiVectorNS < CReaction > & reacs = pModel->getReactions();
	CMatrix < C_FLOAT64 > stoi = pModel->getRedStoi();

	sb_matrix N = sb_createMatrix( stoi.numRows(), stoi.numCols() );

	for  (int i=0; i < N.rows && i < species.size(); ++i)
		if (species[i])
			sb_setRowName(N, i, species[i]->getObjectName().c_str());

	for  (int i=0; i < N.cols && i < reacs.size(); ++i)
		if (reacs[i])
			sb_setColumnName(N, i, reacs[i]->getObjectName().c_str());

	for  (int i=0; i < N.rows; ++i)
		for  (int j=0; j < N.cols; ++j)
			sb_setMatrixValue(N, i, j, stoi(i,j));

	return N;
}

sb_matrix sGetFullStoichiometryMatrix(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
	CCopasiVector< CMetab > & species = pModel->getMetabolites();
	CCopasiVectorNS < CReaction > & reacs = pModel->getReactions();
	CMatrix < C_FLOAT64 > stoi = pModel->getStoi();

	sb_matrix N = sb_createMatrix( stoi.numRows(), stoi.numCols() );

	for  (int i=0; i < N.rows && i < species.size(); ++i)
		if (species[i])
			sb_setRowName(N, i, species[i]->getObjectName().c_str());

	for  (int i=0; i < N.cols && i < reacs.size(); ++i)
		if (reacs[i])
			sb_setColumnName(N, i, reacs[i]->getObjectName().c_str());

	for  (int i=0; i < N.rows; ++i)
		for  (int j=0; j < N.cols; ++j)
			sb_setMatrixValue(N, i, j, stoi(i,j));

	return N;
}

sb_matrix sGetElementaryFluxModes(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);
	sCompileModel(model);
	
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
		return sb_createMatrix(0,0);
	}
	
	const vector< CFluxMode > & fluxModes = pTask->getFluxModes();
	CEFMProblem* pProblem = dynamic_cast<CEFMProblem*>(pTask->getProblem());
	
	if (!pProblem)
		return sb_createMatrix(0,0);

	vector< const CReaction * > & reactions = pProblem->getReorderedReactions();
	sb_matrix M = sb_createMatrix( reactions.size() , fluxModes.size() );
	for (int i=0; i < reactions.size(); ++i)
		sb_setRowName(M, i, reactions[i]->getObjectName().c_str());
	
	for (int i=0; i < fluxModes.size(); ++i)
	{
		CFluxMode::const_iterator itMode = fluxModes[i].begin();
		CFluxMode::const_iterator endMode = fluxModes[i].end();
		for (; itMode != endMode; ++itMode)
			sb_setMatrixValue( M, itMode->first, i, itMode->second);
	}
	return M;
}

/* LIBSTRUCTURAL */

sb_matrix sb_convertFromDoubleMatrix(DoubleMatrix& matrix, vector< string > &rowNames, vector< string > &colNames)
{
	sb_matrix m = sb_createMatrix(matrix.numRows(), matrix.numCols());
	
	for (int i=0; i < m.rows && i < rowNames.size(); ++i)
		sb_setRowName(m, i, rowNames[i].c_str());

	for (int i=0; i < m.cols && i < colNames.size(); ++i)
		sb_setColumnName(m, i, colNames[i].c_str());

	for (int i=0; i < m.rows; ++i)
		for (int j=0; j < m.cols; ++j)
			sb_setMatrixValue(m, i, j, matrix(i,j));
	
	return m;
}

void sb_convertToDoubleMatrix(sb_matrix m, DoubleMatrix & matrix, vector< string > &rowNames, vector< string > &colNames)
{
	matrix.resize(m.rows, m.cols);
	
	rowNames.resize(m.rows);
	colNames.resize(m.cols);
	
	for (int i=0; i < m.rows && i < rowNames.size(); ++i)
		rowNames[i] = string(sb_getRowName(m, i));

	for (int i=0; i < m.cols && i < colNames.size(); ++i)
		colNames[i] = string(sb_getColumnName(m, i));
	
	for (int i=0; i < m.rows; ++i)
		for (int j=0; j < m.cols; ++j)
			matrix(i,j) = sb_getMatrixValue(m, i, j);
}

sb_matrix sGetGammaMatrix(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);

	//get stoichiometry
	sb_matrix sb_N = sGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;	
	DoubleMatrix N;
	
	sb_convertToDoubleMatrix( sb_N , N, rowNames, colNames );
	
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
	sb_matrix m = sb_convertFromDoubleMatrix(*matrix, rowNames, colNames);
	
	//cleanup
	//delete matrix;
	sb_deleteMatrix(sb_N);
	
	return m;
}

sb_matrix sGetKMatrix(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);

	//get stoichiometry
	sb_matrix sb_N = sGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;	
	DoubleMatrix N;
	
	sb_convertToDoubleMatrix( sb_N , N, rowNames, colNames );
	
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
	sb_matrix m = sb_convertFromDoubleMatrix(*matrix, rowNames, colNames);
	
	//cleanup
	//delete matrix;
	sb_deleteMatrix(sb_N);
	
	return m;
}

sb_matrix sGetLinkMatrix(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);

	//get stoichiometry
	sb_matrix sb_N = sGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;	
	DoubleMatrix N;
	
	sb_convertToDoubleMatrix( sb_N , N, rowNames, colNames );
	
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
	sb_matrix m = sb_convertFromDoubleMatrix(*matrix, rowNames, colNames);
	
	//cleanup
	//delete matrix;
	sb_deleteMatrix(sb_N);
	
	return m;
}

sb_matrix sGetK0Matrix(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);

	//get stoichiometry
	sb_matrix sb_N = sGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;	
	DoubleMatrix N;
	
	sb_convertToDoubleMatrix( sb_N , N, rowNames, colNames );
	
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
	sb_matrix m = sb_convertFromDoubleMatrix(*matrix, rowNames, colNames);
	
	//cleanup
	//delete matrix;
	sb_deleteMatrix(sb_N);
	
	return m;
}

sb_matrix sGetL0Matrix(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCopasiDataModel* pDataModel = (CCopasiDataModel*)(model.CopasiDataModelPtr);
	if (!pModel || !pDataModel) return sb_createMatrix(0,0);

	//get stoichiometry
	sb_matrix sb_N = sGetFullStoichiometryMatrix(model);
	vector< string > rowNames, colNames;	
	DoubleMatrix N;
	
	sb_convertToDoubleMatrix( sb_N , N, rowNames, colNames );
	
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
	sb_matrix m = sb_convertFromDoubleMatrix(*matrix, rowNames, colNames);
	
	//cleanup
	//delete matrix;
	sb_deleteMatrix(sb_N);
	
	return m;
}

list<string> sb_splitString(const string& seq, const string& _1cdelim)
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

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT sb_matrix sGetNumberOfReactions (copasiModel model)
{
	  return sb_createMatrix (0, 0);
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT sb_strings sGetReactionNames (copasiModel model)
{
	return sb_createStringsArray(0);
}


// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT int sGetReactionRate(copasiModel model, int index)
{
	return 0;
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT sb_matrix sGetReactionRatesEx(sb_matrix values)
{
  return sb_createMatrix (0, 0);	
}

SBWAPIEXPORT sb_matrix sGetReactionRates(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel) return sb_createMatrix(0,0);

	const CCopasiVectorNS < CReaction > & reactions = pModel->getReactions();

	sb_matrix res  = sb_createMatrix(1, reactions.size());

	for (int i=0; i < reactions.size(); ++i)
		if (reactions[i])
		{
			sb_setColumnName(res, i, reactions[i]->getObjectName().c_str());
			sb_setMatrixValue(res, 0, i, reactions[i]->calculateFlux());
		}
	
	return res;
}

sb_matrix sGetFloatingSpecies(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel) return sb_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	int n = 0;
	for (int i=0; i < species.size(); ++i)
		if (species[i] && 
			(species[i]->getStatus() == CModelEntity::ODE || species[i]->getStatus() == CModelEntity::REACTIONS))
			++n;

	sb_matrix res  = sb_createMatrix(n,1);
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && 
			(species[i]->getStatus() == CModelEntity::ODE || species[i]->getStatus() == CModelEntity::REACTIONS))
		{
			sb_setMatrixValue(res, j, 0, species[i]->getConcentration());
			sb_setRowName(res, j, species[i]->getObjectName().c_str());
			++j;
		}

	return res;
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT sb_strings sGetFloatingSpeciesNames(copasiModel model)
{
	return sb_createStringsArray(0);
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT sb_strings sGetBoundarySpeciesNames(copasiModel model)
{
	return sb_createStringsArray(0);
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT int sGetNumberFloatingSpecies(copasiModel model)
{
	return 0;
}


// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT int sGetNumberOfBoundarySpecies(copasiModel model)
{
  return 0;
}


SBWAPIEXPORT sb_matrix sGetBoundarySpecies(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel) return sb_createMatrix(0,0);
}


// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT void sSetFloatingSpeciesByIndex (copasiModel model, int index)
{
}


// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT double sGetFloatingSpeciesByIndex (copasiModel model, int index)
{
	/*int n = 0;
	for (int i=0; i < species.size(); ++i)
		if (species[i] && 
			(species[i]->getStatus() == CModelEntity::FIXED))
			++n;

	sb_matrix res  = sb_createMatrix(n,1);

	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && 
			(species[i]->getStatus() == CModelEntity::FIXED))
		{
			sb_setMatrixValue(res, j, 0, species[i]->getConcentration());
			sb_setRowName(res, j, species[i]->getObjectName().c_str());
			++j;
		}

  return res;*/
  return 0.0;
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT int sGetNumberOfSpecies(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel) return 0;
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT sb_matrix sGetFloatingSpeciesConcentrations (copasiModel model)
{
	return sb_createMatrix (0, 0);
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT sb_matrix sGetFloatingSpeciesIntitialConcentrations (copasiModel model)
{
  CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	return sb_createMatrix (0, 0);
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT void sSetBoundarySpeciesByIndex (copasiModel model, int index)
{
}


// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT double sGetBoundarySpeciesByIndex (copasiModel model, int index)
{
	return 0.0;
}


SBWAPIEXPORT void sSetFloatingSpeciesIntitialConcentrations (copasiModel model, sb_matrix sp)
{
	if (sp.rows > sp.cols)  //row vector or column vector (lets allow both)
	{
		for (int i=0; i < sp.rows; ++i)
			sSetInitialConcentration(model, sb_getRowName(sp,i), sb_getMatrixValue(sp, i, 0));
	}
	else
	{
		for (int i=0; i < sp.cols; ++i)
			sSetInitialConcentration(model, sb_getColumnName(sp,i), sb_getMatrixValue(sp, 0, i));
	} 
}

SBWAPIEXPORT sb_matrix sGetBoundarySpeciesConcentrations (copasiModel model)
{
	return sb_createMatrix (0, 0);
}


sb_matrix sGetConcentrations(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel) return sb_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	sb_matrix res  = sb_createMatrix(1, species.size());

	for (int i=0; i < species.size(); ++i)
		if (species[i])
		{
			sb_setColumnName(res, i, species[i]->getObjectName().c_str());
			sb_setMatrixValue(res, 0, i, species[i]->getConcentration());
		}
	
	return res;
}

sb_matrix sGetCompartments(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel) return sb_createMatrix(0,0);

	const CCopasiVectorNS< CCompartment > & compartments = pModel->getCompartments();

	sb_matrix res  = sb_createMatrix(1, compartments.size());

	for (int i=0; i < compartments.size(); ++i)
		if (compartments[i])
		{
			sb_setColumnName(res, i, compartments[i]->getObjectName().c_str());
			sb_setMatrixValue(res, 0, i, compartments[i]->getValue());
		}
	
	return res;
}

sb_matrix sGetAllSpecies(copasiModel model)
{
	return sGetConcentrations(model);
}

sb_matrix sGetAmounts(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel) return sb_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	sb_matrix res  = sb_createMatrix(1, species.size());

	for (int i=0; i < species.size(); ++i)
		if (species[i] && species[i]->getCompartment())
		{
			sb_setColumnName(res, i, species[i]->getObjectName().c_str());
			sb_setMatrixValue(res, 0, i, 
					CMetab::convertToNumber( species[i]->getConcentration(), *species[i]->getCompartment(), pModel ));
		}

	return res;
}

double sGetConcentration(copasiModel model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	
	if (!pModel || !sb_contains(hash, s)) return -1.0;
	sb_CopasiPtr p = sb_getHashValue(hash, s);

	if (!p.species || !p.species->getCompartment()) return -1.0;	

	return p.species->getConcentration();
}


double sGetAmount(copasiModel model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s(name);
	
	if (!pModel || !sb_contains(hash, s)) return -1.0;
	sb_CopasiPtr p = sb_getHashValue(hash, s);

	if (!p.species || !p.species->getCompartment()) return -1.0;	

	return CMetab::convertToNumber( p.species->getConcentration(), *p.species->getCompartment(), pModel );
}


sb_matrix sGetRatesOfChange(copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel) return sb_createMatrix(0,0);

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	int n = pModel->getState().getNumVariable();

	if (n != species.size()) return sb_createMatrix(0,0);

	sb_matrix res  = sb_createMatrix(1, n);

	for (int i=0; i < species.size(); ++i)
		if (species[i] && species[i]->getCompartment())
		{
			sb_setColumnName(res, i, species[i]->getObjectName().c_str());
		}
	
	pModel->calculateDerivatives(res.values);
	return res;
}

SBWAPIEXPORT double sGetRateOfChange(copasiModel, int index)
{
	return 0.0;
}

sb_matrix sGetRatesOfChangeEx(copasiModel model, sb_matrix sp)
{
	return sb_createMatrix (0, 0);
}


SBWAPIEXPORT sb_strings sGetRatesOfChangeNames(copasiModel model)
{
	return sb_createStringsArray(0);
}

double sGetFlux(copasiModel model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s (name);	

#ifdef _WIN32
    double NaN = std::numeric_limits<double>::quiet_NaN(); 
#else
	double NaN = 0.0/0.0;

#endif

	if (!pModel || !sb_contains(hash, s)) return NaN;

	sb_CopasiPtr p = sb_getHashValue(hash, s);

	if (!p.reaction) return NaN;

	return p.reaction->calculateFlux();
}

double sGetParticleFlux(copasiModel model, const char * name)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	string s (name);	

#ifdef _WIN32
    double NaN = std::numeric_limits<double>::quiet_NaN(); 
#else
	double NaN = 0.0/0.0;
#endif

	if (!pModel || !sb_contains(hash, s)) return NaN;

	sb_CopasiPtr p = sb_getHashValue(hash, s);

	if (!p.reaction) return NaN;

	return p.reaction->calculateParticleFlux();
}


// ------------------------------------------------------------------
// Parameter Group
// ------------------------------------------------------------------

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT int sGetNumberOfGlobalParameters (copasiModel model)
{
  return 0;
}

SBWAPIEXPORT sb_strings sGetGlobalParameterNames (copasiModel model)
{
	return sb_createStringsArray(0);
}

SBWAPIEXPORT sb_matrix sGetGlobalParameters (copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);

	if (!pModel || !hash) return sb_createMatrix(0,0);

	list<string> names;
	list<double> values;

	for (CCMap::iterator i = hash->begin(); i != hash->end(); i++)
		if ((*i).second.param)
		{
			names.push_back((*i).first);
			values.push_back((*i).second.param->getValue());
		}
	sb_matrix params = sb_createMatrix(names.size(),1);

	int j=0;
	list<string>::iterator i1 = names.begin();
	list<double>::iterator i2 = values.begin();
	for (; i1 != names.end() && i2 != values.end(); i1++, i2++, ++j)
	{
		sb_setRowName(params, j, (*i1).c_str());
		sb_setMatrixValue(params, j, 0, (*i2));
	}
	return params;
}


SBWAPIEXPORT void sSetValues (copasiModel model, sb_matrix gp)
{
	if (gp.rows > gp.cols)  //row vector or column vector (lets allow both)
	{
		for (int i=0; i < gp.rows; ++i)
			sSetValue(model, sb_getRowName(gp,i), sb_getMatrixValue(gp, i, 0));
	}
	else
	{
		for (int i=0; i < gp.cols; ++i)
			sSetValue(model, sb_getColumnName(gp,i), sb_getMatrixValue(gp, 0, i));
	}
}

SBWAPIEXPORT void sSetGlobalParameterValues (copasiModel model, sb_matrix gp)
{
	sSetValues(model, gp);
}


// ------------------------------------------------------------------
// Compartment Group
// ------------------------------------------------------------------

SBWAPIEXPORT int sSetNumberOfCompartments (copasiModel model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel) return 0;
	
	CCopasiVectorNS < CCompartment > & compartments = pModel->getCompartments();
	return compartments.size();

}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT sb_strings sGetCompartmentNames (copasiModel model)
{
  return sb_createStringsArray(0);
}


// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT double sGetCompartmentByIndex (copasiModel model, int index)
{
	 return 0.0;
}

// STUB: NEEDS TO BE IMPLEMENTED
SBWAPIEXPORT void sSetCompartmentByIndex (copasiModel model, int index, double value)
{
	
}

SBWAPIEXPORT void sSetCompartmentVolumes (copasiModel model, sb_matrix v)
{
	sSetValues(model, v);
}
