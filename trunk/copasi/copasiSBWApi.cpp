#include "copasiSBWapi.h"
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
//#define COPASI_MAIN 1
#include "copasi/copasi.h"
#include "copasi/report/CCopasiRootContainer.h"
#include "copasi/CopasiDataModel/CCopasiDataModel.h"
#include "copasi/model/CModel.h"
#include "copasi/model/CCompartment.h"
#include "copasi/model/CMetab.h"
#include "copasi/model/CReaction.h"
#include "copasi/model/CModelValue.h"

char** returnColNamesOnly(c_matrix);

char** returnRowNamesOnly(c_matrix m)
{
	if (m.rows == 1 && m.cols > 1) 
		return returnColNamesOnly(m);

	char ** names = m.rownames.strings;
	m.rownames.strings = 0;
	c_deleteMatrix(m);
	return names;
}

char** returnColNamesOnly(c_matrix m)
{
	if (m.cols == 1 && m.rows > 1) 
		return returnRowNamesOnly(m);

	char ** names = m.colnames.strings;
	m.colnames.strings = 0;
	c_deleteMatrix(m);
	return names;
}

double* returnValuesOnly(c_matrix m)
{
	double* values = m.values;
	m.values = 0;
	c_deleteMatrix(m);
	return values;
}


/*! 
 \brief Read a model from an SBML file
 \param char* file name
 \return copasi_model Copasi model of the SBML file
 \ingroup loadsave
*/
copasi_model readSBMLFile(const char * filename)
{
	return cReadSBMLFile (filename);
}



/*! 
 \brief sGet the number of reactions in the model
 \param copasi_model model
 \return int Returns the number of reactions in the model
 \ingroup reaction
*/
int getNumberOfReactions (copasi_model model)
{
	return cGetNumberOfReactions(model);
}


/*! 
 \brief sGet the list of reaction names
 \param copasi_model model
 \return char** array of char * and length n, where n = number of reactions
 \ingroup reaction
*/
char** getReactionNames (copasi_model model)
{
	c_matrix m = cGetReactionRates(model);
	return returnRowNamesOnly(m);
}


/*! 
 \brief sGet the reaction rate for the ith reaction
 \param copasi_model model
 \param int reactionId
 \return double reaction rate for ith reaction
 \ingroup reaction
*/
double getReactionRate(copasi_model model, int i)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel) return 0.0;

	const CCopasiVectorNS < CReaction > & reactions = pModel->getReactions();
	if (i < reactions.size() && reactions[i])
		return reactions[i]->calculateFlux();

	return 0.0;
}


/*! 
 \brief Returns the vector of current reaction rates
 \param copasi_model model
 \return double array of reaction rates
 \ingroup reaction
*/
double* getReactionRates(copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return 0;

	const CCopasiVectorNS < CReaction > & reactions = pModel->getReactions();

	double * rates = new double[reactions.size()];

	for (int i=0; i < reactions.size(); ++i)
		if (reactions[i])
		{
			rates[i] = reactions[i]->calculateFlux();
		}
	
	return rates;
}


/*! 
 \brief Returns the rates of change given an array of new floating species concentrations
 \param copasi_model model
 \param double array of floating concentrations
 \return double vector of reaction rates
 \ingroup reaction
*/
double* getReactionRatesEx(copasi_model model, double * values)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return 0;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0; i < species.size(); ++i)
		species[i]->setConcentration(values[i]);

	return getReactionRates(model);
}

/*! 
 \brief sGet the number of boundary species
 \param copasi_model model
 \return number of species
 \ingroup boundary
*/
int getNumberOfBoundarySpecies(copasi_model model)
{
	return cGetNumberOfBoundarySpecies(model);
}


/*! 
 \brief sGet a list of boundary species names
 \param copasi_model model
 \return char** array of char * and length n, where n = number of species
 \ingroup boundary
*/
char** getBoundarySpeciesNames(copasi_model model)
{
	c_matrix m = cGetBoundarySpecies(model);
	return returnColNamesOnly(m);
}

/*! 
 \brief Set a boundary species concentration by index
 \param copasi_model model
 \param int index ith boundary species
 \ingroup boundary
*/
void setBoundarySpeciesByIndex (copasi_model model, int index, double value)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::FIXED)
		{
			if (j == index)
			{
				species[i]->setConcentration(value);
				break;
			}
			++j;
		}
}

/*! 
 \brief Set all the boundary species concentration
 \param copasi_model model
 \param sb_matric Vector of boundary species concentrations
 \ingroup boundary
*/
void setBoundarySpeciesConcentrations (copasi_model model, double * d)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::FIXED)
		{
			species[i]->setConcentration(d[j]);
			++j;
		}
}

/*! 
 \brief get all the boundary species concentration
 \param copasi_model model
 \param double * boundary species concentrations
 \ingroup boundary
*/
double * getBoundarySpeciesConcentrations (copasi_model model)
{
	c_matrix m = cGetBoundarySpecies(model);
	return returnValuesOnly(m);
}


/*! 
 \brief sGet a boundary species concentration by index - CURRENTLY NOT IMPLEMENTED
 \param copasi_model model
 \param int index ith boundary species
 \return double concentration of ith boundary species
 \ingroup state
*/
double getBoundarySpeciesByIndex (copasi_model model, int index)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return 0.0;

	double d = 0.0;
	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::FIXED)
		{
			if (j == index)
			{
				d = species[i]->getConcentration();
				break;
			}
			++j;
		}
	return d;
}


/*! 
 \brief sGet the number of floating species
 \param copasi_model model
 \return number of species
 \ingroup floating
*/
int getNumberFloatingSpecies(copasi_model model)
{
	return cGetNumberOfFloatingSpecies(model);
}

/*! 
 \brief sGet a list the floating species names
 \param copasi_model model
 \return c_strings array of char * and length n, where n = number of species
 \ingroup state
*/
char** getFloatingSpeciesNames(copasi_model model)
{
	c_matrix m = cGetFloatingSpeciesConcentrations(model);
	return returnRowNamesOnly(m);
}


/*! 
 \brief Set a floating species concentration by index
 \param copasi_model model
 \param int index ith floating species
 \ingroup state
*/
void setFloatingSpeciesByIndex (copasi_model model, int index, double value)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	if (!pModel) return;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::REACTIONS)
		{
			if (j == index)
			{
				species[i]->setConcentration(value);
				break;
			}
			++j;
		}
}

/*! 
 \brief sGet a floating species concentration by index
 \param copasi_model model
 \param int index ith floating species
 \return double Concentration of ith floating species
 \ingroup state
*/
double getFloatingSpeciesByIndex (copasi_model model, int index)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	if (!pModel) return 0.0;

	double d = 0.0;
	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::REACTIONS)
		{
			if (j == index)
			{
				d = species[i]->getConcentration();
				break;
			}
			++j;
		}
	return d;
}

/*! 
 \brief Set all the floating species concentration
 \param copasi_model model
 \param sb_matric Vector of floating species concentrations
 \ingroup boundary
*/
void setFloatingSpeciesConcentrations (copasi_model model, double * d)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::REACTIONS)
		{
			species[i]->setConcentration(d[j]);
			++j;
		}
}

/*! 
 \brief sGet the initial flreturn NUoating species concentrations  - CURRENTLY NOT IMPLEMENTED
 \param copasi_model model
 \return double * vector of initial floating species concentrations
 \ingroup floating
*/
double* getFloatingSpeciesIntitialConcentrations (copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return 0;

	int n = cGetNumberOfFloatingSpecies(model);
	double * d = new double[n];

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::REACTIONS)
		{
			d[j] = species[i]->getInitialConcentration();
			++j;
		}
	return d;
}


/*! 
 \brief Set the initial floating species concentrations
 \param copasi_model model
 \param c_matrix Vector of initial floating species concentrations
 \ingroup floating
*/
void setFloatingSpeciesIntitialConcentrations (copasi_model model, double * d)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::REACTIONS)
		{
			species[i]->setInitialConcentration(d[j]);
			++j;
		}
}

/*! 
 \brief Set the initial floating species concentration of the ith species
 \param copasi_model model
 \param double value value to set the ith initial floating species concentration
 \ingroup floating
*/
void setFloatingSpeciesIntitialConcentrationByIndex (copasi_model model, int index, double d)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0, j=0; i < species.size(); ++i)
		if (species[i] && species[i]->getStatus() == CModelEntity::REACTIONS)
		{
			if (j == index)
			{
				species[i]->setInitialConcentration(d);
				break;
			}
			++j;
		}
}
/*! 
 \brief sGet the number of global parameters
 \param copasi_model model
 \return int numberOfGlobalParameters. Returns the number of gloabal parameters in the model
 \ingroup parameter
*/
int getNumberOfGlobalParameters (copasi_model model)
{
	return cGetNumberOfGlobalParameters(model);
}

/*! 
 \brief sGet the list of global parameter names
 \param copasi_model model 
 \return char** array of char * and length n, where n = number of global parameters
 \ingroup parameter
*/
char** getGlobalParameterNames (copasi_model model)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return 0;

	c_matrix m = cGetGlobalParameters(model);
	return returnRowNamesOnly(m);	
}

/*! 
 \brief sGet the value of a global parameter by index
 \param copasi_model model 
 \param int index ith global parameter
 \return double returned value of parameter
 \ingroup parameter
*/
double getGlobalParameterByIndex (copasi_model model, int index)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return 0.0;

	double d = 0.0;
	CCopasiVectorN< CModelValue > & params = pModel->getModelValues();
	if (index < params.size() && params[index])
		d = params[index]->getValue();

	return d;
}

/*! 
 \brief Set the value of a global parameter by index
 \param copasi_model model 
 \param int index ith global parameter
 \param double Value to set parameter to
 \ingroup parameter
*/
void setGlobalParameterByIndex (copasi_model model, int index, double value)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return;

	CCopasiVectorN< CModelValue > & params = pModel->getModelValues();
	if (index < params.size() && params[index])
		params[index]->setValue(value);
}

/*! 
 \brief set the vector of global parameters
 \param copasi_model model 
 \paramn c_matrix a vector containing names and values for the global parameters. names must be present
 \ingroup parameter
*/
void setGlobalParameterValues (copasi_model model, double * values)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);	
	if (!pModel) return;

	CCopasiVectorN< CModelValue > & params = pModel->getModelValues();
	for (int i=0; i < params.size(); ++i)
		params[i]->setValue(values[i]);
}

// -----------------------------------------------------------------------
/** \} */
/**
  * @name SBW compartment group
  */
/** \{ */

/*! 
 \brief sGet the number of compartments
 \param copasi_model model
 \return int numberOfCompartments
 \ingroup compartment
*/
int getNumberOfCompartments (copasi_model model)
{
	return cGetNumberOfCompartments(model);
}

/*! 
 \brief sGet the list of compartment names
 \param copasi_model model
 \return char** compartment names
 \ingroup compartment
*/
char** getCompartmentNames (copasi_model model)
{
	c_matrix m = cGetCompartments(model);
	return returnRowNamesOnly(m);
}


/*! 
 \brief sGet the compartment volume by index
 \param copasi_model model
 \param int index ith compartment
 \return double Voluem of compartment 
 \ingroup compartment
*/
double getCompartmentByIndex (copasi_model model, int index)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel) return 0;
	double d = 0.0;
	const CCopasiVectorNS< CCompartment > & compartments = pModel->getCompartments();
	if (index < compartments.size() && compartments[index])
		d = compartments[index]->getValue();
	return d;
}


/*! 
 \brief Set a compartment volume by index
 \param copasi_model model
 \param int index ith compartment
 \param double volume Volume of ith compartment
 \ingroup compartment
*/
void setCompartmentByIndex (copasi_model model, int index, double value)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel) return;

	const CCopasiVectorNS< CCompartment > & compartments = pModel->getCompartments();
	if (index < compartments.size() && compartments[index])
		compartments[index]->setValue(value);
}


/*! 
 \brief Set a compartment volumes using a vector of compartment values
 \param copasi_model model
 \param double volume Vector of compartment volumes
 \ingroup compartment
*/
void setCompartmentVolumes (copasi_model model, double * v)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel) return;

	const CCopasiVectorNS< CCompartment > & compartments = pModel->getCompartments();

	for (int i=0; i < compartments.size(); ++i)
		compartments[i]->setValue(v[i]);
}


// -----------------------------------------------------------------------
/** \} */
/**
  * @name SBW rates of change group
  */
/** \{ */

/*! 
 \brief Compute the current rates of change for all species
 \param copasi_model model
 \return double* array of size = number of species
 \ingroup rateOfChange
*/
double* getRatesOfChange(copasi_model model)
{
	c_matrix m = cGetRatesOfChange(model);
	return returnValuesOnly(m);
}

/*! 
 \brief Compute the current rates of change for the ith species
 \param copasi_model model
 \param ith rate of change to compute
 \return double ith rate of change
 \ingroup rateOfChange
*/
double getRateOfChange(copasi_model model, int index)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	
	if (!pModel)
		return 0.0;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();

	if (species.size() <= index)
		return 0.0;

	species[index]->refreshRate();
	return species[index]->getRate();
}

/*! 
 \brief Returns the names used to represent the rates of change
 \param copasi_model model
 \return char** list of names used to represent the rate of change
 \ingroup rateOfChange
*/
char** getRatesOfChangeNames(copasi_model model)
{
	c_matrix m = cGetRatesOfChange(model);
	return returnColNamesOnly(m);
}

/*! 
 \brief Returns the rates of change given a vector of floating species concentrations
 \param copasi_model model
 \param double* array of floating species concentrations
 \return double* array of rates of change
 \ingroup sbw_rateOfChange
*/
double* getRatesOfChangeEx(copasi_model model, double* sp)
{
	setFloatingSpeciesConcentrations(model, sp);
	return getRatesOfChange(model);
}

// -----------------------------------------------------------------------
/** \} */
/**
  * @name SBW time course simulation
  */
/** \{ */

/*! 
 \brief Carry out a simulation using the Copasi LSODA method
 \param copasiModel model
 \param double start time
 \param double end time
 \param int number of steps in the output
 \return c_matrix matrix of concentration or particles
 \ingroup sbw_rateOfChange
*/
c_matrix simulate(copasi_model model, double timeStart, double timeEnd, int numOfPoints)
{
	return cSimulateDeterministic(model, timeStart, timeEnd, numOfPoints);
}


/*! 
 \brief Simulate the differential equation model over one time step
 \param copasi_model model
 \param double time step
 \return double new time, i.e (current time + timeStep)
 \ingroup simulation
*/
double oneStep(copasi_model model, double timeStep)
{
   return cOneStep (model, timeStep);
}

