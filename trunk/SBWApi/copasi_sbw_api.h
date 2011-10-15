/*! \mainpage SBW API for the COPASI Library
 *
 * \section intro_sec Introduction
 *
 * The developers of COPASI provide COPASI as a reusable library as well as
 * the well known COPASI user interface. The library however has a fairly
 * complex API and can take some time getting used. We have therefore layered 
 * on top of the COPASI library a new C based API that we feel is much simpler
 * to use. For example, to run a simple SBML model and generate time series data
 * we would call:
 *  
 \code
 copasi_model m;
 sb_matrix output;
  
 m = readSBMLFile ("mymodel.xml");
  
 output = simulationDeterministic (m, 0, 10, 100); 
 \endcode
 
 More complex example:
 
 \code
 #include <stdlib.h>
 #include <stdio.h>
 #include "copasi_api.h"

 int main(int nargs, char** argv)
 {
        sb_matrix efm, output, params;
        copasi_model m1, m2;
        
        if (nargs < 2)
        {
            m1 = model1();
        }
        else
        {
            printf("loading model file %s\n", argv[1]);
            m1 = readSBMLFile(argv[1]);        
        }

        writeAntimonyFile(m1, "model.txt");
        
        printf("Antimony file written to model.txt\nSimulating...\n");  
        output = simulateDeterministic(m1, 0, 100, 1000);  //model, start, end, num. points
        printf("output.tab has %i rows and %i columns\n",output.rows, output.cols);
        sb_printMatrixToFile("output.tab", output);
        sb_deleteMatrix(output);
                 
        removeModel(m1);
        copasiEnd();
        return 0;
 }
 \endcode
 * \section install_sec Installation
 *
 * Installation documentation is provided in the main google code page.

 \defgroup loadsave Read and Write models
 \brief Read and write models to files or strings. Support for SBML and Antimony formats.

 \defgroup create Define models
 \brief Create models and set model components using code

 \defgroup state Current state of system
 \brief Compute derivatives, fluxed, and other values of the system at the current state

 \defgroup reaction Reaction group
 \brief Get information about reaction rates
 
 \defgroup rateOfChange Rates of change group
 \brief Get information about rates of change

 \defgroup boundary Boundary species group
 \brief Get information about reaction rates
 
 \defgroup floating Floating species group
 \brief Get information about reaction rates
  
 \defgroup parameters Parameter group
 \brief set and get global and local parameters
 
 \defgroup compartment Compartment group
 \brief set and get information on compartments
 
 \defgroup simulation Time-course simulation
 \brief Deterministic, stochastic, and hybrid simulation algorithms

 \defgroup mca Metabolic Control Analysis
 \brief Calculate control coefficients and sensitivities

 \defgroup matrix Stoichiometry analysis
 \brief Linear algebra based methods for analyzing a reaction network

 \defgroup optim Parameter optimization
 \brief Optimization of parameters to match given data
*/

#ifndef COPASI_SB_C_API
#define COPASI_SB_C_API
 /**
  * @file    copasiSBApi.h
  * @brief   SBW C API for the Copasi C++ library

This is a C API for the COPASI C++ library. Rate equations in COPASI require the "complete name",   
e.g. instead of X, the rate must specify <model.compartment.X>. In this C API, those complete names
are stored in a hash table. The API replaces the simple strings, i.e. "C", with the complete names by
using the hash-table. This is mainly for speed; otherwise, every cSetReactionRate would be searching
through the entire model for each of its variables. The hash-table idea is used for functions such
as cSetValue, which can set the value of a parameter or that of a molecular species. Again, it uses the
hash table to identify what a variable is. 

The C API hides the C++ classes by casting some of the main classes into void pointers inside
C structs. 

std::map is used for performing the hashing (it is not a real hash-table, but close enough).
boost::regex is used for string substitutions.
*/

#include "copasi_api.h"

BEGIN_C_DECLS

// -----------------------------------------------------------------------
/** \} */
/**
  * @name Read and write models
  */
/** \{ */

/*! 
 \brief Create a model from an Antimony, see antimony.sf.net for details of Antimony syntax
 \param char* file name
 \return copasi_model Copasi model of the Antimony file
 \ingroup loadsave
*/
COPASIAPIEXPORT copasi_model readAntimonyFile(const char * filename)
{
	return cReadAntimonyFile(filename);
}


/*! 
 \brief Create a model from an Antimony string
 \param char* Antimony string
 \return copasi_model Copasi model of the Antimony string
 \ingroup loadsave
*/
COPASIAPIEXPORT copasi_model readAntimonyString(const char * antimony);
{
	return cReadAntimonyString(antimony);
}

/*! 
 \brief Create a model from an SBML file
 \param char* file name
 \return copasi_model Copasi model of the SBML file
 \ingroup loadsave
*/
COPASIAPIEXPORT copasi_model readSBMLFile(const char * filename)
{
	return cReadSBMLFile(filename);
}


/*! 
 \brief Create a model from an SBML string
 \param char* SBML string
 \return copasi_model Copasi model of the SBML string
 \ingroup loadsave
*/
COPASIAPIEXPORT copasi_model readSBMLString(const char * sbml)
{
	return cReadSBMLString(filename);
}


/*! 
 \brief Save a model as an SBML file
 \param copasi_model copasi model
 \param char* file name
 \ingroup loadsave
*/
COPASIAPIEXPORT void writeSBMLFile(copasi_model model, const char * filename)
{
	return cWriteSBMLFile(model, filename);
}


/*! 
 \brief Save a model as an Antimony file, see antimony.sf.net for details of Antimony syntax
 \param copasi_model copasi model
 \param char* file name
 \ingroup loadsave
*/
COPASIAPIEXPORT void writeAntimonyFile(copasi_model model, const char * filename)
{
	return cReadAntimonyFile(model, filename);
}


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Create model group
  */
/** \{ */

/*! 
 \brief Create a model
 \param char* model name
 \return copasi_model a new copasi model
 \ingroup create
*/
COPASIAPIEXPORT copasi_model createModel(const char * name)
{
	return cCreateModel(name);
}

/*! 
 \brief This function MUST be called after creating a model or modifying a model (except parameter changes)
             This function was called internally inside every analysis function, but that was inefficient, so it must be
			 called manually. 
			 Note that when models are generated from a file or string (e.g. sbml), they do not need to be compiled again.
 \param copasi_model model
 \ingroup create
*/
COPASIAPIEXPORT void compileModel(copasi_model model)
{
	return cCompileModel(filename);
}

/*! 
 \brief Create compartment
 \param char* compartment name
 \param double volume
 \return copasi_compartment a new compartment
 \ingroup create
*/
COPASIAPIEXPORT copasi_compartment createCompartment(copasi_model model, const char* name, double volume)
{
	return cCreateCompartment(model, name, volume);
}


/*! 
 \brief Set a volume of compartment
 \param copasi_model model
 \param char * compartment name
 \param double volume
 \ingroup create
*/
COPASIAPIEXPORT void setVolume(copasi_model model, const char * compartment, double volume)
{
	return cSetValue(model, compartment, volume);
}


/*! 
 \brief Set the concentration of a species, volume of a compartment, or value of a parameter
      The function will figure out which using the name (fast lookup using hashtables).
      If the name does not exist in the model, a new global parameter will be created.
 \param copasi_model model
 \param char * name
 \param double value
 \return 0 if new variable was created. 1 if existing variable was found
 \ingroup create
*/
COPASIAPIEXPORT int setValue(copasi_model model model, const char * name, double value)
{
	return cSetValue(model, name, value);
}


/*! 
 \brief Add a species to the model
 \param copasi_compartment model
 \param char* species name
 \param double initial value (concentration or count, depending on the model)
 \ingroup create
*/
COPASIAPIEXPORT void createSpecies(copasi_compartment compartment, const char* name, double initialValue)
{
	return cCreateSpecies(compartment, name, initialValue);
}


/*! 
 \brief Set a species as boundary or floating (will remove any assignment rules)
 \param copasi_model model
  \param char * name
 \param int boundary = 1, floating = 0 (default)
 \ingroup create
*/
COPASIAPIEXPORT void setSpeciesType(copasi_model model, const char * species, int isBoundary)
{
	return cSetSpeciesType(model, species, isBoundary);
}


/*! 
 \brief Set a species concentration
 \param copasi_model model
 \param char * species name
 \param double concentration
 \ingroup create
*/
COPASIAPIEXPORT void setConcentration(copasi_model model, const char * species, double conc)
{
	return cSetValue(model, name, initialValue);
}

/*! 
 \brief Set a species amounts
 \param copasi_model model
 \param char * species name
 \param double amount
 \ingroup create
*/
COPASIAPIEXPORT void setSpeciesAmount(copasi_model model, const char * species, double amount)
{
	return cSetSpeciesAmount(model, species, amount);
}


/*! 
 \brief Set the assignment rule for a species (automatically assumes boundary species)
 \param copasi_model model
 \param char * species name
 \param char* formula, use 0 to remove assignment rule
 \return int 0=failed 1=success
 
 \code
 result = cSetAssignmentRule (m, "S1", "sin (time*k1)");
 \endcode
 \ingroup create
*/
COPASIAPIEXPORT int setAssignmentRule(copasi_model model, const char * species, const char * formula)
{
	return cSetAssignmentRule(model, species, formula);
}


/*! 
 \brief Set the value of an existing global parameter or create a new global parameter
 \param copasi_model model
 \param char* parameter name
 \param double value
  \return int 0=new value created 1=found existing value
 \ingroup create
*/
COPASIAPIEXPORT int setGlobalParameter(copasi_model model, const char * name, double value)
{
	return cSetGlobalParameter(model, name, value);
}


/*! 
 \brief Create a new variable that is defined by a formula
 \param copasi_model model
 \param char* name of new variable
 \param char* formula
 \return int 0=failed 1=success
 \ingroup create
*/
COPASIAPIEXPORT int createVariable(copasi_model model, const char * name, const char * formula)
{
	return cCreateVariable(model, name, value);
}


/*! 
 \brief Add a trigger and a response, where the response is defined by a target variable and an assignment formula
 \param copasi_model model
 \param char * event name
 \param char * trigger formula
 \param char * response: name of variable or species
 \param char* response: assignment formula
 \return int 0=failed 1=success
 
 Example Usage. The following code will create an event where the parameter k1 is halved when time > 10.
 \code
 result = cCreateEvent (m, "myEvent", "time > 10", "k1", "k1/2");
 \endcode
 \ingroup create
*/
COPASIAPIEXPORT int createEvent(copasi_model model, const char * name, const char * trigger, const char * variable, const char * formula)
{
	return cCreateEvent(model, name, trigger, variable, formula);
}


/*!
 \brief Create a new reaction with a given name
 \param copasi_model model
 \param char* reaction name
 \return copasi_reaction a new reaction
 
 \code
 r = cCreateReaction (m, "J1")
 \endcode
 \ingroup create
*/
COPASIAPIEXPORT copasi_reaction createReaction(copasi_model model, const char* name)
{
	return cCreateReaction(model, name);
}


/*! 
 \brief Add a reactant to a reaction
 \param copasi_reaction reaction
 \param char * reactant
 \param double stoichiometry
 
 \code
 cCreateReaction (m, "S1", 1);
 \endcode
 \ingroup create
*/
COPASIAPIEXPORT void addReactant(copasi_reaction reaction, const char * species, double stoichiometry)
{
	return cAddReactant(reaction, species, stoichiometry);
}

/*! 
 \brief Add a product to a reaction
 \param copasi_reaction reaction
 \param char * product
 \param double stoichiometry
 
 Create a reaction J1: 2 A -> B + C
 \code
 r = cCreateReaction (m, "J1");
 cAddReactant (r, "A", 2);
 cAddProduct (r, "B", 1);
 cAddProduct (r, "C", 1);
 \endcode

 \ingroup create
*/
COPASIAPIEXPORT void addProduct(copasi_reaction reaction, const char * species, double stoichiometry)
{
	return cAddProduct(reaction, species, stoichiometry);
}

/*! 
 \brief Set reaction rate equation
 \param copasi_reaction reaction
 \param char* custom formula
 \return int success=1 failure=0
 
 \code
 int result;
 result = cSetReactionRate (r, "k1*S1");
 \endcode
 
 \ingroup create
*/
COPASIAPIEXPORT int setReactionRate(copasi_reaction reaction, const char * formula)
{
	return cSetReactionRate(reaction, formula);
}


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Reaction group
  */
/** \{ */

/*! 
 \brief Get the number of reactions in the model
 \param copasi_model model
 \return int Returns the number of reactions in the model
 \ingroup reaction
*/
COPASIAPIEXPORT sb_matrix getNumberOfReactions (copasi_model model)
{
	return cGetNumberOfReactions(model);
}


/*! 
 \brief Get the list of reaction names
 \param copasi_model model
 \return sb_string array of char * and length n, where n = number of reactions
 \ingroup reaction
*/
COPASIAPIEXPORT char** getReactionNames (copasi_model model)
{
	tc_matrix m = cGetReactionRates(model);
	char ** colnames = m.colnames.strings;
	m.colnames.strings = 0; //prevent this from being deleted
	tc_deleteMatrix(m);  //delete everything else in the matrix
	return colnames;
}

/*! 
 \brief Get the reaction rate for the ith reaction
 \param copasi_model model
 \param int reactionId
 \return double reaction rate for ith reaction
 \ingroup reaction
*/
COPASIAPIEXPORT double getReactionRate(copasi_model model, int index)
{
	tc_matrix m = cGetReactionRates(model);
	double rate = tc_getMatrixValue(m, 0, index);  //this will also check bounds of matrix
	tc_deleteMatrix(m); //delete matrix
	return rate;
}

/*! 
 \brief Returns the vector of current reaction rates
 \param copasi_modl model
 \return double array of reaction rates
 \ingroup reaction
*/
COPASIAPIEXPORT double* getReactionRates(copasi_model model)
{
	tc_matrix m = cGetReactionRates(model);
	double * rates = m.values;
	m.values = 0; //prevent this from being deleted
	tc_deleteMatrix(m);  //delete everything else in the matrix
	return rates;
}


/*! 
 \brief Returns the rates of change given an array of new floating species concentrations. 
            Program will crash if the argument is not exaclty the size of getNumberOfSpecies
 \param copasi_model model
 \param double[] Array of floating concentrations
 \return double array of reaction rates
 \ingroup reaction
*/
COPASIAPIEXPORT double[] getReactionRatesEx(copasi_model model, double[] conc)
{
	CModel* pModel = (CModel*)(model.CopasiModelPtr);
	CCMap * hash = (CCMap*)(model.qHash);
	
	if (!pModel) return 0;

	const CCopasiVector< CMetab > & species = pModel->getMetabolites();
	for (int i=0; i < species.size(); ++i)
		if (species[i])
			species[i]->setConcentration(conc[i]); //assume that conc is correct size

	return getReactionRates(model);
}

// -----------------------------------------------------------------------
/** \} */
/**
  * @name Boundary species group
  */
/** \{ */

/*! 
 \brief Get the number of boundary species
 \param copasi_model model
 \return number of species
 \ingroup boundary
*/
COPASIAPIEXPORT int getNumberOfBoundarySpecies(copasi_model model)
{
	return cGetNumberOfBoundarySpecies(model);
}

/*! 
 \brief Get a list of boundary species names
 \param copasi_model model
 \return sb_string array of char * and length n, where n = number of species
 \ingroup boundary
*/
COPASIAPIEXPORT char** getBoundarySpeciesNames(copasi_model model)
{
	tc_matrix m = cGetBoundarySpecies(model);
	char ** colnames = m.colnames.strings;
	m.colnames.strings = 0; //prevent this from being deleted
	tc_deleteMatrix(m);  //delete everything else in the matrix
	return colnames;
}


/*! 
 \brief Set a boundary species concentration by index
 \param copasi_model model
 \param int index ith boundary species
 \ingroup boundary
*/
COPASIAPIEXPORT void setBoundarySpeciesByIndex (copasi_model model, int index)
{
	tc_matrix m = cGetBoundarySpecies(model);
	const char * s = tc_getColumnName(m, index); //this will check bounds
	cSetSpeciesType( model, s, 1);  
	tc_deleteMatrix(m); //delete matrix
}

/*! 
 \brief Set all the boundary species concentration. 
            The argument array must be the same size as getNumberOfBoundarySpecies
 \param copasi_model model
 \param double * 
 \ingroup boundary
*/
COPASIAPIEXPORT void setBoundarySpeciesConcentrations (copasi_model model, double* d)
{
	tc_matrix m = cGetBoundarySpecies(model);
	double * old_d = m.values; //just change the pointers inside the matrix 
	m.values = d;
	cSetBoundarySpecies( model, m );  
	m.values = old_d;
	tc_deleteMatrix(m); //delete matrix
}

/*! 
 \brief Set all the boundary species concentration
 \param copasi_model model
 \param sb_matrix Vector of boundary species concentrations
 \ingroup boundary
*/
COPASIAPIEXPORT double* getBoundarySpeciesConcentrations (copasi_model model)
{
	tc_matrix m = cGetBoundarySpecies(model);
	double * d = m.values;
	m.values = 0;   //prevent this from being deleted
	tc_deleteMatrix(m); //delete matrix
	return d;
}


/*! 
 \brief Get a boundary species concentration by index
 \param copasi_model model
 \param int index ith boundary species
 \return double Concentration of ith boundary species
 \ingroup state
*/
COPASIAPIEXPORT double getBoundarySpeciesByIndex (copasi_model model, int index)
{
	tc_matrix m = cGetBoundarySpecies(model);
	double d = tc_getMatrixValue(m, 0, index); //will check bounds
	tc_deleteMatrix(m); //delete matrix
	return d;
}

// -----------------------------------------------------------------------
/** \} */
/**
  * @name Floating species group
  */
/** \{ */

/*! 
 \brief Get the number of floating species 
 \param copasi_model model
 \return number of species
 \ingroup floating
*/
COPASIAPIEXPORT int getNumberFloatingSpecies(copasi_model model)
{
	return cGetNumberOfFloatingSpecies(model);
}

/*! 
 \brief Get a list the floating species names
 \param copasi_model model
 \return char** array of char * and length n, where n = number of species
 \ingroup state
*/
COPASIAPIEXPORT char** getFloatingSpeciesNames(copasi_model model)
{
	tc_matrix m = cGetFloatingSpecies(model);
	char ** colnames = m.colnames.strings;
	m.colnames.strings = 0; //prevent this from being deleted
	tc_deleteMatrix(m);  //delete everything else in the matrix
	return colnames;
}


/*! 
 \brief Set a floating species concentration by index
 \param copasi_model model
 \param int index ith floating species
 \ingroup state
*/
COPASIAPIEXPORT void setFloatingSpeciesByIndex (copasi_model model, int index)
{
	tc_matrix m = cGetFloatingSpecies(model);
	const char * s = tc_getColumnName(m, index); //this will check bounds
	cSetSpeciesType( model, s, 1);  
	tc_deleteMatrix(m); //delete matrix
}

/*! 
 \brief Get a floating species concentration by index
 \param copasi_model model
 \param int index ith floating species
 \return double Concentration of ith floating species
 \ingroup state
*/
COPASIAPIEXPORT double getFloatingSpeciesByIndex (copasi_model model, int index)
{
	tc_matrix m = cGetFloatingSpecies(model);
	double rate = tc_getMatrixValue(m, 0, index);  //this will also check bounds of matrix
	tc_deleteMatrix(m); //delete matrix
	return rate;
}

/*! 
 \brief Set all the floating species concentration. 
            Argument array must be of same size as getNumberOfFloatingSpecies
 \param copasi_model model
 \param double* vector of floating species concentrations
 \ingroup boundary
*/
COPASIAPIEXPORT void setFloatingSpeciesConcentrations (copasi_model model, sb_matrix sp)
{
	tc_matrix m = cGetFloatingSpecies(model);
	double * old_d = m.values; //just change the pointers inside the matrix 
	m.values = d;
	cSetFloatingSpecies( model, m );  
	m.values = old_d;
	tc_deleteMatrix(m); //delete matrix
}

/*! 
 \brief Set all the floating species concentration 
 \param copasi_model model
 \param double * floating species concentrations
 \ingroup boundary
*/
COPASIAPIEXPORT double* getFloatingSpeciesConcentrations (copasi_model model)
{
	tc_matrix m = cGetFloatingSpecies(model);
	double * d = m.values;
	m.values = 0;   //prevent this from being deleted
	tc_deleteMatrix(m); //delete matrix
	return d;
}

/*! 
 \brief Get the initial floating species concentrations 
 \param copasi_model model
 \return double* initial floating species concentrations
 \ingroup floating
*/
COPASIAPIEXPORT double* getFloatingSpeciesIntitialConcentrations (copasi_model model)
{
	tc_matrix m = cGetFloatingSpeciesIntitialConcentrations(model);
	double * d = m.values;
	m.values = 0;   //prevent this from being deleted
	tc_deleteMatrix(m); //delete matrix
	return d;
}


/*! 
 \brief Set the initial floating species concentrations.
             Array of the same size as numberOfFloatingSpecies
 \param copasi_model model
 \param double * initial floating species concentrations
 \ingroup floating
*/
COPASIAPIEXPORT void setFloatingSpeciesIntitialConcentrations (copasi_model model, double * d)
{
	tc_matrix m = cGetFloatingSpecies(model);
	double * old_d = m.values; //just change the pointers inside the matrix 
	m.values = d;
	cSetFloatingSpeciesIntitialConcentrations( model, m );  
	m.values = old_d;
	tc_deleteMatrix(m); //delete matrix
}

/*! 
 \brief Set the initial floating species concentration of the ith species
 \param copasi_model model
 \param double value value to set the ith initial floating species concentration
 \ingroup floating
*/
COPASIAPIEXPORT void setFloatingSpeciesIntitialConcentrationByIndex (copasi_model model, int index, double sp)
{
	tc_matrix m = cGetFloatingSpecies(model);
	const char * s = tc_getColumnName(m, index);
	cSetFloatingSpeciesIntitialConcentrations( model, s, sp);  
	tc_deleteMatrix(m); //delete matrix
}


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Rates of change group
  */
/** \{ */

/*! 
 \brief Compute the current rates of change for all species
 \param copasi_model model
 \return sb_matrix matrix of with 1 row and n columns, where n = number of species
 \ingroup rateOfChange
*/
COPASIAPIEXPORT sb_matrix getRatesOfChange(copasi_model model);

/*! 
 \brief Compute the current rates of change for the ith species
 \param copasi_model model
 \param ith rate of change to compute
 \return double ith rate of change
 \ingroup rateOfChange
*/
COPASIAPIEXPORT double getRateOfChange(copasi_model model, int index)
{
	tc_matrix m = cGetRatesOfChange(model);
	double d = tc_getMatrixValue(m, index);
	tc_deleteMatrix(m);
	return d;
}

/*! 
 \brief Returns the names used to represent the rates of change
 \param copasi_model model
 \return sb_string List of names used to represent the rate of change
 \ingroup rateOfChange
*/
COPASIAPIEXPORT sb_matrix getRatesOfChangeNames(copasi_model model)
{
	tc_matrix m = cGetRatesOfChange(model);
	double d = tc_getMatrixValue(m, index);
	tc_deleteMatrix(m);
	return d;
}

/*! 
 \brief Returns the rates of change given a vector of floating species concentrations
 \param copasi_model model
 \param sb_matrix vector of floating species concentrations
 \return sb_matrix vector of rates of change
 \ingroup rateOfChange
*/
COPASIAPIEXPORT sb_matrix getRatesOfChangeEx(copasi_model model, sb_matrix sp);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Current state of system
  */
/** \{ */

/*! 
 \brief Get the current concentrations of all species
 \param copasi_model model
 \return sb_matrix matrix of with 1 row and n columns, where n = number of species
 The names of the species are included in the matrix column labels
 \ingroup state
*/
COPASIAPIEXPORT sb_matrix getConcentrations(copasi_model);

/*! 
 \brief Get the current amounts of all species. The amounts are calculated from the concentrations and compartment volume
 \param copasi_model model
 \return sb_matrix matrix of with 1 row and n columns, where n = number of species
 The names of the species are included in the matrix column labels
 \ingroup state
*/
COPASIAPIEXPORT sb_matrix getAmounts(copasi_model);

/*! 
 \brief Get a list of all species names, including floating and boundary species
 \param copasi_model model
 \return sb_string array of char * and length n, where n = number of species
 \ingroup state
*/
COPASIAPIEXPORT sb_strings getAllSpeciesNames(copasi_model model);


/*! 
 \brief Get the current concentration of a species
 \param copasi_model model
 \param string species name
 \return double concentration. -1 indicates that a species by this name was not found
 \ingroup state
*/
COPASIAPIEXPORT double getConcentration(copasi_model, const char * name);

/*! 
 \brief Get the current amount of a species. The amounts are calculated from the concentrations and compartment volume
 \param copasi_model model
 \param string species name
 \return double amount. -1 indicates that a species by this name was not found
 \ingroup state
*/
COPASIAPIEXPORT double getAmount(copasi_model, const char * name);


/*! 
 \brief Compute current flux through the given reactions
 \param copasi_model model
 \param string reaction name, e.g. "J1"
 \return double rate. If reaction by this name does not exist that NaN will be returned
  The names of the fluxes are included in the matrix column labels
 \ingroup state
*/
COPASIAPIEXPORT double getFlux(copasi_model, const char * name);

/*! 
 \brief Compute current flux through the given reactions in terms of particles
 \param copasi_model model
 \param string reaction name, e.g. "J1"
 \return double rate. If reaction by this name does not exist that NaN will be returned
 \ingroup state
*/
COPASIAPIEXPORT double getParticleFlux(copasi_model, const char * name);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Parameter group
  */
/** \{ */

/*! 
 \brief Get the number of global parameters
 \param copasi_model model
 \return int numberOfGlobalParameters. Returns the number of gloabal parameters in the model
 \ingroup parameter
*/
COPASIAPIEXPORT int getNumberOfGlobalParameters (copasi_model);

/*! 
 \brief Get the list of global parameter names
 \param copasi_model model 
 \return sb_string array of char * and length n, where n = number of global parameters
 \ingroup parameter
*/
COPASIAPIEXPORT sbc_string getGlobalParameterNames (copasi_model);

/*! 
 \brief Get the value of a global parameter by index
 \param copasi_model model 
 \param int index ith global parameter
 \return double returned value of parameter
 \ingroup parameter
*/
COPASIAPIEXPORT double getGlobalParameterByIndex (copasi_model, int);

/*! 
 \brief Set the value of a global parameter by index
 \param copasi_model model 
 \param int index ith global parameter
 \param double Value to set parameter to
 \ingroup parameter
*/
COPASIAPIEXPORT void setGlobalParameterByIndex (copasi_model, int, double);

/*! 
 \brief Get a list of the global parameters
 \param copasi_model model 
 \return sb_matrix A vector containing the values for the global parameters.
 \ingroup parameter
*/
COPASIAPIEXPORT sb_matrix getGlobalParameters (copasi_model);

/*! 
 \brief Set the vector of global parameters
 \param copasi_model model 
 \paramn sb_matrix A vector containing the values for the global parameters.
 \ingroup parameter
*/
COPASIAPIEXPORT void setGlobalParameterValues (copasi_model, sb_matrix gp);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Compartment group
  */
/** \{ */

/*! 
 \brief Get the number of compartments
 \param copasi_model model
 \return int numberOfCompartments
 \ingroup compartment
*/
COPASIAPIEXPORT int getNumberOfCompartments (copasi_model);

/*! 
 \brief Get the list of compartment names
 \param copasi_model model
 \return sb_string compartmentNames
 \ingroup compartment
*/
COPASIAPIEXPORT sb_string getCompartmentNames (copasi_model);


/*! 
 \brief Get the compartment volume by index
 \param copasi_model model
 \param int index ith compartment
 \return double Voluem of compartment 
 \ingroup compartment
*/
COPASIAPIEXPORT double getCompartmentByIndex (copasi_model, int);


/*! 
 \brief Set a compartment volume by index
 \param copasi_model model
 \param int index ith compartment
 \param double volume Volume of ith compartment
 \ingroup compartment
*/
COPASIAPIEXPORT void setCompartmentByIndex (copasi_model, int, double);


/*! 
 \brief Set a compartment volumes using a vector of compartment values
 \param copasi_model model
 \param double volume Vector of compartment volumes
 \ingroup compartment
*/
COPASIAPIEXPORT void setCompartmentVolumes (copasi_model, sb_matrix v);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Time course simulation
  */
/** \{ */


/*! 
 \brief Simulate using LSODA numerical integrator
 \param copasi_model model
  \param double start time
 \param double end time
 \param int number of steps in the output
 \return sb_matrix matrix of concentration or particles
 
 \code
 result = cvSimulateDeterministic (m, 0.0, 10.0, 100);
 \endcode
 \ingroup simulation
*/
COPASIAPIEXPORT sb_matrix simulateDeterministic(copasi_model model, double startTime, double endTime, int numSteps);


/*! 
 \brief Simulate the differential equation model over one time step
 \param copasi_model model
 \param double time Step
 \return double New time, i.e t_o + timeStep
 \ingroup simulation
*/
COPASIAPIEXPORT double oneStep(copasi_model model, double timeStep);

/*! 
 \brief Simulate using exact stochastic algorithm
 \param copasi_model model
 \param double start time
 \param double end time
 \param int number of steps in the output
 \return sb_matrix matrix of concentration or particles
 \ingroup simulation
*/
COPASIAPIEXPORT sb_matrix simulateStochastic(copasi_model model, double startTime, double endTime, int numSteps);

/*! 
 \brief Simulate using Hybrid algorithm/deterministic algorithm
 \param copasi_model model
  \param double start time
 \param double end time
 \param int number of steps in the output
 \return sb_matrix matrix of concentration or particles
 \ingroup simulation
*/
COPASIAPIEXPORT sb_matrix simulateHybrid(copasi_model model, double startTime, double endTime, int numSteps);

/*! 
 \brief Simulate using Tau Leap stochastic algorithm
 \param copasi_model model
  \param double start time
 \param double end time
 \param int number of steps in the output
 \return sb_matrix matrix of concentration or particles
 \ingroup simulation
*/
COPASIAPIEXPORT sb_matrix simulateTauLeap(copasi_model model, double startTime, double endTime, int numSteps);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Steady state analysis 
  */
/** \{ */


/*! 
 \brief Bring the system to steady state by solving for the zeros of the ODE's. 
             Performs an initial simulation before solving.
 \param copasi_model model
 \return sb_matrix matrix with 1 row and n columns, where n = number of species
 \ingroup steadystate
*/
COPASIAPIEXPORT sb_matrix getSteadyState(copasi_model model);

/*! 
 \brief Bring the system to steady state by doing repeated simulations.
             Use this is getSteadyState
 \param copasi_model model
 \param int max iterations (each iteration doubles the time duration)
 \return sb_matrix matrix with 1 row and n columns, where n = number of species
 \ingroup steadystate
*/
COPASIAPIEXPORT sb_matrix getSteadyStateUsingSimulation(copasi_model model, int iter);

/*! 
 \brief Get the full Jacobian at the current state
 \param copasi_model model
 \return sb_matrix matrix with n rows and n columns, where n = number of species
 \ingroup steadystate
*/
COPASIAPIEXPORT sb_matrix getJacobian(copasi_model model);
/*! 
 \brief Get the eigenvalues of the Jacobian at the current state
 \param copasi_model model
 \return sb_matrix matrix with 1 row and n columns, each containing an eigenvalue
 \ingroup steadystate
*/
COPASIAPIEXPORT sb_matrix getEigenvalues(copasi_model model);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Metabolic control analysis (MCA)
  */
/** \{ */


/*! 
 \brief Compute the unscaled flux control coefficients
 \param copasi_model model
 \return sb_matrix rows consist of the fluxes that are perturbed, and columns consist 
                             of the fluxes that are affected
 \ingroup mca
*/
COPASIAPIEXPORT sb_matrix getUnscaledFluxControlCoeffs(copasi_model model);


/*! 
 \brief Compute the scaled flux control coefficients
 \param copasi_model model
 \return sb_matrix rows consist of the fluxes that are perturbed, and columns consist 
                             of the fluxes that are affected
 \ingroup mca
*/
COPASIAPIEXPORT sb_matrix getScaledFluxControlCoeffs(copasi_model model);


/*!
 \brief Compute the unscaled concentration control coefficients
 \param copasi_model model
 \return sb_matrix rows consist of the fluxes that are perturbed, and columns consist 
                             of the concentrations that are affected
 \ingroup mca
*/
COPASIAPIEXPORT sb_matrix getUnscaledConcentrationControlCoeffs(copasi_model model);


/*! 
 \brief Compute the scaled concentration control coefficients
 \param copasi_model model
 \return sb_matrix rows consist of the fluxes that are perturbed, and columns consist 
                             of the concentrations that are affected
 \ingroup mca
*/
COPASIAPIEXPORT sb_matrix getScaledConcentrationConcentrationCoeffs(copasi_model model);


/*! 
 \brief Compute the unscaled elasticities
 \param copasi_model model
 \return sb_matrix rows consist of the species that are perturbed, and columns consist 
                             of the reactions that are affected
 \ingroup mca
*/
COPASIAPIEXPORT sb_matrix getUnscaledElasticities(copasi_model model);


/*! 
 \brief Compute the scaled elasticities
 \param copasi_model model
 \return sb_matrix rows consist of the species that are perturbed, and columns consist 
                             of the reactions that are affected
 \ingroup mca
*/
COPASIAPIEXPORT sb_matrix getScaledElasticities(copasi_model model);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Stoichiometry matrix and matrix analysis
  */
/** \{ */


/*! 
 \brief Return the full stoichiometry matrix, N
 \param copasi_model model
 \return sb_matrix rows consist of the species and columns are the reactions
 \ingroup matrix
*/
COPASIAPIEXPORT sb_matrix getFullStoichiometryMatrix(copasi_model model);


/*! 
 \brief Return the reduced stoichiometry matrix, Nr
 \param copasi_model model
 \return sb_matrix rows consist of the species and columns are the reactions
 \ingroup matrix
*/
COPASIAPIEXPORT sb_matrix getReducedStoichiometryMatrix(copasi_model model);


/*! 
 \brief Compute the elementary flux modes
 \param copasi_model model
 \return sb_matrix matrix with reactions as rows and flux modes as columns (no column names)
 \ingroup matrix
*/
COPASIAPIEXPORT sb_matrix getElementaryFluxModes(copasi_model model);


/*! 
 \brief Compute the Gamma matrix (i.e. conservation laws)
 \param copasi_model model
 \return sb_matrix 
 \ingroup matrix
*/
COPASIAPIEXPORT sb_matrix getGammaMatrix(copasi_model model);


/*! 
 \brief Compute the K matrix (right nullspace)
 \param copasi_model model
 \return sb_matrix 
 \ingroup matrix
*/
COPASIAPIEXPORT sb_matrix getKMatrix(copasi_model model);


/*! 
 \brief Compute the K0 matrix
 \param copasi_model model
 \return sb_matrix 
 \ingroup matrix
*/
COPASIAPIEXPORT sb_matrix getK0Matrix(copasi_model model);


/*! 
 \brief Compute the L matrix (link matrix, left nullspace)
 \param copasi_model model
 \return sb_matrix 
 \ingroup matrix
*/
COPASIAPIEXPORT sb_matrix getLinkMatrix(copasi_model model);


/*! 
 \brief Compute the L0 matrix
 \param copasi_model model
 \return sb_matrix 
 \ingroup matrix
*/
COPASIAPIEXPORT sb_matrix getL0Matrix(copasi_model model);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Optimization
  */
/** \{ */


/*! 
 \brief fit the model parameters to time-series data
 \param copasi_model model
 \param char * filename (tab separated)
 \param sb_matrix parameters to optimize. rownames should contain parameter names, column 1 contains parameter min-values, and column 2 contains parameter max values
 \param char * pick method. Use of of the following: "GeneticAlgorithm", "LevenbergMarquardt", "SimulatedAnnealing", "NelderMead", "SRES", "ParticleSwarm", "SteepestDescent", "RandomSearch"
 \ingroup matrix
*/
//COPASIAPIEXPORT void cFitModelToData(copasi_model model, const char * filename, sb_matrix params, const char * method);

/*! 
 \brief use genetic algorithms to generate a distribution of parameter values that satisfy an objective function or fit a data file
 \param copasi_model model
 \param char * objective function or filename
 \param sb_matrix parameter initial values and min and max values (3 columns)
 \return sb_matrix optimized parameters as a column vector
 \ingroup optim
*/
COPASIAPIEXPORT sb_matrix optimize(copasi_model model, const char * objective, sb_matrix input);

/*! 
 \brief set the number of iterations for the genetic algorithm based optimizer
 \param int iterations
 \ingroup optim
*/
COPASIAPIEXPORT void setOptimizerIterations(int);

/*! 
 \brief set the number of random seeds for the genetic algorithm based optimizer
 \param int population size
 \ingroup optim
*/
COPASIAPIEXPORT void setOptimizerSize(int);

/*! 
 \brief set the mutation rate, or step size, for the genetic algorithm based optimizer
 \param double 
 \ingroup optim
*/
COPASIAPIEXPORT void setOptimizerMutationRate(double);

/*! 
 \brief set the probability of crossover for the genetic algorithm based optimizer
 \param double must be between 0 and 1
 \ingroup optim
*/
COPASIAPIEXPORT void setOptimizerCrossoverRate(double);

/*! 
 \brief do not modify assignment rules
             warning: disabling this may cause numerical errors in time-course simulations
 \ingroup cleanup
*/
COPASIAPIEXPORT void disableAssignmentRuleReordering();

/*! 
 \brief modify assignment rules to avoid dependencies between assignment rules (default)
 \ingroup cleanup
*/
COPASIAPIEXPORT void enableAssignmentRuleReordering();

END_C_DECLS
#endif

