/*! \mainpage SBW API for the COPASI Library
 *
 * \section intro_sec Introduction
 *
 * The developers of COPASI provide COPASI as a reusable library as well as
 * the well known COPASI user interface. The library however has a fairly
 * complex API and can take some time sGetting used. We have therefore layered 
 * on top of the COPASI library a new C based API that we feel is much simpler
 * to use. For example, to run a simple SBML model and generate time series data
 * we would call:
 *  
 \code
 copasiModel m;
 sb_matrix output;
  
 m = sReadSBMLFile ("mymodel.xml");
  
 output = simulationDeterministic (m, 0, 10, 100); 
 \endcode
 
 More complex example:
 
 \code
 #include <stdlib.h>
 #include <stdio.h>
 #include "copasiSBApi.h"

 int main(int nargs, char** argv)
 {
        sb_matrix efm, output, params;
        sCcopasiModel m1, m2;
        
        if (nargs < 2)
        {
            m1 = model1();
        }
        else
        {
            printf("loading model file %s\n", argv[1]);
            m1 = sReadSBMLFile(argv[1]);        
        }

        sWriteAntimonyFile(m1, "model.txt");
        
        printf("Antimony file written to model.txt\nSimulating...\n");  
        output = simulateDeterministic(m1, 0, 100, 1000);  //model, start, end, num. points
        printf("output.tab has %i rows and %i columns\n",output.rows, output.cols);
        sb_printMatrixToFile("output.tab", output);
        sb_deleteMatrix(output);
                 
        sRemoveModel(m1);
        sCopasiEnd();
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
 \brief sGet information about reaction rates
 
 \defgroup rateOfChange Rates of change group
 \brief sGet information about rates of change

 \defgroup boundary Boundary species group
 \brief sGet information about reaction rates
 
 \defgroup floating Floating species group
 \brief sGet information about reaction rates
  
 \defgroup parameters Parameter group
 \brief set and sGet global and local parameters
 
 \defgroup compartment Compartment group
 \brief set and sGet information on compartments
 
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

#include "sb_structs.h"

/*!\brief This struct is used to contain a pointer to an instance of a COPASI object*/
typedef struct  
{ 
	void * CopasiModelPtr;
	void * CopasiDataModelPtr;
	void * qHash;
	char * errorMessage;
	char * warningMessage;
} copasiModel;

/*!\brief This struct is used to contain a pointer to an instance of a COPASI reaction object*/
typedef struct  
{
	void * CopasiReactionPtr;
	void * CopasiModelPtr;
	void * qHash; 
} copasiReaction;

/*!\brief This struct is used to contain a pointer to an instance of a COPASI compartment object*/
typedef struct  
{
	void * CopasiCompartmentPtr;
	void * CopasiModelPtr; 
	void * qHash; 
} copasiCompartment;

BEGIN_C_DECLS

// -----------------------------------------------------------------------
/**
  * @name Memory management
  */
/** \{*/

/*! 
 \brief destroy copasi -- MUST BE CALLED at the end of program
 \ingroup memory
*/
SBAPIEXPORT void sCopasiEnd();

/*! 
 \brief remove a model
 \ingroup memory
*/
SBAPIEXPORT void sRemoveModel (copasiModel);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Read and write models
  */
/** \{ */

/*! 
 \brief Create a model from an Antimony, see antimony.sf.net for details of Antimony syntax
 \param char* file name
 \return copasiModel Copasi model of the Antimony file
 \ingroup loadsave
*/
SBAPIEXPORT copasiModel sReadAntimonyFile(const char * filename);


/*! 
 \brief Create a model from an Antimony string
 \param char* Antimony string
 \return copasiModel Copasi model of the Antimony string
 \ingroup loadsave
*/
SBAPIEXPORT copasiModel sRreadAntimonyString(const char * sbml);


/*! 
 \brief Create a model from an SBML file
 \param char* file name
 \return copasiModel Copasi model of the SBML file
 \ingroup loadsave
*/
SBAPIEXPORT copasiModel sReadSBMLFile(const char * filename);


/*! 
 \brief Create a model from an SBML string
 \param char* SBML string
 \return copasiModel Copasi model of the SBML string
 \ingroup loadsave
*/
SBAPIEXPORT copasiModel sRreadSBMLString(const char * sbml);


/*! 
 \brief Save a model as an SBML file
 \param copasiModel copasi model
 \param char* file name
 \ingroup loadsave
*/
SBAPIEXPORT void sWriteSBMLFile(copasiModel model, const char * filename);


/*! 
 \brief Save a model as an Antimony file, see antimony.sf.net for details of Antimony syntax
 \param copasiModel copasi model
 \param char* file name
 \ingroup loadsave
*/
SBAPIEXPORT void sWriteAntimonyFile(copasiModel model, const char * filename);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Create model group
  */
/** \{ */

/*! 
 \brief Create a model
 \param char* model name
 \return copasiModel a new copasi model
 \ingroup create
*/
SBAPIEXPORT copasiModel sCreateModel(const char * name);


/*! 
 \brief This function is only needed for calling COPASI methods not found in this library. This function compiles the COPASI model; it is called internally by the simulate and other anlysis functions. 
 \param copasiModel model
 \ingroup create
*/
SBAPIEXPORT void sCompileModel(copasiModel model);


/*! 
 \brief Create compartment
 \param char* compartment name
 \param double volume
 \return copasiCompartment a new compartment
 \ingroup create
*/
SBAPIEXPORT copasiCompartment sCreateCompartment(copasiModel model, const char* name, double volume);


/*! 
 \brief Set a volume of compartment
 \param copasiModel model
 \param char * compartment name
 \param double volume
 \ingroup create
*/
SBAPIEXPORT void sSetVolume(copasiModel, const char * compartment, double volume);


/*! 
 \brief Set the concentration of a species, volume of a compartment, or value of a parameter
      The function will figure out which using the name (fast lookup using hashtables).
      If the name does not exist in the model, a new global parameter will be created.
 \param copasiModel model
 \param char * name
 \param double value
 \return 0 if new variable was created. 1 if existing variable was found
 \ingroup create
*/
SBAPIEXPORT int sSetValue(copasiModel, const char * name, double value);


/*! 
 \brief Add a species to the model
 \param copasiCompartment model
 \param char* species name
 \param double initial value (concentration or count, depending on the model)
 \ingroup create
*/
SBAPIEXPORT void sCreateSpecies(copasiCompartment compartment, const char* name, double initialValue);


/*! 
 \brief Set a species as boundary or floating (will remove any assignment rules)
 \param copasiModel model
  \param char * name
 \param int boundary = 1, floating = 0 (default)
 \ingroup create
*/
SBAPIEXPORT void sSetSpeciesType(copasiModel model, const char * species, int isBoundary);


/*! 
 \brief Set a species concentration
 \param copasiModel model
 \param char * species name
 \param double concentration
 \ingroup create
*/
SBAPIEXPORT void sSetConcentration(copasiModel, const char * species, double conc);

/*! 
 \brief Set a species amounts
 \param copasiModel model
 \param char * species name
 \param double amount
 \ingroup create
*/
SBAPIEXPORT void sSetAmount(copasiModel, const char * species, double amount);


/*! 
 \brief Set the assignment rule for a species (automatically assumes boundary species)
 \param copasiModel model
 \param char * species name
 \param char* formula, use 0 to remove assignment rule
 \return int 0=failed 1=success
 
 \code
 result = cSetAssignmentRule (m, "S1", "sin (time*k1)");
 \endcode
 \ingroup create
*/
SBAPIEXPORT int sSetAssignmentRule(copasiModel model, const char * species, const char * formula);


/*! 
 \brief Set the value of an existing global parameter or create a new global parameter
 \param copasiModel model
 \param char* parameter name
 \param double value
  \return int 0=new value created 1=found existing value
 \ingroup create
*/
SBAPIEXPORT int sSetGlobalParameter(copasiModel model, const char * name, double value);


/*! 
 \brief Create a new variable that is not a constant by a formula
 \param copasiModel model
 \param char* name of new variable
 \param char* formula
 \return int 0=failed 1=success
 \ingroup create
*/
SBAPIEXPORT int sCreateVariable(copasiModel model, const char * name, const char * formula);


/*! 
 \brief Add a trigger and a response, where the response is defined by a target variable and an assignment formula
 \param copasiModel model
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
SBAPIEXPORT int sCreateEvent(copasiModel model, const char * name, const char * trigger, const char * variable, const char * formula);


/*!
 \brief Create a new reaction with a given name
 \param copasiModel model
 \param char* reaction name
 \return copasiReaction a new reaction
 
 \code
 r = cCreateReaction (m, "J1")
 \endcode
 \ingroup create
*/
SBAPIEXPORT copasiReaction sCreateReaction(copasiModel model, const char* name);


/*! 
 \brief Add a reactant to a reaction
 \param copasiReaction reaction
 \param char * reactant
 \param double stoichiometry
 
 \code
 cCreateReaction (m, "S1", 1);
 \endcode
 \ingroup create
*/
SBAPIEXPORT void sAddReactant(copasiReaction reaction, const char * species, double stoichiometry);

/*! 
 \brief Add a product to a reaction
 \param copasiReaction reaction
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
SBAPIEXPORT void sAddProduct(copasiReaction reaction, const char * species, double stoichiometry);

/*! 
 \brief Set reaction rate equation
 \param copasiReaction reaction
 \param char* custom formula
 \return int success=1 failure=0
 
 \code
 int result;
 result = cSetReactionRate (r, "k1*S1");
 \endcode
 
 \ingroup create
*/
SBAPIEXPORT int sSetReactionRate(copasiReaction reaction, const char * formula);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Reaction group
  */
/** \{ */

/*! 
 \brief sGet the number of reactions in the model
 \param copasiModel model
 \return int Returns the number of reactions in the model
 \ingroup reaction
*/
SBAPIEXPORT sb_matrix sGetNumberOfReactions (copasiModel);


/*! 
 \brief sGet the list of reaction names
 \param copasiModel model
 \return sb_strings array of char * and length n, where n = number of reactions
 \ingroup reaction
*/
SBAPIEXPORT sb_strings sGetReactionNames (copasiModel);


/*! 
 \brief sGet the reaction rate for the ith reaction
 \param copasiModel model
 \param int reactionId
 \return double reaction rate for ith reaction
 \ingroup reaction
*/
SBAPIEXPORT int sGetReactionRate(copasiModel, int);


/*! 
 \brief Returns the vector of current reaction rates
 \param copasiModel model
 \return double array of reaction rates
 \ingroup reaction
*/
SBAPIEXPORT sb_matrix sGetReactionRates(copasiModel);


/*! 
 \brief Returns the rates of change given an array of new floating species concentrations
 \param copasiModel model
 \param sb_matrix Vector of floating concentrations
 \return sb_matrix Vector of reaction rates
 \ingroup reaction
*/
SBAPIEXPORT sb_matrix sGetReactionRatesEx(sb_matrix values);

// -----------------------------------------------------------------------
/** \} */
/**
  * @name Boundary species group
  */
/** \{ */

/*! 
 \brief sGet the number of boundary species - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \return number of species
 \ingroup boundary
*/
SBAPIEXPORT int sGetNumberOfBoundarySpecies(copasiModel model);


/*! 
 \brief sGet a list of boundary species names - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \return sb_strings array of char * and length n, where n = number of species
 \ingroup boundary
*/
SBAPIEXPORT sb_strings sGetBoundarySpeciesNames(copasiModel model);

/*! 
 \brief Set a boundary species concentration by index - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param int index ith boundary species
 \ingroup boundary
*/
SBAPIEXPORT void sSetBoundarySpeciesByIndex (copasiModel model, int index);

/*! 
 \brief Set all the boundary species concentration  - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param sb_matric Vector of boundary species concentrations
 \ingroup boundary
*/
SBAPIEXPORT void sSetBoundarySpeciesConcentrations (copasiModel model, sb_matrix d);

/*! 
 \brief Set all the boundary species concentration  - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param sb_matrix Vector of boundary species concentrations
 \ingroup boundary
*/
SBAPIEXPORT sb_matrix sGetBoundarySpeciesConcentrations (copasiModel model);


/*! 
 \brief sGet a boundary species concentration by index - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param int index ith boundary species
 \return double Concentration of ith boundary species
 \ingroup state
*/
SBAPIEXPORT double sGetBoundarySpeciesByIndex (copasiModel model, int index);

// -----------------------------------------------------------------------
/** \} */
/**
  * @name Floating species group
  */
/** \{ */

/*! 
 \brief sGet the number of floating species - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \return number of species
 \ingroup floating
*/
SBAPIEXPORT int sGetNumberFloatingSpecies(copasiModel model);

/*! 
 \brief sGet a list the floating species names - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \return sb_strings array of char * and length n, where n = number of species
 \ingroup state
*/
SBAPIEXPORT sb_strings sGetFloatingSpeciesNames(copasiModel model);


/*! 
 \brief Set a floating species concentration by index - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param int index ith floating species
 \ingroup state
*/
SBAPIEXPORT void sSetFloatingSpeciesByIndex (copasiModel model, int index);


/*! 
 \brief sGet a floating species concentration by index - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param int index ith floating species
 \return double Concentration of ith floating species
 \ingroup state
*/
SBAPIEXPORT double sGetFloatingSpeciesByIndex (copasiModel model, int index);

/*! 
 \brief Set all the floating species concentration  - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param sb_matric Vector of floating species concentrations
 \ingroup boundary
*/
SBAPIEXPORT void sSetFloatingSpeciesConcentrations (copasiModel model, sb_matrix sp);

/*! 
 \brief Set all the floating species concentration  - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param sb_matrix Vector of floating species concentrations
 \ingroup boundary
*/
SBAPIEXPORT sb_matrix sGetFloatingSpeciesConcentrations (copasiModel model);


/*! 
 \brief sGet the initial floating species concentrations  - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \return sb_matrix Vector of initial floating species concentrations
 \ingroup floating
*/
SBAPIEXPORT sb_matrix sGetFloatingSpeciesIntitialConcentrations (copasiModel model);


/*! 
 \brief Set the initial floating species concentrations  - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param sb_matrix Vector of initial floating species concentrations
 \ingroup floating
*/
SBAPIEXPORT void sSetFloatingSpeciesIntitialConcentrations (copasiModel model, sb_matrix sp);

/*! 
 \brief Set the initial floating species concentration of the ith species  - CURRENTLY NOT IMPLEMENTED
 \param copasiModel model
 \param double value value to set the ith initial floating species concentration
 \ingroup floating
*/
SBAPIEXPORT void sSetFloatingSpeciesIntitialConcentrationByIndex (copasiModel model, int index, double sp);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Rates of change group
  */
/** \{ */

/*! 
 \brief Compute the current rates of change for all species
 \param copasiModel model
 \return sb_matrix matrix of with 1 row and n columns, where n = number of species
 \ingroup rateOfChange
*/
SBAPIEXPORT sb_matrix sGetRatesOfChange(copasiModel);

/*! 
 \brief Compute the current rates of change for the ith species
 \param copasiModel model
 \param ith rate of change to compute
 \return double ith rate of change
 \ingroup rateOfChange
*/
SBAPIEXPORT double sGetRateOfChange(copasiModel, int index);

/*! 
 \brief Returns the names used to represent the rates of change
 \param copasiModel model
 \return sb_strings List of names used to represent the rate of change
 \ingroup rateOfChange
*/
SBAPIEXPORT sb_matrix sGetRatesOfChangeNames(copasiModel);

/*! 
 \brief Returns the rates of change given a vector of floating species concentrations
 \param copasiModel model
 \param sb_matrix vector of floating species concentrations
 \return sb_matrix vector of rates of change
 \ingroup rateOfChange
*/
SBAPIEXPORT sb_matrix sGetRatesOfChangeEx(copasiModel model, sb_matrix sp);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Current state of system
  */
/** \{ */

/*! 
 \brief sGet the current concentrations of all species
 \param copasiModel model
 \return sb_matrix matrix of with 1 row and n columns, where n = number of species
 The names of the species are included in the matrix column labels
 \ingroup state
*/
SBAPIEXPORT sb_matrix sGetConcentrations(copasiModel);

/*! 
 \brief sGet the current amounts of all species. The amounts are calculated from the concentrations and compartment volume
 \param copasiModel model
 \return sb_matrix matrix of with 1 row and n columns, where n = number of species
 The names of the species are included in the matrix column labels
 \ingroup state
*/
SBAPIEXPORT sb_matrix sGetAmounts(copasiModel);

/*! 
 \brief sGet a list of all species names, including floating and boundary species
 \param copasiModel model
 \return sb_strings array of char * and length n, where n = number of species
 \ingroup state
*/
SBAPIEXPORT sb_strings sGetAllSpeciesNames(copasiModel model);


/*! 
 \brief sGet the current concentration of a species
 \param copasiModel model
 \param string species name
 \return double concentration. -1 indicates that a species by this name was not found
 \ingroup state
*/
SBAPIEXPORT double sGetConcentration(copasiModel, const char * name);

/*! 
 \brief sGet the current amount of a species. The amounts are calculated from the concentrations and compartment volume
 \param copasiModel model
 \param string species name
 \return double amount. -1 indicates that a species by this name was not found
 \ingroup state
*/
SBAPIEXPORT double sGgetAmount(copasiModel, const char * name);


/*! 
 \brief Compute current flux through the given reactions
 \param copasiModel model
 \param string reaction name, e.g. "J1"
 \return double rate. If reaction by this name does not exist that NaN will be returned
  The names of the fluxes are included in the matrix column labels
 \ingroup state
*/
SBAPIEXPORT double sGgetFlux(copasiModel, const char * name);

/*! 
 \brief Compute current flux through the given reactions in terms of particles
 \param copasiModel model
 \param string reaction name, e.g. "J1"
 \return double rate. If reaction by this name does not exist that NaN will be returned
 \ingroup state
*/
SBAPIEXPORT double sGgetParticleFlux(copasiModel, const char * name);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Parameter group
  */
/** \{ */

/*! 
 \brief sGet the number of global parameters
 \param copasiModel model
 \return int numberOfGlobalParameters. Returns the number of gloabal parameters in the model
 \ingroup parameter
*/
SBAPIEXPORT int sGetNumberOfGlobalParameters (copasiModel);

/*! 
 \brief sGet the list of global parameter names
 \param copasiModel model 
 \return sb_strings array of char * and length n, where n = number of global parameters
 \ingroup parameter
*/
SBAPIEXPORT sb_strings sGetGlobalParameterNames (copasiModel);

/*! 
 \brief sGet the value of a global parameter by index
 \param copasiModel model 
 \param int index ith global parameter
 \return double returned value of parameter
 \ingroup parameter
*/
SBAPIEXPORT double sGetGlobalParameterByIndex (copasiModel, int);

/*! 
 \brief Set the value of a global parameter by index
 \param copasiModel model 
 \param int index ith global parameter
 \param double Value to set parameter to
 \ingroup parameter
*/
SBAPIEXPORT void setGlobalParameterByIndex (copasiModel, int, double);

/*! 
 \brief sGet a list of the global parameters
 \param copasiModel model 
 \return sb_matrix A vector containing the values for the global parameters.
 \ingroup parameter
*/
SBAPIEXPORT sb_matrix sGetGlobalParameters (copasiModel);

/*! 
 \brief Set the vector of global parameters
 \param copasiModel model 
 \paramn sb_matrix A vector containing the values for the global parameters.
 \ingroup parameter
*/
SBAPIEXPORT void sSetGlobalParameterValues (copasiModel, sb_matrix gp);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Compartment group
  */
/** \{ */

/*! 
 \brief sGet the number of compartments
 \param copasiModel model
 \return int numberOfCompartments
 \ingroup compartment
*/
SBAPIEXPORT int sGetNumberOfCompartments (copasiModel);

/*! 
 \brief sGet the list of compartment names
 \param copasiModel model
 \return sb_strings compartmentNames
 \ingroup compartment
*/
SBAPIEXPORT sb_strings sGetCompartmentNames (copasiModel);


/*! 
 \brief sGet the compartment volume by index
 \param copasiModel model
 \param int index ith compartment
 \return double Voluem of compartment 
 \ingroup compartment
*/
SBAPIEXPORT double sGetCompartmentByIndex (copasiModel, int);


/*! 
 \brief Set a compartment volume by index
 \param copasiModel model
 \param int index ith compartment
 \param double volume Volume of ith compartment
 \ingroup compartment
*/
SBAPIEXPORT void sSetCompartmentByIndex (copasiModel, int, double);


/*! 
 \brief Set a compartment volumes using a vector of compartment values
 \param copasiModel model
 \param double volume Vector of compartment volumes
 \ingroup compartment
*/
SBAPIEXPORT void sSetCompartmentVolumes (copasiModel, sb_matrix v);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Time course simulation
  */
/** \{ */


/*! 
 \brief Simulate using LSODA numerical integrator
 \param copasiModel model
  \param double start time
 \param double end time
 \param int number of steps in the output
 \return sb_matrix matrix of concentration or particles
 
 \code
 result = cvSimulateDeterministic (m, 0.0, 10.0, 100);
 \endcode
 \ingroup simulation
*/
SBAPIEXPORT sb_matrix sSimulateDeterministic (copasiModel model, double startTime, double endTime, int numSteps);


/*! 
 \brief Simulate the differential equation model over one time step
 \param copasiModel model
 \param double time Step
 \return double New time, i.e t_o + timeStep
 \ingroup simulation
*/
SBAPIEXPORT double sOneStep(copasiModel model, double timeStep);

/*! 
 \brief Simulate using exact stochastic algorithm
 \param copasiModel model
 \param double start time
 \param double end time
 \param int number of steps in the output
 \return sb_matrix matrix of concentration or particles
 \ingroup simulation
*/
SBAPIEXPORT sb_matrix sSimulateStochastic(copasiModel model, double startTime, double endTime, int numSteps);

/*! 
 \brief Simulate using Hybrid algorithm/deterministic algorithm
 \param copasiModel model
  \param double start time
 \param double end time
 \param int number of steps in the output
 \return sb_matrix matrix of concentration or particles
 \ingroup simulation
*/
SBAPIEXPORT sb_matrix sSimulateHybrid(copasiModel model, double startTime, double endTime, int numSteps);

/*! 
 \brief Simulate using Tau Leap stochastic algorithm
 \param copasiModel model
  \param double start time
 \param double end time
 \param int number of steps in the output
 \return sb_matrix matrix of concentration or particles
 \ingroup simulation
*/
SBAPIEXPORT sb_matrix sSimulateTauLeap(copasiModel model, double startTime, double endTime, int numSteps);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Steady state analysis 
  */
/** \{ */


/*! 
 \brief Bring the system to steady state by solving for the zeros of the ODE's. 
             Performs an initial simulation before solving.
 \param copasiModel model
 \return sb_matrix matrix with 1 row and n columns, where n = number of species
 \ingroup steadystate
*/
SBAPIEXPORT sb_matrix sGetSteadyState(copasiModel model);

/*! 
 \brief Bring the system to steady state by doing repeated simulations.
             Use this is sGetSteadyState
 \param copasiModel model
 \param int max iterations (each iteration doubles the time duration)
 \return sb_matrix matrix with 1 row and n columns, where n = number of species
 \ingroup steadystate
*/
SBAPIEXPORT sb_matrix sGetSteadyStateUsingSimulation(copasiModel model, int iter);

/*! 
 \brief sGet the full Jacobian at the current state
 \param copasiModel model
 \return sb_matrix matrix with n rows and n columns, where n = number of species
 \ingroup steadystate
*/
SBAPIEXPORT sb_matrix sGetJacobian(copasiModel model);
/*! 
 \brief sGet the eigenvalues of the Jacobian at the current state
 \param copasiModel model
 \return sb_matrix matrix with 1 row and n columns, each containing an eigenvalue
 \ingroup steadystate
*/
SBAPIEXPORT sb_matrix sGetEigenvalues(copasiModel model);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Metabolic control analysis (MCA)
  */
/** \{ */


/*! 
 \brief Compute the unscaled flux control coefficients
 \param copasiModel model
 \return sb_matrix rows consist of the fluxes that are perturbed, and columns consist 
                             of the fluxes that are affected
 \ingroup mca
*/
SBAPIEXPORT sb_matrix sGetUnscaledFluxControlCoeffs(copasiModel model);


/*! 
 \brief Compute the scaled flux control coefficients
 \param copasiModel model
 \return sb_matrix rows consist of the fluxes that are perturbed, and columns consist 
                             of the fluxes that are affected
 \ingroup mca
*/
SBAPIEXPORT sb_matrix sGetScaledFluxControlCoeffs(copasiModel model);


/*!
 \brief Compute the unscaled concentration control coefficients
 \param copasiModel model
 \return sb_matrix rows consist of the fluxes that are perturbed, and columns consist 
                             of the concentrations that are affected
 \ingroup mca
*/
SBAPIEXPORT sb_matrix sGetUnscaledConcentrationControlCoeffs(copasiModel model);


/*! 
 \brief Compute the scaled concentration control coefficients
 \param copasiModel model
 \return sb_matrix rows consist of the fluxes that are perturbed, and columns consist 
                             of the concentrations that are affected
 \ingroup mca
*/
SBAPIEXPORT sb_matrix sGetScaledConcentrationConcentrationCoeffs(copasiModel model);


/*! 
 \brief Compute the unscaled elasticities
 \param copasiModel model
 \return sb_matrix rows consist of the species that are perturbed, and columns consist 
                             of the reactions that are affected
 \ingroup mca
*/
SBAPIEXPORT sb_matrix sGetUnscaledElasticities(copasiModel model);


/*! 
 \brief Compute the scaled elasticities
 \param copasiModel model
 \return sb_matrix rows consist of the species that are perturbed, and columns consist 
                             of the reactions that are affected
 \ingroup mca
*/
SBAPIEXPORT sb_matrix sGetScaledElasticities(copasiModel model);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Stoichiometry matrix and matrix analysis
  */
/** \{ */


/*! 
 \brief Return the full stoichiometry matrix, N
 \param copasiModel model
 \return sb_matrix rows consist of the species and columns are the reactions
 \ingroup matrix
*/
SBAPIEXPORT sb_matrix sGetFullStoichiometryMatrix(copasiModel model);


/*! 
 \brief Return the reduced stoichiometry matrix, Nr
 \param copasiModel model
 \return sb_matrix rows consist of the species and columns are the reactions
 \ingroup matrix
*/
SBAPIEXPORT sb_matrix sGetReducedStoichiometryMatrix(copasiModel model);


/*! 
 \brief Compute the elementary flux modes
 \param copasiModel model
 \return sb_matrix matrix with reactions as rows and flux modes as columns (no column names)
 \ingroup matrix
*/
SBAPIEXPORT sb_matrix sGetElementaryFluxModes(copasiModel model);


/*! 
 \brief Compute the Gamma matrix (i.e. conservation laws)
 \param copasiModel model
 \return sb_matrix 
 \ingroup matrix
*/
SBAPIEXPORT sb_matrix sGetGammaMatrix(copasiModel model);


/*! 
 \brief Compute the K matrix (right nullspace)
 \param copasiModel model
 \return sb_matrix 
 \ingroup matrix
*/
SBAPIEXPORT sb_matrix sGetKMatrix(copasiModel model);


/*! 
 \brief Compute the K0 matrix
 \param copasiModel model
 \return sb_matrix 
 \ingroup matrix
*/
SBAPIEXPORT sb_matrix sGetK0Matrix(copasiModel model);


/*! 
 \brief Compute the L matrix (link matrix, left nullspace)
 \param copasiModel model
 \return sb_matrix 
 \ingroup matrix
*/
SBAPIEXPORT sb_matrix sGetLinkMatrix(copasiModel model);


/*! 
 \brief Compute the L0 matrix
 \param copasiModel model
 \return sb_matrix 
 \ingroup matrix
*/
SBAPIEXPORT sb_matrix sGetL0Matrix(copasiModel model);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name Optimization
  */
/** \{ */


/*! 
 \brief do not modify assignment rules
             warning: disabling this may cause numerical errors in time-course simulations
 \ingroup cleanup
*/
SBAPIEXPORT void sDisableAssignmentRuleReordering();

/*! 
 \brief modify assignment rules to avoid dependencies between assignment rules (default)
 \ingroup cleanup
*/
SBAPIEXPORT void sEnableAssignmentRuleReordering();

END_C_DECLS
#endif

