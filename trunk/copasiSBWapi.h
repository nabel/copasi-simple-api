/*! \mainpage Simplifed API for the COPASI Library
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
 c_matrix output;
  
 m = cReadSBMLFile ("mymodel.xml");
  
 output = cSimulationDeterministic (m, 0, 10, 100); 
 \endcode
 
 More complex example:
 
 \code
 #include <stdlib.h>
 #include <stdio.h>
 #include "copasi_api.h"

 int main(int nargs, char** argv)
 {
        c_matrix efm, output, params;
        copasi_model m1, m2;
        
        if (nargs < 2)
        {
            m1 = model1();
        }
        else
        {
            printf("loading model file %s\n", argv[1]);
            m1 = cReadSBMLFile(argv[1]);        
        }

        cWriteAntimonyFile(m1, "model.txt");
        
        printf("Antimony file written to model.txt\nSimulating...\n");  
        output = cSimulateDeterministic(m1, 0, 100, 1000);  //model, start, end, num. points
        printf("output.tab has %i rows and %i columns\n",output.rows, output.cols);
        c_printMatrixToFile("output.tab", output);
        c_deleteMatrix(output);
                 
        cRemoveModel(m1);
        copasi_end();
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

#ifndef COPASI_SBW_LIKE_C_API
#define COPASI_SBW_LIKE_C_API

 /**
  * @file    copasi_api.h
  * @brief   Simple C API for the Copasi C++ library

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
  * @name Read SBML group
  */
/** \{ */

/*! 
 \brief Read a model from an SBML file
 \param char* file name
 \return copasi_model Copasi model of the SBML file
 \ingroup loadsave
*/
COPASIAPIEXPORT copasi_model readSBMLFile(const char * filename);


// -----------------------------------------------------------------------
/** \} */
/**
  * @name SBW Reaction group
  */
/** \{ */

/*! 
 \brief sGet the number of reactions in the model
 \param copasi_model model
 \return int Returns the number of reactions in the model
 \ingroup reaction
*/
COPASIAPIEXPORT int getNumberOfReactions (copasi_model);


/*! 
 \brief sGet the list of reaction names
 \param copasi_model model
 \return char** array of char * and length n, where n = number of reactions
 \ingroup reaction
*/
COPASIAPIEXPORT char** getReactionNames (copasi_model);


/*! 
 \brief sGet the reaction rate for the ith reaction
 \param copasi_model model
 \param int reactionId
 \return double reaction rate for ith reaction
 \ingroup reaction
*/
COPASIAPIEXPORT double getReactionRate(copasi_model, int);


/*! 
 \brief Returns the vector of current reaction rates
 \param copasi_model model
 \return double array of reaction rates
 \ingroup reaction
*/
COPASIAPIEXPORT double* getReactionRates(copasi_model);


/*! 
 \brief Returns the rates of change given an array of new floating species concentrations
 \param copasi_model model
 \param double array of floating concentrations
 \return double vector of reaction rates
 \ingroup reaction
*/
COPASIAPIEXPORT double* getReactionRatesEx(copasi_model, double * values);

// -----------------------------------------------------------------------
/** \} */
/**
  * @name SBW boundary species group
  */
/** \{ */

/*! 
 \brief sGet the number of boundary species
 \param copasi_model model
 \return number of species
 \ingroup boundary
*/
COPASIAPIEXPORT int getNumberOfBoundarySpecies(copasi_model model);


/*! 
 \brief sGet a list of boundary species names
 \param copasi_model model
 \return char** array of char * and length n, where n = number of species
 \ingroup boundary
*/
COPASIAPIEXPORT char** getBoundarySpeciesNames(copasi_model model);

/*! 
 \brief Set a boundary species concentration by index
 \param copasi_model model
 \param int index ith boundary species
 \ingroup boundary
*/
COPASIAPIEXPORT void setBoundarySpeciesByIndex (copasi_model model, int index, double value);

/*! 
 \brief Set all the boundary species concentration
 \param copasi_model model
 \param sb_matric Vector of boundary species concentrations
 \ingroup boundary
*/
COPASIAPIEXPORT void setBoundarySpeciesConcentrations (copasi_model model, double * d);

/*! 
 \brief Set all the boundary species concentration
 \param copasi_model model
 \param double * boundary species concentrations
 \ingroup boundary
*/
COPASIAPIEXPORT double * getBoundarySpeciesConcentrations (copasi_model model);


/*! 
 \brief sGet a boundary species concentration by index
 \param copasi_model model
 \param int index ith boundary species
 \return double concentration of ith boundary species
 \ingroup state
*/
COPASIAPIEXPORT double getBoundarySpeciesByIndex (copasi_model model, int index);

// -----------------------------------------------------------------------
/** \} */
/**
  * @name SBW floating species group
  */
/** \{ */

/*! 
 \brief Get the number of floating species
 \param copasi_model model
 \return number of species
 \ingroup floating
*/
COPASIAPIEXPORT int getNumberFloatingSpecies(copasi_model model);

/*! 
 \brief Get a list the floating species names
 \param copasi_model model
 \return c_strings array of char * and length n, where n = number of species
 \ingroup state
*/
COPASIAPIEXPORT char** getFloatingSpeciesNames(copasi_model model);


/*! 
 \brief Set a floating species concentration by index
 \param copasi_model model
 \param int index ith floating species
 \ingroup state
*/
COPASIAPIEXPORT void setFloatingSpeciesByIndex (copasi_model model, int index, double value);


/*! 
 \brief sGet a floating species concentration by index
 \param copasi_model model
 \param int index ith floating species
 \return double Concentration of ith floating species
 \ingroup state
*/
COPASIAPIEXPORT double getFloatingSpeciesByIndex (copasi_model model, int index);

/*! 
 \brief Set all the floating species concentration
 \param copasi_model model
 \param sb_matric Vector of floating species concentrations
 \ingroup boundary
*/
COPASIAPIEXPORT void setFloatingSpeciesConcentrations (copasi_model model, double * );

/*! 
 \brief Set all the floating species concentration  - CURRENTLY NOT IMPLEMENTED
 \param copasi_model model
 \param c_matrix Vector of floating species concentrations
 \ingroup boundary
*/
COPASIAPIEXPORT double * getFloatingSpeciesConcentrations (copasi_model model);


/*! 
 \brief Get the initial floating species concentrations
 \param copasi_model model
 \return double * vector of initial floating species concentrations
 \ingroup floating
*/
COPASIAPIEXPORT double* getFloatingSpeciesIntitialConcentrations (copasi_model model);


/*! 
 \brief Set the initial floating species concentrations 
 \param copasi_model model
 \param c_matrix Vector of initial floating species concentrations
 \ingroup floating
*/
COPASIAPIEXPORT void setFloatingSpeciesIntitialConcentrations (copasi_model model, double *);

/*! 
 \brief Set the initial floating species concentration of the ith species 
 \param copasi_model model
 \param double value value to set the ith initial floating species concentration
 \ingroup floating
*/
COPASIAPIEXPORT void setFloatingSpeciesIntitialConcentrationByIndex (copasi_model model, int index, double sp);


/*! 
 \brief sGet the number of global parameters
 \param copasi_model model
 \return int numberOfGlobalParameters. Returns the number of gloabal parameters in the model
 \ingroup parameter
*/
COPASIAPIEXPORT int getNumberOfGlobalParameters (copasi_model);

/*! 
 \brief sGet the list of global parameter names
 \param copasi_model model 
 \return char** array of char * and length n, where n = number of global parameters
 \ingroup parameter
*/
COPASIAPIEXPORT char** getGlobalParameterNames (copasi_model);

/*! 
 \brief sGet the value of a global parameter by index
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
 \brief set the vector of global parameters
 \param copasi_model model 
 \paramn double* a vector containing values for the global parameters in same order as names
 \ingroup parameter
*/
COPASIAPIEXPORT void setGlobalParameterValues (copasi_model, double*);

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
COPASIAPIEXPORT int getNumberOfCompartments (copasi_model);

/*! 
 \brief sGet the list of compartment names
 \param copasi_model model
 \return char** compartment names
 \ingroup compartment
*/
COPASIAPIEXPORT char** getCompartmentNames (copasi_model);


/*! 
 \brief sGet the compartment volume by index
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
COPASIAPIEXPORT void setCompartmentVolumes (copasi_model, c_matrix v);


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
COPASIAPIEXPORT double* getRatesOfChange(copasi_model);

/*! 
 \brief Compute the current rates of change for the ith species
 \param copasi_model model
 \param ith rate of change to compute
 \return double ith rate of change
 \ingroup rateOfChange
*/
COPASIAPIEXPORT double getRateOfChange(copasi_model, int index);

/*! 
 \brief Returns the names used to represent the rates of change
 \param copasi_model model
 \return char** list of names used to represent the rate of change
 \ingroup rateOfChange
*/
COPASIAPIEXPORT char** getRatesOfChangeNames(copasi_model);

/*! 
 \brief Returns the rates of change given a vector of floating species concentrations
 \param copasi_model model
 \param double* array of floating species concentrations
 \return double* array of rates of change
 \ingroup sbw_rateOfChange
*/
COPASIAPIEXPORT double* getRatesOfChangeEx(copasi_model model, double* sp);

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
COPASIAPIEXPORT c_matrix simulate(copasi_model model, double timeStart, double timeEnd, int numOfPoints);

/*! 
 \brief Simulate the differential equation model over one time step
 \param copasi_model model
 \param double time step
 \return double new time, i.e (current time + timeStep)
 \ingroup simulation
*/
COPASIAPIEXPORT double oneStep(copasi_model model, double timeStep);

END_C_DECLS
#endif

