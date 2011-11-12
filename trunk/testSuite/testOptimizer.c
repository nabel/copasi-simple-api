#include <stdlib.h>
#include <stdio.h>
#include "copasi_api.h"

int main(int narg, char** argv)
{
	tc_matrix results;
	copasi_model m;
	
	if (narg < 2)
	{
		printf("Please specify the SBML file\n");
		return 0;
	}

	//cSetSBMLLevelAndVersion(2,3);
	m = cReadSBMLFile(argv[1]);
	cCompileModel(m);

	if (m.errorMessage != NULL)
	{
		printf ("Errors while reading model:\n");
		printf ("%s\n", m.errorMessage);
		getchar();
		return 0;
	}

	cWriteAntimonyFile (m, "antimony.txt");
	printf ("Antimony script written to antimony.txt\n");

	printf("simulating...\n");	
	results = cSimulateDeterministic(m, 0, 20, 100);  //model, start, end, num. points

	printf("results.tab has simulation data\n\n");
	tc_printMatrixToFile("results.tab", results);
	tc_deleteMatrix(results);

	results = cGetReactionRates(m);

	printf("fluxes:\n");
	tc_printOutMatrix(results);
	tc_deleteMatrix(results);

	results = cGetRatesOfChange(m);
	printf("\n\nderivatives:\n");
	tc_printOutMatrix(results);
	tc_deleteMatrix(results);
	
	//printf("%s\n",m1.errorMessage);
	/*
	//the optmization code below works but 
	//has been tested only a couple of times

	//setup for optimization using GA
	params = tc_createMatrix(3,3);
	tc_setRowName(params,0,"k1");
	tc_setRowName(params,1,"k2");
	tc_setRowName(params,2,"k3");

	//intial values
	tc_setMatrixValue(params, 0, 0, 1);
	tc_setMatrixValue(params, 0, 1, 0.0);
	tc_setMatrixValue(params, 0, 2, 5.0);

 	//min values
	tc_setMatrixValue(params, 1, 0, 1);
	tc_setMatrixValue(params, 1, 1, 0.0);
	tc_setMatrixValue(params, 1, 2, 5.0);

	//max values
	tc_setMatrixValue(params, 2, 0, 1);
	tc_setMatrixValue(params, 2, 1, 0.0);
	tc_setMatrixValue(params, 2, 2, 5.0);
	
	cSetValue(m1,"k1",2.0);
	cSetValue(m1,"k2",1.0);
	cSetValue(m1,"k3",1.0);
	
	cSetOptimizerIterations(10);  //set num interations
	results = cOptimize(m1, "output.tab", params); //optimize
	tc_printMatrixToFile("params.out", results);  //optimized parameters
	tc_deleteMatrix(results);
	*/
	
	//cleanup	
	cRemoveModel(m);
	copasi_end();
	printf ("Hit the return key to continue\n");
	getchar();
	return 0;
}

