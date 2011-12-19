#include <stdlib.h>
#include <stdio.h>
#include "copasi_api.h"

int main(int narg, char** argv)
{

	return main1();

	int i;
	double d;
	c_matrix results, params, params2;
	const char * filename = "results.tab";
	copasi_model m;
	
	if (narg < 2)
	{
		printf("Please specify the Antimony file\n");
		return 0;
	}

	//cSetSBMLLevelAndVersion(2,3);
	m = cReadAntimonyFile(argv[1]);

	if (m.errorMessage != NULL)
	{
		printf ("Errors while reading model:\n");
		printf ("%s\n", m.errorMessage);
		getchar();
		return 0;
	}

	printf("simulating...\n");	
	results = cSimulateDeterministic(m, 0, 10, 1000);  //model, start, end, num. points

	printf("%s has original simulation data\n\n",filename);
	c_printMatrixToFile(filename, results);
	c_deleteMatrix(results);

	//jumble the parameters
	params = cGetGlobalParameters(m);
	/*for (i=0; i < params.rows; ++i)
	{
		d = runif(0,10);
		cSetValue(m, c_getRowName(params,i), d);
		c_setMatrixValue(params, i, 0, d);
	}*/

	results = cSimulateDeterministic(m, 0, 10, 1000);  //model, start, end, num. points

	printf("results2.tab has original simulation data\n\n","results2.tab");
	c_printMatrixToFile("results2.tab", results);
	c_deleteMatrix(results);
	c_printOutMatrix(params);
	
	cFitModelToData(m, filename, params, "levenbergmarquardt");
	c_deleteMatrix(params);
	params = cGetGlobalParameters(m);
	c_printOutMatrix(params);
	c_deleteMatrix(params);
	//printf("%s\n",m1.errorMessage);
	/*
	//the optmization code below works but 
	//has been tested only a couple of times

	//setup for optimization using GA
	params = c_createMatrix(3,3);
	c_setRowName(params,0,"k1");
	c_setRowName(params,1,"k2");
	c_setRowName(params,2,"k3");

	//intial values
	c_setMatrixValue(params, 0, 0, 1);
	c_setMatrixValue(params, 0, 1, 0.0);
	c_setMatrixValue(params, 0, 2, 5.0);

 	//min values
	c_setMatrixValue(params, 1, 0, 1);
	c_setMatrixValue(params, 1, 1, 0.0);
	c_setMatrixValue(params, 1, 2, 5.0);

	//max values
	c_setMatrixValue(params, 2, 0, 1);
	c_setMatrixValue(params, 2, 1, 0.0);
	c_setMatrixValue(params, 2, 2, 5.0);
	
	cSetValue(m1,"k1",2.0);
	cSetValue(m1,"k2",1.0);
	cSetValue(m1,"k3",1.0);
	
	cSetOptimizerIterations(10);  //set num interations
	results = cOptimize(m1, "output.tab", params); //optimize
	c_printMatrixToFile("params.out", results);  //optimized parameters
	c_deleteMatrix(results);
	*/
	
	//cleanup	
	cRemoveModel(m);
	copasi_end();
	printf ("Hit the return key to continue\n");
	getchar();
	return 0;
}

