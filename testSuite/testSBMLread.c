#include <stdlib.h>
#include <stdio.h>
#include "copasi_api.h"

int main(int narg, char** argv)
{
	c_matrix results;
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
	c_printMatrixToFile("results.tab", results);
	c_deleteMatrix(results);

	results = cGetReactionRates(m);

	printf("fluxes:\n");
	c_printOutMatrix(results);
	c_deleteMatrix(results);

	results = cGetRatesOfChange(m);
	printf("\n\nderivatives:\n");
	c_printOutMatrix(results);
	c_deleteMatrix(results);
	
	//cleanup	
	cRemoveModel(m);
	copasi_end();
	printf ("Hit the return key to continue\n");
	getchar();
	return 0;
}

