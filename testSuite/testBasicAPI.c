#include <stdlib.h>
#include <stdio.h>
#include "copasi_api.h"

copasi_model model1(); //oscillation
copasi_model model2(); //positive feebdack gene regulation
copasi_model model3();
void eigen(copasi_model, const char*); //compute eigenvalues by changing parameters (similar to root-locus)

int main()
{
	c_matrix results, results2, jac;
	copasi_model m;
	
	printf("creating model...\n");
	m = model1();

    cWriteSBMLFile (m, "model1.xml");
    
	if (m.errorMessage != NULL) {
		printf ("Errors while reading model:\n");
		printf ("%s\n", m.errorMessage);
		getchar();
		return;
	}
	
	cWriteAntimonyFile (m, "antimony.txt");
	printf ("Antimony script written to antimony.txt\n");

	//simulate
	printf("simulating...\n");	
	results = cSimulateDeterministic(m, 0, 20, 100);  //model, start, end, num. points

	//print results to file
	printf("results.tab has simulation data\n");
	c_printMatrixToFile("results.tab", results);
	
	/** perform additional calculations from the simulated data **/
	
	//get derivative from the simulated results
	results2 = cGetDerivativesFromTimeCourse(m,results);
	printf("derivatives.tab has derivatives computed from the simulation data\n");
	c_printMatrixToFile("derivatives.tab", results2); //print to file
	c_deleteMatrix(results2); //delete matrix

	//get reaction rates from the simulated results
	results2 = cGetReactionRatesFromTimeCourse(m,results);
	printf("rates.tab has reaction rates computed from the simulation data\n");
	c_printMatrixToFile("rates.tab", results2); //print to file
	c_deleteMatrix(results2);  //delete matrix

	//get control coeff. from the simulated results
	results2 = cGetCCFromTimeCourse(m,results);
	printf("controlcoeffs.tab has control coefficients computed from the simulation data\n");
	c_printMatrixToFile("controlcoeffs.tab", results2); //print to file
	c_deleteMatrix(results2); //delete matrix

	//get elasticities for reaction rates from the simulated results
	results2 = cGetElasticitiesFromTimeCourse(m,results);
	printf("elasticities.tab has elasticities computed from the simulation data\n");
	c_printMatrixToFile("elasticities.tab", results2); //print to file
	c_deleteMatrix(results2); //delete matrix
	
	/**  cleanup  **/
	c_deleteMatrix(results); //delete simulation results
	cRemoveModel(m);
	copasi_end();
	printf ("\nHit the return key to continue\n");
	getchar();
	return 0;
}

copasi_model model1() //oscillator
{
	//model named M
	copasi_model model = cCreateModel("M");
	copasi_reaction R1, R2, R3;
	copasi_compartment cell;
	
	//species
	cell = cCreateCompartment(model, "cell", 1.0);
	cCreateSpecies(cell, "A", 4);
	cCreateSpecies(cell, "B", 3);
	cCreateSpecies(cell, "C", 2);
	
	//parameters
	cSetValue(model, "k1", 0.2);   //k1
	cSetValue(model, "k2", 0.5);   //k2
	cSetValue(model, "k3", 1);   //k3
	
	//reactions -- make sure all parameters or species are defined BEFORE this step
	R1 = cCreateReaction(model, "R1");  // A+B -> 2B
	cAddReactant(R1, "A", 1.0);
	cAddReactant(R1, "B", 1.0);
	cAddProduct(R1, "B", 2.0);
	cSetReactionRate(R1, "k1*A*B");

	R2 = cCreateReaction(model, "R2");  //B+C -> 2C
	cAddReactant(R2, "B", 1.0);
	cAddReactant(R2, "C", 1.0);
	cAddProduct(R2, "C", 2.0);
	cSetReactionRate(R2, "k2*B*C");

	R3 = cCreateReaction(model, "R3"); //C+A -> 2A
	cAddReactant(R3, "C", 1.0);
	cAddReactant(R3, "A", 1.0);
	cAddProduct(R3, "A", 2.0);
	cSetReactionRate(R3, "k3*C*A");

	cCreateEvent(model, "event1", "time > 10", "k3", "k3/2.0");

	cCompileModel(model); //must be done after creating a model, before analysis

	return model;
}

copasi_model model2() //gene regulation
{
	//model named M
	copasi_model model = cCreateModel("M");
	copasi_compartment cell;
	copasi_reaction R1, R2, R3, R4;
	
	//species
	cell = cCreateCompartment(model, "cell", 1.0);
	cCreateSpecies(cell, "mRNA", 0);
	cCreateSpecies(cell, "Protein", 0);
	
	//parameters	
	cSetValue(model, "d1", 1.0);
	cSetValue(model, "d2", 0.2);  
	cSetValue(model, "k0", 2.0);
	cSetValue(model, "k1", 1.0);
	cSetValue(model, "h", 4.0);  
	cSetValue(model, "Kd", 1.0);
	cSetValue(model, "leak", 0.1);  
	
	//reactions -- make sure all parameters or species are defined BEFORE this step
	R1 = cCreateReaction(model, "R1");  //  mRNA production
	cAddProduct(R1, "mRNA", 1.0);
	cSetReactionRate(R1, "leak + k0 * (Protein^h) / (Kd + (Protein^h))");

	R2 = cCreateReaction(model, "R2");  // Protein production
	cAddProduct(R2, "Protein", 1.0);
	cSetReactionRate(R2, "k1*mRNA");

	R3 = cCreateReaction(model, "R3"); // mRNA degradation
	cAddReactant(R3, "mRNA", 1.0);
	cSetReactionRate(R3, "d1*mRNA");
	
	R4 = cCreateReaction(model, "R4"); // Protein degradation
	cAddReactant(R4, "Protein", 1.0);
	cSetReactionRate(R4, "d2*Protein");
	return model;
}

// eigenvalues
void eigen(copasi_model model, const char* param)
{
	int i, j, k;
	double p;
	FILE * outfile;
	c_matrix ss;
	c_matrix output;
	
	i = 0; j = 0; k = 0;
	//steady states
	
	for (i=0; i < 100; ++i)
	{
		p = (double)(i + 1)/10.0;
		k = cSetValue( model, param, p );
		
		if (k)
			printf("calculating steady state for %s = %lf\n",param, p);
		
		ss = cGetEigenvalues(model);
		//ss = cGetSteadyState(model);

		if (i == 0)
		{
			output = c_createMatrix(100, ss.rows+1);
			c_setColumnName(output, 0, param);
			for (j=0; j < output.cols; ++j)
				c_setColumnName(output, j+1, c_getRowName(ss, j));
		}
		
		c_setMatrixValue(output, i, 0, p);
		for (j=0; j < output.cols; ++j)
			c_setMatrixValue(output, i, j+1, c_getMatrixValue(ss, j, 0));
		
		c_deleteMatrix(ss);
	}
	
	//output
	c_printMatrixToFile("output.tab", output);
	
	printf("\noutput.tab contains the final output\n\n");

	c_deleteMatrix(output);
}

copasi_model model3() //big genetic model
{
	copasi_model model = cCreateModel("M");
	copasi_compartment DefaultCompartment;
	copasi_reaction r0,r1,r2,r3;
	DefaultCompartment = cCreateCompartment(model,"DefaultCompartment",1);
	cCreateSpecies(DefaultCompartment,"dr1_Monomer",0);
	cCreateSpecies(DefaultCompartment,"rs2",1);
	cCreateSpecies(DefaultCompartment,"cod1",1);
	cCreateSpecies(DefaultCompartment,"OUTPUT",5);
	cCreateSpecies(DefaultCompartment,"cod2",1);
	cCreateSpecies(DefaultCompartment,"INPUT",5);
	cCreateSpecies(DefaultCompartment,"as1",1);
	cCreateSpecies(DefaultCompartment,"as2",1);
	cCreateSpecies(DefaultCompartment,"pro1",1);
	cCreateSpecies(DefaultCompartment,"pro2",1);
	cCreateSpecies(DefaultCompartment,"rbs1",1);
	cCreateSpecies(DefaultCompartment,"rbs2",1);
	cCreateSpecies(DefaultCompartment,"ter1",1);
	cCreateSpecies(DefaultCompartment,"ter2",1);
	cSetGlobalParameter(model,"OUTPUT_degradation_rate",0.1);
	cSetGlobalParameter(model,"dr1_degradation_rate",0.1);
	cSetGlobalParameter(model,"dr1_Kd",12);
	cSetGlobalParameter(model,"dr1_h",4);
	cSetGlobalParameter(model,"pro1_strength",5);
	cSetGlobalParameter(model,"pro2_strength",12);
	cSetGlobalParameter(model,"ta1_Kd",5);
	cSetGlobalParameter(model,"ta1_h",4);
	cSetGlobalParameter(model,"ta2_Kd",2);
	cSetGlobalParameter(model,"ta2_h",5);
	cSetAssignmentRule(model, "INPUT","10 * (1 + sin(time * 0.5))");
	cSetAssignmentRule(model, "as1","((1+((INPUT/ta1_Kd)^ta1_h))-1)/((1+((INPUT/ta1_Kd)^ta1_h)))");
	cSetAssignmentRule(model, "as2","((1+((INPUT/ta2_Kd)^ta2_h))-1)/((1+((INPUT/ta2_Kd)^ta2_h)))");
	cSetAssignmentRule(model, "cod1","pro1_strength * (as1)");
	cSetAssignmentRule(model, "cod2","pro2_strength * (( as1 + as2) *(rs2))");
	cSetAssignmentRule(model, "rs2","1/(dr1_Kd+dr1_Monomer^dr1_h)");
	r0 = cCreateReaction(model, "dr1_v1");
	cSetReactionRate(r0,"cod1");
	cAddProduct(r0,"dr1_Monomer",1);
	r1 = cCreateReaction(model, "dr1_v2");
	cSetReactionRate(r1,"dr1_degradation_rate*dr1_Monomer");
	cAddReactant(r1,"dr1_Monomer",1);
	r2 = cCreateReaction(model, "pp1_v1");
	cSetReactionRate(r2,"cod2");
	cAddProduct(r2,"OUTPUT",1);
	r3 = cCreateReaction(model, "pp1_v2");
	cSetReactionRate(r3,"OUTPUT_degradation_rate*OUTPUT");
	cAddReactant(r3,"OUTPUT",1);
	return model;
}




