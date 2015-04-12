#define S_FUNCTION_NAME ode_model


#include "simstruc.h"
#include "stdio.h"
#include "math.h"

/*
 * mdlInitializeSizes - initialize the sizes array
 *
 * The sizes array is used by SIMULINK to determine the S-function block's
 * characteristics (number of inputs, outputs, states, etc.).
 *
 * The direct feedthrough flag can be either 1=yes or 0=no. It should be
 * set to 1 if the input, "u", is used is the mdlOutput function. Setting this
 * to 0 is akin to making a promise that "u" will not be used in the mdlOutput
 * function. If you break the promise, then unpredictable results will occur.
 */
static void mdlInitializeSizes(SimStruct *S)
{

#include "transient/auto_initial_sizes.h"
}

/*
 * mdlInitializeSampleTimes - initialize the sample times array
 *
 * This function is used to specify the sample time(s) for your S-function.
 * You must register the same number of sample times as specified in 
 * ssSetNumSampleTimes. If you specify that you have no sample times, then
 * the S-function is assumed to have one inherited sample time.
 *
 * The sample times are specified period, offset pairs. The valid pairs are:
 *
 *   [CONTINUOUS_SAMPLE_TIME   0     ]  : Continuous sample time.
 *   [period                   offset]  : Discrete sample time where
 *					  period > 0.0 and offset < period or
 *				          equal to 0.0.
 *   [VARIABLE_SAMPLE_TIME     0     ]  : Variable step discrete sample time
 *				          where mdlGetTimeOfNextVarHit is
 *				          called to get the time of the next
 *					  sample hit.
 *
 *  or you can specify that the sample time is inherited from the driving
 *  block in which case the S-function can have only one sample time:
 *    [INHERITED_SAMPLE_TIME    0    ]
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    
    /*
     * SET OTHER SAMPLE TIMES AND OFFSETS HERE
     */
}

/*
 * mdlInitializeConditions - initialize the states
 *
 * In this function, you should initialize the continuous and discrete
 * states for your S-function block.  The initial states are placed
 * in the x0 variable.  You can also perform any other initialization
 * activities that your S-function may require
 */

static void mdlInitializeConditions(double *x0, SimStruct *S)
{
  int i=0;
  while( i < ssGetNumContStates(S)){
    x0[i] =     0.0;
    i = i + 1;
  }
}


/*
 * mdlOutputs - compute the outputs
 *
 * In this function, you compute the outputs of your S-function
 * block.  The outputs are placed in the y variable.
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{

/* begin common block */
#include "transient/auto_common_block.h"
/* end common block */

#include "transient/auto_outputs.h"
/*
*/

}


/*
 * mdlUpdate - perform action at major integration time step
 *
 * This function is called once for every major integration time step.
 * Discrete states are typically updated here, but this function is useful
 * for performing any tasks that should only take place once per integration
 * step.
 */
static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
    /*
     * YOUR CODE GOES HERE
     */
}

/*
 * mdlDerivatives - compute the derivatives
 *
 * In this function, you compute the S-function block's derivatives.
 * The derivatives are placed in the dx variable.
 */

static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)
{

/* begin common block */
#include "transient/auto_common_block.h"
/* end common block */


#include "transient/auto_odes.h"

#include "transient/auto_remap_odes.h"


}

/*
 * mdlTerminate - called when the simulation is terminated.
 *
 * In this function, you should perform any actions that are necessary
 * at the termination of a simulation.  For example, if memory was allocated
 * in mdlInitializeConditions, this is the place to free it.
 */

static void mdlTerminate(SimStruct *S)
{
    /*
     * 
     */
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
