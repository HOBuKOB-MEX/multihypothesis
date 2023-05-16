# multihypothesis
An R script for construction and performance evaluation of sequential   multi-hypothesis tests for the Bernoulli model

This repository contains an R module for  design and performance evaluation of sequential multi-hypothesis tests.

Currently, only the model of binary responses (sampling from a Bernoulli population) is covered.

The detailed description of the method and algorithms can be found in

[*Andrey Novikov. A Numerical Approach to Sequential Multi-Hypothesis Testing for Bernoulli Model, 2023*](https://doi.org/10.1080/07474946.2023.2215825)

## Content description
* The file [multihypothesis.R](multihypothesis.R) contains all the functions providing the  user interface for all the tasks.

The list of functions can be seen below. 

### OptTest

The function for designing an optimal truncated sequential multi-hypothesis test.

Arguments:
* _H_ horizon (maximum number of steps to employ)
* _lam_ matrix of Lagrange multipliers
* _th_ vector of the hypothesized values of the parameter (success probability)
* _gam_ the vector of weights  used to calculate the weighted expected sample size (ESS) 
* _thgam_ parameter values at which ESSs are calculated to  used in the weighted ESS 

Returns:  the designed test which minimizes the weighted ESS

### MSPRT

The function for designing a truncated  MSPRT (see description in [*the cited article*](https://doi.org/10.48550/arXiv.2212.05151)).

Arguments:
* _H_ horizon (maximum number of steps to employ)
* _lA_ matrix of the logarithmic thresholds
* _th_ vector of the hypothesized values of parameter (success probability)

Returns:  the designed truncated MSPRT

### SimplifiedOpt

The function for designing the simplified version of the optimal truncated sequential multi-hypothesis test.
(To be published soon)

Arguments:
* _H_ horizon (maximum number of steps to employ)
* _lam_ matrix of Lagrange multipliers
* _th_ the vector of the hypothesized values of parameter (success probability)
* _gam_ the vector of weights  used to calculate the weighted expected sample size (ESS) 
* _thgam_ parameter values at which ESSs are calculated to  be used in the weighted ESS 

Returns:  the designed test which approximately minimizes the weighted ESS

### PAccept 

The function calculating the probability to accept a hypothesis given a value of the parameter. 

Arguments:
* _test_ test designed by any of the functions _OptTest, SimplifiedOpt, MSPRT_ 
* _th_ parameter value at which the probability to accept the hypothesis is calculated
* _i_ number of the hypothesis whose acceptance probability is evaluated

Returns: the probability to accept

### ESS

The function calculating the expected sample size of a test, given a value of the parameter.

Arguments:
* _test_ test designed by any of the functions _OptTest, SimplifiedOpt, MSPRT_
* _th_ parameter value at which the ESS is calculated

Returns: ESS

### prob_to_stop_after

The function calculating the probability that the test stops after a given stage number.

Arguments:
* _test_ test designed by any of the functions _OptTest, SimplifiedOpt, MSPRT_
* _th_ parameter value at which the probability is calculated
* _k_ the stage number

Returns: the probability

### monte_carlo_simulation

The function for carrying out Monte Carlo simulations

Arguments:

* _test_ test designed by any of the functions _OptTest, SimplifiedOpt, MSPRT_
* _hyp_ parameter value for the simulation
* _K_  number of hypotheses
* _nMC_ number of replications for the simulation

Returns: rates of acceptations and their standard errors, ESS and its standard error
         
### maxNumber

The function to determine the maximum number of stages the test uses

Arguments:

* _test_ test designed by any of the functions _OptTest, SimplifiedOpt, MSPRT_

Returns: maximum number of stages

## Usage example

* The file [usage.R](usage.R) is a usage example of these functions 

