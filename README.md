# Simulation Study

This git repo contains 7 scripts that together form a simulation engine. Here are their contents:
*`00_packages.R`: Script that loads all dependencies, except for the `here` package. 
*`01_population_data_generation.R`: A function that generates data according to the following Data Generating Mechanism. 

- **Z** ~ Normal(0, 1)  
- **U** ~ Normal(0, 1)  
- **Z** independent of **U**

- **X** = b2·Z + b5·U

- Define the centered interaction:
  - **XZc** = X·Z − b2

- **Y** = b1·X + b6·XZc + b3·Z + b4·U

Residual variances are chosen such that each variable has a variance of 1. 

*`02_sampling.R`: Contains two sampling functions. The first pulls a simple random sample from the population data, but excludes variable U. The second pulls a Non-Probability Sample from the population, including all variables, according to the following model: 

Given a population dataset with variables **Z, U, X, Y**, define a selection score for unit *i*:

- **Selection index**:  
  \[
  \eta_i = a_X X_i + a_Y Y_i + a_Z Z_i + a_U U_i
  \]

Convert this to positive sampling weights:

- **Weights** (monotone in \(\eta_i\)):  
  \[
  w_i = \exp(\eta_i)
  \]

Draw an NPS sample of size **n** from the population of size **N**:

- **Sampling rule**: sample **exactly n units without replacement** with inclusion probabilities proportional to \(w_i\):  
  \[
  \Pr(i \text{ selected}) \propto w_i
  \]

So units with larger \(\eta_i\) (higher \(X, Y, Z,\) and/or \(U\) depending on \(a\)) are more likely to be selected, creating selection bias when any \(a\neq 0\).