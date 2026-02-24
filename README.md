# Simulation Study

This git repo contains 7 scripts that together form a simulation engine. Here are their contents:
*`00_packages.R`: Script that loads all dependencies, except for the `here` package. 
*`01_population_data_generation`: A function that generates data according to the following Data Generating Mechanism. 

Exogenous variables:

- \(Z \sim \mathcal{N}(0,1)\)
- \(U \sim \mathcal{N}(0,1)\)
- \(Z \perp U\)

Structural equations:

- \(X = b_2 Z + b_5 U + e_x\), with  
  \(\mathrm{Var}(e_x)= 1-(b_2^2+b_5^2)\)

- Define the **centered interaction** term:
  \[
  XZ_c = XZ - b_2
  \]
  (since \(\mathbb{E}[XZ] = \mathrm{Cov}(X,Z)=b_2\) under this DGM)

- \(Y = b_1 X + b_6 XZ_c + b_3 Z + b_4 U + e_y\)

Error variances (to keep \(\mathrm{Var}(X)=\mathrm{Var}(Y)=1\)):

- \(\mathrm{Var}(e_y)= 1 - \Big(b_1^2 + b_3^2 + b_4^2 + 2b_1b_3b_2 + 2b_1b_4b_5 + b_6^2(1+b_2^2)\Big)\)

Notes:

- Effect modification: the conditional marginal effect of \(X\) on \(Y\) is
  \[
  \frac{\partial Y}{\partial X} = b_1 + b_6 Z
  \]
- Average marginal effect over \(Z\):
  \[
  \mathbb{E}\!\left[\frac{\partial Y}{\partial X}\right]=b_1 \quad \text{because } \mathbb{E}[Z]=0
  \]
- The term \(b_6^2(1+b_2^2)\) comes from \(\mathrm{Var}(XZ_c)=\mathrm{Var}(XZ)=1+b_2^2\) for jointly normal, mean-0, var-1 \(X,Z\) with \(\mathrm{Cov}(X,Z)=b_2\).

