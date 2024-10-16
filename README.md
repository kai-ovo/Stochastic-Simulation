# Stochastic Simulation
This repository contains notebook demonstrations of numerous numerical methods used in stochastic modeling and simulations for scientific problems. Detailed documentations are provided along. Algorithms include 
- Rejection sampling 
- Importance sampling
- Self-normalized importance sampling
- The Cross-Entropy method
- Control variate
- Large-Deviation-Theory-based importance sampling
- Polynomial chaos expansions
- The Nyström method
- Least-squares polynomial approximation
- Pseudo-spectral methods
- Methods for computing Karhunen-Loeve expansions
- Sparse grids
- Sensitivity analysis
- Sobol indices
- ANOVA
- MCMC
- Adaptive Metropolis-Hastings
- Kalman filter
- Ensemble Kalman filter

## Approximate Bayesian Computations
The directory `ABC` contains a notebook demonstration on various advanced algorithms for approximate Bayesian computations (ABC). In particular, the following contents are included in `abc.ipynb`.

The MCMC algorithm and the adaptive Metropolis-Hastings algorithm are applied to solve a representative Bayesian Inverse problem. Mathematical derivations and formulations are walked through in detail. Various evaluation metrics for MCMC algorithms, such as the Integrated Autocorrelation in Time (IAT), are used to judge the algorithmic performances. Numerous common visualization techniques in Bayesina computations are adopted.

The standard Kalman filter is mathematically derived. It is applied on the same Bayesian inverse problem as above. The Ensemble Kalman filter is formulated in detail, and it is tested on a representative data assimilation problem: state estimation in the chaotic Lorentz-63 system.

## Large Deviation Theory
The directory `LDT` contains the notebook `ldt.ipynb` to demonstrate the large-deviation-theory-based adaptive importance sampling method applied on the short column design problem considered in S. Tong, A. Subramanyam, V. Rao, Optimization under rare chance constraints, SIAM Journal on Optimization 32.2:930-958 (2022).

## Forward Uncertain Quantification (Forward UQ)
The notebook `fuq.ipynb` under the directory `FUQ` contains numerous popular advanced algorithms for the forward UQ problem. The forward UQ problem is also known as simluations of stochastic differential equations. 

The notebook includes detailed derivation of the polynomial chaos expansion (PCE)algorithm with Hermite polynomials, which can also be modified into Legendre/Laguerre polynomials readily. The Nyström method is formulated and implemented. Methods for computing the Sobol indices are implemented to conduct detailed sensitivity analysis on using PCE to simulate a stochastic diffusion equation.

## Monte Carlo Sampling and Variance Reduction Schemes
The directory `MC` contains two notebook demonstrations of several basic Monte Carlo schemes along with more advanced variance-reducing techniques. 

In `mc.ipynb`, the rejection sampling and self-normalized are applied to sample heavy-tailed distributions (distributions without higher moments). 

In `vr.ipynb`, the control variate method is implemented and applied on sampling from the Cauchy distirbution to showcase its performance. The cross-entropy method is implemented and applied in the context of conducting rare-event simulations on heavy-tailed distributions





