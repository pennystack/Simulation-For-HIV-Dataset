## Description of the simulation code

 The following code was developed to validate and confirm the accuracy of the estimated model parameters in the non homogeneous semi-Markov model (as described in the *"Parametric and non-homogeneous semi-Markov process for HIV control"* by E.Mathieu, Y.Foucher, P.Dellamonica and JP.Daures) for studying the evolution of the disease in HIV - 1 infected patients. The model is described by the logit-transformed linear transition probabilities (corresponding parameters: $a_{ij}$ and $b_{ij}$) and the weibull duration times (corresponding parameters: $v_{ij}$ and $s_{ij}$).
 

 ###### <small> *Note: The logit transformation was my addition to address a technical issue with the linear transition probabilities. It was not applied in the non-homogeneous semi-Markov model by Mathieu et al.(2007).*
