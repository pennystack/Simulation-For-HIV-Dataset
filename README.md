# 🧬 HIV Dataset Simulation code

<br>

## 📖 Contents
- [📝 Description of the simulation code](#-description-of-the-simulation-code)
- [🗂️ Folder structure](#%EF%B8%8F-folder-structure)
- [💻 How to run the code](#-how-to-run-the-code)
- [🔎 More details about the methodology](#-more-details-about-the-methodology)



<br>

## 📝 Description of the simulation code

 The following code was developed to validate and confirm the accuracy of the estimated model parameters in the non homogeneous semi-Markov model (as described in the *"Parametric and non-homogeneous semi-Markov process for HIV control"* by E.Mathieu, Y.Foucher, P.Dellamonica and JP.Daures) for studying the evolution of the disease in HIV - 1 infected patients. The model is described by the logit-transformed linear transition probabilities (corresponding parameters: $a_{ij}$ and $b_{ij}$) and the Weibull duration times (corresponding parameters: $v_{ij}$ and $s_{ij}$).
 
<br>

> *📌 **Note**: The logit transformation was my addition to address a technical issue with the linear transition probabilities. It was not applied in the non-homogeneous semi-Markov model by Mathieu et al.(2007).*



<br>

## 🗂️ Folder structure

1. Folder **`Parameter estimations`**
   - Files `aij.RData`, `bij.RData`, `sij.RData`, `vij.RData` contain the parameter estimations obtained from the original dataset. These parameters are used for computing the transition probability matrix $P_{ij}$​ and the Weibull duration times, which are used to generate the simulated dataset.

2. Folder **`src`**
   - `05_simulation.R`: Main script where the dataset simulation, parameter estimation on the simulated data and computation of basic statistics for the estimated parameters are performed.
   - `load_functions.R`: Contains functions used for parameter estimation and likelihood maximization.
   - `loglikelihoodnew.cpp`: Contains C++ functions for calculating components used within the likelihood function (e.g., transition probabilities, probability densities, etc.) as well as the likelihood function itself. *Implementing these functions in C++ significantly reduces the runtime of the parameter estimation.*



<br>

## 💻 How to run the code

1. Clone the repository in your computer or download locally the folders **`parameter estimations`**, **`src`** and the file `github upload smm_hiv.Rproj`.
2. Open the `github upload smm_hiv.Rproj` file to load the R environment.
3. Open the R file `05_simulation.R` in the RStudio editor.
4. Run the first 40 lines of the code.
5. When the prompt of defining the numbers of bootstrapping samples that you want to create pops up, you can enter any positive number you wish.
   - For a quick test, you can select to produce 1 or 2 samples *(~ 10-15 minutes runtime)*.
   - For statistically valid results, you should enter 500 or more samples *(Warning! The runtime will be significantly longer)*.
6. Run the code from lines 41-185 to generate the bootstrap samples.
7. Run the remaining code to compute basic statistics about the estimated parameters from the simulated data *(Includes: p-value, t-value, confidence intervals, means)*.

<br>

 > *📌 **Note**: This project was developed using **R version 4.3.1**. Using a different version of R may lead to differences in the estimated parameters from the simulated data, which might not be directly comparable to those obtained from the original data.*



<br>

## 🔎 More details about the methodology

 If you are interested in learning more about this topic, you can find my thesis titled *"Non homogeneous semi-Markov processes with application to HIV"* available in Pergamos, the official unified Institutional Repository/Digital Library of the University of Athens, [here](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://pergamos.lib.uoa.gr/uoa/dl/object/3403042/file.pdf&ved=2ahUKEwjQ7M_MpuSOAxVRIxAIHQVvBBQQFnoECBkQAQ&usg=AOvVaw1tymNuOkbKCGtNwmmVFqkl).
