# Performance-Oriented-DNN-MPC

This repository contains code used in the paper \<insert title\> by \<insert authors\>. The full paper can be found here: \<insert link\>.

## Overview of files
Detailed documentation for the code can be found in the text cells of the Jupyter notebooks. Here, we provide a short overview of the files.

<ul>
  <li> <b>main.ipynb</b>: This is the main file that builds the optimal control problem (OCP) and then solves it in closed-loop, yielding the required model predictive controller (MPC).</li>
 
  <li> <b>case_study </b>: Contains the integrators for the plant and simulation models. </li>
  <li> <b>neural_network.ipynb </b>: Defines the neural_model class which contains all the methods for training the identification residual network. </li>
  
  <li> <b>ResNetExplicit </b>: Uses the neural_network and case_study to generate simulation data based on the HF model and trains the initial predictive model, while saving the initial weights matrix. </li>
 
</ul>
