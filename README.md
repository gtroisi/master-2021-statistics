# Sparse Inference applied to genes for cancer screening.

## Introduction

In this repository, a method for cancer screening was developed. This evaluates 6034 different genes from 102 people. 

## Contents

The code was implemented in R and is made by three steps:
- p values calculation
- MTP for variable reducing (Bonferroni, Sidak, Hochberg et al.)
- Sparse Inference and Lasso Regression for final variable extraction.

## Conclusion
From 6034 variables, only 35 are significant for cancer screening with an accuracy of 97%.
