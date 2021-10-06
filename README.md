# Multi-state models and joint models: a comparison using AIDS data

This project was submitted to University College London in part fullfillment of the requirements for the MSc Data Science degree.

Survival analysis is often used to investigate the association between a time-dependent biomarker and time of death. An example of this is the association between CD4 counts and survival of AIDS patients. The CD4 count is a test that measures the number of CD4 cells in blood.

Multi-state modelling of CD4 count and survival consists of defining several living states corresponding to various levels of CD4 count, and an additional death state. A model for the transitions between the states can then be defined to investigate the joint process of change in CD4 count and survival.

Joint modelling of CD4 count and survival can be defined by linking a model for the change of CD4 count with a model for survival. The link is based on the use of random effects. 

Multi-state models provide a substantial opportunity to gain a deeper understanding of the disease process, and how it can change over time. This project aims to introduce this framework and demonstrate its use by applying it to a dataset. It is also briefly compared with the joint modelling framework. The full report can be found [here](https://github.com/tjgoh/aids-survival-analysis/blob/main/thesis-write-up.pdf).

