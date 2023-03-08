# MarineDNA
A program to analyze eDNA metabarcoding sequence data from multiple trophic levels to identify drivers and patterns of biodiversity 
in marine ecosystems. Created for the 2023 NOAA-NCAR Hackathon for development and optimization. Analysis has three stages: 1) 
modeling sequence read counts as a beta distribution, 2) using random draws from this distrubition to cluster samples and generate Principal Component Analyses, and 3) applying a Bayesian linear regression to identify correlations between two different data sets, generated from two different marker genes targeting different trophic levels.
