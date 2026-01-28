# MAMBO: Metabarcoding Analysis using Modeled Bayesian Occurrences
MAMBO is program to analyze eDNA metabarcoding sequence data from multiple marker genes to identify drivers and patterns of biodiversity in ecosystems. The input consists of two metabarcoding data sets generated from two different marker genes targeting different trophic levels from the same eDNA sample. Analysis has three stages: 1) modeling sequence read counts as a beta distribution and using random draws from this distrubition to generate principal component analyses, 2) applying a Bayesian model to identify correlations between two different data sets, and 3) extracting the features (ASVs) with the greatest effect on the correlations. MAMBO originated from a team of NOAA, academic, and private sector scientists who participated in the 2023 NOAA/NCAR Hackathon in Boulder, Colorado, part of the Open Hackathons program.

####There are two options for installing MAMBO: natively in R, or via Docker. 

## Native installation
In R, type the following command into the terminal:

```
devtools::install_github('nvpatin/mambo/mambo')
```

## Docker installation

We highly recommend installing MAMBO via Docker, which will avoid any potential R library dependencies that may result in conflicts.

These instructions assume you have Docker already installed. Docker Desktop is optional, and instructions for loading the image both with and without the app are below.

### Building the Docker image

1. Clone this repo to your local workspace. 

2. Navigate to the MAMBO/Docker directory and locate the Docker file. It is named "dockerfile", with no extension or capitalizations

3. Build the image in the command line. The dockerfile must be in your working directory. 

```
# -t tags the image, --build-arg provides the build date (optional)
docker build --build-arg WHEN=2025-01-25 --no-cache -t mambo:1.0.3 .
```

### Loading MAMBO from the Docker image
3. Option #1: Run the image from Docker Desktop. In the "Images" menu item, hit the "play" icon for the MAMBO image. Then, expand the "Optional settings" menu. Under "Ports" type "8787" for "Host port." Under "Volumes" provide the local directory path containing the data files you need to run MAMBO (e.g., the "Data" folder in this repo; "Host path") and "/home/rstudio/{optional-name}" in the "Container path," where {optional-name} is a directory name you want for your image working directory. Under "Environment variables" set "Variable" to "PASSWORD" and "Value" to "password" (or your password of choice). Click "run" to start running the image interactively .

4.  Option #2: Run the image from the command line. As with Docker Desktop, set the directory with your test files in /home/rstudio to access them from RStudio.

```
docker run --rm -it -e PASSWORD=password -p 8787:8787 \
-v "/Users/Home/MAMBO/files":/home/rstudio/mambo \
mambotest:1.0.3
```

5. Open a browser and type in "localhost:8787" then log in with the user name "rstudio" and password "password" (or whatever you set as the password variable)
6. Load MAMBO with library(mambo)

## Running MAMBO
To run MAMBO using your own data, you will need two sets each of 1) ASV tables and 2) metadata files. You may also include optional taxonomy tables. The ASV tables should have samples as columns, ASVs as rows, and an empty A1 cell. The cell values should contain raw (untransformed) read counts. The taxonomy tables should have ASVs as rows and taxonomy ranks in subsequent columns, which each rank in its own column. 

See `mamboTutorial()` for examples of running the model with data used in the manuscript. See below for a conceptual overview of the MAMBO workflow.

![alt text](https://github.com/nvpatin/MAMBO/blob/main/FigureS1.png?raw=true)

## Troubleshooting

Both JAGS and its R interface, rjags, should be installed with the 'mambo' library. However, occasionally this installation fails. If you get a warning saying JAGS can't be fund, even after loading the 'rjags' library manually, you may have to install JAGS from source: https://sourceforge.net/projects/mcmc-jags/