# MAMBO: Metabarcoding Analysis using Modeled Bayesian Occurrences
MAMBO is program to analyze eDNA metabarcoding sequence data from multiple marker genes to identify drivers and patterns of biodiversity in ecosystems. The input consists of two metabarcoding data sets generated from two different marker genes targeting different trophic levels from the same eDNA sample. Analysis has three stages: 1) modeling sequence read counts as a beta distribution and using random draws from this distrubition to generate principal component analyses, 2) applying a Bayesian model to identify correlations between two different data sets, and 3) extracting the features (ASVs) with the greatest effect on the correlations. MAMBO originated from a team of NOAA, academica, and private sector scientists who participated in the 2023 NOAA/NCAR Hackathon in Boulder, Colorado, part of the Open Hackathons program.

## Installation
We highly recommend installing MAMBO via Docker. There are several R library dependencies that may otherwise result in conflicts.

These instructions assume you have Docker and Docker Desktop already installed. 

1. Download the Dockerfile provided in this repo. It must be named "dockerfile", with no extension or capitalizations

2. Build the image in the command line. The dockerfile must be in your working directory. 

    ```
# -t tags the image, --build-arg provides the build date (optional)
docker build --build-arg WHEN=2025-01-25 --no-cache -t mambotest:0.1.0 .
```

3. Option #1: Run the image from Docker Desktop. In the "Images" menu item, hit the "play" icon for the MAMBO image. Then, expand the "Optional settings" menu. Under "Ports" type "8787" for "Host port." Under "Volumes" provide the local directory path containing the data files you need to run MAMBO ("Host path") and "/home/rstudio/{optional-name}" in the "Container path," where {optional-name} is a directory name you want for your image working directory. Under "Environment variables" set "Variable" to "PASSWORD" and "Value" to "password" (or your password of choice). Click "run" to start running the image interactively .

4.  Option #2: Run the image from the command line. As with Docker Desktop, set the directory with your test files in /home/rstudio to access them from RStudio.

    ```
docker run --rm -it -e PASSWORD=password -p 8787:8787 \
-v "/Users/Home/MAMBO/test-files":/home/rstudio/mambotest \
mambotest:0.1.0
```

4. Open a browser and type in "localhost:8787" then log in with the user name "rstudio" and password "password" (or whatever you set as the password variable)
5. Load MAMBO with library(mambo)

## Running MAMBO
To run MAMBO, you will need two sets each of 1) ASV tables and 2) taxonomy tables. You may also include an optional metadata file.

The ASV tables should have samples as columns, ASVs as rows, and an empty A1 cell. The cell values should contain raw (untransformed) read counts.

The taxonomy tables should have ASVs as rows and taxonomy ranks in subsequent columns, which each rank in its own column. 

1. In RStudio, import ASV and taxonomy tables 
    ```
asv.16s <- read.csv(file="mambotest/CANON_2018_data/Merged2018_16S_otu_filtered.csv")
asv.18s <- read.csv(file="mambotest/CANON_2018_data/Merged2018_18S_otu_filtered.csv")
tax.16s <- read.csv(file="mambotest/CANON_2018_data/Merged2018_16S_taxa_filtered.csv")
tax.18s <- read.csv(file="mambotest/CANON_2018_data/Merged2018_18S_taxa_filtered.csv")
```

2. Convert ASV tables to data frames

    ```
asv.16s <- data.frame(asv.16s, row.names=1)
asv.18s <- data.frame(asv.18s, row.names=1)
```

3. Run MAMBO

    ```
result <- mambo(
  resp.label = '18s',
  resp.counts = asv.18s,
  pred.label = '16s',
  pred.counts = asv.16s,
  nrep=10,
  chains=3,
  adapt=100,
  burnin=1000,
  total.samples=1000,
  thin=1,
  run.label="mambo",
  output.log=TRUE
)
```

    ```
# to skip the modeling and just run the PCAs, set Bayesian=FALSE
result <- mambo(
  resp.label = '18s',
  resp.counts = asv.18s,
  pred.label = '16s',
  pred.counts = asv.16s,
  nrep=10,
  chains=3,
  adapt=100,
  burnin=1000,
  total.samples=1000,
  thin=1,
  run.label="mambo",
  bayesian=FALSE,
  output.log=TRUE
)
```
4. Summarize results for the first replicate PCA

    ```
# Need those `` symbols
result$reps[[1]]$pca$`16S`
# Individual functions
```

5. Visualize the sign switching

    ```
switchSummary(results = result)
```

6. Plot the PCAs

    ```
plotPCs(result, "16S", pc.x=1, pc.y=2, type="density")
plotPCs(result, "18S", pc.x=1, pc.y=2, type="density")
```

