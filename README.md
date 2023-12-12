# Biostat777 Project 3 Part 1


## Project Description
See: https://www.stephaniehicks.com/jhustatprogramming2023/projects/2023-11-28-project-3/


## Package
STTileMaker


## Author
Caleb Hallinan


## Goal of Package
Create image tiles of each spot in a Spatial Transcriptomics Dataset


## Functions
STTileMaker


## Example

#### Make Tiles of each Spot in a Spatial Transcriptomics Dataset

directory = here("data/f12hr_140_processed/outs/")

tiles = STTileMaker(directory)


## How to Download

#### Please run the following commands in your R terminal to install and use these functions

remotes::install_github(repo = "jhu-statprogramming-fall-2023/biostat777-project3-part1-calebhallinan")

library(STTileMaker)

#### OR

Clone this repository and build the package from within R

## Notes/Thoughts:

- This was a cool project! I am glad I was able to create a function that I can use for my research :)
