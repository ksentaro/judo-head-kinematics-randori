# judo-head-kinematics-randori

Analysis code for:
"Skill-dependent modulation of rotational head velocity
in judo randori"
*Sports Biomechanics* (under review)

## Requirements
- R version 4.4.1
- Packages: glmmTMB, emmeans, ggplot2

## Directory Structure
code/
  01_glmm_analysis.R  # Data loading, AIC selection, model fitting
  02_posthoc.R        # EMM, pairwise comparisons, ICC
  03_figure1.R        # Figure 1 generation
data/                 # Data files (not publicly available)
output/
  models/             # Saved model objects (.rds)
  figures/            # Generated figures

## Usage
Run scripts in order from project root:
  Rscript code/01_glmm_analysis.R
  Rscript code/02_posthoc.R
  Rscript code/03_figure1.R

## Data Availability
Raw data are not publicly available due to ongoing
secondary analyses. Data requests can be directed
to the corresponding author.

## License
MIT License
