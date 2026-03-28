# CHiMiC: Causal High-dimensional Mixed-type Clustering

## Description
This package provides the `CHiMiC` function, which performs clustering on mixed-type clinical data $Z$ (including continuous and categorical variables) guided by treatment effects (TE). It also identifies relevant clinical variables using an MCP (minimax concave penalty)–based variable selection method.

## Installation
You can install the development version of this package from GitHub using:
```R
install.packages("devtools")
devtools::install_github("RuiRen420/CHiMiC")
```

## Usage Example
```R
library(CHiMiC)
set.seed(123)
simdata <- simData(n = 1000, s1 = 20, s2 = 5, ratio = 0.2)

Z <- simdata$Z
G <- simdata$G
q1 <- simdata$q1
q2 <- simdata$q2
cate <- simdata$CATE
w <- 0.5
lambda1 <- 50
lambda2 <- 50

#For demonstration purposes, max_iter and inner_iter are set to small values.
#In practice, larger iterations (e.g., 100) are needed for precise estimation.
chimic_obj <- CHiMiC(Z, q1, q2, TE = cate, w, G, lambda1, lambda2,
                     max_iter = 50, inner_iter = 50, tol = 1e-2)
```

## Data Exploration
To see a detailed example of how to explore the treatment effects and selected high-dimensional features (with complex heatmaps and boxplots) after running the model, please view our vignette:
👉 **[Data Exploration with CHiMiC](vignettes/data-exploration.md)**
