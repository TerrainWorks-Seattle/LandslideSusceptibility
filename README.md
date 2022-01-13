# LandslideSusceptibility

Resources for creating and assessing statistical models for predicting landslide
susceptibility.

## Installation

1. First, make sure you have R installed. Follow instructions [here](https://www.r-project.org/) to install the latest version. 
1. Download this repository to your computer. There are two options 
for how to do this: 
    1. Via git (you must first have git installed and set up on your machine):
        ```
        git clone https://github.com/TerrainWorks-Seattle/LandslideSusceptibility.git
        ```
    1. Download a zip file of all code
        1. Click on the green button "Code"
        1. Select "Download ZIP"
        1. Download to an appropriate location and unzip
1. To run the toolbox in ArcGIS Pro you must install R-ArcGIS bridge (instructions [here](https://github.com/R-ArcGIS/R-Bridge-Tutorial-Notebooks/blob/master/notebooks/01-basics/R-bridge-install-and-setup.ipynb))
1. You will also need to download the DEMutilities scripts to calculate surface metrics used as model inputs. These are part of the ForestedWetlands repository. Follow the instructions in [WetlandTools_for_ArcGIS.pdf](https://github.com/TerrainWorks-Seattle/ForestedWetlands/blob/master/WetlandTools_for_ArcGIS.pdf) to download and run the DEMutilities toolbox in ArcGIS. 


## Tools
### Landslide Susceptibility
Generates an average probability raster of landslide susceptibility. Input a set of landslide initiation points and a set of rasters to be used as explanatory variables. A specified number of models are created -- each trained on an unchanging set of initiation sites and a randomly-generated set of non-initiation points -- and are used to predict a probability raster for the region.

This method begins by creating a set of non-initation points to train a random forest model. To do this, it first filters the region to only include cells which are within a range of values generally associated with landslides, using the supplied set of initiation points. The resulting set of cells are designated as the "analysis region."

A log file is written with details about the model. If you specified multiple iterations, there will be details about each model fit. Interpretation of these outputs requires some knowledge of how [random forest models](https://en.wikipedia.org/wiki/Random_forest) work. Included model statistics include: 
1. Training error rates for [out-of-bag](https://en.wikipedia.org/wiki/Out-of-bag_error) (OOB), initiation, and non-initiation points. 
1. [Confusion matrix](https://en.wikipedia.org/wiki/Confusion_matrix) 
1. Variable importance using Gini index, essentially a score indicating how much each variable increases the fit of the model 

Each model is fit to the data to derive initiation probability for all raster cells within the analysis region. The resulting probability rasters for each iteration are averaged, and the average probability raster is returned.

### Initiation Point Sensitivity
Assesses how much the choice of initiation points affects the prediction skill of landslide susceptibility models, by comparing the prediction skill of a series of models trained on random sets of initiation sites and an unchanging set of non-initiation sites. 

### Cross-Validation
Performs repeated k-fold spatial cross-validation using different sets of generated non-initiation points. Outputs an average landslide susceptibility raster and a log file which records and summarizes model performance. 

## Example workflow
1. Calculate surface metrics for your DEM using the DEMutilities toolbox from the [ForestedWetlands](https://github.com/TerrainWorks-Seattle/ForestedWetlands) resource. These will be the input rasters for the landslide susceptibility tool. You may experiment with using different surface metrics for inputs or even add your own. For your first iteration, try using the surface metrics tool to calculate surface gradient, plan curvature, profile curvature, and elevation deviation at a 15-meter scale. Follow instructions provided in the ForestedWetlands repository for calculating surface metrics. 
1. Run Landslide Susceptibility tool to derive landslide probabilities in all landslide-prone areas. Use your DEM as the reference raster, and the four derived surface metrics as input rasters. You will also need to input a set of landslide initiation points. For now, specify a buffer radius of 15 meters (you may want a larger radius depending on the resolution of your DEM) a initiation range expansion of 150%, and 5 iterations.
1. Find the output and read the log file to assess the model, and visualize the result by adding the probability raster file to your map. 
1. Assess your landslide susceptibility model using initiation point sensitivity and cross-validation tools. 


## Methods

### Generating non-initiation points
A set of landslide initiation points are required as inputs to the model. However, to train the model, we also need a set of points where landslides did not occur. To generate a meaningful set of points, we first identify conditions which are typical of landslides, i.e. with similar slope characteristics to those points where landslides did occur. This is the "analysis region." By only including the analysis region, the model will then be able to characterize landslide initiation points based on more subtle differences between locations. To identify non-initiation points: 
1. A buffer is drawn around all initiation points and the "initiation range" of all input values for all initiation points is calculated. Only cells with values within the initiation range are included in the "analysis region". This analysis region is used for all futre steps.
1. A pre-specified number of points are sampled from the analysis region but outside the initiation buffers, and a buffer drawn around them to create the non-initiation buffers. 
1. Non-initiation points are identified from within non-initiation buffers using one of several metrics: centroid, maximum gradient, or maximum plan curvature. 

### Model calibration
Data is fit to a random forest model using the [randomForest](https://rdrr.io/cran/randomForest/man/randomForest.html) R package. 
A pre-set number of random forest models are trained on the set of initiation and non-initiation points using values from the input rasters. Each iteration generates a raster of probabilities corresponding to predicted landslide susceptibility and logs model metrics including error rate and confusion matrix. The average probabilities for each cell are reported. 


