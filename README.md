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
Generates an average probability raster of landslide susceptibility. A specified number of models are created -- each trained on an unchanging set of initiation sites and a randomly-generated set of non-initiation points -- and are used to predict a probability raster for the region. The probability rasters are averaged together and the results is returned. 

### Initiation Point Sensitivity
Assesses how much the choice of initiation points affects the prediction skill of landslide susceptibility models, by comparing the prediction skill of a series of models trained on random sets of initiation sites and an unchanging set of non-initiation sites. 

### Cross-Validation
Performs repeated k-fold spatial cross-validation using different sets of generated non-initiation points. Outputs an average landslide susceptibility raster and a log file which records and summarizes model performance. 

## Example workflow
1. Calculate surface metrics for your DEM
1. Run Landslide Susceptibility tool to derive landslide probabilities in all landslide-prone areas
1. Generate a proportion raster to identify regions where different proportions of landslide may occur
1. Assess your landslide susceptibility model using initiation point sensitivity and cross-validation


## Methods

### Generating non-initiation points
A set of landslide initiation points are required as inputs to the model. However, to train the model, we also need a set of points where landslides did not occur. To generate a meaningful set of points, we first identify conditions which are typical of landslides, i.e. with similar slope characteristics to those points where landslides did occur. The model will then be able to characterize landslide initiation points based on more subtle differences between locations. To identify non-initiation points: 
1. A buffer is drawn around all initiation points and the range of all input values for all initiation points is calculated. Only cells with values within the resulting range are included in future steps. 
1. A pre-specified number of points are sampled from the resulting set of cells which are within the resulting range of values but outside the initiation buffers, and a buffer drawn around them to create the non-initiation buffers. 
1. Non-initiation points are identified from within non-initiation buffers using one of several metrics: centroid, maximum gradient, or maximum plan curvature. 

### Model calibration
Data is fit to a random forest model using the [randomForest](https://rdrr.io/cran/randomForest/man/randomForest.html) R package. 
A pre-set number of random forest models are trained on the set of initiation and non-initiation points using values from the input rasters. Each iteration generates a raster of probabilities corresponding to predicted landslide susceptibility and logs model metrics including error rate and confusion matrix. The average probabilities for each cell are reported. 


