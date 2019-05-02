# rMSIproc
*rMSIproc* is an open-source R package for mass spectrometry (MS) imaging data pre-processing. The package is a multi-platform tool that has been tested on Linux and Windows systems. It provides optimized routines for a complete MSI pre-processing pipeline including: spectral smoothing, spectral alignment, mass calibration, intensity normalization, peak-picking and peak-binning. *rMSIproc* takes an MSI dataset as input and generates a processed MSI dataset and a peak matrix as output. The supported input formats are: imzML, *rMSI* format (TAR) and XMASS. The output MSI dataset is stored in *rMSI* format (TAR) and the peak matrix is stored in a format readable by *rMSIproc*.
*rMSIproc* complements the previously released *rMSI* package. The *rMSI* package was designed to allow efficient access to large MSI datasets combined with a data visualization graphical user interface (GUI). *rMSIproc* takes advantage of the *rMSI* data handling strategy and adds a full data pre-processing pipeline designed to extract relevant *m/z* features from large datasets. 
*rMSIproc*’s internal processing is implemented in C++ to provide efficient memory management and highly optimized multi-threading execution. However, all the user-relevant methods are exposed as R functions following the classical structure of an R package. The package also provides a graphical user interface (GUI) to facilitate setting up the MS processing parameters.

### Installation
*rMSIproc* needs the *rMSI* package to access MS data. Before installing *rMSIproc*, be sure to have rMSI properly installed. Please, check out *rMSI* github page for instructions of rMSI installation: <https://github.com/prafols/rMSI>.
For Linux platform, you’ll need to install the following libraries using your distribution package manager: fftw and boost. For Windows systems, these libraries are provided as re-distributable binaries inside the *rMSIproc* package, so no extra step is needed. 
To install *rMSIproc* you can use the devtools package to fetch it directly from github:
```R
> devtools::install_github("prafols/rMSIproc", ref = "0.2")
```

### Basic usage
The simplest way to use *rMSIproc* is by using its GUI. Just run the following command to run it:
```R
> rMSIproc::ProcessWizard()
```
A window will appear where you can set up all the processing parameters, the input data and the output directory to store the results. When the processing finishes, use rMSIproc to load a peak matrix in to your R session.
```R
> myPeakMat <- rMSIproc::LoadPeakMatrix("/path/to/your/peakmatrix.zip")
```
A more detailed tutorial that describes basic R, *rMSI* and *rMSIproc* usage is available in pdf at the following link:
<https://github.com/prafols/rMSIproc/blob/master/tutorial/rMSI_rMSIproc_tutorial.pdf>

