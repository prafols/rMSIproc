% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakList2imzML.R
\name{export_imzMLpeakList}
\alias{export_imzMLpeakList}
\title{export_imzMLpeakList.}
\usage{
export_imzMLpeakList(
  peakList,
  posMatrix,
  pixel_size_um,
  filename,
  normalization = rep(1, length(peakList))
)
}
\arguments{
\item{peakList}{an rMSIproc peak list object.}

\item{posMatrix}{an rMSI pos matrix with the pixel coordinates in the image.}

\item{pixel_size_um}{the image pixel size in microns.}

\item{filename}{complete path where the imzML and ibd files will be stored (the .imzML extensions must be omited).}

\item{normalizations}{a numeric vector contaning the normalization value for each pixel (1 by default).}
}
\description{
Export an rMSIproc peak list as an imzML processed dataset.
}
