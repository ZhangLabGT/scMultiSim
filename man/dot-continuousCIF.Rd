% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_main.R
\name{.continuousCIF}
\alias{.continuousCIF}
\title{Generates cifs for cells sampled along the trajectory of cell development}
\usage{
.continuousCIF(
  seed,
  N,
  options,
  ncell_key = "cell",
  is_spatial = FALSE,
  spatial_params = NULL,
  .plot = FALSE,
  .plot.name = "cont_cif.pdf"
)
}
\arguments{
\item{seed}{random seed}

\item{N}{the number list}

\item{options}{the option list}

\item{ncell_key}{the key for the number of cells in N}

\item{is_spatial}{return a list of cifs for spatial}

\item{spatial_params}{the spatial parameters}

\item{.plot}{save the CIF plot}

\item{.plot.name}{plot name}
}
\value{
a list containing the cif and meta data
}
\description{
Generates cifs for cells sampled along the trajectory of cell development
}
