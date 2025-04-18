% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_functions.R
\name{do_efficient_aipw}
\alias{do_efficient_aipw}
\title{Function for efficient AIPW estimator}
\usage{
do_efficient_aipw(
  data,
  models,
  G_name = "G",
  V_name = "V",
  Y_name = "Y",
  return_se = FALSE
)
}
\arguments{
\item{data}{dataset to predict on}

\item{models}{list of pre-fit models needed for estimation}

\item{G_name}{name of growth outcome variable, default G}

\item{V_name}{name of vaccine treatment variable, default V}

\item{Y_name}{name of infection variable, default Y}

\item{return_se}{flag to return standard error, defualt FALSE}
}
\value{
AIPW estimate of growth effect (+ standard error if return_se = TRUE)
}
\description{
Function for efficient AIPW estimator
}
