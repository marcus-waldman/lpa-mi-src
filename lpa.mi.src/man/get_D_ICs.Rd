% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_D_ICs.R
\name{get_D_ICs}
\alias{get_D_ICs}
\title{Get D statistic and D-statistic based difference ICs.}
\usage{
get_D_ICs(LL_df, M, n)
}
\arguments{
\item{LL_df}{(data.frame) with the number of p}

\item{M}{(integer) Number of imputation}

\item{n}{(integer) Sample size}
}
\value{
out_list  (list) with the following elements:
       (A) k_chosen - (data.frame) with the selected model based on ICs with D statistic
       (B) Rel_IC -  (data.frame) with relative information criteria (relative to first class IC values)
       (C) LL - (data.frame) LL_df updated with D statistic, D ICs, and ARIV estimate
}
\description{
Obtain the population parameters given the data
}
\examples{
get_D_ICs(LL_df, M, n)
}
