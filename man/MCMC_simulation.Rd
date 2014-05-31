\name{MCMC_simulation}
\alias{MCMC_simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MCMC_simulation(n_sim, pattern, theta_init, overlap, cluster_coords, p_moves_orig, K, lkhd_z, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_sim}{
%%     ~~Describe \code{n_sim} here~~
}
  \item{pattern}{
%%     ~~Describe \code{pattern} here~~
}
  \item{theta_init}{
%%     ~~Describe \code{theta_init} here~~
}
  \item{overlap}{
%%     ~~Describe \code{overlap} here~~
}
  \item{cluster_coords}{
%%     ~~Describe \code{cluster_coords} here~~
}
  \item{p_moves_orig}{
%%     ~~Describe \code{p_moves_orig} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{lkhd_z}{
%%     ~~Describe \code{lkhd_z} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n_sim, pattern, theta_init, overlap, cluster_coords, 
    p_moves_orig, K, lkhd_z, lambda) 
.Primitive(".Call")(<pointer: 0x10b286b70>, n_sim, pattern, theta_init, 
    overlap, cluster_coords, p_moves_orig, K, lkhd_z, lambda)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
