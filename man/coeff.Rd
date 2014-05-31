\name{coeff}
\alias{coeff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coeff(y_vector, E_vector, a_values, b_values, cluster_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y_vector}{
%%     ~~Describe \code{y_vector} here~~
}
  \item{E_vector}{
%%     ~~Describe \code{E_vector} here~~
}
  \item{a_values}{
%%     ~~Describe \code{a_values} here~~
}
  \item{b_values}{
%%     ~~Describe \code{b_values} here~~
}
  \item{cluster_list}{
%%     ~~Describe \code{cluster_list} here~~
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
function (y_vector, E_vector, a_values, b_values, cluster_list) 
.Primitive(".Call")(<pointer: 0x10b2855f0>, y_vector, E_vector, 
    a_values, b_values, cluster_list)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
