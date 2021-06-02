\description{
\itemize{
\item we base our analysis on \eqn{k} cases, rather than \eqn{k} \emph{other} cases as prescribed in the paper.
\item we do not subtract 1 from the \emph{accumulated numbers of other cases} and \emph{accumulated numbers of others at risk}, as was prescribed in the paper to discount selection bias
\item M is the total number of areas included, not the number of additional areas included.  i.e. \eqn{M} starts at 1, not 0.
\item p-values are not based on the original value of \eqn{k}, rather the actual number of cases observed until we view \eqn{k} or more cases.  Ex:  if \eqn{k = 10}, but as we consider neighbors we encounter 1, 2, 9 then 12 cases, we base our \eqn{p}-values on \eqn{k=12}
\item we do not provide a Monte-Carlo simulated \eqn{R}:  the number of tests that attain significance at a fixed level \eqn{\alpha}
}
The first two and last differences are because we view the testing on an area-by-area level, rather than a case-by-case level.
}
\note{
The \code{clusters} list elements are themselves lists reporting:\cr\cr
\tabular{ll}{
  \code{location.IDs.included} \tab ID's of areas in cluster, in order of distance\cr
  \code{population} \tab population of cluster\cr
  \code{number.of.cases} \tab number of cases in cluster\cr
  \code{expected.cases} \tab expected number of cases in cluster\cr
  \code{SMR} \tab estimated SMR of cluster\cr
  \code{p.value} \tab \eqn{p}-value\cr
}

