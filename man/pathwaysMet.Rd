\name{pathwaysMet}
\alias{pathwaysMet}
\title{
Metabolic Pathways
}

\description{This function finds the metabolite related pathways.}

\usage{
pathwaysMet(doc)
}

\arguments{
\item{doc}{A list of metabolites information produced by \code{\link{readMet}}function.}
}

\value{
A data frame contains the pathways associated with each metabolite.
}


\seealso{
  \code{\link{readMet}},  \code{\link{taxonomyMet}}, \code{\link{enzymesMet}}, \code{\link{diseasesMet}}
   }
   
   
\examples{
\donttest{

ID=c("HMDB0000122","HMDB0000124","HMDB0000243","HMDB0000263")
doc=readMet(ID) 

pat=pathwaysMet(doc)
}
}

