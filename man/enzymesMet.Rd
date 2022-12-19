\name{enzymesMet}

\alias{enzymesMet}

\title{
Metabolite-associated Enzymes
}

\description{This function finds the metabolite related enzymes.}

\usage{
enzymesMet(doc)
}


\arguments{
  \item{doc}{A list of metabolites information produced by \code{\link{readMet}} function.}
}

\value{
A data frame contains the enzymes associated with each metabolite.
}

\seealso{
  \code{\link{pathwaysMet}} ,  \code{\link{taxonomyMet}},  \code{\link{diseasesMet}}
   }


\examples{
\donttest{

ID=c("HMDB0000122","HMDB0000124","HMDB0000243","HMDB0000263")
doc=readMet(ID) 

enz=enzymesMet(doc)
}
}

