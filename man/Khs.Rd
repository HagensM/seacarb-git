\encoding{latin1}
\name{Khs}
\alias{Khs}
\title{Dissociation constant of hydrogen sulfide (mol/kg)}
\description{Dissociation constant of hydrogen sulfide (mol/kg)}
\usage{
Khs(S=35, T=25, P=0, pHscale="T")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{pHscale}{choice of pH scale: "T" for using the total scale, "F" for using the free scale and "SWS" for using the seawater scale, default is total scale}
}

\details{
This formulation is only valid for specific ranges of temperature and salinity: 

\itemize{
\item S ranging between 0 and 45 and T ranging between 0 and 45oC.
}

The pressure correction is applied on the seawater scale. Hence, the values are first transformed from the total scale to the seawater scale, the pressure correction applied as described by Millero (1995), and the value is transformed back to the required scale (T, F or SWS).

Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.
} 

\value{
	\item{Khs}{Dissociation constant of hydrogen sulfide}
}

\references{
Millero F. J., 1995 Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} 59: 661-677.
}
             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl} and Heloise Lavigne
}

%\seealso{}

\examples{
  Khs(S=35,T=25,P=0, pHscale="T")
}

\keyword{utilities}
