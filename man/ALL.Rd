\name{ALL}
\alias{ALL}
\docType{data}
\title{Acute Lymphoblastic Leukemia Data from the Ritz Laboratory}
\description{
 The data consist of microarrays from 128 different individuals with 
 acute lymphoblastic leukemia (ALL). A number of additional covariates
 are available. The data have been normalized (using gcrma) and it is
 the jointly normalized data that are available here. The data are
 presented in the form of an \code{exprSet} object.
}
\usage{data(ALL)}
\format{
 The different covariates are:
\itemize{
 \item \code{cod}: The patient IDs.
  \item \code{diagnosis} The date of diagnosis.
  \item \code{sex}   The sex of the patient, coded as \code{M} and \code{F}.
  \item \code{age}  The age of the patient in years.
  \item \code{BT}  The type and stage of the disease; \code{B} indicates 
                  B-cell ALL while a \code{T} indicates T-cell ALL. 
  \item \code{remission} A factor with two levels, either \code{CR} indicate 
   that remission was achieved or \code{REF} indicating that the patient was
   refractory and remission was not achieved. 
  \item \code{CR}  As above, but with extra information on whether the patient
    died while remission was being attempted.
  \item \code{date.cr} The date on which remission was achieved.
  \item \code{t(4;11)}  A logical vector indicating whether a t(4;11) 
      translocation was detected.
  \item \code{t(9;22)}  A logical vector indicating whether a t(9;22) 
    translocation was detected.
  \item \code{cyto.normal} A logical vector indicating whether the 
    cytogenetics were normal.
  \item \code{citog}  A vector indicating the various cytogenetic abnormalities
    that were detected.
  \item \code{mol.biol} The assigned molecular biology of the cancer (mainly 
    for those with B-cell ALL), BCR/ABL, ALL/AF4, E2APBX etc. 
  \item \code{fusion protein} For those with BCR/ABL which of the fusion 
     proteins was detected,  \code{p190}, \code{p190/p210}, \code{p210}.
  \ite \code{mdr}   The patients response to multidrug resistance, either
   \code{NEG}, or \code{POS}.
  \item \code{kinet}  Not sure what this is.
  \item \code{ccr}   What does it really mean
 \item \code{relapse}  Definition
 \item \code{transplant}  Did the patient receive a bone marrow transplant or
    not.
 \item \code{f.u} List the possible values and what they mean "BMT / DEATH IN CR" "rel" "rel" "rel " ...
 \item { date last seen} Date the patient was last seen.
}
}
\details{
 Put some detail here
}
\source{
 Put a reference here
}
\examples{
data(ALL)
}
\keyword{datasets}
