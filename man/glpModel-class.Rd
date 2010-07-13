\name{glpModel-class}
\Rdversion{1.1}
\docType{class}
\alias{glpModel-class}
\title{Class "glpModel" of General Linear Prediction Models}
\description{
  The class \code{"glpModel"} conceptionally contains a very large class
  of \emph{\dQuote{General Linear Prediction Models}}.

  Its \code{resp} slot (of class \code{"\linkS4class{respModule}"}) may
  model linear, non-linear, generalized linear and non-linear
  generalized response models.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("glpModel", ...)},
  but typically rather are returned by our modeling functions, e.g., the
  (experimental, hence currently hidden) \code{glm4()}.
}
\section{Slots}{
  \describe{
    \item{\code{resp}:}{a \code{"\linkS4class{respModule}"} object.}
    \item{\code{pred}:}{a \code{"\linkS4class{predModule}"} object.}
  }
}
\section{Methods}{
 There are no methods, currently, but that will change.
}
%\author{Douglas Bates}
\seealso{
  The constituents of this class are
  \code{\linkS4class{respModule}} and
  \code{\linkS4class{predModule}}, both of which have several sub
  classes.
}
\examples{
showClass("glpModel")

### All the following is very experimental -- and probably will change: -------

data(CO2, package="datasets")
## dense linear model
str(Matrix:::glm4(uptake ~ 0 + Type*Treatment, data=CO2, doFit = FALSE), 4)
## sparse linear model
str(Matrix:::glm4(uptake ~ 0 + Type*Treatment, data=CO2, doFit = FALSE,
                  sparse = TRUE), 4)

## From example(glm): -----------------

## Dobson (1990) Page 93: Randomized Controlled Trial :
str(trial <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                        outcome=gl(3,1,9,labels=LETTERS[1:3]),
                        treatment=gl(3,3,labels=letters[1:3])))
glm.D93 <- glm(counts ~ outcome + treatment, family=poisson, data=trial)
summary(glm.D93)
glmM <- Matrix:::glm4(counts ~ outcome + treatment, family = poisson, data=trial)
stopifnot( all.equal(glmM@pred@coef, unname(coef(glm.D93))) )

## Watch the iterations --- and use no intercept --> more sparse X
## 1) dense generalized linear model
glmM <- Matrix:::glm4(counts ~ 0+outcome + treatment, poisson, trial,
                      verbose = TRUE)
## 2) sparse generalized linear model
glmS <- Matrix:::glm4(counts ~ 0+outcome + treatment, poisson, trial,
                      verbose = TRUE, sparse = TRUE)
str(glmS, max.lev = 4)
stopifnot( all.equal(glmM@pred@coef, glmS@pred@coef),
           all.equal(glmM@pred@Vtr,  glmS@pred@Vtr) )


## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(u = c(5,10,15,20,30,40,60,80,100),
                       lot1 = c(118,58,42,35,27,25,21,19,18),
                       lot2 = c(69,35,26,21,18,16,13,12,12))
str(gMN <- Matrix:::glm4(lot1 ~ log(u), data=clotting, family=Gamma, verbose=TRUE))
glm. <- glm(lot1 ~ log(u), data=clotting, family=Gamma)
stopifnot( all.equal(gMN@pred@coef, unname(coef(glm.))) )
}
\keyword{classes}
