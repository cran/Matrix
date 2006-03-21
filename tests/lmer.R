library(Matrix)
library(lattice)
options(show.signif.stars = FALSE)

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm1a <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, method = "ML"))
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))

## should produce a warning but fit by REML
(fm1b <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, method = "AGQ"))

## transformed vars [failed in 0.995-1]
(fm2l <- lmer(log(Reaction) ~ log(Days+1) + (log(Days+1)|Subject),
              data = sleepstudy, method = "ML"))

## generalized linear mixed model
(fm3 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, family = poisson(), method = "PQL"))

## PQL is used per default:
fm3. <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, family = poisson())
fm3.@call <- fm3@call
stopifnot(all.equal(fm3, fm3., tol = 0))

## Laplace approximation {takes time}
(fm4 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
              OrchardSprays, poisson(), method = "Laplace"))

## Simple example by Andrew Gelman (2006-01-10) ----
n.groups <- 10 ; n.reps <- 2
n <- length(group.id <- gl(n.groups, n.reps))
## simulate the varying parameters and the data:
set.seed(0)
a.group <- rnorm(n.groups, 1, 2)
y <- rnorm (n, a.group[group.id], 1)
## fit and summarize the model
fit.1 <- lmer (y ~ 1 + (1 | group.id))
coef (fit.1)# failed in Matrix 0.99-6
(sf1 <- summary(fit.1)) # show() is as without summary()
## ranef and coef
rr <- ranef(fm1)
stopifnot(is.list(rr), length(rr) == 1, class(rr[[1]]) == "data.frame")
print(plot(rr))
cc <- coef(fm1)
stopifnot(is.list(cc), length(cc) == 1, class(cc[[1]]) == "data.frame")
print(plot(cc))
rr <- ranef(fm2)
stopifnot(is.list(rr), length(rr) == 2,
          all((sapply(rr, class) == "data.frame")))
print(plot(rr))
cc <- coef(fm2)
stopifnot(is.list(cc), length(cc) == 2,
          all((sapply(cc, class) == "data.frame")))
print(plot(cc))


## Many family = binomial cases
if (isTRUE(try(data(Contraception, package = 'mlmRev')) == 'Contraception')) {
    print(fm.1 <- lmer(use ~ urban + age + livch + (1 | district),
                       Contraception, binomial))
    print(system.time(fm1 <- lmer(use ~ urban + age + livch + (1 | district),
                                  Contraception, binomial), gc = TRUE))
    print(fm.2 <- lmer(use ~ urban + age + livch + (1 | district),
                       Contraception, binomial, method = 'Laplace'))
    print(system.time(lmer(use ~ urban + age + livch + (1 | district),
                           Contraception, binomial, method = 'Laplace'),
                      gc = TRUE))
##     print(fm.2a <- lmer(use ~ urban + age + livch + (1 | district),
##                         Contraception, binomial, method = 'AGQ'))
##     print(system.time(lmer(use ~ urban + age + livch + (1 | district),
##                            Contraception, binomial, method = 'AGQ'),
##                       gc = TRUE))
    print(fm.3 <- lmer(use ~ urban + age + livch + (urban | district),
                       Contraception, binomial))
    print(fm.4 <- lmer(use ~ urban + age + livch + (urban | district),
                       Contraception, binomial, method = 'Laplace'))
}

if (require('MASS', quietly = TRUE)) {
    bacteria$wk2 <- bacteria$week > 2
    contrasts(bacteria$trt) <-
        structure(contr.sdif(3),
                  dimnames = list(NULL, c("diag", "encourage")))
    print(fm5 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial))
    print(system.time(fm5 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial),
                      gc = TRUE))
    print(fm6 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial,
                      method = 'Laplace'))
    print(system.time(lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial,
                           method = 'Laplace'), gc = TRUE))
##     print(fm6a <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial,
##                        method = 'AGQ'))
##     print(system.time(lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial,
##                            method = 'AGQ'), gc = TRUE))
}


### mcmcsamp() :
## From: Andrew Gelman <gelman@stat.columbia.edu>
## Date: Wed, 18 Jan 2006 22:00:53 -0500

## Very simple example
y <- 1:10
group <- gl(2,5)
(M1 <- lmer (y ~ 1 + (1 | group))) # works fine
(r1 <- mcmcsamp (M1))              # dito
r2 <- mcmcsamp (M1, saveb = TRUE)  # gave error in 0.99-* and 0.995-[12]
(r10 <- mcmcsamp (M1, n = 10, saveb = TRUE))

## another one, still simple
y <- (1:20)*pi
x <- (1:20)^2
group <- gl(2,10)
M1 <- lmer (y ~ 1 + (1 | group)) # << MM: why is the "1 + " needed ?
mcmcsamp (M1, n = 2, saveb=TRUE) # fine

M2 <- lmer (y ~ 1 + x + (1 + x | group)) # false convergence
## should be identical (and is)
M2 <- lmer (y ~ x + ( x | group))#  false convergence -> simulation doesn't work:
if(FALSE) ## try(..) fails here (in R CMD check) [[why ??]]
    mcmcsamp (M2, saveb=TRUE)
## Error: inconsistent degrees of freedom and dimension ...

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

