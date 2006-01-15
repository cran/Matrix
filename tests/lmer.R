library(Matrix)
options(show.signif.stars = FALSE)

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, method = "ML"))

## should produce a warning but fit by REML
(fm1b <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, method = "AGQ"))

## generalized linear mixed model
(fm3 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, family = poisson(), method = "PQL"))

## should use PQL
(fm3 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, family = poisson()))

## Laplace approximation {takes time}
## (fm4 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
##              OrchardSprays, poisson(), method = "Laplace"))

## Simple example by Andrew Gelman (2006-01-10) ----
n.groups <- 10 ; n.reps <- 2
n <- length(group.id <- gl(n.groups, n.reps))
## simulate the varying parameters and the data:
set.seed(0)
a.group <- rnorm(n.groups, 1, 2)
y <- rnorm (n, a.group[group.id], 1)
## fit and summarize the model
fit.1 <- lmer (y ~ 1 + (1 | group.id))
coef (fit.1)# failed in Matrix 0.99-6 -- FIXME: should get a show() method
try(summary(fit.1))

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

proc.time() # for ``statistical reasons''
