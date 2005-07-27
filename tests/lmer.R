library(Matrix)
options(show.signif.stars = FALSE)

(fm1 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays))
(fm2 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, method = "ML"))

## should produce a warning but fit by REML
(fm1 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, method = "AGQ"))

## generalized linear mixed model
(fm3 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, poisson(), method = "PQL"))

## should use PQL
(fm3 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, poisson()))

## Laplace approximation
#(fm4 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
#             OrchardSprays, poisson(), method = "Laplace"))

if (isTRUE(try(data(Contraception, package = 'mlmRev')) ==
           'Contraception') && exists("nlminb", mode = "function")) {
    print(fm1 <- lmer(use ~ urban + age + livch + (1 | district),
                      Contraception, binomial))
    print(system.time(fm1 <- lmer(use ~ urban + age + livch + (1 | district),
                                  Contraception, binomial), gc = TRUE))
    print(fm2 <- lmer(use ~ urban + age + livch + (1 | district),
                      Contraception, binomial, method = 'Laplace'))
    print(system.time(fm2 <- lmer(use ~ urban + age + livch + (1 | district),
                                  Contraception, binomial, method = 'Laplace'),
                      gc = TRUE))
    print(fm2a <- lmer(use ~ urban + age + livch + (1 | district),
                       Contraception, binomial, method = 'AGQ'))
    print(system.time(fm2a <- lmer(use ~ urban + age + livch + (1 | district),
                                   Contraception, binomial, method = 'AGQ'),
                      gc = TRUE))
    print(fm3 <- lmer(use ~ urban + age + livch + (urban | district),
                 Contraception, binomial))
    print(fm4 <- lmer(use ~ urban + age + livch + (urban | district),
                 Contraception, binomial, method = 'Laplace'))
}
if (require('MASS', quietly = TRUE) && exists("nlminb", mode = "function")) {
    bacteria$wk2 <- bacteria$week > 2
    contrasts(bacteria$trt) <-
        structure(contr.sdif(3),
                  dimnames = list(NULL, c("diag", "encourage")))
    print(fm5 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial))
    print(system.time(fm5 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial),
                      gc = TRUE))
    print(fm6 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial,
                 method = 'Laplace'))
    print(system.time(fm6 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial,
                                  method = 'Laplace'), gc = TRUE))
    print(fm6a <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial,
                  method = 'AGQ'))
    print(system.time(fm6a <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial,
                                   method = 'AGQ'), gc = TRUE))
}
q('no')

