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
(fm4 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays, poisson(), method = "Laplace"))

