rgcm.example <- data.frame(colour=c(rnorm(14, 0.2, 0.09),rnorm(12, 0.8, 0.09), rnorm(14, 0.5, 0.15)),
                           size=c(rnorm(40, 10, 1)), 
                           shape=factor(c(sample(c("a", "b", "c"), 26, replace=T, prob=c(0.8,0.1,0.1)), sample(c("a", "b", "c"), 14, replace=T, prob=c(0.2,0.3,0.5)))),
                           symmetry=factor(c(sample(c("yes", "no"), 40, replace=T, prob=c(0.6,0.4)))), 
                           category=factor(c(rep("X", 26), rep("Y", 14))))


