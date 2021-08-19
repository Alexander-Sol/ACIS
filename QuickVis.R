# Quick Vis
library(ggplot2)

MenonPR.results <- readRDS("Results/MenonPR/MenonPR.rds")
mpr.ari <- lapply(MenonPR.results, pluck, "ARI") %>% lapply('[', 1:5) %>% unlist()
mrp.algo <- rep(names(MenonPR.results), each = 5)
mrp.run <- rep(1:5, 7)
mrp.logRT <- lapply(MenonPR.results, pluck, "runtime") %>% lapply('[', 1:5) %>% lapply('*', 60) %>% unlist() %>% log()
mrp.df <- data.frame(Algorithm = mrp.algo,
                     Run = mrp.run,
                     ARI = mpr.ari,
                     `Log Runtime` = mrp.logRT)
plot.df <- mrp.df[-1*(1:5), ] 
# plot.df <- plot.df[-1*((0:5)*5+1), ]

plot <- ggplot(data = plot.df, aes(x = Log.Runtime, y = ARI, color = Algorithm)) + geom_point()
