  
## See also file plot_sim_example.R

library(ggplot2)
## plot of results

# Load file with results from "simulation_results.R"
# load("..")

ggdata <- data.frame(threshold = c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10),
                     value = as.numeric(rbind(eksp1, eksp2, eksp3, eksp4, eksp5)),
                     experiment = rep(as.factor(1:5), each = 7),
                     id = rep(factor(c("sensitivity", "false positive rate", "fdr", "fdr, unadjusted")), each = 35))

p <- ggplot(aes(x = threshold, y = value, group = id, col = experiment), data = ggdata) + 
  geom_point(shape = 4) + facet_wrap(~ id, scales = "free") +  xlab("threshold") + ylab("")



#pdf("...", width = 8, height = 6)
p + geom_abline(aes(slope = 1, intercept = 0), data = data.frame(id = factor(c("fdr","fdr, unadjusted"))), linetype="dotted") +
  geom_abline(aes(slope = 0.717, intercept = 0), data = data.frame(id = factor(c("fdr","fdr, unadjusted"))), linetype="dotted")
#dev.off()

