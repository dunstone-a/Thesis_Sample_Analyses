# Linear model section
# Amelia Dunstone
# 2024-10-11

ps <- c(0.02, 0.48)
vars <- 2*ps*(1-ps)
vars

p_j <- seq(0, 1, 0.01)
var_p_j <- 2*p_j*(1-p_j)
df <- data.frame(p_j, var_p_j)

plot(p_j, var_p_j, type = "l", xlab = expression(p[j]), ylab = expression(Var(p[j])))



library(ggplot2)

# Create figure 2.1 ------------------------------------------------------------


x1 <- 0.05
y1 <- 2*x1*(1-x1)

colours <- c("#315990", "#5AE3EA", "#EA9B60")

p1 <- ggplot(data.frame(x = c(0, 1)), aes(x = x)) +
    stat_function(fun = function(x) 2 * x * (1 - x), colour = colours[1], linewidth = 1) +  # Main function
    annotate("segment", x = x1, xend = x1, y = 0, yend = y1, linetype = 2, colour = colours[3]) +  # Vertical dashed line
    annotate("segment", x = 0, xend = x1, y = y1, yend = y1, linetype = 2, colour = colours[3])+  # Horizontal dashed line
    geom_area(stat = "function", fun = function(x) 2 * x * (1 - x), xlim = c(0, x1), fill = colours[3], alpha = 0.2) +  # Shaded area
    labs(x = expression(p[j]), y = expression(Var(g[ij]))) +  # Axis labels with math notation
    theme_minimal(base_family = "serif")  # Use serif font
p1

png("figures/genotype_variance_plot.png", width = 800, height = 600, res = 300)
p1
dev.off()



