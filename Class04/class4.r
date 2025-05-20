#' ---
#' title: "Intro to R"
#' author: "Grace Wang"
#' ---

#Class 04 R script

x <- 1:50

plot(x, sin(x))
plot(x, sin(x), typ = "l", col = "turquoise", lwd = 3, 
     xlab = "Silly x axis", ylab = "Sensible y axis")
