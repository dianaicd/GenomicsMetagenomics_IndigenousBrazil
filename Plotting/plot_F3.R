# Functions to make plot related to F3-statistics


ciempies_plot <- function(f3, main = NULL, 
                          doLayout = T, doBar = T,
                          cex.axis = 0.5){
  f3 <- f3[order(f3$f_3, decreasing = F),]
  rownames(f3) <- 1:nrow(f3)
  
  color=sequential_hcl(length(f3$f_3), rev = T)[rank(f3$Z)]
  if(doLayout){
    layout(matrix(seq(1, 2), byrow = F, nrow = 1, ncol = 2),
           widths = c(3,1), heights = c(6))
  }
  # Plot with values
  plot(y = rownames(f3), x = f3$f_3, bty = "n", axes = F,
       ylab = NA, xlab = expression(F[3]),
       col = color, pch = 16, main = main)
  segments(y0 = as.integer(rownames(f3)), 
           x0 = f3$f_3 - qnorm(0.995)*f3$std.err,
           x1 = f3$f_3 + qnorm(0.995)*f3$std.err,
           col = color)
  axis(side = 2, 
       at =  1:nrow(f3) , labels = f3$Source2, las = 2, cex.axis = cex.axis)
  axis(side=1)
  grid()
  
  if(doBar){
    # Color bar for Z-score
    par(font.main = 1)
    barplot(cbind(seq(min(f3$Z), max(f3$Z), 
                      length.out = 100)),
            border = NA,
            col = sequential_hcl(100, rev = T), 
            yaxt="n", horiz = F, main = "Z-score")
    axis(side=4, at = seq(min(f3$Z), 
                          max(cumsum(seq(min(f3$Z), 
                                         max(f3$Z)))),
                          length.out = 6),
         labels = round(seq(min(f3$Z), max(f3$Z), 
                            length.out = 6), 1),
         las = 2)
  }
}
