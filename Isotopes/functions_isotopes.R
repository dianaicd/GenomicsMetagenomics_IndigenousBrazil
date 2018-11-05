
cross_isotope <- function(v, tef, sour, color){ # Name of the sources, data frames
  N <- tef[tef$Source == v, 2] + sour[sour$Source == v, 2]
  N_err <- sqrt(tef[tef$Source == v, 3]**2 + sour[sour$Source == v, 3]**2)
  C <- tef[tef$Sour ==v, 4] + sour[sour$Source == v, 4]
  C_err <- sqrt(tef[tef$Source == v, 5]**2 + sour[sour$Source == v, 5]**2)
  
  p <- c(N, N_err, C, C_err)
  names(p) <- c("N", "N_err", "C", "C_err")
  rect(xleft =  N - N_err*2,
      ybottom = C - C_err*2,
      xright =  N + N_err*2,
      ytop = C + C_err*2, 
      col= alpha(color, 0.2),
      border = NA)
  
  segments(x0 = N - N_err*2, x1 = N + N_err*2,
           y0 = C, col = color, lty = 5)
  segments(y0 = C - C_err*2, y1 = C + C_err*2,
           x0 = N, col = color, lty = 5)
  points(x = N, y = C, col = color, pch = 15) 
  #return(p)
}
