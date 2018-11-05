#------------------------------------------------------------------------------
# Script to plot D-statistics
plot_dstat <- function(directory = "~/Projects/Botocudos/Files/Dstat/fromFlorian/",
                       fileName = "ABBAresult.Observed.txt",
                       h3 = "Log02WithError",
                       h2 =  "French", start = 0, end = 0.1, pval = 0.05){

setwd(directory)
abba <- read.table(fileName, header = T,
                   stringsAsFactors = F)
#head(abba)

# Dstat of the form (h1, h2, target, outgroup)
# usually h2 and target are fixed

#------------------------------------------------------------------------------
# Select trees where H3=target is tested
# Switch the H2 from H1 to H2 when necessary
abba <- abba[abba$H3 == h3 & (abba$H1 == h2 |abba$H2 == h2),]
index <- abba$H1 == h2
h1 <- abba[index, "H2"]
abba[index, "H1"] <- h1
abba[index, "H2"] <- h2

#------------------------------------------------------------------------------
# Order data frame by Dstat
abba <- abba[order(abba$JK.D),]
abba$H1 <- factor(abba$H1, levels = unique(abba$H1), ordered = T)

#------------------------------------------------------------------------------
# Plot
require(ggplot2)
low <- col_dstat(10)[1]
high <- col_dstat(10)[10]
limits <- c(start, end)
breaks <- c(start, pval, end)
if(!is.character(title)){
  title <- paste("D(H1, ", h2, "; ", h3, ", outgroup)", sep = "")
}
p <- ggplot(abba, aes(x = JK.D, y = H1, 
                 xmin = JK.D - V.JK.D.,
                 xmax = JK.D + V.JK.D., color = pvalue)) +
  geom_point(position = position_dodge(0.8)) +
  geom_errorbarh(position = position_dodge(0.8),
                height = 0.2) +
  theme_classic() +
  labs(x = "Dstat", y = "H1", title = title) +
  scale_color_gradient2(low = low, high = high,
                       limits = limits, breaks = breaks, mid = "white",
                       midpoint = pval)

return(p)

}

# Palette for Dstat

col_dstat <- function (n, h = -158, c. = c(100, 0), l = c(54, 0), power = 3, 
                       fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
  if (!is.null(gamma)) 
    warning("'gamma' is deprecated and has no effect")
  if (n < 1L) 
    return(character(0L))
  c <- rep(c., length.out = 2L)
  l <- rep(l, length.out = 2L)
  power <- rep(power, length.out = 2L)
  rval <- seq(1, 0, length = n)
  rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
                       C = c[2L] - diff(c) * rval^power[1L], H = h[1L]), fixup = fixup, 
              ...)
  if (!missing(alpha)) {
    alpha <- pmax(pmin(alpha, 1), 0)
    alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                    width = 2L, upper.case = TRUE)
    rval <- paste(rval, alpha, sep = "")
  }
  return(rval)
}

