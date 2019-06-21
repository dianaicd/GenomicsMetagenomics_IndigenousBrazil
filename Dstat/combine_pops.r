library(plyr)
panelName <- "~/archive/Panels/Magic.panel"
indName <- "h.B.P.ASM.VMAnc.modified.ind"
outgroup <- "Yoruba"

ind <- read.table(indName)
panel <- read.table(panelName, header = T)
colnames(ind) <- c("id", "U", "population")

miniPanel <- unique(panel[,c("population", "region")])
extendedInd <- join(ind, miniPanel, by = "population")

popInH1 <- 
popOfInterest <- unique(extendedInd$population[extendedInd$region %in% c("Botocudos","Americas")])
other <- unique(extendedInd$population[!(extendedInd$region %in% c("Botocudos","Americas"))])
other <- other[-which(other == outgroup)]

combPops <- t(combn(popOfInterest, 2))
remain <- length(other)

dstat <- data.frame(h1 = rep(combPops[,1], each = remain), 
                    h2 = rep(combPops[,2], each = remain), 
                    h3 = rep(other, nrow(combPops)), 
                    h4 = outgroup)


write.table(dstat, sub(".modified.ind", ".list_qpDstat", indName),
            col.names = F, row.names = F, sep = " ", quote = F)