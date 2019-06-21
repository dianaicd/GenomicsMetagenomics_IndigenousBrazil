ma2mn <- function(ma){
  db <- read.csv("~/Projects/Botocudos/Files/database_names.csv", header = T)
  mn <- db$MN[which(ma %in% db$MA)]
  return(mn)
}

mn2ma <- function(mn){
  db <- read.csv("~/Projects/Botocudos/Files/database_names.csv", header = T)
  ma <- db$MA[which(mn %in% db$MN)]
  return(ma)
}
