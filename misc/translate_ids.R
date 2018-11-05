ma2mn <- function(ma){
  db <- read.csv("~/Projects/Botocudos/Files/database_names.csv", header = T)
  mn <- db$MN[db$MA==ma]
  return(mn)
}

mn2ma <- function(mn){
  db <- read.csv("~/Projects/Botocudos/Files/database_names.csv", header = T)
  ma <- db$MA[db$MN==mn]
  return(ma)
}
