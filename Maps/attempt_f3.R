library(ggplot2)
library(maptools )
library(stringi)
library(descr)
library(rgeos)
library(mapdata)

m_coord <- read.table("~/Projects/Botocudos/Files/Maps/Maanasa/MaanasaCoordsForDiana.txt")
colnames(m_coord) <- c("Pop", "y", "x")
m <- read.table("~/Projects/Botocudos/Files/Panels/Maanasa_0.01_minind.panel", header = T)
m <- unique(m[,c("Pop", "Region")])
m_coord <- join(m_coord, m, by = "Pop")
m_coord <- m_coord[complete.cases(m_coord),]
m_coord <- m_coord[m_coord$Region == "Americas",]
xl <- range(m_coord$x)
yl <- range(m_coord$y)
# png("~/Projects/Botocudos/Plots/Map/F3_map",
#     width = 7, height = 7, res = 350, units = "in")


map("worldHires","Mexico", col="cornsilk2", fill=T,
    bg = "azure3",
    xlim = c(-120,-10), ylim = c(-170, 150), mar = c(0,0,0,0), border = NA
    #xlim = c(-52, -38), ylim = c(-25, -13)
)

points(x = m_coord$x,
         y = m_coord$y,
         col = alpha("white",0.7),pch = 16)
  