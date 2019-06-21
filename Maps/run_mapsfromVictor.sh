##Parameter order
#1. title
#2. color breaks
#3. infile. It only requires the first 3 columns. The first two, would be the coordinates, the way they are defined in google maps. and the 3rd would be the f3 value. 
#4. center (use 0 to get default parameters. otherwise, you invoke a routine for moving the map that is not efficient at all)
#5. decimal places to have in the color scale
#6. cex for the circles
#7. side for color legend (L or anything else) (L only works for zoom-ins apparently)
#outfile
#8. region specification. 2 options. see belo 
##ALL for the whole map. If you specify this, then you are done
##x1, x2, y1, y2, pdf_height, pdf_width for a zoomin. You have to specify the 6 of them
##there's some other options for just plotting a specific country/region by name, but has not been tested. For example: -177 -19 -60 90 corresponds to two points that encompass the americas. In google maps: -60,-177 is the bottom-left corner of the map and 90,-19 is the top-right. 

#map for whole world
Rscript WorldHeatmap_0_6.R sisi 20 BC25_f3_coords.txt 0 5 5 L out ALL
#map for the americas
Rscript WorldHeatmap_0_6.R sisi 20 BC25_f3_coords.txt 0 5 5 R sisi -177 -19 -60 90 6 5.5
#map for the americas with legend on the left
Rscript WorldHeatmap_0_6.R sisi 20 BC25_f3_coords.txt 0 5 5 L sisi -177 -19 -60 90 6 5.5

