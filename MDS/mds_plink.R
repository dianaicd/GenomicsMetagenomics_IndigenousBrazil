d_tab <- read.table("~/Projects/Botocudos/Files/MDS/2018_12_21/dits_maanasa_boto_victor.dist")
d_mat <- as.matrix(d_tab)
mean.d_mat = mean(d_mat,na.rm = T)
d_mat[which(is.na(d_mat))] <- mean.d_mat
d <- as.dist(d_tab)

x <- cmdscale(d=d_mat, eig = T, k = 6)

nb.individual <- nrow(d_mat)
see.where.na <- matrix(rep(0,nb.individual^2),nrow = nb.individual)
see.where.na[which(is.na(d_mat))] <- mean.d_mat
image(see.where.na)

id <- read.table(file = "~/Projects/Botocudos/Files/MDS/2018_12_21/dits_maanasa_boto_victor.dist.id")
id$col <- match(id$V1, "Botocudo")
id$col[is.na(id$col)] <- rep(2, nb.individual)[is.na(id$col)]


plot(x$points[,1], x$points[,2], col = id$col)
