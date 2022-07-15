# --------------------------------------------------------------------------------------------------
args <- commandArgs(TRUE)

# ## Help section

if("--help" %in% args) {
    cat("

Arguments:
    --str_basename          
    --prefix_out                              
    --models
    --poplabels
    --n_bad_entries
    --do_everything
    --numRuns                                   
    --get_best_run
    --print_sfs
    --print_worst_fit
    --make_pdf
    --make_png

    --help                                - print this text

Example:

Rscript make_different_plots.r                                                          \\
    --str_basename='sim_{model}/{model}_rep{i}/{model}_rep{i}'                          \\
    --prefix_out=summary                                                                \\
    --models='c('Botocudos_Adm')'                                                       \\
    --numRuns=100                                                                       \\                      
    --poplabels='list(Botocudos_Adm = c(\"Botocudos\", \"Karitiana\", \"Surui\"))'      \\
    --n_bad_entries=30                                                                  \\
    --do_everything=T                                                                   \\
    --get_best_run=T                                                                    \\
    --print_sfs=T                                                                       \\
    --print_worst_fit=T                                                                 \\
    --make_pdf=T                                                                        \\
    --make_png=T                                                                
\n\n"
         
         )

    q(save="no")
}

source("parse_arguments.R")
# --------------------------------------------------------------------------------------------------
#Specify how many runs were done per model (here same number is assumed for all models)
# Code modified by Diana I. Cruz Davalos 
# get arguments from the command line
prefix_out = get_args(name  = "prefix_out", default = "summary")

#numRuns=50 ##TO BE MODIFIED
numRuns = get_args(name = "numRuns", mandatory = TRUE)

#--- Define which models to analyse here ---------------------
curModels = get_args("models", "c('Botocudos_Adm')", mandatory = TRUE, eval_string = TRUE )

#--- Specify which pop labels are valid for the different models 
pop.Labels = get_args(
    "poplabels", 
    'list(Botocudos_Adm = c("Botocudos", "Karitiana", "Surui"))', 
    mandatory = TRUE, eval_string = TRUE
    ) 

numBadEntriesToPLot = get_args("n_bad_entries", 30) #For worse fitted SFS entries in each model
doEverything = get_args(name = "do_everything", default = FALSE)

if (doEverything) {
   extractBestRun        = TRUE
   #drawModel             = TRUE
   drawMarginalSFS       = TRUE
   printPDF              = TRUE
   printPNG              = TRUE
   #summarizeRresults     = TRUE
   plotWorstFitSFS       = TRUE
} else {  #Select what to do
    extractBestRun        = get_args(name = "get_best_run", default = TRUE)
   #drawModel             = get_args(name = "draw_model", default = FALSE)
   drawMarginalSFS       = get_args(name = "print_sfs", default = TRUE)
   printPDF              = get_args(name = "make_pdf", default = TRUE)
   printPNG              = get_args(name = "make_png", default = TRUE)
   #summarizeRresults     = get_args(name = "make_summary", default = TRUE)
   plotWorstFitSFS       = get_args(name = "print_worst_fit", default = TRUE)
}


MatLhoods=c()


str_basename                = get_args(
    name = "str_basename",
    default = "sim_{model}/{model}_rep{i}/{model}_rep{i}", 
    mandatory = TRUE
    )

#maxLpar_path                = str_glue("{str_basename}_maxL.par")



# --------------------------------------------------------------------------------------------------
# Functions to extract and visualize the maximum likelihoods of the runs
plot_maxEstlhood <- function(MatBrent, model, numRuns, bestRun, BEST){
  x_value = seq(1,nrow(MatBrent))
  plot(
    x_value, MatBrent[,1],
    ylim = c(1.002*max(MatBrent), max(MatBrent)), 
    type = 'l', lwd = 1, xlab = "Update", ylab = "", cex.main = 0.7,
    main = str_glue("MaxEstLhood for {model}"), las = 2
  )
  
  for (i in c(1:numRuns)){
    lines(x_value, MatBrent[,i], type='l', lwd=1)
  }
  
  lines(x_value, BEST, type='l', lwd=1, col="red")
  legend("bottomright",lty=1, col="red",
         legend=str_glue("Best run: {bestRun}"), 
         inset=0.01, box.col = NA
  )
}

extract_best_lhoods <- function(numRuns, model,
                                str_basename, 
                                bestlhoods_path, 
                                all_lhoods_path, 
                                brentlhoods_path,
                                plot_maxEstlhood_path,
                                printPDF,
                                printPNG){
  
  cat("Extracting best lhood files\n")
  all_runs <- do.call(
    rbind, 
    lapply(1:numRuns, function(i) read.table( str_glue(bestlhoods_path), header = T  ))
    )
  all_runs$run <- 1:numRuns
  # Save summary of all runs
  write.table(all_runs, all_lhoods_path, col.names = T, row.names = F, sep = "\t", quote = F)
  # Get best run
  bestRun <- which(all_runs$MaxEstLhood == max(all_runs$MaxEstLhood))[1]   
  str_bestmodel = gsub("{i}", "{bestRun}", str_basename, fixed = T)
  # Copy files to output directory
  system(paste0(" cp -r ",  str_glue(str_bestmodel), "* ", out_dir))
  cat("done!\n\n")
  cat("Extracting and plotting lhood\n")
  
  MatBrent <- do.call(cbind, 
                      lapply(1:numRuns, 
                             function(i){
                               BRENT=read.delim(str_glue(brentlhoods_path))
                               BRENT=BRENT[-nrow(BRENT),"MaxEstLhood"] # remove line with ---
                             }
                      )
  )
  
  BEST=MatBrent[,bestRun]
  
  if (printPDF) {
    pdfFileName=paste0(plot_maxEstlhood_path, ".pdf")
    cat(paste("   Plotting likelihoods to file", pdfFileName), sep="\n")
    pdf(pdfFileName, height=10, width=8)
    plot_maxEstlhood(MatBrent, model, numRuns, bestRun, BEST)
    dev.off()
  }
  
  if (printPNG) {
    pngFileName=paste0(plot_maxEstlhood_path, ".png")
    cat(paste("   Plotting likelihoods to file", pngFileName, "\n"))
    png(pngFileName, height=3000, width=2000, res=300)
    par(mar=c(5,5,2,1))
    plot_maxEstlhood(MatBrent, model, numRuns, bestRun, BEST)
    dev.off()
  }
  cat("done!\n\n")
  return(bestRun)
}

# --------------------------------------------------------------------------------------------------
# Functions to get and plot the marginal SFS

prepare_sfs <- function(sfsFile, transpose = T){
  dim.sfs=read.table(file=sfsFile, skip=1, header=F,nrows=1)
  numDims=as.numeric(dim.sfs[1])
  dim.sfs=dim.sfs[-1]
  #Inverting number of dimension to reflect the difference in SFS ordering and normal R matrix ordering
  dim.sfs=as.numeric(dim.sfs+1, nrow=1)[numDims:1]
  
  #---- Reading obs or exp multidimensional sfs ---------------
  sfs=as.matrix(read.table(file=sfsFile, skip=2, header=F, nrows=1))
  if(transpose){
    #Extract number of dimensions
    sfs=array(data=sfs, dim=dim.sfs)
    sfs[1]=0
    # Transposition of the matrix entries
    sfs=aperm(sfs, numDims:1)
  }
  return(sfs)
}

plot.marginal1D.sfs=function(obs.sfs, exp.sfs, mar, title, log) { 
   #cols=c("orange", "cornflowerblue")
   cols = c("firebrick2", "black")
   #Remove observed monomorphic sites
   obs.sfs[1] = 0
   #Get marginal SFS
   obs=apply(obs.sfs, mar, sum)
   exp=apply(exp.sfs, mar, sum)
   tot.num.sites=sum(obs)
   exp.abs=exp*tot.num.sites
   
   miny=min(c(obs, exp.abs))
   maxy=max(c(obs, exp.abs))
   if(log) logy="y" else logy=""
   plot(
    x=0:(length(obs)-1), y = obs, 
    type="b", col=cols[1], lwd=2, xlab="i", ylab="sfs[i]", pch=16, cex=1.2, 
        main = title, log=logy, ylim=c(miny, maxy))
   lines(
    x=0:(length(obs)-1),y = exp.abs, 
    col=cols[2], lwd=2, type="b" , pch=16, cex=1.2
    )
   legend("topright", legend=c("observed", "expected"), col=cols, lwd=2)
}


get_marginal_SFS <- function(bestrun_bestlhoods_path, model,
                             sfsFile.obs, sfsFile.exp,
                             pop.Labels, curcase, printPDF, printPNG,
                            plot_basename
                             ){
  
  cat("Drawing marginal SFS ...\n")
  
  curBestParams = read.table(bestrun_bestlhoods_path, header=T)
  maxEstLhood   = curBestParams[["MaxEstLhood"]]
  maxObsLhood   = curBestParams[["MaxObsLhood"]]
  
  cat(paste("Printing sketch of demographic scenario.....\n")) 
  
  
  #---- Reading number of dimensions and sample sizes ---------------
  dim.sfs=read.table(file=sfsFile.obs, skip=1, header=F,nrows=1)
  numDims=as.numeric(dim.sfs[1])
  sfs.obs <- prepare_sfs(sfsFile = sfsFile.obs)
  sfs.exp <- prepare_sfs(sfsFile = sfsFile.exp)
  
  plot_sfs <- function(do_log){
    plot.new()
    mtext(paste(model,"_sfs_fit\n","EstLhood=", maxEstLhood, "   ObsLhood=", maxObsLhood, sep=""), side=1, cex=0.8)
    par(mar=c(4,4,2,1))
    for (i in 1:numDims) {
      plot.marginal1D.sfs(obs.sfs=sfs.obs, exp.sfs=sfs.exp, mar=i, pop.Labels[[curcase]][i], log=do_log)
    }
  }
  
  if (printPDF) {
    cat(paste("   Printing obs and exp sfs to file", plot_basename), sep="\n")
    pdf(paste0(plot_basename, ".pdf"), height=10, width=8)
    mat.layout=rbind(c(1,1), c(2,3), c(4,5), c(6,7))
    layout(mat.layout, height=c(0.15,1,1,1))
    par(mar=c(1,1,1,1))
    plot_sfs(do_log = F)
    dev.off()
  }
  
  if (printPNG) {

    cat(paste("   Printing obs and exp sfs to file", plot_basename, "\n"))
    png(paste0(plot_basename, ".png"), height=3000, width=2000, res=300)
    mat.layout=rbind(c(1,1), c(2,3), c(4,5), c(6,7))
    layout(mat.layout, height=c(0.15,1,1,1))
    par(mar=c(1,0,0,0))
    plot_sfs(do_log=F)
    dev.off()
    

    png(paste0(plot_basename, "_log_scale.png"), height=3000, width=2000, res=300)
    layout(mat.layout, height=c(0.15,1,1,1))
    par(mar=c(1,0,0,0))
    plot_sfs(do_log=T)
    dev.off()
  }
  cat("done!\n", sep="\n")
}

# --------------------------------------------------------------------------------------------------
# worst fitted entries


plot_worst_entries <- function(sfsFile.obs, sfsFile.exp, printPDF, printPNG,
                               plot_basename){
    do_barplot <- function(){
        par(mfrow=c(2,1))
        par(mar=c(7,4,2,1))
        barplot(dif.lhoods[ordered.indices][1:numBadEntriesToPLotBIS], 
                names.arg = worse.sfs.labels,las=2,
                ylab="logLhood difference (max-est logLhoods)", main=paste(numBadEntriesToPLotBIS, "worse fitted SFS entries"), 
                cex.names=0.8)
        barplot(rbind(sfs.obs.pol,sfs.exp.pol)[,ordered.indices][,1:numBadEntriesToPLotBIS], col=c("black","grey"), 
                names.arg = worse.sfs.labels, beside = TRUE, las=2, ylab= "number of sites", cex.lab=1, 
                main=paste(numBadEntriesToPLotBIS, "worse fitted SFS entries"), 
                cex.names=0.8)
        legend("top", legend =c("Observed", "Expected"), pch=15,  col=c("black","grey"), cex=0.8)
        par(mfrow=c(1,1))
    }
  
    #---- Reading obs and exp multidimensional sfs ---------------
    sfs.obs = prepare_sfs(sfsFile.obs, transpose = F)

    cat(paste("Plotting worse fitted SFS entries", curcase), sep="\n")
    cat("   Reading obs and exp SFS ...\n")
    #---  Extract number of polymorphic sites
    mono.pos=c(1,length(sfs.obs))
    sfs.obs.pol=sfs.obs[-mono.pos]
    num.polym.sites=sum(sfs.obs.pol)
    #---- Reading exp multidimensional sfs ---------------
    sfs.exp = prepare_sfs(sfsFile.exp, transpose = F)
    penalty=1e-10
    #replace all zero entries in exp sfs by a small arbitrary value
    sfs.exp[sfs.exp==0]=penalty
    #--- Getting expected SFS in absolute numbers
    sfs.exp.pol=sfs.exp[-mono.pos]*num.polym.sites
    cat("   done!\n")

    #--- Computing logLhood individual entries
    lhoods.obs=sfs.obs.pol*log10(sfs.exp[-mono.pos])
    lhood.obs.total=sum(lhoods.obs)
    #--- computing max possible lhoods
    sfs.exp.max=sfs.obs.pol/num.polym.sites
    #Set zero obs entries to 1
    sfs.exp.max[sfs.obs.pol==0]=1
    lhoods.exp.max=sfs.obs.pol*log10(sfs.exp.max)
    lhoods.exp.max.total=sum(lhoods.exp.max)
    #--- Computing difference in lhoods and order by decreasing difference
    dif.lhoods=lhoods.exp.max-lhoods.obs
    ordered.indices=order(abs(dif.lhoods),decreasing =T)
    #--- extract SFS entries labels
    popSizes=read.table(file=sfsFile.obs, skip=1, header=F, nrows=1)[-1]
    sfs.labels=getEntriesLabels(popSizes)[-mono.pos]

    #--- Plot worse SFS entries
    numBadEntriesToPLotBIS=min(length(sfs.labels[ordered.indices]),numBadEntriesToPLot)
    worse.sfs.labels=sfs.labels[ordered.indices][1:numBadEntriesToPLotBIS]
    if (printPDF) {
    cat(paste("   Plotting worse fitted sfs entries to file", plot_basename, "\n"))
    pdf(paste0(plot_basename, ".pdf"), height=10, width=8)
    do_barplot()
    dev.off()
    }
    if (printPNG) {
    cat(paste("   Plotting worse fitted sfs entries to file", plot_basename, "\n"))
    png(paste0(plot_basename, ".png"), height=3000, width=2000, res=300)
    do_barplot()
    dev.off()
    }
    cat("done!\n\n")
}

getEntriesLabels <- function(pop.sizes) {
   if(prod(pop.sizes) >= 1e6) {
      stop("Too many entries (> 1e6).")
   }
   inv.pop.sizes=pop.sizes[length(pop.sizes):1]
   li <- lapply(inv.pop.sizes, seq, from = 0, by = 1)
   entr <- expand.grid(li)
   entr <- rev(entr)
   entries.lab <- apply(entr, 1, function(x) {
      paste0("(", paste(x, collapse = ","), ")")
   })
   return(entries.lab)
} 
   
# --------------------------------------------------------------------------------------------------
# Call the functions for every model

# Strings to replace once you know the best run


bestEstLhoods=1:length(curModels)*0
bestObsLhoods=1:length(curModels)*0
bestParams=vector("list", length=length(curModels))
   
curcase=0

for (model in curModels) {
    curcase=curcase +1 
    #curDir=paste(baseDir,curDirs[curcase], sep="")
    #setwd(curDir)
    cat("----------------------------------------------------------------------------------------\n")
    cat(paste("Working on model ", model, "\n", sep=""))
    cat("----------------------------------------------------------------------------------------\n")

    out_dir                     = str_glue("{prefix_out}_{model}")
    all_lhoods_path             = str_glue("{out_dir}/{model}.lhhoods")
    plot_maxEstlhood_path       = str_glue("{out_dir}/{model}_MaxEstLhoodCurve")

    bestlhoods_path             = str_glue("{str_basename}.bestlhoods")
    brentlhoods_path            = str_glue("{str_basename}.brent_lhoods")

    popLabels=pop.Labels[[curcase]]

    system( str_glue("mkdir -p {out_dir} ") ) 


    bestRun = extract_best_lhoods(numRuns = numRuns, model = model,
                            str_basename = str_basename, 
                            bestlhoods_path = bestlhoods_path, 
                            all_lhoods_path = all_lhoods_path, 
                            brentlhoods_path = brentlhoods_path,
                            plot_maxEstlhood_path = plot_maxEstlhood_path,
                            printPDF = printPDF,
                            printPNG = printPNG)


    files_basename              = str_glue(gsub("{i}", "{bestRun}", basename(str_basename), fixed = T))
    bestrun_bestlhoods_path     = str_glue("{out_dir}/{files_basename}.bestlhoods")
    sfsFile.obs                 = str_glue("{model}_DSFS.obs")
    sfsFile.exp                 = sub(".bestlhoods", "_DSFS.txt", bestrun_bestlhoods_path)
    
    if(drawMarginalSFS){

        sfs_sit_path                = str_glue("{out_dir}/{model}_sfs_fit")
        sfs_worse                   = str_glue("{out_dir}/{model}_worst_sfs_entries")
        plot_basename               = str_glue("{out_dir}/{model}_{bestRun}_sfs_fit")
        get_marginal_SFS(
            bestrun_bestlhoods_path = bestrun_bestlhoods_path, 
            model = model,
            sfsFile.obs = sfsFile.obs, 
            sfsFile.exp = sfsFile.exp,
            pop.Labels = pop.Labels, 
            curcase = curcase, 
            printPDF = printPDF, printPNG = printPNG,
            plot_basename = plot_basename
                             )

    }

    if(plotWorstFitSFS){
        # {model}_{n}_worse_sfs_entries.{pdf,png}
        plot_basename               = str_glue("{out_dir}/{model}_{numBadEntriesToPLot}_worse_sfs_entries")
        plot_worst_entries(
            sfsFile.obs = sfsFile.obs,
            sfsFile.exp = sfsFile.exp,
            printPDF = printPDF, printPNG = printPNG,
            plot_basename = plot_basename
            )
    }
}