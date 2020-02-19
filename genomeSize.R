#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
          args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      genomeSize.R - Estimate genome size based on Jellyfish histogram data
 
      Arguments:
      --in	- input file
      --help    - print this Help

      Example:
      ./genomeSize.R --in=\"input1.histo\" 
      
      Daniel Guariz Pinheiro
      FCAV/UNESP - Univ Estadual Paulista
      dgpinheiro@gmail.com
      \n\n")

  q(save="no")
}


find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
     pks <- unlist(pks)
     pks
}

#https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data


## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL[['in']])) {
        sink(stderr())
        cat("\nERROR: Missing input file (--in) !\n\n")
        sink()
        q(save="no")
}

k.df <- read.table(argsL[['in']]) #load the data into dataframe

peak.n <- find_peaks(k.df[1:dim(k.df)[1],'V2'],m=1)[1]

valley.first.n <- find_peaks( k.df[peak.n,'V2']-k.df[1:dim(k.df)[1],'V2'] , m=1)[1]

valley.second.n <- find_peaks(k.df[peak.n:dim(k.df)[1],'V2']-k.df[peak.n,'V2'], m=1)[1]

#plot(k.df[,c('V1','V2')],type="l", xlim=c(valley.first.n,valley.second.n), ylim=c(0,max(k.df[valley.first.n:valley.second.n,'V2'])) )
#points(k.df[peak.n,c('V1','V2')],col="red")

cat(paste("First valley: ", valley.first.n, " (",k.df[valley.first.n,'V2'],")\n", sep=""))
cat(paste("Peak: ", peak.n, " (",k.df[peak.n,'V2'],")\n", sep=""))
cat(paste("Second valley: ", valley.second.n, " (",k.df[valley.second.n,'V2'],")\n", sep=""))

genomesize <- sum(as.numeric(k.df[valley.first.n:dim(k.df)[1],'V1']*k.df[valley.first.n:dim(k.df)[1],'V2']))/peak.n
#genomesize <- sum(as.numeric(k.df[2:dim(k.df)[1],1]*k.df[2:dim(k.df)[1],2]))/peak.n
#genomesize <- sum(as.numeric(k.df[valley.first.n:valley.second.n,1]*k.df[valley.first.n:valley.second.n,2]))/peak.n

cat(paste("Genome size: ", genomesize/1000000, " Mb\n\n",sep=""))
