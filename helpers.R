get.plot.number <- function(name="Rplot"){
  dir.list = dir(path="./",pattern=paste("^", name, sep=''))  #all dirs & files starting with "Rplot"
  if (length(dir.list)==0){return(1)} else{
    dir.last = tail(dir.list, 1)  #the last item of the list
    num.plot = strsplit(dir.last, ".", fixed=TRUE)[[1]][1]
    return(as.numeric(gsub(name, "", num.plot, fixed=TRUE))+1)
  }
}