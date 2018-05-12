get.plot.number <- function(){
  dir.list = dir(path="./",pattern='^Rplot')  #all dirs & files starting with "Rplot"
  if (length(dir.list)==0){return(1)} else{
    dir.last = tail(dir.list, 1)  #the last item of the list
    num.plot = strsplit(dir.last, "-")[[1]][1]
    return(as.numeric(gsub("Rplot", "", num.plot))+1)
  }
}