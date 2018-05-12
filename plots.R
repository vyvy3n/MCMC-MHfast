source("helpers.R")

plt.legend <- function(acc.list, loc="topright"){
  color = c("red", "darkorange", "blue", "darkseagreen3")
  legend(loc, col = color, lty = c(1,1,1,1), bty = "n",
         legend = c(paste(format(acc.list[[1]], digits=3), "ULA"),
                    paste(format(acc.list[[2]], digits=3), "ULA large"),
                    paste(format(acc.list[[3]], digits=3), "MALA"),
                    paste(format(acc.list[[4]], digits=3), "MRW")))
}

# ########### 
# #  Error  #
# ###########__________________________________________________

cal.err<- function(out, iterSeq) {
  err = matrix(nrow=length(iterSeq), ncol=1)
  count = 0
  for (j in iterSeq) {
    count = count + 1
    theta = 0
    for (i in 1:n.exps) {
      theta = theta + out[[i]][, j]
    }
    err[count] = mean(abs(theta/n.exps - Theta))
  }
  return(err)
}

plt.err <- function(chain.list, acc.list, params, 
                    save=TRUE, num=0, logScale=TRUE, maxiter=FALSE) {  
  if (maxiter!=FALSE) {
    params$n.iters = min(params$n.iters, maxiter)
  }
  iterSeq = seq(1, params$n.iters, by=5)#length.out = 200)
  color = c("red", "darkorange", "blue", "darkseagreen3")
  par(mfrow = c(1, 1))
  if (save==TRUE){
    num = get.plot.number("err")
    dev.copy(png, paste("err", num, ".png", sep=""), width =500, height = 400)
  }
  if (log==TRUE){
    plot( iterSeq, log(cal.err(chain.list[[1]], iterSeq)), col=color[1],
          type='l', ann=FALSE, ylim=c(-1, 0))
    lines(iterSeq, log(cal.err(chain.list[[2]], iterSeq)), col=color[2])
    lines(iterSeq, log(cal.err(chain.list[[3]], iterSeq)), col=color[3])
    lines(iterSeq, log(cal.err(chain.list[[4]], iterSeq)), col=color[4])
  } else {
    plot( iterSeq, cal.err(chain.list[[1]], iterSeq), col=color[1],
          type='l', ann=FALSE, ylim=c(0, 1))
    lines(iterSeq, cal.err(chain.list[[2]], iterSeq), col=color[2])
    lines(iterSeq, cal.err(chain.list[[3]], iterSeq), col=color[3])
    lines(iterSeq, cal.err(chain.list[[4]], iterSeq), col=color[4])
  }
  title(paste("n.samples=", params$n.samples, ", n.exps=", params$n.exps, 
              ", pre-cond=", params$pre.cond, ", alpha=", params$alpha, 
              ", lazy=", lazy.step, sep=''),
        xlab="iterations", ylab="L1 error")
  plt.legend(acc.list)
  
  if (save==TRUE){dev.off()}
}
