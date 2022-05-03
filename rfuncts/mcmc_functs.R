# trace loop
trace_loop <- function(theta.samps, j, type_j, K, B){
  samps <- theta.samps[[j]]
  tps <- list()
  if (type_j == "c"){
    for (h in 1:K) {
      var <- samps[h, 1, -(1:B)]
      R <- length(var)
      g <- ggplot(data.frame(var = var, time = 1:R), aes(x = time, y = var)) + geom_point() + xlab("iters") + ylab("mean") +
        labs(title = paste0("h=", toString(h)))
      tps[[h]] <- g
    }
  }
  else if (type_j == "b") {
    # histogram loop
    for (h in 1:K) {
      var <- samps[-(1:B), h]
      R <- length(p)
      g <- ggplot(data.frame(var = var, time = 1:R), aes(x = time, y = var)) + geom_point() + xlab("iters") + ylab("mean") +
        labs(title = paste0("h=", toString(h)))
      tps[[h]] <- g
    }
  }
  return(tps)
}

hist_loop <- function(theta.samps, j, type_j, K, B) {
  samps <- theta.samps[[j]]
  hists <- list()
  class_names <- c("Alpha Class", "Gamma Class", "Delta Class", "New Class")
  if (type_j == "c"){
    m <- mean(samps[, 1, -(1:B)])
    sd_m <- sd(samps[,1,-(1:B)])
    for (h in 1:K) {
      var <- samps[h, 1, -(1:B)]
      g <- ggplot(data.frame(var = var), aes(x = var)) + geom_histogram(color = "blue", fill = "skyblue1") +
        xlab(class_names[h]) + xlim(m - 3 * sd_m, m + 3*sd_m)
      hists[[h]] <- g
    }
  }
  else if (type_j == "b") {
    # histogram loop
    for (h in 1:K) {
      var <- samps[-(1:B), h]
      g <- ggplot(data.frame(var = var), aes(x = var)) + geom_histogram()
      hists[[h]] <- g
    }
  }
  return(hists)
}