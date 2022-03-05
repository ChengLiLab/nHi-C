.get.init <- function (data, nstate = 2, iter.max = 20, nstart = 3) 
{
  cluster <- kmeans(data, nstate, iter.max = iter.max, nstart = nstart)
  means <- cluster$centers
  vars <- cluster$withinss/(cluster$size - 1)
  ret <- list()
  ret[["mean"]] <- means
  ret[["var"]] <- vars
  ret[["cluster"]] <- cluster$cluster
  ret[["size"]] <- cluster$size
  ret
}

hmm.setup.our <- function (data, state = c("enriched", "non-enriched"), probe.region = 35, 
                           frag.size = 1000, pos.state = 1, em.type = "tDist", max.prob = 1, 
                           df = 9) 
{
  comp <- numeric(length(state)) + 1
  names(comp) <- state
  if (is(data, "list")) 
    data <- c(data, recursive = TRUE)
  #cl <- .get.init(data, nstate = c(min(data), max(data)))
  cl <- .get.init(data, nstate = 2)
  if(cl$mean[1] > cl$mean[2]){
    mean.pos <- cl$mean[1]
    mean.neg <- cl$mean[2]
    var.pos <- cl$var[1]
    var.neg <- cl$var[2]
  }else{
    mean.pos <- cl$mean[2]
    mean.neg <- cl$mean[1]
    var.pos <- cl$var[2]
    var.neg <- cl$var[1]
  }
  
  weight.pos <- 1
  weight.neg <- 1
  emission.comp <- list()
  emission <- list()
  emission.comp[[state[1]]] <- cbind(weight.pos, mean.pos, 
                                     var.pos)
  emission.comp[[state[2]]] <- cbind(weight.neg, mean.neg, 
                                     var.neg)
  for (i in 1:length(state)) {
    if (em.type == "tDist") {
      emission[[state[i]]] <- new("tDist", mean = emission.comp[[state[i]]][, 
                                                                            2], var = emission.comp[[state[i]]][, 3], df = df[(i%%length(df)) + 
                                                                                                                                1])
    }
    else {
      stop(paste("Requested emission distribution of type", 
                 em.type, "is not supported."))
    }
  }
  if(cl$mean[1] > cl$mean[2]){
    pos.index <- cl$cluster == 1
  }else{
    pos.index <- cl$cluster == 2
  }
  
  neg.index <- !pos.index
  freq.pos <- sum(pos.index)/length(data)
  freq.pos <- freq.pos * (probe.region/frag.size)
  if (pos.state == 1) {
    prob.aa <- min(1 - probe.region/frag.size, max.prob)
    prob.bb <- min(1 - freq.pos, max.prob)
  }
  else {
    prob.aa <- min(1 - freq.pos, max.prob)
    prob.bb <- min(1 - probe.region/frag.size, max.prob)
  }
  neg.state <- 2
  if (pos.state == 2) 
    neg.state <- 1
  transition <- list()
  transition[[state[1]]] <- new("discDist", alpha = state, 
                                prob = c(prob.aa, 1 - prob.aa))
  transition[[state[2]]] <- new("discDist", alpha = state, 
                                prob = c(1 - prob.bb, prob.bb))
  a <- transition[[state[1]]][2]
  b <- transition[[state[2]]][1]
  pi <- new("discDist", alpha = state, prob = c(b/(a + b), 
                                                a/(a + b)))
  hmm.init <- new("contHMM", transition = transition, emission = emission, 
                  init = pi)
  return(hmm.init)
}
getNAD <- function(mat){
  ## step1 get observe sequence
  obs <- mat$ratio
  ## step2 generate hmm model 
  hmm.init <- hmm.setup.our(obs, state = c("enriched","non-enriched"), pos.state = 1, probe.region = 1000, frag.size = 3000)
  show(hmm.init)
  plot(hmm.init)
  
  ## step2 parameter optimisation 
  max.gap <- 2000000
  obs.gap <- diff(mat[['position']]) > max.gap
  gap.idx <- which(obs.gap)
  if(length(gap.idx) !=0){
    start <- c(1, gap.idx + 1)
    end <- c(gap.idx, length(obs))
    obs.lst <- mapply(function(s, e, data) data[s:e], start, end,
                      MoreArgs = list(obs))
    hmm.opt <- viterbiEM(hmm.init, obs.lst, df = 9)
    mat.lst <- split(mat, cumsum(1:nrow(mat) %in% (end+1)))
    #mat.lst <- mapply(function(s,e,data) data[s:e,], start, end, MoreArgs = list(mat) )
    #print(hmm.opt)
  }else{
    obs.lst <- list(obs)
    hmm.opt <- hmm.init
    mat.lst <- list(mat)
    print(hmm.opt)
  }
  #hmm.opt <- hmm.init
  #plot(hmm.opt)
  
  
  ## step3 calling probes
  post.lst <- lapply(obs.lst, posterior, hmm.opt)
  state.seq.lst <- lapply(post.lst, apply, 2, which.max)
  regions <- data.frame(chr = as.character(), start = as.numeric(), end = as.numeric())
  callregion <- function(state.seq, post, mat){
    state.seq <- states(hmm.opt)[state.seq]
    regions.idx <- region.position(state.seq, region = 'enriched')
    regions.pos <- matrix(mat[regions.idx, 2],
                          nrow = 2, ncol = dim(regions.idx)[2])
    post.enriched <- exp(post)
    post.enriched <- post.enriched[seq(2,2*nrow(mat),2)]
    region.score <- apply(regions.idx, 2,
                          function(reg, post) mean(post[reg[1]:reg[2]]), post.enriched)
    region.len <- apply(regions.pos, 2, diff)
    regions.clean <- remove.short(regions.idx, post.enriched,
                                  mat[ , 1:2], min.length = 1000, min.score = 0.8)
    regions.pos.clean <- data.frame(chr=mat[regions.clean[1,], "chromosome"], 
                                    start=mat[regions.clean[1,], "position"],
                                    end = mat[regions.clean[2,], "position"])
    return(regions.pos.clean)
  }
  
  
  for(i in 1:length(state.seq.lst)){
    regions <- rbind(regions, callregion(state.seq.lst[[i]], post.lst[[i]], mat.lst[[i]]))
  }
  print(dim(regions))
  return(regions)
}

