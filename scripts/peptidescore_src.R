do.hash.table <- function(mat = div.mat,
                          grid = som.grid,
                          params = simhash_params) {
  do.som.hash <- function(mat,k,B,sgrid) {
    .train.mat <- mat[sample(B),]
    .model <- som(.train.mat, 
                  grid = sgrid)
    .hs <- predict(object = .model,
                   newdata = mat)[["unit.classif"]]
    return(.hs)
  }
  
  .hash.table <- lapply(seq(params$L),function(i){
    ## compute first layer
    htable <- cbind(rep(0,params$N),
                    matrix(do.som.hash(mat,
                                       params$k,
                                       params$B,
                                       grid),
                           ncol=1))
    hs <- rowSums(sweep(htable,2,(1+params$k^2)^seq(0,ncol(htable)-1),'*'))
    hss <- table(hs)
    
    ## compute extra layers until all contain less than B cells
    while(any(hss>params$B)) {
      for(i in as.numeric(names(which(hss>params$B)))) {
        htable[hs==i,1] <- do.som.hash(mat[hs==i,],
                                       params$k,
                                       params$B,
                                       grid)
      }
      hs <- rowSums(sweep(htable,2,(1+params$k^2)^seq(0,ncol(htable)-1),'*'))
      hss <- table(hs)
      htable <- cbind(rep(0,params$N),
                      htable)
    }
    return(hs)
  })
  .hash.table <- do.call(cbind,.hash.table)
  rownames(.hash.table) <- rownames(mat)
  .hash.table
}

do.hash.dist <- function(mat = div.mat,
                         htable = hash.table,
                         params = simhash_params) {
  
  .div.dist <- lapply(seq(params$L),function(cur.hash) {
    cur.buckets <- htable[,cur.hash]
    bks <- sort(unique(cur.buckets))
    cur.neighbors <- lapply(bks,function(bk) {
      # cat('\tbucket #',bk,'\n')
      .ms <- mat[which(cur.buckets==bk),]
      .ms <- t(.ms)
      .res <- crossprod(.ms)/(sqrt(tcrossprod(colSums(.ms^2))))
      
      diag(.res) <- 0 # avoid self loops
      
      .res.min <- apply(.res,1,rank,
                        ties.method = "first")
      .res.min <- nrow(.res)+1-.res.min
      .res.min <- which(.res.min<=params$K,
                        arr.ind = TRUE)
      .res.min <- data.frame(from=dimnames(.res)[[2]][.res.min[,2]],
                             to=dimnames(.res)[[1]][.res.min[,1]],
                             weight=.res[.res.min])
      
      return(.res.min)
    })
    cur.neighbors <- bind_rows(cur.neighbors) %>%
      mutate(Hash=cur.hash)
    return(cur.neighbors)
  })
  .div.dist <- bind_rows(.div.dist)
  
  .div.dist %>%
    group_by(from,to) %>%
    summarize(weight=max(weight)) %>% 
    group_by(from) %>% 
    mutate(rweight = length(to)+1-rank(weight,ties.method = "first")) %>% 
    ungroup()
}