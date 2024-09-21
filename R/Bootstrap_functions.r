
#' @export

bootstrapS <- function (B, meta.info, group.name ,paired)
{
  groups <- meta.info[,group.name]
  bootsamples <- matrix(nrow = B, ncol = length(groups))
  for (i in seq_len(B)) {
    for (g in unique(groups)) {
      g.names <- row.names(meta.info)[which(groups == g)]
      bootsamples[i, which(groups == g)] <- sample(g.names, length(g.names), replace = TRUE)
    }
  }
  if (paired) {
    for (i in 1:B) {
        g.names1 <-  bootsamples[i, which(groups == unique(groups)[1])]

        g.names2 <- match(g.names1 , row.names(meta.info)) + length(g.names1)

        bootsamples[i, which(groups == unique(groups)[2])] <- row.names(meta.info)[g.names2]

    }
  }
  return(bootsamples)
}


#' @export
permutatedS <- function (meta.info, B)
{
  persamples <- matrix(nrow = B, ncol = nrow(meta.info))
  for (i in seq_len(B)) {
    persamples[i, ] <- sample(row.names(meta.info))
  }
  return(persamples)
}



#' @export
bootstrapSamples.limRots <- function (data, B, meta.info ,group.name)
{
  labels <- as.numeric( meta.info[,group.name] )
  samples <- matrix(nrow = B, ncol = length(labels))
  for (i in 1:B) {
    for (label in unique(labels)) {
      pos <- which(labels == label)
      meta.info.pos <- meta.info[meta.info[,group.name] == label,]
      meta.info.factors <- c()
      for (j in 1:ncol(meta.info)){
        if(is.factor(meta.info.pos[,j])){
          meta.info.factors <- c(meta.info.factors, colnames(meta.info.pos)[j])
        }
      }
      meta.info.factors <- meta.info.factors[meta.info.factors != group.name]
      meta.info.pos$stratum <- interaction(meta.info.pos[,meta.info.factors])
      stratum_sizes <- table(meta.info.pos$stratum)
      stratum_samples <- round(length(pos) * prop.table(stratum_sizes))
      sampled_indices <- unlist(lapply(names(stratum_samples), function(stratum) {
        stratum_indices <- row.names(meta.info.pos)[which(meta.info.pos$stratum ==  stratum )]
        sample(stratum_indices, stratum_samples[stratum], replace = TRUE)
      }))
      samples[i, pos] <- sampled_indices
    }
  }
  return(samples)
}


