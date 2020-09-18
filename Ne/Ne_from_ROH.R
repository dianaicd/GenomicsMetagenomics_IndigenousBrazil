obs_sharing <- function(sample, roh, min_l_Mb = 1.6,
                        genome_size_cM = 3545.83, haploid_inds = 2){
  sharing <- sum( roh$KB[ roh$IID == sample & roh$KB >= (min_l_Mb*1000) ]/1000 ) / ( genome_size_cM * choose(n = haploid_inds, 2))
  sharing
}

Ne_estimate <- function(sharing, min_l_Mb = 1.6){
  Ne <- (50 * (1 - sharing + sqrt( 1 - sharing )) ) / ( min_l_Mb * sharing)
  Ne
}

get_Ne <- function(sample, roh){
  sample_sharing <- obs_sharing(roh = roh, 
                                sample = sample)
  Ne <- Ne_estimate(sample_sharing, 1.6)
  return(Ne)
}
