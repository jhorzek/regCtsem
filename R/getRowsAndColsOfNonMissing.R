getRowsAndColsOfNonMissing <- function(uniqueMissingPattern){
  nonMissing <- !uniqueMissingPattern
  indicators <- c()
  for(i in seq_len(length(nonMissing))){
    if(nonMissing[i]){
      indicators <- c(indicators, i)
    }
  }
  # for cpp:
  indicators = indicators -1
  return(indicators)
}

