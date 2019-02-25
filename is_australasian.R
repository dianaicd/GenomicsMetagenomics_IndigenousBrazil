is_australasian <- function(pop){
  australasians <- c("Andaman", "Papuan", "WCD05", "WCD02", "WCD03", "WCD08", "WCD01", "WCD13", "WCD04")
  if(pop %in% australasians){
    return("Australasian")
  }else{
    return("Other")
  }
}