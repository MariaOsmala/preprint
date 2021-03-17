distToTss <-function(peak, tss){
  
  peak <- resize(peak, width=1, fix="center")
  idx <- nearest(peak, tss) #(x, subject) returns indexes of the nearest subject from x
  naind=which(is.na(idx))
  notnaind=which(!is.na(idx))
  if(length(naind)!=0){
    idx=idx[-naind]
  }
  sgn <- as.numeric(ifelse(strand(tss)[idx] == "+", 1, -1))
  distances=data.frame(rep(5000000, length(peak)))
  distances[notnaind,]=(as.data.frame(ranges(peak))$start[notnaind] - as.data.frame(ranges(tss))$start[idx]) * sgn
  distances
}
