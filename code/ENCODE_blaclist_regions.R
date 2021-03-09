ENCODE_blaclist_regions <- function(path_to_dir) {
  # regions in the human genome that have anomalous, unstructured, high signal/read counts in NGS experiments independent of cell line and type of experiment
  # empirical from many ENCODE cell lines (not RNA-seq)
  
  path <- paste0(path_to_dir, "/blacklists/")
  DAC <- rtracklayer::import(paste0(path, "wgEncodeDacMapabilityConsensusExcludable.bed.gz"))
  Duke <- rtracklayer::import(paste0(path, "wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))

  blacklist = union(DAC, Duke)
  strand(blacklist) <- "+"

  blacklist
}
