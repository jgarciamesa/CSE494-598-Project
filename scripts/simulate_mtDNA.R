

main = function(n) {
  # read mtDNA haplogroup database
  dataset = read.table("data/mtDNA/mtDNA_dataset.tsv", header = FALSE, stringsAsFactors = FALSE)
  
  # vector of random numbers
  random = runif(n*2,1,nrow(dataset))
  
  # synthetic fasta and record of "known" profile files
  fasta = paste0("results/fasta/mixed_",n,".fasta")
  known_profile = paste0("results/fasta/mixed_",n,".txt")
  
  # create empty synthetic fasta file
  file.create(fasta)
  
  # add 2 mtDNA haplogroups per individual (2*n)
  for(i in 1:(2*n)) {
    download.file(paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",dataset$V1[random[i]],"&rettype=fasta&retmode=text"), destfile=fasta, mode = "a")
  }
  
  # create empty known profile file
  file.create(known_profile)
  
  # save two haplogroups
  writeLines("mtDNA Haplogroups", known_profile)
  cat(dataset[random[1],]$V2, file=known_profile, append=TRUE, sep = "\n")
  cat(dataset[random[2],]$V2, file=known_profile, append=TRUE, sep = "\n")
}

if(!interactive()) {
  ARGS = commandArgs(trailing=TRUE)
  main(ARGS[1])
}