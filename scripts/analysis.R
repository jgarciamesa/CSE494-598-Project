

# x = variants from known profile
# y = variants found by haplogrep

# score = which(x=y)/x

# Usage: Rscript --vanilla --quiet scripts/analysis.R $(find -name haplogrep_n.txt)

main = function(haplogrep_files) {
  known_files = gsub("haplogrep", "fasta/mixed", haplogrep_files)
  score = vector()
  
  for(i in 1:length(haplogrep_files)) {
    
    known_variants_list = read.csv(known_files[i], skip = 2, header = FALSE, nrows = 1, stringsAsFactors = FALSE)
    known_variants = unlist(known_variants_list, use.names = FALSE)
    
    haplogrep_result_dataframe = read.csv(haplogrep_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    haplogrep_result = strsplit(haplogrep_result_dataframe$Found_Polys, " ")[[1]]
    
    count = 0
    
    for(j in known_variants) {
      if(trimws(j) %in% haplogrep_result) {
        count = count + 1
      }
    }
    
    score = rbind(score, count / length(known_variants))
    
  }
  #print(sum(score)/length(haplogrep_files))
  #print(gsub("[A-z.]",'', basename(haplogrep_files[1])))
  cat(paste0("mtDNA identification score with ",gsub("[A-z.]",'', basename(haplogrep_files[1]))," individuals and ",length(haplogrep_files)," runs is ",sum(score)/length(haplogrep_files)), file="analysis.txt", append = TRUE, sep = "\n")
  
}

if(!interactive()) {
  ARGS = commandArgs(trailing=TRUE)
  main(ARGS)
}