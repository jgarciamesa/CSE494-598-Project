
# read mixed_n.txt

# read haplogrep_n.result

# x = variants from known profile
# y = variants found by haplogrep

# score = which(x=y)/x

known_variants_list = read.csv("results/fasta/mixed_5.txt", skip = 2, header = FALSE, nrows = 1, stringsAsFactors = FALSE)
known_variants = unlist(known_variants_list, use.names = FALSE)

haplogrep_result_dataframe = read.csv("results/haplogrep_5.txt", header = TRUE, sep = "\t")
haplogrep_result = strsplit(haplogrep_result_dataframe$Found_Polys, " ")[[1]]

count = 0

for(i in known_variants) {
  if(trimws(i) %in% haplogrep_result) {
    count = count + 1
  }
}

score = count / length(known_variants)

writeLines()
