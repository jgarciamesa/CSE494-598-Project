# Usage: Rscript mtdna_tally.R

main = function() {
	library(dplyr)

	all_haplogrep_files = list.files(".","haplogrep_[0-9]+.txt", recursive = TRUE) %>% basename()
	N = gsub("[A-z.]", "", all_haplogrep_files) %>% unique()

	if(file.exists("simulations/mtdna_tally.txt")) {
		file.remove("simulations/mtdna_tally.txt")
	}

	for(n in N) {
		haplogrep_files = list.files(".", paste0("haplogrep_",n,".txt"), recursive = TRUE)
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
		cat(paste0("mtDNA identification score with ",gsub("[A-z.]",'', basename(haplogrep_files[1]))," individuals and ",length(haplogrep_files)," runs is ",sum(score)/length(haplogrep_files)), file="simulations/mtdna_tally.txt", append = TRUE, sep = "\n")
	}

}

if(!interactive()) {
	main()
}
