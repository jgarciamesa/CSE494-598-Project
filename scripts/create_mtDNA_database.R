library(rvest)
library(rentrez)

# vector with all the haplotypes
haplogroups = c("A","B","C","D","E","F","G","H","HV","I","J","K","L0","L1","L2","L3","L4","L5","L6","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")

# data frame to contain the GenBank ID and haplogroup for sequences from mitomap
all_haplotypes = vector()

for(group in haplogroups) {
  tmp = paste0("https://www.mitomap.org/cgi-bin/haplo_distribution?hg=",group,"&hpp=10000") %>% read_html() %>% html_nodes(xpath = '//*[@id="haplo-list"]/table') %>% html_table()
  all_haplotypes = rbind(all_haplotypes, tmp[[1]][,c(2,3,4)])
}

all_haplotypes$Haplogroup = gsub(" ","-",all_haplotypes$Haplogroup)

write.table(all_haplotypes, file = "data/mtDNA/mtDNA_dataset.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)#, quote = FALSE)
