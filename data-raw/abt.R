## code to prepare `abt` dataset goes here
#load dataset
load('abt.rda')

#Subset the data to patient F
abt = subset_samples(ps,ind == "F")

#find lowest known level of classification for each taxa
tax = tax_table(abt)
for(i in 1:nrow(tax)){
  name = tax[i,1]
  for(j in 2:ncol(tax)){
    if(tax[i,j] %in% c("","uncultured","Incertae Sedis")){
      tax[i,j] = name
    } else {
      name = tax[i,j]
    }
  }
  print(i/nrow(tax))
}
colnames(tax)
tax_table(abt) = tax

#merge reads at a higher level
ps_merge  = tax_glom(abt,taxrank = "Taxon_8",NArm = TRUE)

#filter for taxa which appear at least 100 times across all 54 samples
ps = filter_taxa(ps_merge,function(x) sum(x) > 100, TRUE)


usethis::use_data(abt, overwrite = TRUE)
