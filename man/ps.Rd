% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ps}
\alias{ps}
\title{The microbiome data from Dethlefsen and Relman 2011.}
\format{
A phyloseq object with an otu table, sample data, taxonomic table, and a phylogenetic tree.  
\describe{
  \item{otu_table()}{A count matrix detailing the counts of 44 unique taxa across 54 samples.}
  \item{sample_data()}{For each sample, records the individual, treatment condition, and timepoint.}
  \item{tax_table()}{For each sequencing read, records the Domain, Phylum, and Taxon_3 through Taxon_8 to the finest known level.}
  \item{phy_tree()}{Contains an edge matrix encoding the phylogenetic tree, the number of internal nodes, leaf labels, edge lengths, and node labels.}
}
}
\source{
\url{https://github.com/krisrs1128/microbiome_plvm/tree/master/data/antibiotics-study}
}
\usage{
ps
}
\description{
A dataset containing microbiome data from the experiment detailed in Dethlefsen and Relman 2011, where three patients were given antibiotic regimens and the microbiome response was monitored.  The data has been filtered to only include Patient F, taxa were merged at their lowest known level of classification, and taxa were pruned if did not appear at least 100 times across all 54 samples.
}
\references{
Les Dethlefsen and David A. Relman.  Incomplete recovery and individualized responses of the human distal gut microbiota to repeated antibiotic perturbation.  Proceedings of the National Academy of the Sciences of the United States of America.  18(Supplement 1): 4554-4561, 2011.
}
\keyword{datasets}
