% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coexpression_ppi_auc.R
\name{coexpression_ppi_auc}
\alias{coexpression_ppi_auc}
\title{AUC of a co-expression networks using STRING protein-protein interactions (PPIs).}
\usage{
coexpression_ppi_auc(
  net,
  gene_annot = NULL,
  symbol_col = "gene_name",
  directed = FALSE,
  string_version = "10",
  species = 9606,
  string_dir = "",
  max_homology_bitscore = Inf,
  benchmark_pathway = NULL
)
}
\arguments{
\item{net}{matrix. A gene x gene matrix where each entry represents the weight of the edge between corresponding genes. For a directed network, rows correspond to source genes and columns correspond to target genes. 
Row and column names of \code{net} should be HGNC gene symbols or gene ids. If they are not HGNC gene symbols, gene annotations (\code{gene_annot}) must be provided.}

\item{gene_annot}{data.frame. Gene annotations, where the row names correspond to gene ids and one of it columns (\code{symbol_col}) contains HGNC gene symbols.}

\item{symbol_col}{character or integer. Name or index of the column containg HGNC gene symbol in \code{gene_annot}.}

\item{directed}{logical. Directed interactions? If TRUE, it is assumed that the g1_col gene regulates the g2_col gene.}

\item{string_version}{character. Version of STRING database.}

\item{species}{numeric. Species code for STRING (9606 for homo sapiens).}

\item{string_dir}{character. Local path where STRING database will be stored locally so that it can be used offline. If \code{string_dir = ""}, a temporary directory is used.}

\item{max_homology_bitscore}{numeric. Maximum homology score. \code{0} to filter all interactions between homologous proteins, \code{Inf} not to filter any interaction.}

\item{benchmark_pathway}{character. Pathway type to benchmark. \code{NULL} to use all pathway types. Examples: \code{"KEGG" / "REACTOME" / "BIOCARTA" / "Disease" / "Pfam" / "NCI" / "ECOCYC"}}
}
\value{
Returns a list object with the following items.
\item{pr}{Precision curve object.}
\item{roc}{Roc curve object.}
}
\description{
This function computes area under the curve (AUC) of a co-expression network using protein-protein interactions (PPIs) from STRING.
}
\details{
This function computes the area under the ROC curve and the precision-recall curve. The area is computed using the \code{PRROC} package using weighted PPI scores. 
Note: higher values in \code{net} should correspond to higher confidence in the edge. If the sign of values mean positive or negative association between genes, you probably should provide absolute values. 

This function does not yet handle directed networks.
}
\examples{
require('reshape2')
genes = c('TP53', 'RBM3', 'SF3', 'LIM12', 'ATM', 'TMEM160', 'BCL2L1', 'MDM2', 'PDR', 'MEG3', 'EGFR', 'CD96', 'KEAP1', 'SRSF1', 'TSEN2')
dummy_net = matrix(rnorm(length(genes)^2), nrow = length(genes), dimnames = list(genes, genes))
dummy_net = (dummy_net + t(dummy_net))/2 # symmetric undirected nework
auc_res = coexpression_ppi_auc(net = dummy_net)
print(sprintf('AUC under ROC curve: \%s', auc_res$roc$auc))
print(sprintf('AUC under PR curve: \%s', auc_res$pr$auc.integral))
plot(auc_res$roc)
plot(auc_res$pr)
}
