### sample code
library('ioutil')
library('genomicsutil')
library('miscutil')
library('mappabilityutil')

### files are on marcc
expr_fn = "/work-zfs/abattle4/ashis/progdata/misc/cross_mappability/gtex_v7/processed_expression/Whole_Blood.v7.normalized_expression.txt"
tpm_fn = "/work-zfs/abattle4/lab_data/GTEx_v7/rna_seq_from_portal/gene_tpm/WholeBlood.txt"
count_fn = "/work-zfs/abattle4/lab_data/GTEx_v7/rna_seq_from_portal/gene_count/WholeBlood.txt"
gene_annot_fn = "/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"
mappability_fn = "/work-zfs/abattle4/lab_data/annotation/mappability_hg19_gencode19/avg_mappability_Exon_UTR.txt"

# show time required to read data
system.time({expr_df = read.table(expr_fn, sep = '\t', header = T, row.names = 1, stringsAsFactors=F, check.names=F, quote="", comment.char="")})
system.time({expr_df = read_df(expr_fn)})


# read data
expr_df = read_df(expr_fn)
tpm_df = read_df(tpm_fn)
count_df = read_df(count_fn)
gene_annot_df = read_df(gene_annot_fn) #, row.names = T, col.names = T)
mappability_df = read_df(mappability_fn, header = F, row.names = T)

# change sample names in tpm and count data
colnames(tpm_df) = sapply(colnames(tpm_df), function(s) paste(parse_delimitted_string(s, delim = "-")[1:2], collapse = '-') )
colnames(count_df) = sapply(colnames(count_df), function(s) paste(parse_delimitted_string(s, delim = "-")[1:2], collapse = '-') )

# show head of expr_df
head_df(expr_df)
dim(expr_df)

# take only autosomal genes
expr_df = filter_expr_by_chr(expr_df, annot.gene = gene_annot_df, chr.include = c(1:22))
dim(expr_df)

# take only protein-coding genes
expr_df = filter_expr_by_gene_type(expr_df, annot.gene = gene_annot_df, type.col = 'gene_type', type.values = 'protein_coding')
dim(expr_df)

# take genes with mappability >= 0.8
expr_df = filter_expr_by_mappability(expr_df, annot.mappability = mappability_df, min.mappability = 0.8, mappability.col = 1)
dim(expr_df)

# filter on min tpm and min count in min samples
expr_df = filter_expr_by_tpm_read(expr_df, tpm.df = tpm_df, count.df = count_df, min.tpm = 1, min.count = 10, min.samples = 20)
dim(expr_df)

# filter on gene coefficient of variation
expr_df = filter_expr_by_coeff_of_variation(expr_df, raw.df = tpm_df, n = 1000)
dim(expr_df)

