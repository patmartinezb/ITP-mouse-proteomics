
# LOAD LIBRARIES ---------------------------

lib <-  c("readxl", "dplyr", "openxlsx", "goseq", "Mus.musculus",
          "org.Mm.eg.db", "msigdbr")

invisible(lapply(lib, library, character.only = TRUE))

# LOAD DATA ---------------------------

rm(list = ls())
module <- read_excel("ResultsWGCNA.xlsx", sheet = "brown") # Re-run script for every color/module


all <- read_excel("ResultsWGCNA.xlsx") # Universe


# ENRICHMENT ANALYSIS ---------------------------

# Vectors
query <- as.character(module$Protein.IDs) # Significant proteins, included present/not present proteins
query.all <- as.character(all$Protein.IDs) # Universe


# Annotate with ENTREZID and GO
txdb <- Mus.musculus

res.sig <- AnnotationDbi::select(txdb, keys = query, columns = c("ENTREZID", "GO", "SYMBOL"), keytype = "UNIPROT")

res.all <- AnnotationDbi::select(txdb, keys = query.all, columns = c("ENTREZID", "GO", "SYMBOL"), keytype = "UNIPROT")



df.sig <- all[all$Protein.IDs %in% module$Protein.IDs,] # Keep proteins that are differentialy expressed after limma


# We need mm10 mappings:
# ls("package:org.Mm.eg.db")
# columns(org.Mm.eg.db)

# Get the ENTREZIDs and the KEGG pathways ("PATH")
entrezs <- AnnotationDbi::select(org.Mm.eg.db, keys = keys(org.Mm.eg.db, keytype="UNIPROT"), 
                                 columns = c("ENTREZID","PATH", "SYMBOL"), keytype="UNIPROT")



df.all <- all[all$Protein.IDs %in% res.all$UNIPROT,] # Select proteins from df pre-limma that later are gonna be differentialy expressed - we use this dataset cause it has all the intensities, to calculate de means

df.all.means <- rowMeans(df.all[,c(2:20)], na.rm = TRUE) # Calculate the row means
names(df.all.means) <- df.all$Protein.IDs

res.all$means <- df.all.means[match(res.all$UNIPROT, names(df.all.means))] # Assign the means to every hit of the annotation, since most are repeated

res.all.uni <- res.all[!(duplicated(res.all$UNIPROT)),]
res.all.unique <- res.all.uni[!(duplicated(na.omit(res.all.uni$ENTREZID))),] # Get rid of duplicates




genes <- rep(0,length(unique(res.all.unique$ENTREZID))) # Create vector with as many 0 as genes in res.all with an ENTREZID 
names(genes) <- unique(res.all.unique$ENTREZID) # Give names to the vector
genes[names(genes) %in% res.sig$ENTREZID ] <- 1 # 1 for those that are find in the significant group

table(names(genes) %in% res.all.unique$ENTREZID)
genes.ordered <- genes[as.character(res.all.unique$ENTREZID)]

# Remove NA gene
genes.ordered <- genes.ordered[!is.na(names(genes.ordered))]
res.all.unique <- res.all.unique[!is.na(res.all.unique$ENTREZID),]          


pwf <- nullp(genes.ordered, bias.data = (as.numeric(res.all.unique$means) - min(as.numeric(res.all.unique$means))))

# Run the enrichment analysis

GO.goseq <- goseq(pwf,gene2cat = entrezs, use_genes_without_cat = TRUE) # Take bias into account
GO.goseq$over_represented_pvalue_adj <-  p.adjust(GO.goseq$over_represented_pvalue, method = "BH")


# Retrieve MSigDB categories and subcategories

m.H <- as.data.frame(msigdbr(species = "Mus musculus", category = "H"))

m.C1 <- as.data.frame(msigdbr(species = "Mus musculus", category = "C1"))

m.C2 <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2"))

m.C2.cgp <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP"))

m.C2.cp <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP"))

m.C2.bio <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:BIOCARTA"))

m.C2.kegg <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG"))

m.C2.reac <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"))

m.C2.pid <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:PID"))

m.C3 <- as.data.frame(msigdbr(species = "Mus musculus", category = "C3"))


m.C4 <- as.data.frame(msigdbr(species = "Mus musculus", category = "C4"))
m.C4.cgn <- as.data.frame(msigdbr(species = "Mus musculus", category = "C4", subcategory = "CGN"))
m.C4.cm <- as.data.frame(msigdbr(species = "Mus musculus", category = "C4", subcategory = "CM"))

m.C5 <- as.data.frame(msigdbr(species = "Mus musculus", category = "C5"))

m.C5.bp <- as.data.frame(msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP"))

m.C5.cc <- as.data.frame(msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC"))

m.C5.mf <- as.data.frame(msigdbr(species = "Mus musculus", category = "C5", subcategory = "MF"))

m.C6 <- as.data.frame(msigdbr(species = "Mus musculus", category = "C6"))

m.C7 <- as.data.frame(msigdbr(species = "Mus musculus", category = "C7"))                    

# Get the datasets into goseq() format
# Get gene2cat annotations for pathway enrichment
m.H.cat <- m.H[,c("entrez_gene","gs_name", "gs_exact_source")]
colnames(m.H.cat) <- c("ENTREZID","PATH", "ID")

m.C1.cat <- m.C1[,c("entrez_gene","gs_name")]
colnames(m.C1.cat) <- c("ENTREZID","PATH")

m.C2.cat <- m.C2[,c("entrez_gene","gs_name")]
colnames(m.C2.cat) <- c("ENTREZID","PATH")

m.C2.bio.cat <- m.C2.bio[,c("entrez_gene","gs_name")]
colnames(m.C2.bio.cat) <- c("ENTREZID","PATH")

m.C2.cgp.cat <- m.C2.cgp[,c("entrez_gene","gs_name")]
colnames(m.C2.cgp.cat) <- c("ENTREZID","PATH")

m.C2.cp.cat <- m.C2.cp[,c("entrez_gene","gs_name")]
colnames(m.C2.cp.cat) <- c("ENTREZID","PATH")

m.C2.kegg.cat <- m.C2.kegg[,c("entrez_gene","gs_name", "gs_exact_source")]
colnames(m.C2.kegg.cat) <- c("ENTREZID","PATH", "ID")

m.C2.reac.cat <- m.C2.reac[,c("entrez_gene","gs_name")]
colnames(m.C2.reac.cat) <- c("ENTREZID","PATH")

m.C2.pid.cat <- m.C2.pid[,c("entrez_gene","gs_name")]
colnames(m.C2.pid.cat) <- c("ENTREZID","PATH")

m.C3.cat <- m.C3[,c("entrez_gene","gs_name")]
colnames(m.C3.cat) <- c("ENTREZID","PATH")


m.C4.cat <- m.C4[,c("entrez_gene","gs_name")]
colnames(m.C4.cat) <- c("ENTREZID","PATH")

m.C4.cgn.cat <- m.C4.cgn[,c("entrez_gene","gs_name")]
colnames(m.C4.cgn.cat) <- c("ENTREZID","PATH")

m.C4.cm.cat <- m.C4.cm[,c("entrez_gene","gs_name")]
colnames(m.C4.cm.cat) <- c("ENTREZID","PATH")

m.C5.cat <- m.C5[,c("entrez_gene","gs_name")]
colnames(m.C5.cat) <- c("ENTREZID","PATH")
m.C5.bp.cat <- m.C5.bp[,c("entrez_gene","gs_name", "gs_exact_source")]
colnames(m.C5.bp.cat) <- c("ENTREZID","PATH", "ID")
m.C5.cc.cat <- m.C5.cc[,c("entrez_gene","gs_name", "gs_exact_source")]
colnames(m.C5.cc.cat) <- c("ENTREZID","PATH", "ID")
m.C5.mf.cat <- m.C5.mf[,c("entrez_gene","gs_name", "gs_exact_source")]
colnames(m.C5.mf.cat) <- c("ENTREZID","PATH", "ID")

m.C6.cat <- m.C6[,c("entrez_gene","gs_name")]
colnames(m.C6.cat) <- c("ENTREZID","PATH")

m.C7.cat <- m.C7[,c("entrez_gene","gs_name")]
colnames(m.C7.cat) <- c("ENTREZID","PATH")


gene2cats <- list(m.H.cat,m.C1.cat,m.C2.cat,m.C2.bio.cat,m.C2.cgp.cat,m.C2.cp.cat,
                  m.C2.kegg.cat,m.C2.reac.cat,m.C2.pid.cat,m.C3.cat,
                  m.C4.cat,m.C4.cgn.cat,m.C4.cm.cat,m.C5.cat,m.C5.bp.cat,m.C5.cc.cat,
                  m.C5.mf.cat,m.C6.cat,m.C7.cat)
names(gene2cats) <- c("H","C1","C2all","C2bio","C2cgp","C2cp","C2kegg","C2reac","C2pid",
                      "C3all","C4all","C4cgn","C4cm","C5all",
                      "C5bp","C5cc","C5mf","C6","C7")


listy <- list()
for (i in 1:length(gene2cats)){
  listy[[i]] <- goseq(pwf,gene2cat = gene2cats[[i]], use_genes_without_cat = TRUE)  # REMEMBER: CORRECTING BIAS ! (use "Hypergeometric" to avoid correction)
  listy[[i]]$over_represented_pvalue_adj <- p.adjust(listy[[i]]$over_represented_pvalue, method = "BH")
}
names(listy) <- names(gene2cats)


# Obtain significant terms
goterms_bp <- listy[["C5bp"]]$category[listy[["C5bp"]]$over_represented_pvalue_adj < 0.05]
goterms_cc <- listy[["C5cc"]]$category[listy[["C5cc"]]$over_represented_pvalue_adj < 0.05]
goterms_mf <- listy[["C5mf"]]$category[listy[["C5mf"]]$over_represented_pvalue_adj < 0.05]

keggterms <- listy[["C2kegg"]]$category[listy[["C2kegg"]]$over_represented_pvalue_adj < 0.05]

c5bp <- gene2cats[["C5bp"]]
c5cc <- gene2cats[["C5cc"]]
c5mf <- gene2cats[["C5mf"]]

kegg <- gene2cats[["C2kegg"]]

# Run function to get the symbol ids of the significant genes that accompany those terms
getGeneLists_GO <- function(df, goterms){
  cat2gene <- unstack(df, ENTREZID ~ PATH)
  
  out <- list()
  for(term in goterms){
    tmp <- rownames(pwf)[(pwf$DEgenes == 1) & (rownames(pwf) %in% cat2gene[[term]])]
    out[[term]] <- tmp
  }
  out
}


getGeneLists_KEGG <- function(df, keggerms){
  cat2gene <- unstack(df, ENTREZID ~ PATH)
  
  out <- list()
  for(term in keggterms){
    tmp <- rownames(pwf)[(pwf$DEgenes == 1) & (rownames(pwf) %in% cat2gene[[term]])]
    out[[term]] <- tmp
  }
  out
}


getGeneLists_ID <- function(df, terms){
  cat2gene <- unique(df[-1])
  
  out <- list()
  for(i in terms){
    tmp <- cat2gene$ID[which(cat2gene$PATH == i)]
    out[[i]] <- tmp
  }
  out
}

# Ids
GO_bp_terms <- getGeneLists_ID(c5bp, goterms_bp)
GO_cc_terms <- getGeneLists_ID(c5cc, goterms_cc)
GO_mf_terms <- getGeneLists_ID(c5mf, goterms_mf)

KEGG_terms <- getGeneLists_ID(kegg, keggterms)

# ENTREZIDs
GO_bp_ids <- getGeneLists_GO(c5bp, goterms_bp)
GO_cc_ids <- getGeneLists_GO(c5cc, goterms_cc)
GO_mf_ids <- getGeneLists_GO(c5mf, goterms_mf)

KEGG_ids <- getGeneLists_KEGG(kegg, keggterms)

# Attach 

listy[["C5bp"]]$IDs <- sapply(listy[["C5bp"]]$category, function(x) paste0(GO_bp_terms[[x]], collapse = ", "))
listy[["C5cc"]]$IDs <- sapply(listy[["C5cc"]]$category, function(x) paste0(GO_cc_terms[[x]], collapse = ", "))
listy[["C5mf"]]$IDs <- sapply(listy[["C5mf"]]$category, function(x) paste0(GO_mf_terms[[x]], collapse = ", "))

listy[["C2kegg"]]$IDs <- sapply(listy[["C2kegg"]]$category, function(x) paste0(KEGG_terms[[x]], collapse = ", "))

listy[["C2kegg"]]$ENTREZID <- sapply(listy[["C2kegg"]]$category, function(x) paste0(KEGG_ids[[x]], collapse = ", "))

listy[["C5bp"]]$ENTREZID <- sapply(listy[["C5bp"]]$category, function(x) paste0(GO_bp_ids[[x]], collapse = ", "))
listy[["C5cc"]]$ENTREZID <- sapply(listy[["C5cc"]]$category, function(x) paste0(GO_cc_ids[[x]], collapse = ", "))
listy[["C5mf"]]$ENTREZID <- sapply(listy[["C5mf"]]$category, function(x) paste0(GO_mf_ids[[x]], collapse = ", "))


# Write file

write.xlsx(listy, file = "brown_goids.xlsx", row.names = FALSE, col.names = TRUE)
