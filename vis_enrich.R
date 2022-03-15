
# LOAD LIBRARIES ---------------------------

lib <-  c("readxl", "dplyr", "ggplot2", "stringr", "wesanderson")

invisible(lapply(lib, library, character.only = TRUE))

# LOAD DATA ---------------------------

# Iterate for each module/color
rm(list = ls())

module_bp <- read_excel("pink_goids.xlsx", sheet = "C5bp")

module_cc <- read_excel("pink_goids.xlsx", sheet = "C5cc")

module_mf <- read_excel("pink_goids.xlsx", sheet = "C5mf")

# KEGG and Reactome
kegg <- read_excel("pink_goids.xlsx", sheet = "C2kegg") %>%
  select(-IDs,
         -ENTREZID)

reactome <- read_excel("pink_goids.xlsx", sheet = "C2reac")


# Name

module_bp$cat <- c("BP")
module_cc$cat <- c("CC")
module_mf$cat <- c("MF")

plot_enrich <- rbind(module_bp, module_cc, module_mf, deparse.level = 1)

plot_enrich <- plot_enrich[which(plot_enrich$over_represented_pvalue_adj < 0.05),]

plot_enrich$cat <- as.factor(plot_enrich$cat)

plot_enrich <- plot_enrich %>%
  select(cat,
         IDs,
         category,
         ENTREZID,
         over_represented_pvalue_adj,
         numDEInCat) %>%
  rename(Category = "cat",
         ID = "IDs",
         Term = "category",
         Genes = "ENTREZID",
         adj_pval = "over_represented_pvalue_adj",
         Count = "numDEInCat")


# Clean GO terms

plot_enrich$Term <- gsub("GO_", "", plot_enrich$Term)

plot_enrich$Term <- gsub("_", " ", plot_enrich$Term)

plot_enrich$Term <- str_to_sentence(plot_enrich$Term)


# Choose palette

pal <- wes_palette("Zissou1", 5, type = "continuous")
pal2 <- c("#d2fbd4","#a5dbc2","#7bbcb0","#559c9e","#3a7c89","#235d72","#123f5a")


# Name
kegg$cat <- c("Kegg")
reactome$cat <- c("Reactome")


plot_path <- rbind(kegg, reactome, deparse.level = 1)

plot_path <- plot_path[which(plot_path$over_represented_pvalue_adj < 0.05),]

plot_path$cat <- as.factor(plot_path$cat)

plot_path <- plot_path %>%
  select(cat,
         category,
         over_represented_pvalue_adj,
         numDEInCat) %>%
  rename(Category = "cat",
         Term = "category",
         adj_pval = "over_represented_pvalue_adj",
         Count = "numDEInCat")

# Clean terms
plot_path$Term <- gsub("KEGG_|REACTOME_", "", plot_path$Term)

plot_path$Term <- gsub("_", " ", plot_path$Term)

plot_path$Term <- str_to_sentence(plot_path$Term)



# Joint visualization
plot_enrich <- plot_enrich[,which(colnames(plot_enrich) %in% colnames(plot_path))]

plot_vis <- rbind(plot_enrich, plot_path)

ggplot(plot_vis, aes(Count, Term, fill = adj_pval)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  facet_grid(Category ~ ., scales = "free", space = "free") +
  scale_fill_gradientn(colours = rev(pal2), name = "Adjusted p-value") +
  xlab("Gene count") +
  ggtitle("MEPink") +
  theme(plot.title = element_text(hjust = 0.5))

