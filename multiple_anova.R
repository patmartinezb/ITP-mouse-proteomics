
# LOAD DATA ---------------------------

rm(list = ls())
df_init <- readRDS("df_mbr_v1617.rds")

# LOAD LIBRARIES ---------------------------

lib <-  c("dplyr", "tidyr", "preprocessCore", "factoextra", "ggplot2", "stringr",
          "PerformanceAnalytics", "limma", "ggsci", "openxlsx", "imputeLCMD", 
          "PCAtools", "RColorBrewer", "pheatmap")

invisible(lapply(lib, library, character.only = TRUE))


# DATA PRE-PROCESSING ---------------------------

# See performance and then remove q_value variable
summary(as.numeric(df_init$q_value))

df_init <- select(df_init,
                  -q_value)

# Check if any gene name is duplicated (make sure to count out the NAs)
df_init$gene_names %>% 
  na.omit() %>% 
  duplicated() %>% 
  any()

table(duplicated(na.omit(df_init$gene_names)))

na.omit(df_init) %>% 
  group_by(gene_names) %>% 
  summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% 
  filter(frequency > 1)

# Transform missing values into NAs
df_init[df_init == 0] <- NA

sapply(df_init, summary)

# Log2 transform
df_init[,c(3:ncol(df_init))] <- log2(df_init[,c(3:ncol(df_init))])
sapply(df_init, summary)


# MISSING VALUES ---------------------------

# Count
missing.values <- df_init %>%
  gather(key = "key", value = "val") %>%
  mutate(is.missing = is.na(val)) %>%
  group_by(key, is.missing) %>%
  summarise(num.missing = n()) %>%
  filter(is.missing == T) %>%
  select(-is.missing) %>%
  arrange(desc(num.missing)) 

missing.values


# Boxplots for higher number of NAs in treated
plotmis <- missing.values %>% 
  filter(str_detect(key, "i_baq"))

ivig <- grep("cc", plotmis$key, value = TRUE, perl = TRUE)
tto <- grep("tto", plotmis$key, value = TRUE, perl = TRUE)
cc_p <- grep("([1-3])._*", plotmis$key, value = TRUE, perl = TRUE)
d7 <- grep("([4-6]).*_7$", plotmis$key, value = TRUE, perl = TRUE)
d3 <- grep("([4-6]).*_3$", plotmis$key, value = TRUE, perl = TRUE)

plotmis <- plotmis[-c(1:3),]

# png("supp_fig_1A2.png", 2000, 2000, res=300)
plotmis %>% 
  mutate(group = case_when(key %in% ivig ~ "IVIg",
                           key %in% tto ~ "Active ITP",
                           key %in% cc_p ~ "Control",
                           key %in% d7 ~ "Day 7 - Passive ITP",
                           key %in% d3 ~ "Day 3 - Passive ITP")) %>% 
  ggplot(aes(group, num.missing, fill = group)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_locuszoom(alpha = 0.6) +
  theme(legend.position = "none") +
  xlab("Group") +
  ylab("Number of missing values")

# dev.off()


# Plot bars with percentages
missing.values <- df_init %>%
  gather(key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, isna) %>%
  summarise(num.isna = n()) %>%
  mutate(pct = num.isna / total * 100)

levels <- (missing.values %>% filter(isna == T) %>% filter(str_detect(key, "i_baq")) %>% arrange(desc(pct)))$key

percentage.plot.NA <- missing.values %>%
  filter(str_detect(key, "i_baq")) %>%
  ggplot() +
  geom_bar(aes(x = reorder(key, dplyr::desc(pct)), 
               y = pct, fill=isna), 
           stat = 'identity', alpha=0.8) +
  scale_x_discrete(labels = c("i_baq_tto_pti_4" = "Active ITP 4",
                              "i_baq_tto_pti_1" = "Active ITP 1",
                              "i_baq_1_dag_3" = "Day 3 - Control 1",
                              "i_baq_tto_pti_5" = "Active ITP 5",
                              "i_baq_6_dag_3" = "Day 3 - Passive ITP 6",
                              "i_baq_4_dag_3" = "Day 3 - Passive ITP 4",
                              "i_baq_5_dag_3" = "Day 3 - Passive ITP 5",
                              "i_baq_tto_pti_2" = "Active ITP 2",
                              "i_baq_cc_pti_10" = "IVIg 10",
                              "i_baq_6_dag_7" = "Day 7 - Passive ITP 6",
                              "i_baq_5_dag_7" = "Day 7 - Passive ITP 5",
                              "i_baq_tto_pti_3" = "Active ITP 3",
                              "i_baq_4_dag_7" = "Day 7 - Passive ITP 4",
                              "i_baq_3_dag_3" = "Day 3 - Control 3",
                              "i_baq_2_dag_7" = "Day 7 - Control 2",
                              "i_baq_cc_pti_9" = "IVIg 9",
                              "i_baq_cc_pti_6" = "IVIg 6",
                              "i_baq_3_dag_7" = "Day 7 - Control 3",
                              "i_baq_cc_pti_7" = "IVIg 7",
                              "i_baq_1_dag_7" = "Day 7 - Control 1",
                              "i_baq_2_dag_3" = "Day 3 - Control 2",
                              "i_baq_cc_pti_8" = "IVIg 8"), limits = levels) +
  scale_fill_manual(name = "", 
                    values = c('gray43', 'mediumpurple1'), labels = c("Present", "Missing")) +
  coord_flip() +
  labs(x = 'Samples', y = "% of missing values") +
  theme_minimal() +
  geom_hline(yintercept=35, linetype="dashed", color = "black")

# png("supp_fig_1A.png", 3500, 2000, res=300)

percentage.plot.NA

# dev.off()

# Select only iBAQs
dfNA <- select(df_init,
               protein_i_ds,
               gene_names,
               starts_with("i_baq"))

# Remove i_baq_ from the names
colnames(dfNA) <- gsub("i_baq_", "", colnames(dfNA))


# Make groups out of intensity_names
grouping <- function(x, pattern){
  for (i in 1:ncol(x)){
    new_vector <- grep(pattern, names(x), value = TRUE, perl = TRUE)
  }
  return(new_vector)
}

cc <- grouping(dfNA, "cc")
tto <- grouping(dfNA, "tto")
cc_p <- grouping(dfNA, "(^[1-3])._*")
d7 <- grouping(dfNA, "(^[4-6]).*_7$")
d3 <- grouping(dfNA, "(^[4-6]).*_3$")
intensity_names <- c(cc_p, d3, d7, tto, cc)

dfNA <- select(dfNA,
               protein_i_ds,
               gene_names,
               intensity_names)

# Multi-scatter plot
# chart.Correlation(dfNA[,-c(1,2)], histogram = TRUE, pch = 19)

# Histograms
# multi.hist(dfNA[cc], density = FALSE)
# multi.hist(dfNA[tto], density = FALSE)
# multi.hist(dfNA[cc_p], density = FALSE)
# multi.hist(dfNA[d3], density = FALSE)
# multi.hist(dfNA[d7], density = FALSE)



# Remove samples with too many NAs
dfNA <- dfNA[,-c(3, 15, 18)]

cc <- grouping(dfNA, "cc")
tto <- grouping(dfNA, "tto")
cc_p <- grouping(dfNA, "(^[1-3])._*")
d7 <- grouping(dfNA, "(^[4-6]).*_7$")
d3 <- grouping(dfNA, "(^[4-6]).*_3$")
intensity_names <- c(cc_p, d3, d7, tto, cc)

dfNA <- select(dfNA,
               protein_i_ds,
               gene_names,
               intensity_names)


# FILTERING ---------------------------
# Filter out those rows with all NAs
filt.0 <- apply(dfNA[intensity_names], 1, function(x){
  sum(is.na(x)) == length(x) | sum(is.na(x))== length(x) - 1
})
dfNA_filt <- dfNA[!filt.0,]

# Keep those cases with all valid values in at least one group
filt3 <- apply(dfNA_filt, 1, function(x){
  sum(is.na(x[cc])) == 0 | sum(is.na(x[tto])) == 0 | sum(is.na(x[cc_p])) == 0 | sum(is.na(x[d3])) == 0 | sum(is.na(x[d7])) == 0
})
dfNA_filt_3 <- dfNA_filt[filt3,] # This the universe


# FILTERING ---------------------------

## Pre-normalization
boxplot(dfNA_filt_3[intensity_names], outline = FALSE)

boxplot(list(passive_model = as.numeric(unlist(dfNA_filt_3[,3:13])), active_model = as.numeric(unlist(dfNA_filt_3[,14:21]))), outline = FALSE, main = "Pre-normalization") # pre-normalization

# png("supp_fig_1B_up.png", 3500, 2000, res=300)

dfNA_filt_3[,2:21] %>% 
  pivot_longer(!gene_names, names_to = "group", values_to = "counts") %>%
  mutate(Model = ifelse((group %in% cc) | (group %in% tto), "Active", "Passive")) %>% 
  ggplot(aes(x=reorder(group,counts,na.rm = TRUE), y=counts, fill = Model)) +
  geom_boxplot() +
  labs(y="Intensity", x="Samples") +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, vjust=0.5, size=12)) +
  scale_x_discrete(labels = c("cc_pti_10" = "IVIg 10",
                              "tto_pti_5" = "Active ITP 5",
                              "cc_pti_6" = "IVIg 6",
                              "tto_pti_2" = "Active ITP 2",
                              "6_dag_7" = "Day 7 - Passive ITP 6",
                              "6_dag_3" = "Day 3 - Passive ITP 6",
                              "4_dag_3" = "Day 3 - Passive ITP 4",
                              "5_dag_3" = "Day 3 - Passive ITP 5",
                              "tto_pti_2" = "Active ITP 2",
                              "cc_pti_10" = "IVIg 10",
                              "6_dag_7" = "Day 7 - Passive ITP 6",
                              "5_dag_7" = "Day 7 - Passive ITP 5",
                              "tto_pti_3" = "Active ITP 3",
                              "4_dag_7" = "Day 7 - Passive ITP 4",
                              "3_dag_3" = "Day 3 - Control 3",
                              "2_dag_7" = "Day 7 - Control 2",
                              "cc_pti_9" = "IVIg 9",
                              "cc_pti_6" = "IVIg 6",
                              "3_dag_7" = "Day 7 - Control 3",
                              "cc_pti_7" = "IVIg 7",
                              "1_dag_7" = "Day 7 - Control 1",
                              "2_dag_3" = "Day 3 - Control 2",
                              "cc_pti_8" = "IVIg 8"))

# dev.off()


# NORMALIZATION ---------------------------

# Quantile regression
set.seed(42)
df_norm <- normalize.quantiles(as.matrix(dfNA_filt_3[intensity_names]))
colnames(df_norm) <- intensity_names

boxplot(list(passive_model = as.numeric(unlist(df_norm[,1:11])), active_model = as.numeric(unlist(df_norm[,12:19]))), outline = FALSE, main = "Quantile normalized")

dfn <- as.data.frame(df_norm)
dfn$ids <- dfNA_filt_3$gene_names
dfn %>% 
  pivot_longer(!ids, names_to = "group", values_to = "counts") %>%
  mutate(Model = ifelse((group %in% cc) | (group %in% tto), "Active", "Passive")) %>% 
  ggplot(aes(x=reorder(group,counts,na.rm = TRUE), y=counts, fill = Model)) +
  geom_boxplot() +
  labs(y="Intensity", x="Samples", 
       subtitle="Quantile normalization (no NAs)") +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, vjust=0.5, size=12))


dfNA.norm.uni <- as.data.frame(df_norm)
rownames(dfNA.norm.uni) <- dfNA_filt_3$protein_i_ds


# IMPUTATION ---------------------------
#Imputation of left-censored data using QRILC
set.seed(42)

dfNAimp <- impute.QRILC(dfNA.norm.uni[intensity_names])
dfNAimp <- dfNAimp[[1]] 

dfNAimp$Protein.IDs <- rownames(dfNAimp)
dfNAimp$gene_names <- dfNA_filt_3$gene_names


dev.off()


# After imputation
boxplot(list(passive_model = as.numeric(unlist(dfNAimp[,1:11])), active_model = as.numeric(unlist(dfNAimp[,12:19]))), outline = FALSE, main = "Quantile normalized and imputed")

# png("supp_fig_1B_down.png", 3500, 2000, res=300)

dfNAimp %>%
  select(-Protein.IDs) %>% 
  pivot_longer(!gene_names, names_to = "group", values_to = "counts") %>%
  mutate(Model = ifelse((group %in% cc) | (group %in% tto), "Active", "Passive")) %>% 
  ggplot(aes(x=reorder(group,counts,na.rm = TRUE), y=counts, fill = Model)) +
  geom_boxplot() +
  labs(y="Intensity", x="Samples") +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, vjust=0.5, size=12))  +
  scale_x_discrete(labels = c("cc_pti_10" = "IVIg 10",
                              "tto_pti_5" = "Active ITP 5",
                              "cc_pti_6" = "IVIg 6",
                              "tto_pti_2" = "Active ITP 2",
                              "6_dag_7" = "Day 7 - Passive ITP 6",
                              "6_dag_3" = "Day 3 - Passive ITP 6",
                              "4_dag_3" = "Day 3 - Passive ITP 4",
                              "5_dag_3" = "Day 3 - Passive ITP 5",
                              "tto_pti_2" = "Active ITP 2",
                              "6_dag_7" = "Day 7 - Passive ITP 6",
                              "5_dag_7" = "Day 7 - Passive ITP 5",
                              "tto_pti_3" = "Active ITP 3",
                              "4_dag_7" = "Day 7 - Passive ITP 4",
                              "3_dag_3" = "Day 3 - Control 3",
                              "2_dag_7" = "Day 7 - Control 2",
                              "cc_pti_9" = "IVIg 9",
                              "3_dag_7" = "Day 7 - Control 3",
                              "cc_pti_7" = "IVIg 7",
                              "1_dag_7" = "Day 7 - Control 1",
                              "2_dag_3" = "Day 3 - Control 2",
                              "cc_pti_8" = "IVIg 8"))

# dev.off()


# PLOTS PRE-DE ---------------------------
p <- pca(dfNAimp[,-c(20:21)], removeVar = 0.1)

screeplot(p,
          components = getComponents(p, 1:6),
          hline = 80, vline = 14) + 
  geom_text(aes(5, 80, label = '80% explained variation', vjust = -1))

# Bi-plot
biplot(p, 
       hline = 0, vline = c(-25, 0, 25),
       gridlines.major = FALSE, gridlines.minor = FALSE,
       drawConnectors = FALSE, lab = intensity_names,
       legendPosition = 'right')



# Other PCA
tt <- as.data.frame(t(dfNAimp[,-c(20:21)]))
tt <- tt %>% 
  mutate(group = c(c("Control day 7 - Passive ITP"),
                   c("Control day 3 - Passive ITP"),
                   c("Control day 7 - Passive ITP"),
                   c("Control day 3 - Passive ITP"),
                   c("Control day 7 - Passive ITP"),
                   rep("Day 3 - Passive ITP", times=3),
                   rep("Day 7 - Passive ITP", times=3),
                   rep("Active ITP", times=3),
                   rep("IVIg", times=5)),
         model = c(rep("Passive", times = 11),
                   rep("Active", times = 8))) %>% 
  select(group, everything())

rownames(tt) <- colnames(dfNAimp[,c(1:19)])

res.pca <- prcomp(select_if(tt, is.numeric),  scale = FALSE)

# png("fig_1C_pre.png", 3500, 2700, res=300)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (but not "text")
             habillage = tt$model, # color by groups
             legend.title = "", title="",
             invisible="quali", 
             pointsize = 4,
             pointshape = 19) +
  theme_classic()+
  # scale_color_brewer(palette="Dark2")+
  theme(legend.title = element_blank())+
  scale_color_npg()

# dev.off()


# Heatmap
data_sig2 <- dfNAimp[intensity_names]

#tiff("PCA.tiff", width = 14, height = 10, units = 'in', res = 300)
pheatmap(data_sig2, fontsize_row = 4, clustering_distance_cols = "correlation")
#graphics.off()


# REMOVE BATCH EFFECT ---------------------------

pal <- brewer.pal(5, "Dark2")

quan2 <- dfNAimp[,-c(20:21)]
pheno <- data.frame(group = c("c7","c3","c7","c3","c7",
                              "t3","t3","t3","t7","t7","t7",
                              "pt","pt","pt",
                              "pc","pc","pc","pc","pc"),
model <- c(rep("model1",11),rep("model2",8)),
row.names  <- colnames(quan2), stringsAsFactors = F)

set.seed(42)
quan2.b <- removeBatchEffect(quan2, batch = pheno$model)

layout(matrix(1:4,2,2))

plotMDS(quan2, gene.selection = "common", top = nrow(quan2), col = pal[factor(pheno$group)], main = "With batch - by group")
plotMDS(quan2.b, gene.selection = "common", top = nrow(quan2.b), col = pal[factor(pheno$group)], main = "Without batch - by group")

plotMDS(quan2, gene.selection = "common", top = nrow(quan2), col = pal[factor(pheno$model)], main = "With batch - by model")
plotMDS(quan2.b, gene.selection = "common", top = nrow(quan2.b), col = pal[factor(pheno$model)], main = "Without batch - by model")


tt <- as.data.frame(t(quan2.b))
tt <- tt %>% 
  mutate(group = c(c("Control day 7 - Passive ITP"),
                   c("Control day 3 - Passive ITP"),
                   c("Control day 7 - Passive ITP"),
                   c("Control day 3 - Passive ITP"),
                   c("Control day 7 - Passive ITP"),
                   rep("Day 3 - Passive ITP", times=3),
                   rep("Day 7 - Passive ITP", times=3),
                   rep("Active ITP", times=3),
                   rep("IVIg", times=5)),
         model = c(rep("Passive", times = 11),
                   rep("Active", times = 8))) %>% 
  select(group, model, everything())

rownames(tt) <- colnames(quan2.b)

res.pca <- prcomp(select_if(tt, is.numeric),  scale = FALSE)


# png("fig_1B_all.png", 3500, 2700, res=300)
# png("fig_1C_post.png", 3500, 2700, res=300)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             habillage = tt$group, # color by groups
             legend.title = "Groups", title="",
             invisible="quali", 
             pointsize = 4,
             pointshape = 19) +
  theme_classic()+
  # scale_color_brewer(palette="Dark2")+
  theme(legend.title = element_blank())+
  scale_color_npg()

# dev.off()

dfNAimp <- as.data.frame(quan2.b)
rownames(dfNAimp) <- dfNA_filt_3$protein_i_ds


dfNAimp %>%
  mutate(Protein.IDs = rownames(dfNAimp)) %>% 
  pivot_longer(!Protein.IDs, names_to = "group", values_to = "counts") %>%
  mutate(Model = ifelse((group %in% cc) | (group %in% tto), "Active", "Passive")) %>% 
  ggplot(aes(x=reorder(group,counts,na.rm = TRUE), y=counts, fill = Model)) +
  geom_boxplot() +
  labs(y="Intensity", x="Samples") +
  scale_fill_npg() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, vjust=0.5, size=12))  +
  scale_x_discrete(labels = c("cc_pti_10" = "IVIg 10",
                              "tto_pti_5" = "Active ITP 5",
                              "cc_pti_6" = "IVIg 6",
                              "tto_pti_2" = "Active ITP 2",
                              "6_dag_7" = "Day 7 - Passive ITP 6",
                              "6_dag_3" = "Day 3 - Passive ITP 6",
                              "4_dag_3" = "Day 3 - Passive ITP 4",
                              "5_dag_3" = "Day 3 - Passive ITP 5",
                              "tto_pti_2" = "Active ITP 2",
                              "6_dag_7" = "Day 7 - Passive ITP 6",
                              "5_dag_7" = "Day 7 - Passive ITP 5",
                              "tto_pti_3" = "Active ITP 3",
                              "4_dag_7" = "Day 7 - Passive ITP 4",
                              "3_dag_3" = "Day 3 - Control 3",
                              "2_dag_7" = "Day 7 - Control 2",
                              "cc_pti_9" = "IVIg 9",
                              "3_dag_7" = "Day 7 - Control 3",
                              "cc_pti_7" = "IVIg 7",
                              "1_dag_7" = "Day 7 - Control 1",
                              "2_dag_3" = "Day 3 - Control 2",
                              "cc_pti_8" = "IVIg 8"))




# universe <- dfNAimp %>% 
#   mutate(Protein.IDs = rownames(dfNAimp)) %>% 
#   dplyr::select(Protein.IDs, everything())
# 
# write.xlsx(universe, file = "universe.xlsx", row.names=FALSE, col.names=TRUE)


# DIFFERENTIAL EXPRESSION ANALYSIS LIMMA ---------------------------

anno <- rownames(dfNAimp)

# Model
mod <- model.matrix(~0+as.factor(group), data = pheno)
table(rownames(mod) == colnames(dfNAimp))
colnames(mod) <- c("c3","c7","pc","pt","t3","t7")


# We use limma::lmFit to fit the model
# Fit the model
fit <- lmFit(dfNAimp, mod)

# Comparisons
contrast.matrix <- makeContrasts(Control_passive_D3_passive = t3 - ((c3+c7)/2),
                                 Control_passive_D7_passive = t7 - ((c3+c7)/2),
                                 Control_passive_Exp_active = pt - ((c3+c7)/2),
                                 Control_passive_Control_active = pc - ((c3+c7)/2),
                                 D3_passive_D7_passive = t3 - t7,
                                 D3_passive_Exp_active = t3 - pt,
                                 D3_passive_Control_active = t3 - pc,
                                 D7_passive_Exp_active = t7 - pt,
                                 D7_passive_Control_active = t7 - pc,
                                 Exp_active_Control_active = pt - pc,
                                 levels = mod)


contrast.matrix

# Fit the comparisons
fit2 <- contrasts.fit(fit, contrast.matrix)

head(fit2$coefficients) # Now our model coefficient is the DIFFERENCE IN MEAN control - exp

# Test these coefficients being != 0 with limma::eBayes
fit2 <- eBayes(fit2)

# Check the results for coefficients != 0 
summary(decideTests(fit2))


# all proteins DE at least in one comparison
all <- topTable(fit2, number = Inf, coef = 1:10, genelist = as.data.frame(anno), p.value = 0.05)

df_all <- dfNAimp[which(rownames(dfNAimp) %in% all$anno),]

# write.xlsx(df_all, "sig_prot_for_wgcna_NAremove_nobatch.xlsx", col.names = TRUE, row.names=TRUE)

# Heatmap
pheatmap(df_all)


# PCA
tt3 <- as.data.frame(t(df_all))
tt3 <- tt3 %>% 
  mutate(group = c(c("Control day 7 - Passive ITP"),
                   c("Control day 3 - Passive ITP"),
                   c("Control day 7 - Passive ITP"),
                   c("Control day 3 - Passive ITP"),
                   c("Control day 7 - Passive ITP"),
                   rep("Day 3 - Passive ITP", times=3),
                   rep("Day 7 - Passive ITP", times=3),
                   rep("Active ITP", times=3),
                   rep("IVIg", times=5))) %>% 
  select(group, everything())

res.pca3 <- prcomp(select_if(tt3, is.numeric),  scale = FALSE)


# png("fig_1B_DE.png", 3500, 2700, res=300)
fviz_pca_ind(res.pca3,
             geom.ind = "point", # show points only (nbut not "text")
             habillage = tt3$group, # color by groups
             legend.title = "Groups", title="",
             invisible="quali", 
             pointsize = 4,
             pointshape = 19) +
  theme_classic()+
  # scale_color_brewer(palette="Dark2")+
  theme(legend.title = element_blank())+
  scale_color_npg()

# dev.off()

# SAVE DATA FOR WGCNA ---------------------------

saveRDS(df_all, file = "df_wgcna.rds")

