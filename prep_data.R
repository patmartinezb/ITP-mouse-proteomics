# PREP DATA ---------------------------

# Load libraries ---------------------------
lib <-  c("readr", "janitor", "dplyr")

invisible(lapply(lib, library, character.only = TRUE))

# Load data ---------------------------
df_mbr_v1617 <- read_delim("proteinGroups_mbr_v1.6.17.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


# Explore ---------------------------
colnames(df_mbr_v1617)

# Select variables and clean data ---------------------------
df_mbr_v1617 <- clean_names(df_mbr_v1617)

df_mbr_v1617 <- select(df_mbr_v1617, 
                       protein_i_ds,
                       q_value,
                       gene_names,
                       starts_with("lfq"),
                       starts_with("i_baq_"),
                       only_identified_by_site, 
                       reverse,
                       potential_contaminant,
                       -i_baq_peptides)


# Data clean-up of contaminants, reverse proteins and proteins only identified by site
df_mbr_v1617 <- filter(df_mbr_v1617, 
                       is.na(potential_contaminant),
                       is.na(reverse),
                       is.na(only_identified_by_site))


# Choose first protein id for further identification and annotation
df_mbr_v1617$protein_i_ds <- sub(";.*", "", df_mbr_v1617$protein_i_ds)


# Remove columns with contaminant info
df_mbr_v1617 <- select(df_mbr_v1617,
                       -only_identified_by_site,
                       -reverse,
                       -potential_contaminant)



# Save data for further analysis ---------------------------
saveRDS(df_mbr_v1617, file = "df_mbr_v1617.rds")