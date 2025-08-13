library(dplyr); library(ggplot2); library(ggforce)
setwd("/home/mukesh/binukumar/Project_ASD/asd_maf/maf_pop/gnomAD_download/")


#------------ Data Modification for Analysis ------------# 

# Read the data
gnomad <- read.table("./all_variant.txt", sep = "\t", header = TRUE)
cols_to_keep <- grepl("variant_id|_ac|_an", colnames(gnomad), ignore.case = TRUE) &
  !grepl("_hom|_hemi", colnames(gnomad), ignore.case = TRUE)

gnomad1 <- gnomad[,cols_to_keep]

# Modify column names
old <- c("variant_id",  "joint_", "ac", "an", "afr", "amr", "asj", "eas", "nfe", "fin", "mid", "sas", "ami", "remaining")
new <- c("VCF", "", "AC", "AN", "AFR", "AMR", "ASJ", "EAS", "NFE", "FIN", "MID", "SAS", "AMI", "OTH")

colnames(gnomad1) <- sapply(colnames(gnomad1), function(x) {
  for (i in seq_along(old)) x <- gsub(old[i], new[i], x, ignore.case = FALSE)
  x
})
gnomad1$VCF <- paste0("chr", gnomad1$VCF)

# select variants of interest
indigen <- read.table("./var_indigen.txt", sep = "\t", header = T)
filtered_data <- merge(x = gnomad1, y = indigen[, c("VCF", "IndiGen_AC", "IndiGen_AN")], by = "VCF", all.y = TRUE)

# prepare AN absent columns
add_subtracted_columns <- function(data) {
  start_indices <- seq(3, 25, by = 2)
  for (i in start_indices) {
    if (i <= ncol(data)) {
      new_col_name <- paste0(colnames(data)[i], "_ab")
      data[[new_col_name]] <- data[[i]] - data[[i - 1]]
    }
  }
  return(data)
}

# Apply the function to filtered_data
ft_data <- add_subtracted_columns(filtered_data)

# Sort columns based on column name & Reorder
cat("c(", paste(sprintf('"%s"', colnames(ft_data)), collapse = ", "), ")", sep = "")
cols <- c("VCF", "IndiGen_AC", "IndiGen_AN", "IndiGen_AN_ab",
          "AC", "AN", "AN_ab","AFR_AC", "AFR_AN", "AFR_AN_ab",
          "AMR_AC", "AMR_AN", "AMR_AN_ab", "ASJ_AC", "ASJ_AN", "ASJ_AN_ab",
          "EAS_AC", "EAS_AN", "EAS_AN_ab", "FIN_AC", "FIN_AN", "FIN_AN_ab",
          "MID_AC", "MID_AN", "MID_AN_ab", "NFE_AC", "NFE_AN", "NFE_AN_ab",
          "AMI_AC", "AMI_AN", "AMI_AN_ab", "SAS_AC", "SAS_AN", "SAS_AN_ab",
          "OTH_AC", "OTH_AN", "OTH_AN_ab")

# Get the indices and reordering of the desired column names in ft_data
indices <- match(cols, colnames(ft_data))
ft_data_reordered <- ft_data[, indices]
df <- ft_data_reordered[, c( -3, -6, -9, -12, -15, -18, -21, -24, -27, -30, -33, -36)]
#write.table(df, "./gnomad_ranbp2.txt", sep = "\t", row.names = F)


#------------ Statistics: Fisher-exact/Chi Square Test ------------# 

# df data to f variable keep original for trobuleshoot
f = df[, 2:25]
rownames(f) = df$VCF

# Define the column pairs for each group
column_pairs = list(IndiGen = c(1, 2, 1, 2), ALL = c(1, 2, 3, 4), AFR = c(1, 2, 5, 6), 
                    AMR = c(1, 2, 7, 8), ASJ = c(1, 2, 9, 10), EAS = c(1, 2, 11, 12), 
                    FIN = c(1, 2, 13, 14), MID = c(1, 2, 15, 16), NFE = c(1, 2, 17, 18), 
                    AMI = c(1, 2, 19, 20), SAS = c(1, 2, 21, 22), OTH = c(1, 2, 23, 24))

# Perform stat and store p_corrected or adjusted values
for (group in names(column_pairs)) {
  f[[paste0("gnomad_", group, "_bonf")]] <- sapply(1:nrow(f), function(i) {
    vals <- as.numeric(f[i, column_pairs[[group]]])
    if (length(vals) == 4 && all(!is.na(vals)) && all(vals >= 0)) {
      mat <- matrix(vals, nrow = 2)
      if (all(rowSums(mat) > 0) && all(colSums(mat) > 0)) {
        return(tryCatch(fisher.test(mat)$p.value, error = function(e) NA))
      }
    }
    NA
  }) |> p.adjust(method = "bonferroni")
}


#------------  Graph Preparation ------------#

# ---> Step-1. Arrange cols for MAF Preparation
column_pairs <- list(
  c("IndiGen_AC", "IndiGen_AN", "IndiGen_MAF"), c("AC", "AN", "ALL_MAF"), c("AFR_AC", "AFR_AN", "AFR_MAF"),  
  c("AMR_AC", "AMR_AN", "AMR_MAF"), c("ASJ_AC", "ASJ_AN", "ASJ_MAF"),  c("EAS_AC", "EAS_AN", "EAS_MAF"),  
  c("FIN_AC", "FIN_AN", "FIN_MAF"), c("MID_AC", "MID_AN", "MID_MAF"),  c("NFE_AC", "NFE_AN", "NFE_MAF"),  
  c("AMI_AC", "AMI_AN", "AMI_MAF"),  c("SAS_AC", "SAS_AN", "SAS_MAF"),  c("OTH_AC", "OTH_AN", "OTH_MAF"))
  
# Separate the first column
first_col <- filtered_data[, 1]
data_rest <- filtered_data[, -1]
# Apply NA replacement logic to columns 2–25
for (i in seq_along(data_rest)) {
  data_rest[[i]][is.na(data_rest[[i]])] <- if (i %% 2 == 1) 0 else max(data_rest[[i]], na.rm = TRUE)
}
# Combine the first column back
filtered_data_cleaned <- cbind(first_col, data_rest)
colnames(filtered_data_cleaned)[1] <- colnames(filtered_data)[1]  # retain original column name
# trim bigger legth for VCF
filtered_data_cleaned$VCF <- sapply(filtered_data_cleaned$VCF, function(x) {
  len <- nchar(x)
  if (len > 32) {
    trimmed <- substr(x, 1, 32)
    paste0(trimmed, "_", len - 32, "nts")
  } else {
    x
  }
})

# ---> Step-2. Calculate the new MAF columns
for (pair in column_pairs) {
  filtered_data_cleaned[[pair[3]]] <- filtered_data_cleaned[[pair[1]]] / filtered_data_cleaned[[pair[2]]]
}

# ---> Step-3. Extract relevant columns for long-data preparation
vcf_ids <- filtered_data_cleaned$VCF
maf_cols <- grep("_MAF$", names(filtered_data_cleaned), value = TRUE)
bonf_cols <- grep("gnomad_.*_bonf$", names(f), value = TRUE)

# Create long-format MAF data
maf_long <- do.call(rbind, lapply(maf_cols, function(col) {
  data.frame(VCF = vcf_ids, Population = col, MAF = filtered_data_cleaned[[col]])
}))

# Create long-format Bonferroni p-value data
bonf_long <- do.call(rbind, lapply(bonf_cols, function(col) {
  pop <- sub("^gnomad_", "", sub("_bonf$", "_MAF", col))  # Match MAF naming
  data.frame(VCF = vcf_ids, Population = pop, f_pvalue = f[[col]])
}))

# Merge both on VCF and Population
final_data <- merge(maf_long, bonf_long, by = c("VCF", "Population"))

# ---> Step-4. Add ACMG parameters to data 
master <- read.table("./var_indigen.txt", header = T, sep = "\t")
#trim bigger legth for VCF=VCF(filtered_data_cleaned)
master$VCF <- sapply(master$VCF, function(x) {
  len <- nchar(x)
  if (len > 32) {
    trimmed <- substr(x, 1, 32)
    paste0(trimmed, "_", len - 32, "nts")
  } else {
    x
  }
})

acmg <- arrange(merge(x=master[,c("VCF", "ACMG", "gene")], y= final_data, by = "VCF", all.y = T))

# Specify the order of the Population levels
level <- c("IndiGen_MAF", "ALL_MAF", "AFR_MAF","AMR_MAF","ASJ_MAF", "EAS_MAF","FIN_MAF","MID_MAF",
           "NFE_MAF", "AMI_MAF", "SAS_MAF","OTH_MAF")
acmg$Population <- factor(acmg$Population, levels = level)


# <---------------------- PLOT ----------------------> #

# --- Prepare gene-variant data ---
uv <- unique(data.frame(VCF=acmg$VCF, gene=acmg$gene))
uv_plot <- uv %>%
  distinct(VCF, gene) %>% mutate(VCF = factor(VCF, levels = unique(VCF)),
                                 global_x = as.numeric(VCF))
vcf_levels <- levels(uv_plot$VCF)

gene_blocks <- uv_plot %>% group_by(gene) %>%  summarise(
  xmin = min(global_x) - 0.5,  xmax = max(global_x) + 0.5,  xmid = mean(global_x),
  .groups = "drop")

# --- Main ACMG plot (top) ---
acmg$VCF <- factor(acmg$VCF, levels = vcf_levels)

main_plot <- ggplot(acmg, aes(x = VCF, y = Population, colour = ACMG, size = MAF)) + 
  geom_point() + 
  geom_point(data = subset(acmg, f_pvalue < 0.05), shape = 1, colour = "#D36D68", size = 6) +
  theme_minimal() +
  theme(panel.grid = element_line(),  axis.line = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90,  hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

# --- Bottom gene/variant axis plot ---
bottom_axis_plot <- ggplot() +
  geom_rect(data = gene_blocks, aes(xmin = xmin, xmax = xmax, ymin = 0.7, ymax = 0.77, fill = gene), color = "black") +
  geom_text(data = gene_blocks, aes(x = xmid, y = 0.735, label = gene), angle = 90, size = 2.8, color = "white") +
  scale_x_continuous(limits = c(0.5, length(vcf_levels) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.68, 0.79), expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(3, 10, 25, 10)) +
  coord_cartesian(clip = "off")

# --- Combine both plots tightly ---

png("./Figure7.png", height = 2800, width = 4000,  units = "px", res = 300, pointsize = 12) 
library(patchwork)

bottom_axis_plot / main_plot + 
  plot_layout(heights = c(0.5, 3)) +
  plot_annotation(
    title = "Figure 7",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0),  # hjust = 0 aligns to the left
      plot.margin = margin(10, 10, 10, 10)
    )
  )


dev.off()
