# Load required libraries
library(ggplot2)
library(cluster)
library(dplyr); library(tidyr)
library(ggrepel); library(treemap); library(ggpubr); library(plotly)


setwd("/home/mukesh/binukumar/Project_ASD/plots/")
missense <- read.table("./revel0.75.tsv", sep = "\t",  header = T)
mis <- select(missense, c(6, 7, 8, 9, 10, 13, 11,  20, 21, 30,64, 112, 174, 183, 187, 197, 215, 289))

loftee <- read.table("./loftee.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", comment.char = "")
lof <- loftee[, c("Variant", "HET_AC", "HOM_AC", "AC_total", "Func.refGeneWithVer", "ExonicFunc.refGeneWithVer", 
                  "Gene.refGeneWithVer", "avsnp151", "SIFT_score", "Polyphen2_HVAR_score", "REVEL_score",
                  "CADD_phred", "CLNSIG", "IndiGen_AF", "X1000g2015aug_all", "gnomad41_exome_AF", "gnomad41_genome_AF", "ACMG")]

# Figure 2: Pie Chart Overview of LoF+Missense
td <- rbind(mis, lof)
td[td[[6]] == ".", 6] <- td[td[[6]] == ".", 5]
# pie chart
td1 <- as.data.frame(table(td$ExonicFunc.refGeneWithVer))
colnames(td1) <- c("Variant_class", "Count")

# Prepare data
td1 <- td1 %>%
  arrange(desc(Count)) %>%       # Optional: Largest slices first
  mutate(
    Percent = round(Count / sum(Count) * 100, 1),
    Label = paste0(Variant_class, "\n", Count, " (", Percent, "%)"),
    ypos = cumsum(Count) - 0.5 * Count
  )

# Pie chart with ggrepel
# Prepare legend labels with counts and percentages
counts <- td1$Count
percentages <- round(100 * counts / sum(counts), 1)
labels <- paste0(td1$Variant_class, " (", counts, ")")


png("./Figure2.png", pointsize = 10, res = 350, width=3100, height=2500)
# Set margins: less space around plot, more on right for legend
par(mar = c(1, 1, 3, 16))  # bottom, left, top, right

# Plot pie chart
par(cex.main = 1.3) 
pie(td1$Count,
    labels = paste0(round(100 * td1$Count / sum(td1$Count), 1), "%"),
    col = c("#e31a1c", "#1f78b4", "#ff7f00", "#33a02c", "#6a3d9a", "#b15928"), 
    border = "black",
    cex = 1.3, clockwise = T, init.angle = 150)

mtext(expression(bold("Figure 2")), side = 3, adj = 0, line = 0.5, cex = 1.3)

# Add legend outside plot area
legend("right",
       inset = c(-0.42, 0),  # move legend outside plot
       legend = labels,
       fill = c("#e31a1c", "#1f78b4", "#ff7f00", "#33a02c", "#6a3d9a", "#b15928"),
       title = "Variant Class",
       bty = "n",
       xpd = TRUE,  # allow drawing outside plot region
       cex = 1.5)

dev.off()




# Figure 3: Create heatmap
plot_data1 <- mis %>%
  group_by(Gene = .[[7]], ACMG = .[[18]]) %>%
  summarise(Variant_count = n(), .groups = 'drop')

m1 <- plot_data1[1:70,]
m1c <- m1 %>% tidyr::complete(Gene, ACMG, fill = list(Variant_count = 0))

m2 <- plot_data1[71:141,]
m2c <-  m2 %>% tidyr::complete(Gene, ACMG, fill = list(Variant_count = 0))
  
p1 <- ggplot(m1c, aes(x = Gene, y = ACMG, fill = Variant_count)) +
  geom_tile(color = "grey", size = 0.1, width = 0.98, height = 0.99) +
  scale_fill_gradient(low = "#f3eded", high = "red", na.value = "#f3eded") +
  geom_text(data = dplyr::filter(m1c, Variant_count != 0),
            aes(label = Variant_count), color = "black", size = 3) +   # Only non-zero labels
  labs(title = "A) Set-1",
       x = "Gene (n=64)", y = "ACMG Classification") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.9, size = 5.5),
        plot.title = element_text(face = "bold"))

p2 <- ggplot(m2c, aes(x = Gene, y = ACMG, fill = Variant_count)) +
  geom_tile(color = "grey", size = 0.1, width = 0.98, height = 0.99) +
  scale_fill_gradient(low = "#f3eded", high = "red", na.value = "#f3eded") +
  geom_text(data = dplyr::filter(m2c, Variant_count != 0),
            aes(label = Variant_count), color = "black", size = 3) +   # Only non-zero labels
  labs(title = "B) Set-2",
       x = "Gene (n=66)", y = "ACMG Classification") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.9, size = 5.5),
        plot.title = element_text(face = "bold"))

# Save Panel
png("./Figure3.png", pointsize = 8, res = 500, width=5000, height=2500)
annotate_figure(
  ggarrange(p1, p2, ncol = 1, nrow = 2, label.x = 0.18, label.y = 0.95),
  top = text_grob("Figure 3", face = "bold", size = 14, hjust = 0, x = 0 ))

dev.off()
  


# Figure 4:  Stackbar for LoF: ACMG==LP/P
sb <- dplyr::filter(lof, ACMG=="Pathogenic")
gn <- unique(sb$Gene.refGeneWithVer)
lof_ft <- lof[lof$Gene.refGeneWithVer %in% gn, ]

# Step 1: Group and count
lof_bar_data <- lof_ft %>%
  count(Gene.refGeneWithVer, ACMG, name = "Count")
# Step 2: Plot stacked bar chart
png("./Figure4.png", pointsize = 8, res = 500, width=4500, height=5000)

ggplot(lof_bar_data, aes(x = Count, y = Gene.refGeneWithVer, fill = ACMG)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("Pathogenic" = "#71221d","VUS" = "grey")) +
  scale_x_continuous(labels = number_format(accuracy = 1)) +                # decimal remove
  labs(title = "Figure 4", x = "Variant Count", y = "Gene")+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1), 
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 14))

dev.off()


# Figure 5: heatmap

# ---> step 1
library(tidyverse)
# Read your data
acmg_atr <- read.table("acmg_atrributes.txt", sep = "\t", header = TRUE)
# Get only the ACMG attribute columns (exclude "VCF" and "ACMG")
long_acmg <- acmg_atr %>%
  pivot_longer(
    cols = -c(VCF, ACMG),               # pivot everything except VCF and ACMG
    names_to = "ACMG_attribute",
    values_to = "value"
  )
# ---> step 2
#png("./Figure4.png", pointsize = 8, res = 500, width=4500, height=5000)
long_acmg$value <- as.character(long_acmg$value)
long_acmg$VCF <- sapply(long_acmg$VCF, function(x) {
  len <- nchar(x)
  if (len > 32) {
    trimmed <- substr(x, 1, 32)
    paste0(trimmed, "_", len - 32, "nts")
  } else {
    x
  }
})

# save
png("./Figure5.png", pointsize = 8, res = 500, width=6000, height=4000)
ggplot(long_acmg, aes(x = VCF, y = ACMG_attribute, fill = ACMG)) +
  geom_tile(color = "black") +
  geom_text(aes(label = value), size = 3, color = "black", na.rm = TRUE) +
  scale_fill_manual(values = c(P = "#008080", LP = "#F4A384"), na.value = "gray") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Figure 5", x = "VCF", y = "ACMG-Attribute", fill = "ACMG") 

dev.off()

# Figure 6: Genotype count data
a <- read.table("./f6_var.txt", header = T, sep = "\t")
a$variant <- sapply(a$variant, function(x) {
  len <- nchar(x)
  if (len > 32) {
    trimmed <- substr(x, 1, 32)
    paste0(trimmed, "_", len - 32, "nts")
  } else {
    x
  }
})

b <- read.table("./f6_gene_hethom.txt", header = T, sep = "\t")
c <- read.table("./f6_gene_uhet_hom.txt", header = T, sep = "\t")

# a
a1 <- ggplot(a, aes(x = zygosity, y = variant, fill = Individuals_count)) +
  geom_tile(color = "grey", size = 0.1, width = 0.98, height = 0.99) +
  scale_fill_gradient(low = "#f3eded", high = "red", na.value = "#f3eded") +
  geom_text(data = dplyr::filter(a, Individuals_count != 0),
            aes(label = Individuals_count), color = "black", size = 3) +   # Only non-zero labels
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.9, size = 5.5),
        plot.title = element_text(face = "bold"))

# b 
b1 <- ggplot(b, aes(x = zygosity, y = gene, fill = Individuals_count)) +
  geom_tile(color = "grey", size = 0.1, width = 0.98, height = 0.99) +
  scale_fill_gradient(low = "#f3eded", high = "red", na.value = "#f3eded") +
  geom_text(data = dplyr::filter(b, Individuals_count != 0),
            aes(label = Individuals_count), color = "black", size = 3) +   # Only non-zero labels
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.9, size = 5.5),
        plot.title = element_text(face = "bold"))



#c 
c1 <- ggplot(c, aes(x = zygosity, y = gene, fill = Individuals_count)) +
  geom_tile(color = "grey", size = 0.1, width = 0.98, height = 0.99) +
  scale_fill_gradient(low = "#f3eded", high = "red", na.value = "#f3eded") +
  geom_text(data = dplyr::filter(c, Individuals_count != 0),
            aes(label = Individuals_count), color = "black", size = 3) +   # Only non-zero labels
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.9, size = 5.5),
        plot.title = element_text(face = "bold"))


library(ggpubr)
library(grid)  # for unit()

# Arrange with subplot labels
combined <- ggarrange(a1, ggarrange(b1, c1, ncol = 2, nrow = 1, labels = c("B", "C"), font.label = list(size = 14)),
  ncol = 2, widths = c(2, 2), labels = c("A", ""),      # Label for a1, leave nested blank
  font.label = list(size = 14))

# Add "Figure 6" at upper left
png("./Figure6.png", pointsize = 18, res = 300, width = 6000, height = 4000)
annotate_figure(combined, 
                top = text_grob("Figure 6", x = 0, hjust = 0, face = "bold", size = 16) # left align
)

dev.off()
