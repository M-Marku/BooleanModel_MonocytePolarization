# Load the libraries

library(BoolNet)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(proxy)
library(plot.matrix)
library(clValid)
library(pheatmap)
library(digest)
library(cluster)
library(pca3d)
library(xlsx)

# 1. Load the Boolean model from your file
mono_rules <- loadNetwork("rules_extended.txt")
saveNetwork(mono_rules,file = "mono_Bmodel.net") #save the Boolean netwrok in a bnet file
#plotNetworkWiring(mono_rules)  # to plot the network

# 2. Get the attractors
# Method 1: get all the possible attractors, starting from every possible initial condition
# attr <- getAttractors(mono_rules, type = "synchronous", returnTable = T) # choose between synchronous or asynchronous

# For more than 29 genes, choose sat.restricted instead
attr <- getAttractors(mono_rules, method = "sat.restricted",maxAttractorLength = 1)

Attr_Matrix <- plotAttractors(attr, borderColor = NA, main = "Fixed point attractors")


#################### Knock-outs and perturbations
# Only NLC signals ON
Sig_NLC <- getAttractors(mono_rules, method = "sat.restricted",maxAttractorLength = 1,
                         genesOFF = c("IFNG", "GMCSF", "IL1", "LPS", "IC", "IL4", "IL13", "IL10"))
Sig_NLC_Matrix <- plotAttractors(Sig_NLC, borderColor = NA, main = "Fixed point attractors")

# Only M1 stimuli ON
Sig_M1 <- getAttractors(mono_rules, method = "sat.restricted",maxAttractorLength = 1,
                        genesOFF = c("IC", "IL4", "IL13", "IL10", "MCSF", "HMGB1"))
Sig_M1_Matrix <- plotAttractors(Sig_M1, borderColor = NA, main = "Fixed point attractors")

# Only M2 stimuli ON
Sig_M2 <- getAttractors(mono_rules, method = "sat.restricted",maxAttractorLength = 1,
                        genesOFF = c("IFNG", "GMCSF", "IL1", "LPS","MCSF", "HMGB1"))
Sig_M2_Matrix <- plotAttractors(Sig_M2, borderColor = NA, main = "Fixed point attractors")

# Knock-outs
KO <- getAttractors(mono_rules, method = "sat.restricted",maxAttractorLength = 1,
                    genesOFF = c("IRF4")) # here STAT3 is chosen to be turned OFF
KO_Matrix <- plotAttractors(KO, borderColor = NA, main = "Fixed point attractors")

# Overexpression
OE <- getAttractors(mono_rules, method = "sat.restricted",maxAttractorLength = 1,
                    genesON = c("STAT3")) # here STAT3 is chosen to be turned OFF
OE_Matrix <- plotAttractors(OE, borderColor = NA, main = "Fixed point attractors")

#########################

# Remove the inputs from attractors
inputnames <- c("IFNG", "GMCSF", "IL1", "LPS", "IC", "IL4", "IL10", "MCSF", "IL13", "HMGB1")
Rem_inputs <- Attr_Matrix$`1`[!rownames(Attr_Matrix$`1`) %in% inputnames, ]
outputnames <- c("M1", "M2", "TAM")
Rem_outs <- Rem_inputs[!rownames(Rem_inputs) %in% outputnames, ]

# Detect and remove the dublicates
# The analysis will be based on the expression of intracellular components, 

Rem_inputs_df <- as.data.frame(Rem_inputs)
Rem_dubs <- Rem_inputs_df[!duplicated(lapply(Rem_inputs_df, digest))]
Attr_sorted <- Rem_dubs[ order(row.names(Rem_dubs)), ]
Rem_dubs_mat <- as.matrix(Attr_sorted)
plot(Rem_dubs_mat, xlab = NA, ylab = NA, main = "Fixed point attractors", las = 2)


# 2.1.Organize the attractors' space
attr_frame <- as_tibble(Attr_sorted, rownames = "gene") #creates the attractors' frame
attr_frame_t <- attr_frame %>% rownames_to_column() %>% gather(var, value, -rowname) %>% spread(rowname, value)
colnames(attr_frame_t) <- attr_frame_t[nrow(attr_frame_t),]
attr_frame_t <- attr_frame_t[-nrow(attr_frame_t),]
colnames(attr_frame_t)[1] <- "Attr"

attr_frame_t$Attr <- sapply(attr_frame_t$Attr, function(attr) {
  str_split(attr,'\\.')[[1]][[1]] })

attr_frame_t_numeric <- attr_frame_t %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
attr_frame_t_bool <- attr_frame_t_numeric %>% mutate_all(funs(as.logical(.)))
attr_frame_t_bool$Attr <- attr_frame_t$Attr
attr_frame_t_bool %<>% select(Attr,everything())

# 2.2.Detect phenotype categories
# M0 category
m_0 <- attr_frame_t_bool %>% filter(!xor(M1,M2))
m0_num <- m_0 %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
m0_bio <- m0_num[- c(32, 33),]
rownames(m0_bio) <- m_0$Attr
m0_num_order <- m0_bio[order(names(m0_bio))]
m0_num_order$M1 <- NULL
m0_num_order$M2 <- NULL
m0_num_order$TAM <- NULL
m0_num_order$phenotype <- NULL
m0_mat <- as.matrix(m0_num_order)
rownames(m0_mat) <- m_0$Attr
write.xlsx(t(m0_mat), "M0_category.xlsx", row.names = T, col.names = T)
plot(m0_mat, las=2, xlab=NA, ylab=NA, main = "M0 attractors")
m0_means <- as.matrix(colMeans(m0_mat))
write.xlsx(m0_means, "M0_category.xlsx", row.names = T, col.names = T)
barplot(t(m0_means), las = 2, main = "Averaged expression profile of M0 attractors", col = "#67B0D0")
m0_sub_cat <- pheatmap(m0_mat, border = NA, color = c("#F94545", "#67C170"), cutree_rows = 2,
                       show_colnames = T,show_rownames = F, fontsize = 8, main = "M0 subcategories")

# M1 category
m_1 <- attr_frame_t_bool %>% filter(M1)
m1_num <- m_1 %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
rownames(m1_num) <- m_1$Attr
m1_num_order <- m1_num[order(names(m1_num))]
m1_num_order$M1 <- NULL
m1_num_order$M2 <- NULL
m1_num_order$TAM <- NULL
m1_num_order$phenotype <- NULL
m1_mat <- as.matrix(m1_num_order)
rownames(m1_mat) <- m_1$Attr
write.xlsx(m1_mat, "M1_category.xlsx", row.names = T, col.names = T)
plot(t(m1_mat), las=2, xlab=NA, ylab=NA, main = "M1 attractors")
m1_means <- as.matrix(colMeans(m1_mat))
write.xlsx(m1_means, "M1_category.xlsx", row.names = T, col.names = T)
barplot(t(m1_means), las = 2, main = "Averaged expression profile of M1 attractors", col = "#E9BD65")

# M2 category
m_2 <- attr_frame_t_bool %>% filter(M2)
m2_num <- m_2 %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
m2_num_order <- m2_num[order(names(m2_num))]
m2_num_order$M1 <- NULL
m2_num_order$M2 <- NULL
m2_num_order$TAM <- NULL
m2_num_order$phenotype <- NULL
m2_mat <- as.matrix(m2_num_order)
rownames(m2_mat) <- m_2$Attr
plot(m2_mat, las=2, xlab=NA, ylab=NA, main = "M2 attractors")
m2_means <- as.matrix(colMeans(m2_mat))
write.xlsx(m2_means, "M2_category.xlsx", row.names = T, col.names = T)
barplot(t(m2_means), las = 2, main = "Averaged expression profile of M2 attractors", col = "#C5ABE4")

# Search for M2 subcategories
m2_sub_cat <- pheatmap(m2_mat, border = NA, color = c("#F94545", "#67C170"), 
              cutree_rows = 4, show_colnames = T,show_rownames = F, fontsize = 8, main = "M2 subcategories")
m2_cutree <- cutree(m2_sub_cat$tree_row, k = 4)
m2_subclass <- as.data.frame(m2_mat)
all(names(m2_cutree) == rownames(m2_subclass))
m2_subclass$Clusters <- m2_cutree

M2a <- m2_subclass %>% filter(Clusters == 1)
M2a$Clusters <- NULL
m2a_means <- as.matrix(colMeans(M2a))
barplot(t(m2a_means), las = 2, main = "Averaged expression profile of M2 subcategory 1", col = "#C5ABE4")

M2b <- m2_subclass %>% filter(Clusters == 2)
M2b$Clusters <- NULL
m2b_means <- as.matrix(colMeans(M2b))
barplot(t(m2b_means), las = 2, main = "Averaged expression profile of M2 subcategory 2", col = "#C5ABE4")

M2c <- m2_subclass %>% filter(Clusters == 3)
M2c$Clusters <- NULL
m2c_means <- as.matrix(colMeans(M2c))
barplot(t(m2c_means), las = 2, main = "Averaged expression profile of M2 subcategory 3", col = "#C5ABE4")

M2d <- m2_subclass %>% filter(Clusters == 4)
M2d$Clusters <- NULL
m2d_means <- as.matrix(colMeans(M2d))
barplot(t(m2d_means), las = 2, main = "Averaged expression profile of M2 subcategory 4", col = "#C5ABE4")

# TAM/NLC category
tam <- attr_frame_t_bool %>% filter(TAM)
tam_num <- tam %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
tam_num_order <- tam_num[order(names(tam_num))]
tam_num_order$M1 <- NULL
tam_num_order$M2 <- NULL
tam_num_order$TAM <- NULL
tam_num_order$phenotype <- NULL
tam_mat <- as.matrix(tam_num_order)
rownames(tam_mat) <- tam$Attr
write.xlsx(t(tam_mat), "NLC_category.xlsx", row.names = T, col.names = T)
plot(tam_mat, las=2, xlab=NA, ylab=NA, main = "TAM attractors")
tam_means <- as.matrix(colMeans(tam_mat))
write.xlsx(tam_means, "NLC_category.xlsx", row.names = T, col.names = T)
barplot(t(tam_means), las = 2, main = "Averaged expression profile of NLC attractors", col = "#A1AB70")


# 2.3. Plot the heatmap

attr_frame_t_bool[attr_frame_t_bool$Attr %in% m_0$Attr, 'phenotype'] <- "M0"
attr_frame_t_bool[attr_frame_t_bool$Attr %in% m_1$Attr, 'phenotype'] <- "M1"
attr_frame_t_bool[attr_frame_t_bool$Attr %in% m_2$Attr, 'phenotype'] <- "M2"
attr_frame_t_bool[attr_frame_t_bool$Attr %in% tam$Attr, 'phenotype'] <- "NLC"
colors <- c("#67B0D0", "#E9BD65", "#C5ABE4", "#A1AB70")
phenotypes <- as.factor(attr_frame_t_bool$phenotype)

attr_frame_t_numeric <- attr_frame_t_bool %>% select(-c(Attr, phenotype)) %>% mutate_all(funs(as.numeric(.)))
attr_frame_t_numeric$M1 <- NULL
attr_frame_t_numeric$M2 <- NULL
attr_frame_t_numeric$TAM <- NULL

attr_frame_matrix <- as.matrix(attr_frame_t_numeric)
rownames(attr_frame_matrix) <- attr_frame_t$Attr
emr_col <- as.matrix(rownames(attr_frame_matrix))
phenotype_colors <- list(phenotype = c(M0 = "#67B0D0", M1 = "#E9BD65", M2 = "#C5ABE4", NLC = "#A1AB70"))
phenotypes_frame <- data.frame(phenotype = phenotypes)
rownames(phenotypes_frame) <- rownames(attr_frame_matrix)
pheatmap(t(attr_frame_matrix), border_color = NA, color = c("#F94545", "#67C170"), 
         legend_labels = c("OFF", "ON"), legend_breaks = 0:1, 
         annotation_col = phenotypes_frame, annotation_colors = phenotype_colors, 
         annotation_names_col = F, fontsize = 6, rownames = F, colnames = F)


#######################

# 2.4.Conctruct the pca

pca <- prcomp(attr_frame_matrix)
ind <- get_pca_ind(pca) 
deep_ana <- pca2d(pca, group = phenotypes, palette = colors, show.centroids = T, show.group.labels = T)



########################

# Usupervized categorization: clustering
# 3. Calculate the Jaccard distance between attractors
jac_dist <- dist(attr_frame_matrix, attr_frame_matrix, method = "jaccard")
m <- matrix(jac_dist, nrow = length(Rem_dubs))
plot(m, breaks = NULL, col = hcl.colors(20, palette = "PRGn"), border = NA, key = list(side = 4, cex.axis = 0.75), las = 2, 
     main = "Jaccard distance matrix", xlab = "Attractors", ylab = "Attractors")
matrica <- as.matrix(m)
colnames(matrica) <- t(rownames(attr_frame_matrix))
rownames(matrica) <- rownames(attr_frame_matrix)

h <- pheatmap(matrica, border = NA, col = hcl.colors(20, palette = "PRGn"), 
              cutree_rows = 5, cutree_cols = 5, show_colnames = T,show_rownames = F, fontsize = 3)


# 3.1. Plot the k clusters
S <- sort(cutree(h$tree_row, k = 5))
Clust_frame <- melt(lapply(split(S, names(S)), unname))
colnames(Clust_frame) <- c("Cluster", "Attr")
Clust_frame_attr <- merge(Clust_frame, attr_frame_t, by = intersect(names(Clust_frame), names(attr_frame_t)))
Clust_order <- Clust_frame_attr[order(names(Clust_frame_attr))]
Clust_order$M1 <- NULL
Clust_order$M2 <- NULL
Clust_order$TAM <- NULL

Clust1 <- Clust_order %>% filter(Cluster == 1)
Clust1_mat <- as.matrix(Clust1[3:32])
rownames(Clust1_mat) <- Clust1[,1]
plot(Clust1_mat, main = "Cluster 1", xlab = "Genes", ylab = "Attractors", border = "grey", cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75, las = 2)
cl1_num <- Clust1 %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
cl1_num$Cluster <- NULL
cl1_mat <- as.matrix(cl1_num)
cl1_means <- as.matrix(colMeans(cl1_mat))
write.xlsx(cl1_means, "Cluster1.xlsx", row.names = T, col.names = T)
barplot(t(cl1_means), las = 2, main = "Averaged expression profile of Cluster1 attractors")

Clust2 <- Clust_order %>% filter(Cluster == 2)
Clust2_mat <- as.matrix(Clust2[3:32])
rownames(Clust2_mat) <- Clust2[,1]
plot(Clust2_mat, main = "Cluster 2", xlab = "Genes", ylab = "Attractors", border = "grey", cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75, las = 2)
cl2_num <- Clust2 %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
cl2_num$Cluster <- NULL
cl2_mat <- as.matrix(cl2_num)
cl2_means <- as.matrix(colMeans(cl2_mat))
write.xlsx(cl2_means, "Cluster2.xlsx", row.names = T, col.names = T)
barplot(t(cl2_means), las = 2, main = "Averaged expression profile of Cluster2 attractors_biased")

Clust3 <- Clust_order %>% filter(Cluster == 3)
Clust3_mat <- as.matrix(Clust3[3:32])
rownames(Clust3_mat) <- Clust3[,1]
plot(Clust3_mat, main = "Cluster 3", xlab = "Genes", ylab = "Attractors", border = "grey", cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75, las = 2)
cl3_num <- Clust3 %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
cl3_num$Cluster <- NULL
cl3_mat <- as.matrix(cl3_num)
cl3_means <- as.matrix(colMeans(cl3_mat))
write.xlsx(cl3_means, "Cluster3.xlsx", row.names = T, col.names = T)
barplot(t(cl3_means), las = 2, main = "Averaged expression profile of Cluster3 attractors_biased")

Clust4 <- Clust_order %>% filter(Cluster == 4)
Clust4_mat <- as.matrix(Clust4[3:32])
rownames(Clust4_mat) <- Clust4[,1]
plot(Clust4_mat, main = "Cluster 4", xlab = "Genes", ylab = "Attractors", border = "grey", cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75, las = 2)
cl4_num <- Clust4 %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
cl4_num$Cluster <- NULL
cl4_mat <- as.matrix(cl4_num)
cl4_means <- as.matrix(colMeans(cl4_mat))
write.xlsx(cl4_means, "Cluster4.xlsx", row.names = T, col.names = T)
barplot(t(cl4_means), las = 2, main = "Averaged expression profile of Cluster4 attractors")

Clust5 <- Clust_order %>% filter(Cluster == 5)
Clust5_mat <- as.matrix(Clust5[3:32])
rownames(Clust5_mat) <- Clust5[,1]
plot(Clust5_mat, main = "Cluster 4", xlab = "Genes", ylab = "Attractors", border = "grey", cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75, las = 2)
cl5_num <- Clust5 %>% select(-1) %>% mutate_all(funs(as.numeric(.)))
cl5_num$Cluster <- NULL
cl5_mat <- as.matrix(cl5_num)
cl5_means <- as.matrix(colMeans(cl5_mat))
write.xlsx(cl5_means, "Cluster5.xlsx", row.names = T, col.names = T)
barplot(t(cl5_means), las = 2, main = "Averaged expression profile of Cluster5 attractors_biased")

