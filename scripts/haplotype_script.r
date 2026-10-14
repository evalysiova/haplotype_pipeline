#Set coordinates of interest
# region 1 => 7:114362480-114650436 / region 2 => 7:147017367-147323381 / region 3 => 7:148057724-148333434
# +-20kb / region 1 => 7:114342480-114670436 / region 2 => 7:146997367-147343381 / region 3 => 7:148037724-148353434

# Download raw file from repository
dir.create(file.path("data", "raw"), recursive = TRUE)
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/contr_phased.vcf.gz"
download.file(url, destfile = file.path("data", "raw", "contr_phased.vcf.gz"))

output_dir <- file.path("output_haplo", "figures")
dir.create(output_dir, recursive = TRUE)

library(vcfR); library(ggplot2); library(reshape2); library(tidyr); library(tidyverse)
library(randomcoloR); library(adegenet); library(ade4); library(hierfstat); library(dplyr); library(stringdist); library(ape); library(plotly); library(GWLD)

vcf <- read.vcfR(file.path("data", "raw", "contr_phased.vcf.gz"))

#check for unique IDs
with(list(ids = getID(vcf)), length(unique(ids, incomparables = NA)) == length(ids)) 

genotypes <- extract.gt(vcf, element = "GT")
fix <- getFIX(vcf)
haplo <- matrix(NA, nrow = nrow(genotypes), ncol = ncol(genotypes)*2)
haplo1 <- apply(genotypes, 2, function(x) as.numeric(substr(x, 1, 1)))
haplo2 <- apply(genotypes, 2, function(x) as.numeric(substr(x, 3, 3)))
haplo[, seq(1, ncol(haplo), 2)] <- haplo1
haplo[, seq(2, ncol(haplo), 2)] <- haplo2

rownames(haplo) <- rownames(genotypes)
colnames(haplo) <- paste0(rep(colnames(genotypes), each = 2), c("_1", "_2"))

all_names <- data.frame(var = vcf@fix[, "ID"], pos = vcf@fix[, "POS"])

#subset for 3 different regions
region1_names <- all_names[all_names[,2] >= 114342480 & all_names[,2] <= 114670436, ]
region2_names <- all_names[all_names[, 2] >= 146997367 & all_names[, 2] <= 147343381, ]
region3_names <- all_names[all_names[, 2] >= 148037724 & all_names[, 2] <= 148353434, ]

regions_names <- list(region1 = region1_names, region2 = region2_names, region3 = region3_names)
#insert only IDs for contributing variants
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/contrib_ids.txt"
download.file(url, destfile = file.path("data", "raw", "contrib_ids.txt"))
contrib_ids <- readLines(file.path("data", "raw", "contrib_ids.txt"))

for(i in seq_along(regions_names))
    {region <- names(regions_names)[i]
     ids <- regions_names[[i]][,1]
     contr_haplo <- haplo[rownames(haplo) %in% ids, ]
     assign(paste0("contr_haplo_", region), contr_haplo)
     }

haplo_all_regions <- list( region1 = contr_haplo_region1, region2 = contr_haplo_region2, region3 = contr_haplo_region3)

# Population specific
url <- "https://github.com/evalysiova/FOXP2_CNTNAP2_Population_Analysis/releases/download/v1.0/pop.csv"
download.file(url, destfile = file.path("data", "raw", "pop.csv"))
pop <- read.table(file.path("data", "raw", "pop.csv"), header = FALSE, sep = ",")
pop <- pop[rep(1:nrow(pop), each = 2), ]
suf <- rep(c("_1", "_2"), times = nrow(pop)/2)
pop[, 2] <- paste0(pop[, 2], suf)
pop_unique <- unique(pop[,3])
for(i in seq_along(pop_unique))
    {pop_s <- pop_unique[[i]]
     samples <- pop[pop[, 3] == pop_unique[[i]], ]
     ids <- samples[, 2]
     for(y in seq_along(haplo_all_regions))
            {b <- haplo_all_regions[[y]]
             reg <- names(haplo_all_regions)[y]
             haplo_p <- b[, colnames(b) %in% ids]
             assign(paste0("contr_haplo_", reg, "_", pop_s), haplo_p)
            }
    }

reg_names <- ls(pattern = "^contr_haplo_region")
contr_reg <- mget(reg_names)
names(contr_reg) <- gsub("contr_haplo_", "", names(contr_reg))

for(i in seq_along(contr_reg))
    {hap <- contr_reg[[i]]
     n <- names(contr_reg)[i]
     string <- apply(hap, 2, paste0, collapse = "")
     n_string <- length(string)
     string <- table(string)/n_string
     string <- data.frame(haplo = names(string), freq = as.numeric(string))
     assign(paste0("string_", n), string)
     }

#create total matrices for 3 regions
strs_names <- ls(pattern = "^string_")
strs <- mget(strs_names)
names(strs) <- gsub("string_", "", names(strs))
strs_r <- list(region1 = string_region1, region2 = string_region2, region3 = string_region3)

for(i in seq_along(strs_r))
    {x <- names(strs_r)[i]
     df <- data.frame()
     z <- strs_r[[i]]
     if(grepl("region1*", x))
        {haplos <- unique(strs_r[[i]][, 1])
         df <- data.frame(haplotype = haplos)
         assign(paste0(x, "_df"), df)}
     else if (grepl("region2*", x))
        {haplos <- unique(strs_r[[i]][, 1])
         df <- data.frame(haplotype = haplos)
         assign(paste0(x, "_df"), df)}
     else if (grepl("region3*", x))
        {haplos <- unique(strs_r[[i]][, 1])
         df <- data.frame(haplotype = haplos)
         assign(paste0(x, "_df"), df)
        }
    }


for(i in seq_along(strs))
    {v <- strs[[i]]
     d <- names(strs)[i]
     if(grepl("region1", d))
        {region1_df[[d]] <- v[,2][ match(region1_df[,1], v[, 1])]}
    else if (grepl("region2", d))
        {region2_df[[d]] <- v[,2][ match(region2_df[,1], v[, 1])]}
    else if (grepl("region3", d))
        {region3_df[[d]] <- v[,2][ match(region3_df[,1], v[, 1])]}
    }



superpop_map <- c(
  ACB="AFR", ASW="AFR", ESN="AFR", GWD="AFR", LWK="AFR", MSL="AFR", YRI="AFR",
  CEU="EUR", FIN="EUR", GBR="EUR", IBS="EUR", TSI="EUR",
  CDX="EAS", CHB="EAS", CHS="EAS", JPT="EAS", KHV="EAS",
  CLM="AMR", MXL="AMR", PEL="AMR", PUR="AMR",
  BEB="SAS", GIH="SAS", ITU="SAS", PJL="SAS", STU="SAS")


dfs <- list(region1 = region1_df, region2 = region2_df, region3 = region3_df)
hap_order_per_region <- list()
for(i in seq_along(dfs))
{region_df <- dfs[[i]]
reg <- names(dfs)[i]

plot_df <- region_df %>%
  select(haplotype, matches("_")) %>%
  pivot_longer(cols = where(is.numeric), names_to = "population", values_to = "freq") %>%
  mutate(population = gsub("^.*_", "", population), superpop   = superpop_map[population])

global_freq <- reg

#freq
global_hap_order <- region_df %>% arrange(desc(.data[[global_freq]])) %>% pull(haplotype)
hap_labels <- setNames(as.character(seq_along(global_hap_order)), global_hap_order)
hap_order_per_region[[i]] <- global_hap_order
plot_df <- plot_df %>% mutate(hap_label = hap_labels[as.character(haplotype)])

all_hap_labels <- as.character(seq_along(global_hap_order))
n_haps_global  <- length(all_hap_labels)

set.seed(42)
global_cols <- distinctColorPalette(n_haps_global)
names(global_cols) <- all_hap_labels

plot_df <- plot_df %>%
  mutate(superpop = factor(superpop, levels = c("AFR", "EUR", "EAS", "AMR", "SAS")))

superpops <- c("AFR", "EUR", "EAS", "AMR", "SAS")
plot_list <- list()

for (sp in superpops) 
  {sp_df <- plot_df %>% filter(superpop == sp)
  
  hap_order_sp <- sp_df %>% group_by(hap_label) %>% summarise(mean_freq = mean(freq, na.rm = TRUE)) %>%
    arrange(mean_freq) %>% pull(hap_label)
  
  sp_df <- sp_df %>% mutate(hap_label = factor(hap_label, levels = hap_order_sp))

  sp_df <- sp_df %>% arrange(population, desc(hap_label)) %>% group_by(population) %>%
    mutate(cum_freq  = cumsum(freq), label_pos = cum_freq - freq / 2, label     = ifelse(freq >= 0.05, as.character(hap_label), NA)) %>% ungroup()
  
  p <- ggplot(sp_df, aes(x = population, y = freq, fill = hap_label)) +
    geom_bar(stat = "identity", width = 0.85, colour = "white", linewidth = 0.2) +
    geom_text(aes(y = label_pos, label = label),
      size = 3, fontface = "bold", colour = "grey15", na.rm = TRUE) +
    scale_fill_manual(values = global_cols, name = "Haplotype") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0), limits = c(0, 1.01)) +
    labs(x = "Population", y = "Haplotype Frequency", title = paste0(sp)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), legend.position = "none", plot.title = element_text(face = "bold"))
  
  plot_list[[sp]] <- p
  
  ggsave(filename = file.path("output_haplo", "figures", paste0(reg, sp, "_haplotypes.png")), plot = p, width = 6, height = 5)
 }}
 
 pop_clean <- read.table(file.path("data", "raw", "pop.csv"), header = FALSE, sep = ",")
total_df <- data.frame(row.names = pop_clean[,2])

pdf(file.path("output_haplo", "figures", "haplo_analysis.pdf"))
#haplo_all_regions list
for(i in seq_along(haplo_all_regions))
  {contr_haplo <- haplo_all_regions[[i]]
   reg <- names(haplo_all_regions)[i]
   hap <- t(contr_haplo)
   hap1 <- hap[grepl("_1$", rownames(hap)), ]
   rownames(hap1) <- gsub("_1$", "", rownames(hap1))
   hap2 <- hap[grepl("_2$", rownames(hap)), ]
   rownames(hap2) <- gsub("_2$", "", rownames(hap2))
   hap1_string <- apply(hap1, 1, paste, collapse = "")
   hap2_string <- apply(hap2, 1, paste, collapse = "")
   genot_df <- data.frame(hap = paste(hap1_string, hap2_string, sep = "/"), row.names = rownames(hap1))
   total_df <- cbind(total_df, genot_df)
   genind <- df2genind(X = genot_df, pop = pop_clean[, 3], sep = "/")
   assign(paste0("genind_", reg), genind)
   genpop <- genind2genpop(genind)
   assign(paste0("genpop_", reg), genpop)}

colnames(total_df) <- c("hap1", "hap2", "hap3")
genind_total <- df2genind(X = total_df, pop = pop_clean[, 3], sep = "/")

X <- tab(genind_total, freq=TRUE, NA.method = "mean")
pca1 <- dudi.pca(genind_total, scannf = FALSE, nf = 3)
s.class(pca1$li, fac = pop(genind_total), xax = 1, yax = 2, col = transp(Col1, 0.6), axesell = FALSE, cstar = 0, cpoint = 3, grid = FALSE)

#Filter outliers
apply(pca1$li, 2, range)
#       Axis1     Axis2      Axis3
#[1,] -3.196912 -3.989578 -28.919199
#[2,] 60.921483 79.421274   8.557886
rownames(pca1$li)[which(pca1$li[,1] > 20)]
#[1] "HG02944"
rownames(pca1$li)[which(pca1$li[,2] > 20)]
#[1] "HG02944"
rownames(pca1$li)[which(pca1$li[,3] < -10)]
#[1] "NA19144"
outliers <- c("HG02944", "NA19144")
genind_clean <- genind_total[!indNames(genind_total) %in% outliers, ]
pca2 <- dudi.pca(genind_clean, scannf = FALSE, nf = 3)
s.class(pca2$li, fac = pop(genind_clean), xax = 1, yax = 2, col = transp(Col1, 0.6), axesell = FALSE, cstar = 0, cpoint = 3, grid = FALSE)

#Second round of filtering
apply(pca2$li, 2, range)
#Axis1     Axis2     Axis3
#[1,] -3.936060 -8.297135 -4.763799
#[2,]  6.139918 25.062148 75.828089
rownames(pca2$li)[which(pca2$li[,2] > 20)]
#[1] "HG02282"
rownames(pca2$li)[which(pca2$li[,3] > 30)]
#[1] "HG02282"
outliers2 <- c("HG02282")
genind_clean2 <- genind_clean[!indNames(genind_clean) %in% outliers2, ]
pca3 <- dudi.pca(genind_clean2, scannf = FALSE, nf = 3)
s.class(pca3$li, pop(genind_clean2),xax=1,yax=2, col=Col1)
title("PCA of diversity \ axes 1-2")
s.class(pca3$li, fac = pop(genind_clean2), xax=1,yax=2, col=transp(col,.6), axesell=FALSE,cstar=0, cpoint=3, grid=FALSE)
colorplot(pca3$li[c(1,2)], pca3$li, transp=TRUE, cex=1.5, xlab="PC 1", ylab="PC 2")
title("PCA of human genes \ axes 1-2")
abline(v=0,h=0,col="grey", lty=2)


genpop_clean2 <- genind2genpop(genind_clean2)
ca <- dudi.coa(tab(genpop_clean2), scannf = FALSE, nf = 2)

barplot(ca$eig, main = "CA eigenvalues ", col = heat.colors(length(ca$eig)))
s.label(ca$li, sub = "CA 1,2", csub = 2)
add.scatter.eig(ca$eig, nf = 3, xax = 1, yax = 3, posi = "bottomright")

grp <- find.clusters(genind_clean2, max.n.clust = 100, n.pca = 150, n.clust = 5)
#150 PCs
#5 clusters
dapc1 <- dapc(genind_clean2, grp$grp, , scannf = FALSE, n.pca = 100, n.da = 4)
#100 PCs - 4 discr

myCol <- c("orange", "red", "blue", "green", "pink")
myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  par(mgp = c(1.2, 0.5, 0))
  plot(temp, col=rep(c("black", "lightgray"),
       c(dapc1$n.pca,1000)), ylim=c(0,100), xlab="PCA axis", ylab="Cumulated variance",
       cex=0.5, pch=10, type="h", lwd=2, cex.lab = 0.7, cex.axis = 0.7) }

scatter(dapc1,xax=1,yax=3, ratio.pca=0.3, bg="white", pch=20, cell=0, cstar=0, col=Col1,solid=.4, cex=1.5, clab=0, mstree=TRUE, scree.da=FALSE
        , posi.pca="topright",leg=TRUE,txt.leg=paste("Cluster",1:5))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,3], pch=4, cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,3], pch=4, cex=3, lwd=2, col=myCol)
add.scatter(myInset(), posi="topleft", inset=c(-0.05,-0.1), ratio=.10, bg=transp("white"))
scatter(dapc1,3,3, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

dapc2 <- dapc(genind_clean2, n.da = 4, n.pca = 150)
#150 PCs

optim_pcs <- optim.a.score(dapc2)

dapc3 <- dapc(genind_clean2, n.pca = 36, n.da = 4)
scatter(dapc3,xax=1,yax=2, col=transp(myPal(8)), scree.da=FALSE, cell=1.5, cex=1, bg="white",cstar=0, legend=TRUE)
compoplot(dapc3, cleg = .5, posi = list(x = 0, y = 1.2, lab = pop(genind_clean2)))

contrib <- loadingplot(dapc3$var.contr, axis=2, thres=0.04, lab.jitter=1)
contr_haplo <- contrib$var.names
loadings <- dapc3$var.contr[, 2]
loadingplot(dapc3$var.contr, axis = 2, thres = 0.04, lab.jitter = 1, lab = " ")

df2 <- data.frame(haplotype = names(loadings), loading = loadings, index = seq_along(loadings))

haplo_palette <- c(
  "hap1.0011011"= "seagreen", 
  "hap1.0111011" = "palegreen3",
  "hap1.0000000" = "aquamarine",
  "hap1.0100000" = "darkolivegreen1",
  "hap2.00000000010000000000000000000000100" = "magenta4",
  "hap2.00000000000000000000000000000000000" = "magenta2",
  "hap2.10101111101111111111111111111111011" = "purple3",
  "hap2.10000100010000000000000000000000100" = "lightpink3",
  "hap3.0000000000000000000000000000000001" = "royalblue4", 
  "hap3.0000000000000000000000000000000010" = "skyblue", 
  "hap3.0000000000000000000000000000000000" = "royalblue"
)


alleles_df <- data.frame()
vcf <- as.data.frame(getFIX(vcf))
al <- data.frame(VAR = vcf$ID, POS = vcf$POS, REF = vcf$REF, ALT = vcf$ALT, stringsAsFactors = FALSE)
alleles_df <- rbind(alleles_df, al)
alleles_df <- alleles_df[order(alleles_df[, 2]), ]

al_con <- alleles_df[alleles_df$VAR %in% contrib_ids, ]
r1_df <- al_con[al_con[, 2] >= 114342480 & al_con[,2] <= 114670436, ]
r1_df <- separate(r1_df, ALT, into = c("ALT", "ALT2"), sep = ",", fill = "right")
r2_df <- al_con[al_con[, 2] >= 146997367 & al_con[, 2] <= 147343381, ]
r2_df <- separate(r2_df, ALT, into = c("ALT", "ALT2"), sep = ",", fill = "right")
r3_df <- al_con[al_con[, 2] >= 148037724 & al_con[, 2] <= 148353434, ]
r3_df <- separate(r3_df, ALT, into = c("ALT", "ALT2"), sep = ",", fill = "right")
r1_df$ALT <- c("C", "G----", "G-", "C", "G-", "TA", "C")
r2_df$REF[1] <- "A--"; r2_df$ALT[1] <- "AT-"
r2_df$ALT2[is.na(r2_df$ALT2)] <- "-"
r3_df$ALT[11] <- "T--"; r3_df$ALT[14] <- "T---"
r3_df$ALT2[is.na(r3_df$ALT2)] <- "-"; r3_df$ALT[29] <- "A-"; r3_df$ALT[30] <- "C--"; r3_df$ALT[31] <- "A--"
r3_df$REF[32] <- "A-"; r3_df$ALT[33] <- "C--"; r3_df$ALT[34] <- "G-"

r_dfs <- list(hap1 = r1_df, hap2 = r2_df, hap3 = r3_df)
contr_haplo_list <- split(gsub(".*\\.", "", contr_haplo), gsub("\\..*", "", contr_haplo))
total_haplo_seq <- list()

for(r in seq_along(r_dfs))
  {r_df <- r_dfs[[r]]
   haplos_i <- contr_haplo_list[[r]]
   total_haplo_i <- c()
   for(l in seq_along(haplos_i))
      {haplo0 <- strsplit((haplos_i), "")[[l]]
       mat <- matrix(ncol = length(haplo0), nrow = 2)
       mat[1, ] <- r_df[,1]; mat[2, ] <- haplo0
       for(y in 1:ncol(mat))
          {id <- mat[1, y]
           ref <- r_df$REF[r_df$VAR == id][1]
           alt <- r_df$ALT[r_df$VAR == id][1]
           alt2 <- r_df$ALT2[r_df$VAR == id][1]
           if(mat[2, y] == "0")    {mat[2, y] <- ref}
           else if(mat[2, y] == "1")    {mat[2, y] <- alt}
           else if(mat[2, y] == "2")    {mat[2, y] <- alt2}
          }
    temp_matrix <- mat[2, ]
    temp_string <- paste(temp_matrix, collapse = "")
    total_haplo_i <- c(total_haplo_i, temp_string)
    }
  total_haplo_seq[[r]] <- total_haplo_i
} 

haplo_seq_df <- data.frame(region = rep(seq_along(total_haplo_seq), lengths(total_haplo_seq)), haplotype = unlist(total_haplo_seq), row.names = NULL)
#edit for region1 (coexisting variants rs71999274 & rs7799652)

haplo_seq_df$haplotype[1] <- "TGTTTAG-CTAC"
haplo_seq_df$haplotype[2] <- "TG----G-CTAC"
haplo_seq_df$haplotype[3] <- "TGTTTAGATGAA"
haplo_seq_df$haplotype[4] <- "TG----GATGAA"

haplo_long <- haplo_long <- haplo_seq_df %>% rowwise() %>% mutate(pos  = list(seq_len(nchar(haplotype))),
         char = list(strsplit(haplotype, "")[[1]])) %>% unnest(c(pos, char)) %>% ungroup()

haplo_palette2 <- setNames(c("seagreen", "palegreen3", "aquamarine", "darkolivegreen1", "magenta4",
                   "magenta2", "purple3", "lightpink3", "royalblue4", "skyblue", "royalblue"), haplo_seq_df$haplotype)

ggplot(haplo_long, aes(x = pos, y = haplotype, fill = haplotype)) +
  geom_tile(colour = "grey8", linewidth = 0) +                    
  geom_text(aes(label = char), colour = "grey8", family = "Courier", fontface = "bold", size = 3) +
  scale_fill_manual(values = haplo_palette2) +
  facet_wrap(~region, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    strip.text = element_blank(), legend.position = "none")

# PCoA, Hamming and visualization for all haplotypes
all_haplo_list <- genind_clean2@all.names
all_haplo_list <- lapply(seq_along(all_haplo_list), function(r)
  {all_haplo_list[[r]][match(hap_order_per_region[[r]], all_haplo_list[[r]])]})
all_haplo_seq <- list()

for(r in seq_along(r_dfs))
  {r_df <- r_dfs[[r]]
   haplos_i <- all_haplo_list[[r]]
   total_haplo_i <- c()
   for(l in seq_along(haplos_i))
      {haplo0 <- strsplit((haplos_i), "")[[l]]
       mat <- matrix(ncol = length(haplo0), nrow = 2)
       mat[1, ] <- r_df[,1]
       mat[2, ] <- haplo0
       for(y in 1:ncol(mat))
        {id <- mat[1, y]
         ref <- r_df$REF[r_df$VAR == id][1]
         alt <- r_df$ALT[r_df$VAR == id][1]
         alt2 <- r_df$ALT2[r_df$VAR == id][1]
         if(mat[2, y] == "0")
            {mat[2, y] <- ref}
         else if(mat[2, y] == "1")
            {mat[2, y] <- alt}
         else if(mat[2, y] == "2")
            {mat[2, y] <- alt2}
    }
    temp_matrix <- mat[2, ]
    temp_string <- paste(temp_matrix, collapse = "")
    total_haplo_i <- c(total_haplo_i, temp_string)
}
  all_haplo_seq[[r]] <- total_haplo_i
} 

all_haplo_seq_df <- data.frame(region = rep(seq_along(all_haplo_seq), lengths(all_haplo_seq)), haplotype = unlist(all_haplo_seq), row.names = NULL)

# Hamming Distances

hamming2 <- lapply(seq_along(all_haplo_seq), function(i)
  {seq <- all_haplo_seq[[i]]
  mat_i <- stringdistmatrix(seq, seq, method = "hamming")
  rownames(mat_i) <- colnames(mat_i) <- paste0("hap", seq_along(seq))
  mat_i})


pcoa_per_region <- lapply(hamming2, function(mat) {
  pcoa_res <- ape::pcoa(mat)
  pcoa_res
})

coords_per_region <- lapply(pcoa_per_region, function(p) p$vectors[, 1:3])

var_per_region <- lapply(pcoa_per_region, function(p) {
  eigenvalues <- p$values$Eigenvalues
  eigenvalues <- pmax(eigenvalues, 0)
  (eigenvalues / sum(eigenvalues)) * 100
})

contributing_per_region <- lapply(1:3, function(i) 
 {contrib_seqs <- haplo_seq_df$haplotype[haplo_seq_df$region == i]
  contrib_idx <- which(all_haplo_seq[[i]] %in% contrib_seqs)
  paste0("hap", contrib_idx)
 })

colors_l1.2 <- setNames(c( "hap1" = "seagreen", "hap6" = "palegreen3", "hap9" = "aquamarine", "hap23" = "darkolivegreen1"), contributing_per_region[[1]])
colors_l2.2 <- setNames(c("hap7" = "magenta4", "hap9" = "magenta2", "hap81" = "purple3", "hap109" = "lightpink3"), contributing_per_region[[2]])
colors_l3.2 <- setNames(c("hap1" = "royalblue4", "hap2" = "skyblue", "hap3" = "royalblue"), contributing_per_region[[3]])

color_groups <- list(colors_l1.2, colors_l2.2, colors_l3.2)

# Contributing haplotypes per region
for (i in 1:3) 
{coords_i <- coords_per_region[[i]]
 var_i <- var_per_region[[i]]
 hap_names <- paste0("hap", seq_len(nrow(coords_i)))
 contributing <- contributing_per_region[[i]]
 contrib_colors <- color_groups[[i]]
 point_colors <- rep("gray40", length(hap_names))
 names(point_colors) <- hap_names
 point_colors[contributing] <- contrib_colors[contributing]
 point_colors_hex <- sapply(point_colors, col2hex)
  
 fig <- plot_ly()
  for (j in seq_len(nrow(coords_i))) 
  {fig <- fig %>% add_trace(type = "scatter3d", mode = "markers", x = coords_i[j, 1], 
           y = coords_i[j, 2], z = coords_i[j, 3], text = hap_names[j], hoverinfo = "text",
           marker = list(color = point_colors_hex[j], 
           size = ifelse(hap_names[j] %in% contributing, 7, 4)), showlegend = FALSE)
  }

  contrib_idx <- which(hap_names %in% contributing)
  for (j in contrib_idx) 
  {col_rgb <- col2rgb(point_colors[j])
    fill_color <- sprintf("rgba(%d,%d,%d,0.7)", col_rgb[1], col_rgb[2], col_rgb[3])
    fig <- fig %>% add_trace(type = "scatter3d", mode = "markers",
      x = coords_i[j, 1], y = coords_i[j, 2], z = coords_i[j, 3], text = hap_names[j], hoverinfo = "text", 
      marker = list(color = fill_color, size = 12, width = 0), showlegend = FALSE)
  } 

  fig <- fig %>% layout(
    title = paste("Locus", i),
    scene = list(xaxis = list(title = paste0("PCoA1 (", round(var_i[1], 1), "%)")),
                 yaxis = list(title = paste0("PCoA2 (", round(var_i[2], 1), "%)")),
                 zaxis = list(title = paste0("PCoA3 (", round(var_i[3], 1), "%)"))))
  print(fig)

}

for (i in 1:3) 
{coords_i <- coords_per_region[[i]]
  var_i <- var_per_region[[i]]
  hap_names <- paste0("hap", seq_len(nrow(coords_i)))
  contributing <- contributing_per_region[[i]]
  contrib_colors <- color_groups[[i]]
  
  point_colors <- rep("gray40", length(hap_names))
  names(point_colors) <- hap_names
  point_colors[contributing] <- contrib_colors[contributing]
  point_colors_hex <- sapply(point_colors, col2hex)
  
  plot(coords_i[, 1], coords_i[, 2],
       xlab = paste0("PCoA1 (", round(var_i[1], 1), "%)"),
       ylab = paste0("PCoA2 (", round(var_i[2], 1), "%)"),
       main = paste0("Locus ", i, " (", round(sum(var_i[1:2]), 1), "% variance)"),
       pch = 19, col = point_colors_hex, cex = 0.8,
       xlim = range(coords_i[, 1]) * 1.3,
       ylim = range(coords_i[, 2]) * 1.3)
  
  contrib_idx <- which(hap_names %in% contributing)
  points(coords_i[contrib_idx, 1], coords_i[contrib_idx, 2],
         pch = 19, col = point_colors_hex[contrib_idx], cex = 1.5)

  points(coords_i[contrib_idx, 1], coords_i[contrib_idx, 2],
         pch = 21, col = point_colors_hex[contrib_idx],
         bg = NA, cex = 2.5, lwd = 1.5)
  abline(h = 0, v = 0, lty = 2, col = "grey70")
}


# PCoA, Hamming and visualization for all haplotypes
all_haplo_list <- genind_clean2@all.names
all_haplo_list <- lapply(seq_along(all_haplo_list), function(r)
  {all_haplo_list[[r]][match(hap_order_per_region[[r]], all_haplo_list[[r]])]})
all_haplo_seq <- list()

for(r in seq_along(r_dfs))
  {r_df <- r_dfs[[r]]
   haplos_i <- all_haplo_list[[r]]
   total_haplo_i <- c()
   for(l in seq_along(haplos_i))
      {haplo0 <- strsplit((haplos_i), "")[[l]]
       mat <- matrix(ncol = length(haplo0), nrow = 2)
       mat[1, ] <- r_df[,1]
       mat[2, ] <- haplo0
       for(y in 1:ncol(mat))
        {id <- mat[1, y]
         ref <- r_df$REF[r_df$VAR == id][1]
         alt <- r_df$ALT[r_df$VAR == id][1]
         alt2 <- r_df$ALT2[r_df$VAR == id][1]
         if(mat[2, y] == "0")
            {mat[2, y] <- ref}
         else if(mat[2, y] == "1")
            {mat[2, y] <- alt}
         else if(mat[2, y] == "2")
            {mat[2, y] <- alt2}
    }
    temp_matrix <- mat[2, ]
    temp_string <- paste(temp_matrix, collapse = "")
    total_haplo_i <- c(total_haplo_i, temp_string)
}
  all_haplo_seq[[r]] <- total_haplo_i
} 

all_haplo_seq_df <- data.frame(region = rep(seq_along(all_haplo_seq), lengths(all_haplo_seq)), haplotype = unlist(all_haplo_seq), row.names = NULL)

# Hamming Distances

hamming2 <- lapply(seq_along(all_haplo_seq), function(i)
  {seq <- all_haplo_seq[[i]]
  mat_i <- stringdistmatrix(seq, seq, method = "hamming")
  rownames(mat_i) <- colnames(mat_i) <- paste0("hap", seq_along(seq))
  mat_i})


pcoa_per_region <- lapply(hamming2, function(mat) {
  pcoa_res <- ape::pcoa(mat)
  pcoa_res
})

coords_per_region <- lapply(pcoa_per_region, function(p) p$vectors[, 1:3])

var_per_region <- lapply(pcoa_per_region, function(p) {
  eigenvalues <- p$values$Eigenvalues
  eigenvalues <- pmax(eigenvalues, 0)
  (eigenvalues / sum(eigenvalues)) * 100
})

contributing_per_region <- lapply(1:3, function(i) 
 {contrib_seqs <- haplo_seq_df$haplotype[haplo_seq_df$region == i]
  contrib_idx <- which(all_haplo_seq[[i]] %in% contrib_seqs)
  paste0("hap", contrib_idx)
 })

colors_l1.2 <- setNames(c( "hap1" = "seagreen", "hap6" = "palegreen3", "hap9" = "aquamarine", "hap23" = "darkolivegreen1"), contributing_per_region[[1]])
colors_l2.2 <- setNames(c("hap7" = "magenta4", "hap9" = "magenta2", "hap81" = "purple3", "hap109" = "lightpink3"), contributing_per_region[[2]])
colors_l3.2 <- setNames(c("hap1" = "royalblue4", "hap2" = "skyblue", "hap3" = "royalblue"), contributing_per_region[[3]])

color_groups <- list(colors_l1.2, colors_l2.2, colors_l3.2)

# Contributing haplotypes per region

col2hex <- function(col) {
  rgb_vals <- col2rgb(col)
  rgb(rgb_vals[1, ], rgb_vals[2, ], rgb_vals[3, ], maxColorValue = 255)
} 

for (i in 1:3) 
{coords_i <- coords_per_region[[i]]
 var_i <- var_per_region[[i]]
 hap_names <- paste0("hap", seq_len(nrow(coords_i)))
 contributing <- contributing_per_region[[i]]
 contrib_colors <- color_groups[[i]]
 point_colors <- rep("gray40", length(hap_names))
 names(point_colors) <- hap_names
 point_colors[contributing] <- contrib_colors[contributing]
 point_colors_hex <- sapply(point_colors, col2hex)
  
 fig <- plot_ly()
  for (j in seq_len(nrow(coords_i))) 
  {fig <- fig %>% add_trace(type = "scatter3d", mode = "markers", x = coords_i[j, 1], 
           y = coords_i[j, 2], z = coords_i[j, 3], text = hap_names[j], hoverinfo = "text",
           marker = list(color = point_colors_hex[j], 
           size = ifelse(hap_names[j] %in% contributing, 7, 4)), showlegend = FALSE)
  }

  contrib_idx <- which(hap_names %in% contributing)
  for (j in contrib_idx) 
  {col_rgb <- col2rgb(point_colors[j])
    fill_color <- sprintf("rgba(%d,%d,%d,0.7)", col_rgb[1], col_rgb[2], col_rgb[3])
    fig <- fig %>% add_trace(type = "scatter3d", mode = "markers",
      x = coords_i[j, 1], y = coords_i[j, 2], z = coords_i[j, 3], text = hap_names[j], hoverinfo = "text", 
      marker = list(color = fill_color, size = 12, width = 0), showlegend = FALSE)
  } 

  fig <- fig %>% layout(
    title = paste("Locus", i),
    scene = list(xaxis = list(title = paste0("PCoA1 (", round(var_i[1], 1), "%)")),
                 yaxis = list(title = paste0("PCoA2 (", round(var_i[2], 1), "%)")),
                 zaxis = list(title = paste0("PCoA3 (", round(var_i[3], 1), "%)"))))
  print(fig)

}

for (i in 1:3) 
{coords_i <- coords_per_region[[i]]
  var_i <- var_per_region[[i]]
  hap_names <- paste0("hap", seq_len(nrow(coords_i)))
  contributing <- contributing_per_region[[i]]
  contrib_colors <- color_groups[[i]]
  
  point_colors <- rep("gray40", length(hap_names))
  names(point_colors) <- hap_names
  point_colors[contributing] <- contrib_colors[contributing]
  point_colors_hex <- sapply(point_colors, col2hex)
  
  plot(coords_i[, 1], coords_i[, 2],
       xlab = paste0("PCoA1 (", round(var_i[1], 1), "%)"),
       ylab = paste0("PCoA2 (", round(var_i[2], 1), "%)"),
       main = paste0("Locus ", i, " (", round(sum(var_i[1:2]), 1), "% variance)"),
       pch = 19, col = point_colors_hex, cex = 0.8,
       xlim = range(coords_i[, 1]) * 1.3,
       ylim = range(coords_i[, 2]) * 1.3)
  
  contrib_idx <- which(hap_names %in% contributing)
  points(coords_i[contrib_idx, 1], coords_i[contrib_idx, 2],
         pch = 19, col = point_colors_hex[contrib_idx], cex = 1.5)

  points(coords_i[contrib_idx, 1], coords_i[contrib_idx, 2],
         pch = 21, col = point_colors_hex[contrib_idx],
         bg = NA, cex = 2.5, lwd = 1.5)
  abline(h = 0, v = 0, lty = 2, col = "grey70")
}

library(pheatmap) 

reg <- names(haplo_all_regions)[1]
gi  <- get(paste0("genind_", reg))
gp  <- get(paste0("genpop_", reg))

hapstr <- alleles(gi)[[1]]
counts <- t(gp@tab)


tot    <- colSums(counts)

keep   <- rowSums(counts) >= 1
hapstr <- hapstr[keep]; counts <- counts[keep, , drop = FALSE]

freq <- sweep(counts, 2, tot, "/")

# Order haplotypes per Hamming distance
bin <- matrix(as.integer(do.call(rbind, strsplit(hapstr, ""))), nrow = length(hapstr))
hc  <- hclust(dist(bin, method = "manhattan"), method = "average")

sp_map <- c(ACB="AFR", ASW="AFR", ESN="AFR", GWD="AFR", LWK="AFR", MSL="AFR", YRI="AFR",
            CLM="AMR", MXL="AMR", PEL="AMR", PUR="AMR",
            CDX="EAS", CHB="EAS", CHS="EAS", JPT="EAS", KHV="EAS",
            CEU="EUR", FIN="EUR", GBR="EUR", IBS="EUR", TSI="EUR",
            BEB="SAS", GIH="SAS", ITU="SAS", PJL="SAS", STU="SAS")
ann <- data.frame(superpop = unname(sp_map[colnames(freq)]), row.names = colnames(freq))

pheatmap(freq, cluster_rows = hc, cluster_cols = TRUE, annotation_col = ann, annotation_names_col = FALSE, show_rownames = FALSE, border_color = NA)

# Co-occurence of variants (pairs)
snps <- c("rs7782412", "rs7799652")

vcf <- read.vcfR(file.path("data", "raw", "contr_phased.vcf.gz"))
gt  <- extract.gt(vcf, element = "GT")[snps, ]

h1 <- sub("\\|.*", "", gt)
h2 <- sub(".*\\|", "", gt)

pop_vec <- pop(gen)[match(colnames(gt), indNames(gen))]

hap_freq <- data.frame(
  pop       = rep(pop_vec, 2),
  haplotype = c(paste0(h1[snps[1], ], h1[snps[2], ]),
                paste0(h2[snps[1], ], h2[snps[2], ]))
) |>
  dplyr::count(pop, haplotype) |> group_by(pop) |> mutate(frq = n / sum(n)) |>
  ungroup()

pop_order <- names(sp_map)
hap_freq$pop <- factor(hap_freq$pop, levels = pop_order)

ggplot(hap_freq, aes(x = pop, y = frq, fill = haplotype)) +
  geom_col() +
  labs(x = NULL, y = "Frequency", fill = "Haplotype") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))     

# Tree
library(ape)

genind_var <- vcfR2genind(vcf)
pop(genind_var) <- pop_clean[, 3]
genpop_var <- genind2genpop(genind_var)

distance <- dist.genpop(genpop_var, method = 1)

plot(hclust(distance), type = "unrooted")
