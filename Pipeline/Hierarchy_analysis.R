# Load Hirnet
source('http://www.dartmouth.edu/~chaocheng/software/hirnet/HirNet_functions.R')
# Load network and node attributes
dip <- read.csv('Dataset/DIP/dip.csv', sep='\t', header = TRUE, as.is=TRUE)
mynet <- dip[, c(1,2)]
myres <- cal_hier_score26(mynet, kmax=10000, ptim=100, anneal.coeff=1e-6,  myoutf = "Hierarchical_DIP2.txt") 
mint <- read.csv('~/Dropbox/CBB752_Final_Amy_Hussein_Frank/Dataset/MINT/mint.csv', sep='\t', header = TRUE, as.is=TRUE)
mynet <- mint[, c(1,2)]
myres <- cal_hier_score26(mynet, kmax=10000, ptim=100, anneal.coeff=1e-6,  myoutf = "Hierarchical_mint.txt") 
dip_cyto <- read.csv('~/Dropbox/CBB752_Final_Amy_Hussein_Frank/Cytoscape/dip.csv', sep=',', header=TRUE, as.is=TRUE)
mint_cyto <- read.csv('~/Dropbox/CBB752_Final_Amy_Hussein_Frank/Cytoscape/mint.csv', sep=',', header = TRUE, as.is=TRUE)

# Get DIP hierarchy
hier_dip2 <- scan('Hierarchical_DIP2.txt', what = character(), sep='\n')
hier_dip2 <- scan('Hierarchical_DIP2.txt', what = character(), sep='\n', skip = grep('Lev=6', hier_dip2)+1)
hier_dip2 <- t(sapply(hier_dip2, function(x) strsplit(x, split = '\t')[[1]], simplify = TRUE))
rownames(hier_dip2) <- hier_dip2[, 1]
# Determine hierarchy
hierarchy <- apply(hier_dip2, 1, function(x) which.max(x[2:7]))
# Get the hierarchy for genes with SNP or without SNP
dip_SNP <- dip_cyto$SNP == 'Yes'
names(dip_SNP) <- dip_cyto$name
dip_with_SNP <- hierarchy[names(which(dip_SNP))]
dip_without_SNP <- hierarchy[names(which(!dip_SNP))]
# Do fisher exact test for each hierarchy to test the enrichment of SNP
DIP_result <- c()
for(i in 1:6){
  conting <- matrix(c(table(dip_with_SNP)[i], table(dip_without_SNP)[i], 
                      sum(table(dip_with_SNP)[-i]), sum(table(dip_without_SNP)[-i])), ncol=2)
  DIP_result[i] <- fisher.test(conting, alternative = 'greater')$p.value
}

# Get MINT hierarchy

hier_mint <- scan('Hierarchical_MINT.txt', what = character(), sep='\n')
hier_mint <- scan('Hierarchical_MINT.txt', what = character(), sep='\n', skip = grep('Lev=6', hier_mint)+1)
hier_mint <- t(sapply(hier_mint, function(x) strsplit(x, split = '\t')[[1]], simplify = TRUE))
rownames(hier_mint) <- hier_mint[, 1]
# Determine hierarchy
hierarchy <- apply(hier_mint, 1, function(x) which.max(x[2:7]))
# Get the hierarchy for genes with SNP or without SNP
mint_SNP <- mint_cyto$SNP == 'Yes'
names(mint_SNP) <- mint_cyto$name
mint_with_SNP <- hierarchy[names(which(mint_SNP))]
mint_without_SNP <- hierarchy[names(which(!mint_SNP))]
# Do fisher exact test for each hierarchy to test the enrichment of SNP
mint_result <- c()
for(i in 1:6){
  conting <- matrix(c(table(mint_with_SNP)[i], table(mint_without_SNP)[i], 
                      sum(table(mint_with_SNP)[-i]), sum(table(mint_without_SNP)[-i])), ncol=2)
  mint_result[i] <- fisher.test(conting, alternative = 'greater')$p.value
}

Result <- cbind(1:6, DIP_result, mint_result)
colnames(Result) <- c('Hierarchy', 'DIP-p.value', 'MINT-p.value')
write.table(Result, file = 'Hierarchy_result.csv',
            sep=',', row.names = FALSE)


