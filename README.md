# GBM: Gene-based Method using summary data
# Introduction
Genome-wide association studies (GWAS) have been widely used for identifying common variants associated with complex diseases. Due to the small effect sizes of common variants, the power to detect individual risk variants is generally low. Complementary to SNP-level analysis, a variety of gene-based association tests have been proposed. However, the power of existing gene-based tests is often dependent on the underlying genetic models, and it is not known a priori which test is optimal.  Here we proposed Gene-based Method multiple existing gene-based tests, including BT, SKAT, SKATO, ACAT, TPM, FCP, Simes, GM, ExtSimes, RTP, SimpleM, ARTP, DOT, TQ, GATES, VEGAS, aSPUs, HYST, COMBAT, ART, ART.A, ACATO, BJ, GBJ, HC, GHC, MLR, FLM and PCA. The algorithm behind thess methods are described in Shao et al.
# Example
library(data.table)
input.data = data.frame(fread("input.data.txt"))
ref.data = data.frame(fread("ref.data.txt"))
P.BT = GBM(input.data,ref.data,weight.matrix = NULL, method = "BT")
