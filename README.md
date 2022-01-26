# MCA: Multilocus association methods
===============================================================================
## Background
As an effective supplementary strategy of single-marker analysis, multilocus methods have been increasingly applied. Multilocus analysis often jointly examines a set of SNPs that are pre-defined within a functional unit such as gene to evaluate the overall association evidence at the gene level; it is thus also referred to as SNP-set or gene-based approach. Due to the usefulness, distinct SNP-set methods have been recently developed, many of which can be implemented with only GWAS summary association statistics, greatly generalizing their applicability due to the widespread availability of summary-level data. With distinct SNP-set approaches for multilocus studies, one naturally wonders which one should be chosen in practice. Moreover, existing SNP-set methods are not used without deficiencies, potential limitations include insufficient power, inability to provide statistically valid tests under certain parameter settings, and reliance on permutation sampling. Unfortunately, despite the importance of multilocus analysis in GWAS and the vast number of SNP-set methods, few comprehensive comparison studies have been performed to evaluate their effectiveness. Subsequently, due to the lack of consensus on the most suitable SNP-set method, the realization of the above advantages and benefits is to some extent currently hindered.

In the present work we sought to fill this knowledge gap by conducting a comprehensive comparison for 22 commonly-used summary-statistics based SNP-set methods in the hope that our results could serve as an important guidance for practitioners on how to choose appropriate methods for SNP-set analysis in post-GWAS era. Including: MLR: multiple linear regression; FLM: functional multiple linear regression model; HC: higher criticism test; GHC: generalized higher criticism test; BJ: Berk-Jones test; GBJ: generalized Berk-Jones test; DOT: decorrelation by orthogonal transformation; BT: burden test; SKATO: optimal sequence kernel association test; SKAT: sequence kernel association test; Simes: Simes’s test; FCP: fisher combined probability; TPM: truncated product method; RTP: rank truncated product; ART: augmented rank truncation; ART-A: adaptive augmented rank truncation; GM: gamma method; GATES: gene-based association test that uses extended Simes procedure; HMP: The harmonic mean P value test; ACAT: aggregated Cauchy association test.

**[MCA](https://github.com/biostatpzeng/MCA)** is implemented in R statistical environment.

## Example
```ruby
library(data.table)
library(dplyr)
source("main function.R")

INPUT = data.frame(fread("input.data.txt"))
REF = data.frame(fread("ref.data.txt"))

P.BT = MCA(input.data=INPUT,ref.data=REF,weight.matrix = NULL, method = "BT")
```

## Cite
Zhonghe Shao, Ting Wang, Shuiping Huang and [Ping Zeng](https://github.com/biostatpzeng) (2021). A comprehensive comparison of multilocus association methods with summary statistics in genome-wide association studies.

## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.

## Update
2022-01-25 MCA version 1.0.
