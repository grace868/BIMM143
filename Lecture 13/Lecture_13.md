---
title: "Lecture 13"
author: "Grace Lee"
date: "May 15, 2018"
output: 
  html_document:
    keep_md: yes
---



# Section 1: 1000 Genomes Data Analysis
Downloaded genotype data from Ensemble for Asthma associated SNP rs8067378 
See: (MXL) Mecican Ancestry in Los Angeles, California

## Question 5
Read CSV file

```r
genotype <- read.csv("Lecture 13 data 373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head(genotype)
```

```
##   Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s.
## 1                  NA19648 (F)                       A|A ALL, AMR, MXL
## 2                  NA19649 (M)                       G|G ALL, AMR, MXL
## 3                  NA19651 (F)                       A|A ALL, AMR, MXL
## 4                  NA19652 (M)                       G|G ALL, AMR, MXL
## 5                  NA19654 (F)                       G|G ALL, AMR, MXL
## 6                  NA19655 (M)                       A|G ALL, AMR, MXL
##   Father Mother
## 1      -      -
## 2      -      -
## 3      -      -
## 4      -      -
## 5      -      -
## 6      -      -
```

Find how many of each genotype there are in this dataset

```r
genotype.summary <- table(genotype[,2])
genotype.summary
```

```
## 
## A|A A|G G|A G|G 
##  22  21  12   9
```

Proportion of each genotype

```r
genotype.summary/nrow(genotype) * 100
```

```
## 
##     A|A     A|G     G|A     G|G 
## 34.3750 32.8125 18.7500 14.0625
```

**Q5** What proportion of the Mexican Ancestry in Los Angeles sample population (MXL)
are homozygous for the asthma associated SNP (G|G)?
Ans: 14%

# Section 2: RNA - Seq
Look more closely at sample HG00109 with G|G genotype for this SNP (from GBR population)
Use http://149.165.168.199/

From Lecture 13 Slide #58
Side-note on Fastq quality scores

```r
#install.packages("seqinr")
#install.packages("gtools")
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred
```

```
##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31
```


```r
prob <- 10**(-phred/10)
prob
```

```
##            D            D            D            D            C 
## 0.0003162278 0.0003162278 0.0003162278 0.0003162278 0.0003981072 
##            D            E            D            C            D 
## 0.0003162278 0.0002511886 0.0003162278 0.0003981072 0.0003162278 
##            D            D            D            B            B 
## 0.0003162278 0.0003162278 0.0003162278 0.0005011872 0.0005011872 
##            D            D            D            C            C 
## 0.0003162278 0.0003162278 0.0003162278 0.0003981072 0.0003981072 
##            @ 
## 0.0007943282
```

# Section 4

```r
expr <- read.table("rs8067378_ENSG00000172057.6.txt")
head(expr)
```

```
##    sample geno      exp
## 1 HG00367  A/G 28.96038
## 2 NA20768  A/G 20.24449
## 3 HG00361  A/A 31.32628
## 4 HG00135  A/A 34.11169
## 5 NA18870  G/G 18.25141
## 6 NA11993  A/A 32.89721
```

**Q14** Sample size for each genotype

```r
summary(expr)
```

```
##      sample     geno          exp        
##  HG00096:  1   A/A:108   Min.   : 6.675  
##  HG00097:  1   A/G:233   1st Qu.:20.004  
##  HG00099:  1   G/G:121   Median :25.116  
##  HG00100:  1             Mean   :25.640  
##  HG00101:  1             3rd Qu.:30.779  
##  HG00102:  1             Max.   :51.518  
##  (Other):456
```

```r
table(expr$geno)
```

```
## 
## A/A A/G G/G 
## 108 233 121
```

How many samples?

```r
nrow(expr)
```

```
## [1] 462
```


```r
summary( expr$exp[ (expr$geno == "G/G")])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   6.675  16.903  20.074  20.594  24.457  33.956
```

```r
summary( expr$exp[ (expr$geno == "A/A")])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   11.40   27.02   31.25   31.82   35.92   51.52
```

```r
summary( expr$exp[ (expr$geno == "A/G")])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   7.075  20.626  25.065  25.397  30.552  48.034
```

Boxplot (genotype vs expression levelss)

```r
boxplot(exp ~ geno, data = expr, notch = TRUE)
```

![](Lecture_13_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Use ggplot

```r
library(ggplot2)
```


```r
ggplot(expr, aes(geno, exp)) + geom_boxplot()
```

![](Lecture_13_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

Histogram of the exp column with ggplot2

```r
ggplot(expr, aes(exp, fill = geno)) + geom_density()
```

![](Lecture_13_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


```r
# Boxplot with the data shown
ggplot(expr, aes(geno, exp, fill=geno)) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.4)
```

![](Lecture_13_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

