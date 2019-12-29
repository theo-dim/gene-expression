# Gene Expression
Analyzing Gene Expression Data in R.


#General Information

The central dogma of molecular biology is DNA $\rightarrow$ RNA $\rightarrow$ Protein. This means that for a given gene, DNA is transcribed into RNA and then RNA is translated into proteins. One way to characterize the level to which genes are active -- or "expressed" -- is to quantify their RNA abundances in a sample of cells. This provides a snapshot of which biological processes are active.  We can do this simultaneously for nearly every gene in the genome, which is called *genome-wide gene expression profiling*.

In this project, I am working with this type of data measured on tumor biopsies from many individuals diagnosed with different types of cancer. The specific genome-wide gene expression profiling technique considered in this project is called RNA-Seq.

The data in this project was obtained from this paper: <http://www.nature.com/nbt/journal/v33/n3/full/nbt.3080.html>

The raw gene expression measurements were transformed into a measure called "RPKM": <http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/>


```{r global_options, include=FALSE}
knitr::opts_chunk$set(warning=FALSE, message=FALSE,
                                     cache=TRUE, echo=TRUE, results='markup',
                                     fig.height=8, fig.width=12, out.width="6in",
                                     out.height="4in", dev="pdf")
```


#Part 1: Data Wrangling

```{r setup}
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(broom)
library(magrittr)
library(reshape2)
library(lattice)
library(caret)
```

```{r read_files}

glioma_melanoma <- read.table("/Users/Theo/Desktop/SML_201/Projects/project_4/glioma_melanoma.txt", sep="\t", header=TRUE, quote="")
head(glioma_melanoma)
tail(glioma_melanoma)

with_covariates <- read.table("/Users/Theo/Desktop/SML_201/Projects/project_4/with_covariates.txt", sep="\t",
                              header=TRUE, quote="")
head(with_covariates)
tail(with_covariates)

gene_ids <- read.table("gene_ids.txt", sep="\t",
                              header=TRUE, quote="")
head(gene_ids)
tail(gene_ids)

design <- read.table("design.txt", sep="\t",
                              header=TRUE, quote="")
head(design)
tail(design)
names(design) <- c("sample","cell_line","organ","disease","age","sex","ethnicity")
```

`glioma_melanoma.txt` contains just gene expression in these two cancers (glioma and melanoma), while `with_covariates.txt` contains additional cancer types which have been filtered to contain only observations where both age and sex recorded. These are two subsets of the above cited very large study. The entire dataset is difficult to fit into memory on ordinary computers, so these subsets were sectioned out of the original dataset.

```{r merge_data}

glioma_melanoma_gene_ids <- inner_join(glioma_melanoma,gene_ids, by="gene_id")
glioma_melanoma_tidy <- inner_join(glioma_melanoma_gene_ids,design,by="sample")

```

There are a lot of missing values in `glioma_melanoma.txt`. A small test illustrates the details.

```{r na_manipulation}

not_available <- glioma_melanoma_tidy == "not available"
glioma_melanoma_tidy[not_available] <- NA

```
The way R deals with missing values is with NA. Using "not available" does not really help with the various built-in functions that R has, as those recognize NA. Also, "not available" would be a string or factor, while NA is not a string or a numeric value, but a flag that indicates a missing value. This is helpful in creating vectors and manipulating data accurately. Finally, NA cannot be used in comparisons, unlike other statistical languages that assign a crazy numeric value to a missing value (looking at you SAS!), which might lead to potential errors in our data manipulation.

```{r sample_names}

glioma_melanoma_tidy <- select(glioma_melanoma_tidy,-sample)
glioma_melanoma_tidy <- glioma_melanoma_tidy %>% group_by(organ) %>% mutate(sample = paste(organ,1:n(),sep = ""), transformed_rpkm = log10(rpkm+0.5))
glioma_melanoma_tidy <- glioma_melanoma_tidy %>% select(gene_id,sample,rpkm,gene_name,cell_line,organ,disease,age,sex,ethnicity,transformed_rpkm)

```


#Part 2: Becoming familiar with the dataset.

**Genes that have gene expression measurements available in `glioma_melanoma.txt` are shown here:**
```{r}
length(unique(glioma_melanoma_tidy$gene_id))

```


**This is a histogram of the recorded ages of the individuals in `with_covariates.txt`.**

```{r}

ggplot()+geom_histogram(data= with_covariates, mapping = aes(age))

```

The RPKM gene expression data in `glioma_melanoma.txt` are gathered in a single vector and plotted on a histogram.

```{r}
RPKM <- glioma_melanoma$rpkm
hist(RPKM, breaks=10)

```
There is one tall bar on the left because the data is right-skewed.

**Locations where the data does not look Normal:**
```{r}
RPKM_log1 <- log10(RPKM+0.5)
hist(RPKM_log1, breaks=10)


```
The data does not look normal for the negative numbers. There are a lot of data between 0 and 5 and a lot especially between 0 and 1, giving the data a right skew. 

```{r}
RPKM_2 <- with_covariates$rpkm
hist(RPKM_2, breaks=10)
RPKM_log2 <- log10(RPKM_2+0.5)
hist(RPKM_log2, breaks=10)

with_covariates$transformed_RPKM <- log10(with_covariates$rpkm+0.5)
```

The with_covariates data has even more data between 0 and 1 and has a lot of extremely small numbers close to 0, which makes the data right skewed, and thus harder to tranform.


**Here, I perform a hypothesis test of whether there a mean difference in gene expression between males and females in `glioma_melanoma.txt`.**
```{r}

male <- subset(glioma_melanoma_tidy, subset = glioma_melanoma_tidy$sex =="male")
female <- subset(glioma_melanoma_tidy, subset = glioma_melanoma_tidy$sex =="female")

femlae_RPKM <- female$transformed_rpkm
male_RPKM <- male$transformed_rpkm

t.test(male_RPKM, femlae_RPKM)

var(femlae_RPKM)
var(male_RPKM)
```
Assumptions:
* The variances of the two populations are equal.
* The data are normally distributed.
* The data are independent and continuous.

(The assumptions were minimized by looking at the variance of the two samples to see if they are equal. The transformation applied to RPKM above also made the data normally distributed.)

#Part 3: Differences Between Diseases

Glioma is a cancer of the brain, while melanoma is a cancer of the skin. It is interesting to examine differences in gene expression between these two very different diseases to formulate an understanding of what is biologically different between them.

**For each gene individually, I perform a hypothesis test of whether there is a population mean difference in expression between the two cancer types.**

```{r mean_difference}

gene_p_values <- glioma_melanoma_tidy %>% group_by(gene_name) %>% do(t = t.test(.$transformed_rpkm~.$disease)$p.value)
names(gene_p_values) <- c("gene_name","p_values")
gene_p_values <- transform(gene_p_values, p_values = as.numeric(p_values))
sapply(gene_p_values, mode)
sapply(gene_p_values, class)

gene_t_values <- glioma_melanoma_tidy %>% group_by(gene_name) %>% do(t = t.test(.$transformed_rpkm~.$disease)$statistic)
names(gene_t_values) <- c("gene_name","t_statistics")
gene_t_values <- transform(gene_t_values, t_statistics = as.numeric(t_statistics))
sapply(gene_t_values, mode)
sapply(gene_t_values, class)

gene_estimates <- glioma_melanoma_tidy %>% group_by(gene_name) %>% do(t = t.test(.$transformed_rpkm~.$disease)$estimate)
names(gene_estimates) <- c("gene_name","estimates")
gene_estimates[,"lower_bound_est"] <- NA
gene_estimates[,"upper_bound_est"] <- NA

#extract bounds:
for (i in 1:nrow(gene_estimates)) {
  gene_estimates$lower_bound_est[i] <- gene_estimates[[2]][[i]][[1]]
  gene_estimates$upper_bound_est[i] <- gene_estimates[[2]][[i]][[2]]
}

#get rid of estimate intervals column: 
gene_estimates <- select(gene_estimates,-estimates)
#add effect size column:
gene_estimates <- gene_estimates %>% mutate(effect_size = lower_bound_est-upper_bound_est)
sapply(gene_estimates, mode)
sapply(gene_estimates, class)


gene_ci <- glioma_melanoma_tidy %>% group_by(gene_name) %>% do(t = t.test(.$transformed_rpkm~.$disease)$conf.int)
names(gene_ci) <- c("gene_name","confidence_intervals")

gene_statistics <- inner_join(gene_p_values,gene_t_values,by="gene_name")
gene_statistics <- inner_join(gene_statistics,gene_ci,by="gene_name")
gene_statistics <- inner_join(gene_statistics,gene_estimates,by="gene_name")

#split confidence intervals into two columns:
gene_statistics[,"lower_bound_ci"] <- NA
gene_statistics[,"upper_bound_ci"] <- NA

#extract bounds:
for (i in 1:nrow(gene_statistics)) {
  gene_statistics$lower_bound_ci[i] <- gene_statistics[[4]][[i]][[1]]
  gene_statistics$upper_bound_ci[i] <- gene_statistics[[4]][[i]][[2]]
}

#get rid of confidence_intervals column: 
gene_statistics <- select(gene_statistics,-confidence_intervals)

```

Some thoughts on scale:
For this part I used the original scale of the data. I did this to avoid writing more code that would discern the two types of cancer for each gene and executing the t.test function. Since each data point has the same meaning as the original, there is no difference in the interpretation of the t.test above.

**This is a histogram of the resulting p-values:**

```{r p_value_histogram}

#The default one:
ggplot(data=gene_p_values, aes(gene_p_values$p_values)) + 
  geom_histogram() +
  xlab("P-Values") +
  ylab("Frequency") +
  ggtitle("P-Value distribution") +
  theme_bw()

```

The distribution of the p-values seems to be greatly skewed to the right, meaning that they are for the most part extremely small, reinforcing the idea that there is a mean difference in the expression of the various genes i.e.null hypothesis is false.

**The two genes below are vastly different in how they are expressed among the two cancer types. This boxplot illustrates their value distributions:**

```{r two_genes}

#which ones are the most significant?
most_significant <- head(gene_statistics[order(gene_p_values[,2]),])
most_significant

#which ones are the least significant?
least_significant <- head(gene_statistics[order(-gene_p_values[,2]),])
least_significant

#extract data for genes of interest
important_genes <- glioma_melanoma_tidy %>% filter(gene_name == "CCL28" | gene_name == "PRAME")
important_genes <- important_genes %>% group_by(gene_name)

#assemble plot
ggplot(important_genes, aes(gene_name,transformed_rpkm)) +
  geom_boxplot() +
  xlab("Gene Name") +
  ylab("RPKM Values") +
  ggtitle("RPKM value distribution for most significant genes") +
  geom_jitter(width = 0.15, alpha = 0.2) +
  theme_bw()

```

Testing for problematic values:
```{r problematic_value}

which(is.infinite(gene_p_values$p_values))
which(is.infinite(gene_t_values$t_statistics))
which(is.nan(gene_p_values$p_values))
which(is.nan(gene_t_values$t_statistics))

```

**The following plot shows the 95% confidence intervals for the mean difference in the five most significant genes and the five least significant genes.**

```{r ci_plot}
#which ones are the most significant?
most_significant
#which ones are the least significant?
least_significant
library(ggplot2)
#plot of the t_statistics and the corresponding error bars:
ggplot(data = NULL, mapping = aes(x = gene_name, y = effect_size)) +
  geom_point(data = most_significant, color = "green") +
  geom_errorbar(data=most_significant, aes(ymin = lower_bound_ci, ymax = upper_bound_ci)) +
  geom_point(data = least_significant, color = "red") +
  geom_errorbar(data=least_significant, aes(ymin = lower_bound_ci, ymax = upper_bound_ci)) +
  ggtitle("95% CI's for the effect size of 5 most/least important genes") +
  xlab("Gene Name") +
  ylab("Effect Size") +
  theme_bw()

```
It is worthwhile noting that all of the estimates fall within their correspodnign confidence intervals. Also, all of the least significant genes have a very low effect size (makes sense) and the most significant ones are scattered on both sides.

Securing the super-senior risk of pure causality between the selected genes:
```{r p_value_theoretical}

theoretical_cutoff <- nrow(gene_p_values)*0.05
theoretical_cutoff

```

In the case where all null hypotheses are true, then everything would be purely chance and 5% of the p-values observed will be less than 0.05.

**I confirm that there is mean difference in the selected expressions by comparing the theoretical cutoff to the number of p-values less than 0.05 observed in the data:**

```{r p_value_quantity}

actual_cutoff <- gene_p_values %>% filter(p_values < 0.05)
actual_cutoff <- nrow(actual_cutoff)
actual_cutoff

```
The above quantity is the total p-values < 0.05 calculated originally. If all of the p-values observed were all larger than 0.05 then the null hypothesis could not be rejected. More values reject the null hypothesis than confirm it (8072 p-values < 0.05), leading us to the potential conclusion that there might actually be mean difference in the expression of genes among different cancer types.

This is a volcano plot, with observed mean difference (the effect size) on the x-axis and the $-{\rm log}_10$ transformed p-values on the y-axis. Link](http://www.r-bloggers.com/using-volcano-plots-in-r-to-visualize-microarray-and-rna-seq-results/).

```{r volcano_plot}

#transform t_statistics
var_glioma <- var(filter(glioma_melanoma_tidy, disease == "glioma")$rpkm)
var_glioma
var_melanoma <- var(filter(glioma_melanoma_tidy, disease == "melanoma")$rpkm)
var_melanoma
size_glioma <- nrow(filter(glioma_melanoma_tidy, disease == "glioma"))
size_glioma
size_melanoma <- nrow(filter(glioma_melanoma_tidy, disease == "melanoma"))
size_melanoma
gene_statistics <- gene_statistics %>% mutate(mean_difference = t_statistics*sqrt((var_glioma^2/size_glioma)+(var_melanoma^2/size_melanoma)))



#a simple plot that accounts for the overplotting that occurs:
ggplot(data = gene_statistics, aes(x = effect_size, y = -log10(p_values))) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  ggtitle("Volcano Plot") +
  xlab("Observed Mean Difference (effect size)") +
  ylab("-log10(P-Values)")

```
The regions that are statistically and practically significant lie on the edges of the two "tails" that were generated and, ideally, as far away from the main shape as possible. The points on either the left- or right-hand side exhibit a large magnitude fold change (putting them to the left- or right- of center) as well as high statistical significance (putting them closer to the top of the graph).

#Part 4: Modelling Gene Expression

**I begin by fitting a model regressing RPKM on organ affected and age (long execution times).**
```{r}
library(broom)

fit_model <- function(t) {
  m = lm(lm(transformed_RPKM~organ+age, t))
  return(tidy(m))
}

m <- with_covariates %>%
  group_by(gene_name) %>%
  do(fit_model(.)) 

head(m)

tail(m)
```

The intercept is the expected value of Y when Organ is at the baseline category, brain, and age is held to 0. The coefficient on organcolon, organovary, organpancreas, and organstomach all represent the change in Y with respect to the baseline category, brain, holding all else constant. The coefficient on the age term represents the expected change in Y for a 1 unit increase in age, holding all else constant. 

Sample model fit:
```{r}
ACTB <- subset(with_covariates, subset = with_covariates$gene_name=="ACTB")

model2 <- lm(transformed_RPKM~organ+age, data=ACTB)

```
- All values reported for `(Intercept)`, `age`, and `organpancreas`.  This includes `Estimate`, `Std. Error`, `t value`, and `Pr(>|t|)`
- `Multiple R-squared`
```{r}
summary(model2)$coefficients[1,]
summary(model2)$coefficients[6,]
summary(model2)$coefficients[4,]

summary(model2)$r.squared


```

The intercept represents the value of transformed RPKM for organBrain, which is 3.331520. This estimate of the transformed RPKM has a standard error of 7.258113e-02. The t value tests the hypotheses that the corresponding population parameters are 0. A large t-value, 4.590063e+01, shows that the parameters are not 0. The p-value is the probability of obtaining a t-value greater than the absolute value of the t-value we obtained. The probability of that for organBrain is 9.487170e-60, or insignificant.

The coefficient on the age variable represents the value of transformed RPKM for a 1 unit increase in age, which is 0.0002983308 added to the intercept, holding all else constant. This estimate of the transformed RPKM has a standard error of 0.0014867244. The t value tests the hypotheses that the corresponding population parameters are 0. A small t-value, 0.2006631239, shows that the parameters could be 0. The p-value is the probability of obtaining a t-value greater than the absolute value of the t-value we obtained. The probability of that for age is 0.8414649901, very big. Meaning that the coefficient on age is not statistically discernable from 0. 

The coefficient on the variable organPancreas represents the value of trasnformed RPKM compared to the baseline category organBrain which has a value of 3.331520. The organPacreas variable represents an increase of 0.07272315 in transformed RPKM compared to the intercept. This estimate of the RPKM has a standard error of 0.10115147. The t value tests the hypotheses that the corresponding population parameters are 0. A small t-value, 0.71895300, shows that the parameter could be 0. The p-value is the probability of obtaining a t-value greater than the absolute value of the t-value we obtained. The probability of that for organPancreas is 0.47423890, large. This means that the coefficient on organPancreas is not statistically discernable from 0. 

The r-squared value can be thought of as the percent of the variance in the data explained by our model. The model above has a r-squared value of  0.08310237, meaning that the model only explains 8.3% of the variance in the data.

Key assumptions on inference:
```{r}
plot(model2, which=1)
plot(model2, which=2)
```
The residuals and fitted values do appear to have some trends with respect to each other. There are 3 clusters apparent. 

The residuals do appear to be distributed normally.  

Plotting against age:
```{r}

#HEAD
ggplot(ACTB)+geom_point(mapping=aes(age, rpkm))
plot(ACTB$organ, ACTB$rpkm, xlab="Organ", ylab="RPKM")
#####################
ggplot(ACTB)+geom_point(mapping=aes(age, transformed_RPKM))

plot(ACTB$organ, ACTB$transformed_RPKM, xlab="Organ", ylab="RPKM")
#c77c2f809e696e695db121ab1f418396ae39d043
points(ACTB$organ, model2$fitted.values, col="blue", pch=20, cex=2)

```

Sample gene with large effect size:
```{r}

ATP <- subset(with_covariates, subset= with_covariates$gene_name=="ATP8")
ggplot(ATP, mapping=aes(age, transformed_RPKM))+geom_point(mapping=aes(age, transformed_RPKM))+geom_smooth(method="lm", se=FALSE, formula=y~x)

```
In the data for ATP8, RPKM increases with age. There are some irregularities, e.g. a point at 70 years of age that looks to be an outlier and there are multiple other points after the age of 50 that heavily influence the fit of the line. 

Potential variables for the linear model:
```{r}
fit_model <- function(t) {
  m = lm(lm(transformed_RPKM~organ+age+sex, t))
  return(tidy(m))
}
model.1 <- with_covariates %>%
  group_by(gene_name) %>%
  do(fit_model(.)) 
fit_model <- function(t) {
  m = lm(lm(transformed_RPKM~organ+age+sex+age:sex, t))
  return(tidy(m))
}
model.2<-with_covariates %>%
  group_by(gene_name) %>%
  do(fit_model(.)) 


```

Model results:
```{r}
TRNE <- subset(with_covariates, subset = gene_name == "TRNE")
KRT18 <- subset(with_covariates, subset = gene_name == "KRT18")
model.trne <- lm(transformed_RPKM~organ+age+sex, data=TRNE)
model.KRT18 <- lm(transformed_RPKM~organ+age+sex+age:sex, data = KRT18)

plot(TRNE$organ, TRNE$transformed_RPKM, xlab="Organ", ylab="RPKM")
points(TRNE$organ, model.trne$fitted.values, col="blue", pch=20, cex=2)

plot(KRT18$organ, KRT18$transformed_RPKM, xlab="Organ", ylab="RPKM")
points(KRT18$organ, model.KRT18$fitted.values, col="blue", pch=20, cex=2)


```

* TRNE is notable for its large coefficient for organOvary that was very significant and had the largest effect for a variable in model1. 

* KRT18 also had a large coefficient for organStomach and organPancreas that was very statistically significant and had the largest effect for a variable in model2.

**To test the two models, I use ANOVA:**
```{r}
LOC <- subset(with_covariates, subset=with_covariates$gene_name=="KRT18")
model.loc <- lm(transformed_RPKM~organ+age+sex, data=LOC)
model.loc2 <- lm(transformed_RPKM~organ+age+sex+age:sex, data=LOC)

anova(model.loc, model.loc2)

```
With a p=value of 0.3442 the null hypothesis i.e. the additional terms have coefficients equal to zero is retained. The first model is not significantly different from the second model at a 0.05 significance level. 

#Part 5: Prediction

p53 (TP53) is a well-known oncogene (gene associated with the development of cancer).

**Testing prediction accuracy:**
```{r prediction_accuracy}

set.seed(201)
TP53_data <- with_covariates %>% filter(gene_name == "TP53")
TP53_data <- TP53_data %>% mutate(transformed_rpkm = log10(rpkm+0.5))
TP53_Train <- createDataPartition(TP53_data$transformed_rpkm, p=0.6, list=FALSE)
training <-TP53_data[TP53_Train,]
testing <-TP53_data[-TP53_Train,]
TP53_model_fit_lm <- train(transformed_rpkm~age, data=training, method="lm")
TP53_model_fit_lm

#now we calculate the predictions:
predictions <- predict(TP53_model_fit_lm, newdata = testing)
predictions <- as.data.frame(predictions)
names(predictions) <- c("predicted_rpkm")

#create matrix of actual and predicted values:
actual <- testing %>% select(transformed_rpkm,age)
actual_predicted <- cbind(predictions,actual)

```

**Prediction model:**
```{r prediction_visualization}

#main plot:
ggplot(data = actual_predicted) +
  geom_point(aes(x = age, y = transformed_rpkm, shape = "b"),color = "red") +
  geom_point(aes(x = age, y = predicted_rpkm, shape = "c"),color = "green") +
  scale_shape(solid = FALSE) +
  ggtitle("Predicted and Actual RPKM Values vs. Age") +
  xlab("Age") +
  ylab("RPKM Values") +
  theme_bw() +
  theme(legend.position="none")

```

Limitations:
This model utilizes the lm method within train() does not use the testing set in any way, and is therefore a biased measure of error. lm() objects or a different package (such as MLR that has more functionalities than caret) could provide an unbiased measure of error. 

...
#Part 6: PCA and Clustering

##PCA
```{r}
glioma_melanoma$transformed_RPKM <- log10(glioma_melanoma$rpkm+0.5)

pca_dat_raw = glioma_melanoma[,c("sample","transformed_RPKM", "gene_id")]
pca_dat = acast(pca_dat_raw, gene_id ~ sample, value.var = "transformed_RPKM")

pca <- function(x, space=c("rows", "columns"),
 center=TRUE, scale=FALSE) {
 space <- match.arg(space)
 if(space=="columns") {x <- t(x)}
 x <- t(scale(t(x), center=center, scale=scale))
 s <- svd(x)
 loading <- s$u
 colnames(loading) <- paste0("Loading", 1:ncol(loading))
 rownames(loading) <- rownames(x)
 pc <- diag(s$d) %*% t(s$v)
 rownames(pc) <- paste0("PC", 1:nrow(pc))
 colnames(pc) <- colnames(x)
 pve <- s$d^2 / sum(s$d^2)
 if(space=="columns") {pc <- t(pc); loading <- t(loading)}
 return(list(pc=pc, loading=loading, pve=pve))
 }

mypca <- pca(pca_dat, space="rows")
```
PCA relies on the assumption of normally distributed data, and the transformed RPKM data is suitable for PCA. 

**Proportion of variance:**
```{r}
data.frame(Component=1:length(mypca$pve), PVE=mypca$pve) %>%
 ggplot() + geom_point(aes(x=Component, y=PVE), size=2)

sum(mypca$pve[1:45])
```
45 PC's are needed to explain 90% of the variance. 

**Sample PC-PC plotting:**
```{r}
data.frame(PC1=mypca$pc[1,], PC2=mypca$pc[2,]) %>%
 ggplot() + geom_point(aes(x=PC1, y=PC2), size=2)

data.frame(PC1=mypca$pc[2,], PC2=mypca$pc[3,]) %>%
 ggplot() + geom_point(aes(x=PC1, y=PC2), size=2)

data.frame(PC1=mypca$pc[3,], PC2=mypca$pc[4,]) %>%
     ggplot() + geom_point(aes(x=PC1, y=PC2), size=2)
```

Plotting PC1 vs PC2 shows a parabolic relationship between the two PC's and two relatively undefined clusters.

Plotting PC2 vs PC3 shows a quasi linear relationship between the two PC's and two loose clusters.

Plotting PC3 vs PC4 we see that there is a single cluster around 0,0 with some outliers and no real defined relationship.

This means that the points that are close together correspond to observations that have similar scores on the components displayed in the plot, which is mainly determined by the RPKM values that are similar for genes that come from the same sample. 

##Hierarchical Clustering

***Perform hierarchical clustering on the gene expression data in `with_covariates.txt`. Explain what you observe from this clustering.

**Hierarchical clustering:**
```{r}
with_covariates<- distinct(with_covariates, gene_name, sample)
cluster_dat_raw = with_covariates[,c("sample","rpkm", "gene_name","transformed_RPKM")]
cluster_dat = acast(cluster_dat_raw, sample~gene_name, value.var = "rpkm")

mydist <- dist(cluster_dat, method = "euclidean")
myhclust <- hclust(mydist, method="complete")
plot(myhclust)
plot(as.dendrogram(myhclust))
```

Samples that are close together are more closesly related in the expression of RPKM for a gene than samples that are further away. 

# Session Information  
  
Session information always included for reproducibility!

```{r sessionInformation, cache=FALSE}
sessionInfo()
```
