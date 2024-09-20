## ----results="hide", message=FALSE, warning=FALSE-----------------------------
# Load necessary packages
library(LimROTS, warn.conflicts = FALSE, quietly = TRUE)
library(parallel, warn.conflicts = FALSE, quietly = TRUE)



## ----results="hide", message=FALSE, warning=FALSE, eval=FALSE-----------------
## 
## # Load the dataset
## data("UPS1.Case0")
## print(UPS1.Case0)
## 
## # Set random seed for reproducibility
## set.seed(1234, kind = "default")
## 
## # Set metadata and formula for LimROTS analysis
## meta.info <- "Conc."
## B <- 1000  # Number of bootstrap samples
## K <- dim(UPS1.Case0)[1]/4  # Set the value for K based on the data size
## K <- floor(K)
## num_cores <- 10  # Number of cores for parallel processing
## cluster <- makeCluster(num_cores)  # Create a parallel cluster
## group.name <- "Conc."
## formula.str <- "~ 0 + Conc."  # Formula for group comparison
## 
## # Run LimROTS analysis with trend and robust settings enabled
## limrots.result <- LimROTS(data.exp = UPS1.Case0,
##                           B = B, K = K, meta.info = meta.info,
##                           cluster = cluster, group.name = group.name,
##                           formula.str = formula.str, trend = TRUE, robust = TRUE)


## ----results="hide", message=FALSE, warning=FALSE, echo=FALSE-----------------

file_path <- system.file("extdata", "UPS1.Case0.LimROTS.rds", package = "LimROTS")
limrots.result <- readRDS(file_path)


## ----results="hide" , message=FALSE, warning=FALSE----------------------------
## Quality Control Plots

# Plot of q-values
plot(limrots.result$q_values, main = "Q-values", xlab = "Index", ylab = "Q-value")

# Histogram of q-values
hist(limrots.result$q_values, main = "Q-value Distribution", xlab = "Q-value", col = "lightgreen", border = "white")

# Summary of q-values
summary(limrots.result$q_values)


## -----------------------------------------------------------------------------
# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf,
                                q.value = limrots.result$q_values$qvalues)

# Filter for significant proteins (FDR < 0.01)
limrots.result.df <- limrots.result.df[limrots.result.df$q.value < 0.01,]

# Mark proteins as true positives (HUMAN UPS1 proteins)
limrots.result.df$TP <- ifelse(grepl("HUMAN", limrots.result.df$proteins), "HUMAN_TP", "ECOLI.FP")

# Count the number of true positives
table(limrots.result.df$TP)


## -----------------------------------------------------------------------------
# Load ggplot2 for visualization
library(ggplot2)

# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf,
                                q.value = limrots.result$q_values$qvalues)

# Mark proteins as true positives (HUMAN UPS1 proteins)
limrots.result.df$TP <- ifelse(grepl("HUMAN", limrots.result.df$proteins), "HUMAN_TP", "ECOLI_FP")

# Create a volcano plot
ggplot(limrots.result.df, aes(x = LimROTS.FC, y = -log10(q.value), color = factor(TP))) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10 q.value", color = "True Positive") +
  scale_color_manual(values = c("grey", "red"))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "blue")+
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "black")


## ----results="hide", message=FALSE, warning=FALSE, eval=FALSE-----------------
## # Load the dataset
## data("UPS1.Case1")
## print(UPS1.Case1)
## 
## # Set random seed for reproducibility
## set.seed(1234, kind = "default")
## 
## # Set metadata and formula for LimROTS analysis
## meta.info <- c("Conc." , "inj")
## B <- 1000  # Number of bootstrap samples
## K <- dim(UPS1.Case1)[1]/4 # Set the value for K based on the data size
## K <- floor(K)
## num_cores <- 10  # Number of cores for parallel processing
## cluster <- makeCluster(num_cores)  # Create a parallel cluster
## group.name <- "Conc."
## formula.str <- "~ 0 + Conc.+ inj + Conc.*inj"  # Formula for group comparison + injections
## 
## # Run LimROTS analysis with trend and robust settings enabled
## limrots.result <- LimROTS(data.exp = UPS1.Case1,
##                           B = B, K = K, meta.info = meta.info,
##                           cluster = cluster, group.name = group.name,
##                           formula.str = formula.str, trend = TRUE, robust = TRUE)


## ----results="hide", message=FALSE, warning=FALSE, echo=FALSE-----------------

file_path <- system.file("extdata", "UPS1.Case1.LimROTS.rds", package = "LimROTS")
limrots.result <- readRDS(file_path)


## -----------------------------------------------------------------------------
# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf,
                                q.value = limrots.result$q_values$qvalues)

# Filter for significant proteins (FDR < 0.01)
limrots.result.df <- limrots.result.df[limrots.result.df$q.value < 0.01,]

dim(limrots.result.df)


## -----------------------------------------------------------------------------
# Load ggplot2 for visualization
library(ggplot2)

# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf,
                                q.value = limrots.result$q_values$qvalues)

# Create a volcano plot
ggplot(limrots.result.df, aes(x = LimROTS.FC, y = -log10(q.value))) +
  geom_point(alpha = 0.6, color = "black") +
  theme_bw() +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10 q.value") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "blue")+
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "black")


## ----results="hide", message=FALSE, warning=FALSE, eval=FALSE-----------------
## 
## # Load the dataset
## data("UPS1.Case2")
## print(UPS1.Case2)
## 
## # Set random seed for reproducibility
## set.seed(1234, kind = "default")
## 
## # Set metadata and formula for LimROTS analysis
## meta.info <- c("Conc." , "tool")
## B <- 1000  # Number of bootstrap samples
## K <- dim(UPS1.Case2)[1]/4 # Set the value for K based on the data size
## K <- floor(K)
## num_cores <- 10  # Number of cores for parallel processing
## cluster <- makeCluster(num_cores)  # Create a parallel cluster
## group.name <- "Conc."
## formula.str <- "~ 0 + Conc. + tool"  # Formula for group comparison
## 
## # Run LimROTS analysis with trend and robust settings enabled
## limrots.result <- LimROTS(data.exp = UPS1.Case2,
##                           B = B, K = K, meta.info = meta.info,
##                           cluster = cluster, group.name = group.name,
##                           formula.str = formula.str, trend = TRUE, robust = TRUE)


## ----results="hide", message=FALSE, warning=FALSE, echo=FALSE-----------------

file_path <- system.file("extdata", "UPS1.Case2.LimROTS.rds", package = "LimROTS")
limrots.result <- readRDS(file_path)


## -----------------------------------------------------------------------------
# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf,
                                q.value = limrots.result$q_values$qvalues, bh = limrots.result$BH.pvalue, p.value = limrots.result$pvalue)

# Filter for significant proteins (FDR < 0.01)
limrots.result.df <- limrots.result.df[limrots.result.df$q.value < 0.01,]

dim(limrots.result.df)


## ----results="hide", message=FALSE, warning=FALSE, eval=FALSE-----------------
## 
## # Load the dataset
## data("UPS1.Case3")
## print(UPS1.Case3)
## 
## # Set random seed for reproducibility
## set.seed(1234, kind = "default")
## 
## # Set metadata and formula for LimROTS analysis
## meta.info <- c("Conc." , "tool" , "inj")
## B <- 1000  # Number of bootstrap samples
## K <- dim(UPS1.Case3)[1]/4 # Set the value for K based on the data size
## K <- floor(K)
## num_cores <- 10  # Number of cores for parallel processing
## cluster <- makeCluster(num_cores)  # Create a parallel cluster
## group.name <- "tool"
## formula.str <- "~ 0 + tool + Conc. + inj + tool:Conc. + tool:inj"  # Formula for group comparison
## 
## # Run LimROTS analysis with trend and robust settings enabled
## limrots.result <- LimROTS(data.exp = UPS1.Case3,
##                           B = B, K = K , meta.info = meta.info,
##                           cluster = cluster, group.name = group.name,
##                           formula.str = formula.str, trend = TRUE, robust = TRUE)


## ----results="hide", message=FALSE, warning=FALSE, echo=FALSE-----------------

file_path <- system.file("extdata", "UPS1.Case1.LimROTS.rds", package = "LimROTS")
limrots.result <- readRDS(file_path)


## -----------------------------------------------------------------------------
# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf,
                                q.value = limrots.result$q_values$qvalues, bh = limrots.result$BH.pvalue, p.value = limrots.result$pvalue)

# Filter for significant proteins (FDR < 0.01)
limrots.result.df <- limrots.result.df[limrots.result.df$q.value < 0.01,]

dim(limrots.result.df)


## ----results="hide", message=FALSE, warning=FALSE, eval=FALSE-----------------
## 
## # Load the dataset
## data("UPS1.Case4")
## print(UPS1.Case4)
## 
## # Set random seed for reproducibility
## set.seed(1234, kind = "default")
## 
## # Set metadata and formula for LimROTS analysis
## meta.info <- c("Conc." , "tool" , "fake.batch")
## B <- 1000  # Number of bootstrap samples
## K <- dim(UPS1.Case4)[1]/4 # Set the value for K based on the data size
## K <- floor(K)
## num_cores <- 10  # Number of cores for parallel processing
## cluster <- makeCluster(num_cores)  # Create a parallel cluster
## group.name <- "Conc."
## formula.str <- "~ 0 + Conc. + tool + fake.batch"  # Formula for group comparison
## 
## # Run LimROTS analysis with trend and robust settings enabled
## limrots.result <- LimROTS(data.exp = UPS1.Case4,
##                           B = B, K = K, meta.info = meta.info,
##                           cluster = cluster, group.name = group.name,
##                           formula.str = formula.str, trend = TRUE, robust = TRUE)


## ----results="hide", message=FALSE, warning=FALSE, echo=FALSE-----------------

file_path <- system.file("extdata", "UPS1.Case1.LimROTS.rds", package = "LimROTS")
limrots.result <- readRDS(file_path)



## -----------------------------------------------------------------------------
# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf,
                                q.value = limrots.result$q_values$qvalues)

# Filter for significant proteins (FDR < 0.01)
limrots.result.df <- limrots.result.df[limrots.result.df$q.value < 0.01,]

# Mark proteins as true positives (HUMAN UPS1 proteins)
limrots.result.df$TP <- ifelse(grepl("HUMAN", limrots.result.df$proteins), "HUMAN_TP", "ECOLI.FP")

# Count the number of true positives
table(limrots.result.df$TP)


## -----------------------------------------------------------------------------
# Load ggplot2 for visualization
library(ggplot2)

# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf, q.value = limrots.result$q_values$qvalues)
# Mark proteins as true positives (HUMAN UPS1 proteins)
limrots.result.df$TP <- ifelse(grepl("HUMAN", limrots.result.df$proteins), "HUMAN_TP", "ECOLI_FP")

# Create a volcano plot
ggplot(limrots.result.df, aes(x = LimROTS.FC, y = -log10(q.value), color = factor(TP))) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10 q.value", color = "True Positive") +
  scale_color_manual(values = c("grey", "red"))+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "blue")+
  geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "black")


## ----results="hide", message=FALSE, warning=FALSE, eval=FALSE-----------------
## 
## # Load the dataset
## data(BM21.Case1)
## 
## print(BM21.Case1)
## 
## # Set random seed for reproducibility
## set.seed(1234, kind = "default")
## 
## # Set metadata and formula for LimROTS analysis
## meta.info <- "mixing.ratio"
## B <- 1000  # Number of bootstrap samples
## K <- dim(BM21.Case1)[1]/4  # Set the value for K based on the data size
## K <- floor(K)
## num_cores <- 30   # Number of cores for parallel processing
## cluster <- makeCluster(num_cores)  # Create a parallel cluster
## group.name <- "mixing.ratio"
## formula.str <- "~ 0 + mixing.ratio"  # Formula for group comparison
## 
## 
## # Run LimROTS analysis with trend and robust settings enabled
## limrots.result <- LimROTS(data.exp = BM21.Case1,
##                           B = B, K = K, meta.info = meta.info,
##                           cluster = cluster, group.name = group.name,
##                           formula.str = formula.str, trend = TRUE, robust = TRUE)
## 
## 


## ----results="hide", message=FALSE, warning=FALSE, echo=FALSE-----------------

file_path <- system.file("extdata", "UPS1.Case1.LimROTS.rds", package = "LimROTS")
limrots.result <- readRDS(file_path)


## -----------------------------------------------------------------------------
# Create a data frame from the LimROTS results
limrots.result.df <- data.frame(proteins = row.names(limrots.result$data),
                                LimROTS.FC = limrots.result$corrected.logf,
                                q.value = limrots.result$q_values$qvalues, bh = limrots.result$BH.pvalue, p.value = limrots.result$pvalue)

# Filter for significant proteins (FDR < 0.01)
limrots.result.df <- limrots.result.df[limrots.result.df$q.value < 0.01,]

dim(limrots.result.df)

