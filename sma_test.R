library(reticulate)
library(SJD)
library(tidyverse)
use_condaenv("base")

source_python("/Users/Vermouth/Desktop/Thesis_Project/sma_code/sma_single.py")
print("Start")
filelist <- list.files(pattern = "\\.txt$")
inVitro_bulk <- read.table(filelist[[2]], stringsAsFactors = FALSE, header = TRUE) %>% select(-1) %>% as.matrix
inVitro_sc <- read.table(filelist[[3]], stringsAsFactors = FALSE, header = TRUE) %>%  select(-1) %>% as.matrix
inVivo_bulk <- read.table(filelist[[5]], stringsAsFactors = FALSE, header = TRUE) %>% select(-1) %>% as.matrix
inVivo_sc <- read.table(filelist[[7]], stringsAsFactors = FALSE, header = TRUE) %>% select(-1) %>% as.matrix

group = list(c(1, 2, 3, 4))
comp_num = c(40)

SMA <- function(dataset, group, comp_num, encoder_config, decoder_config, train_epoch = 1000){
  
  ## Obtain names for dataset, gene and samples
  dataset_name = datasetNameExtractor(dataset)
  gene_name = geneNameExtractor(dataset)
  sample_name = sampleNameExtractor(dataset)
  group_name = groupNameExtractor(group)
  
  dataset = frameToMatrix(dataset)

  ## Parameters to be initialized
  N = length(dataset)
  K = length(group)
  M = sum(comp_num)
  p = nrow(dataset[[1]])
  N_dataset = unlist(lapply(dataset, ncol))
  ## Output the component and scores
  list_component = list()
  list_score = list()
  removed_common = list()
  for(j in 1 : N){
    list_score[[j]] = list()
  }
  
  for(i in 1 : K){
    list_component[[i]] = matrix(0, nrow = p, ncol = comp_num[i])
    for(j in 1 : N){
      list_score[[j]][[i]] = matrix(0, nrow = comp_num[i], ncol = N_dataset[j])
    }
  }
  
  ## SMA called from Python by reticulate
  out = StructuredMaskedAutoencoder(dataset, group, comp_num, N, p, c(p, encoder_config, sum(comp_num)), c(sum(comp_num), decoder_config, p), train_epoch)
  output_embedding = out[[1]]
  recovered_dataset = out[[2]]
  for (j in 1:N){
    removed_common[[j]] = dataset[[j]] - t(recovered_dataset[[j]])
  }
  ## compute the score for each dataset from SMA
  index = 1
  for(i in 1 : K){
    for(j in 1 : N){
      if (j %in% group[[i]]){
        list_score[[j]][[i]] = t(output_embedding[[j]][, index: (index + comp_num[i] - 1)])
      }
    }
    index = index + comp_num[i]
  }
  return(list(score_list = list_score, recovered_dataset = recovered_dataset, removed_common = removed_common))
}
# split dataset
training_frac <- 0.75
train_val_split <- function(X, training_frac){
  train_size <- floor(training_frac * ncol(X))
  set.seed(1)
  train_idx <- sample(seq_len(ncol(X)), size = train_size)
  X_train <- X[, train_idx]
  X_val <- X[ , -train_idx]
  return(list(X_train = X_train, X_val = X_val))
}
# normalization
normalization <- function(data_list, train_mean){
  normalized_data_list <- lapply(data_list, function(x) {(x - train_mean)})
  return(normalized_data_list)
}
inVitro_sc_split <- train_val_split(inVitro_sc, training_frac)
train_mean <- colMeans(t(inVitro_sc_split$X_train))
dataset <- normalization(list(inVitro_sc_split$X_train, inVitro_sc_split$X_val, inVitro_sc, inVivo_sc), train_mean)
print("normalization complete")
vitro <- SMA(dataset = dataset, group = group, comp_num = comp_num, encoder_config = c(512), decoder_config = c(512), train_epoch = 2000)
print("Training complete")