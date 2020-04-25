rm(list = ls())

require(data.table)
require(caret)
require(doParallel)
require(lars)
require(DMwR)
# require(earth)
# require(pscl)
# require(stats)

######## Define Functions ########
AutoTune <- function(training, mtry_range, min_node_range = 2:8) {
  
  # Set Control
  control <- trainControl(method = 'cv', number = 10, classProbs = TRUE)
  
  ### Random Forest
  grid_ranger <- expand.grid(mtry = mtry_range, splitrule = 'gini', min.node.size = min_node_range)
  model_ranger <- train(make.names(count) ~ ., data = training, method = 'ranger', trControl = control, tuneGrid = grid_ranger)
  
  return(model_ranger)
}

WillSample <- function(count_data, nonCount_data, ZipRange = 30, nTime = 1, poissonModifier = 1) {
  newdata_nonLoc <- setNames(data.frame(matrix(ncol = ncol(nonCount_data), nrow = 0)), colnames(nonCount_data))
  newdata_nonLoc <- data.table(newdata_nonLoc)
  
  # For each Zip location that has count, randomly find nTime*poissonModifier number of non_count locations
  for (i in 1:nrow(count_data)) {
    zip_code <- as.numeric(count_data[i, 'zip'])
    temp_zipRange <- zip_code-ZipRange:zip_code+ZipRange
    temp_subset <- nonCount_data[(zip %in% temp_zipRange) & !(zip %in% newdata_nonLoc$zip)]
    sampled_nonLoc <- temp_subset[sample(1:nrow(temp_subset), nTime*poissonModifier)]
    
    newdata_nonLoc <- rbind(newdata_nonLoc, sampled_nonLoc)
  }
  
  # Smote Resample
  if (nTime == 1) {
    new_data <- count_data[, -1]
  } else {
    temp_all <- rbind(count_data[,-1], nonCount_data[,-1])
    new_data <- SMOTE(count~., data = temp_all, perc.over = (nTime-1)*100, perc.under = 0)
  }
  
  # Re-assemble
  # new_data <- data.table()
  
  # for (i in 1:nTime) {
  #   new_data <- rbind(new_data, count_data)
  # }
  new_data <- rbind(new_data, newdata_nonLoc[,-1])
  
  # Return WillSampled Data
  return(new_data)
}

# FindFeature = function(train_data, threshold = 10^(-5)){
#   X = as.matrix(train_data[ , -"count"])
#   Y = as.matrix(train_data[ , "count"])
#   lar = lars(X, Y, type = "lasso")
#   stepNum = unname(which.min(lar$Cp))
#   coef = coef.lars(lar,mode="step",s=stepNum)
#   coef = coef[abs(coef)>threshold]
#   return(coef)
# }

####### Loading Data Directory ----
data_folder <- './Mosaic Household/'
data_dir <- dir(data_folder)
roc_zip <- c(14602:14627, 14638, 14639, 14642:14644, 14646, 14647, 14649:14653, 14692, 14694)

# Set Seed 323
set.seed(323)

# Set Time Upsample
# In case for need of upsampling minority class (usually datapoint where count of store > 0).
nTime = 1
# Inc ase want to maintain a poisson balanced dataset, 
# if believe poisson distribution of count of stores would work better in model.
poissonModifier = 1

# Creat a empty place holder to store result, contains Rochester Area Zip code that are included in Mosaic Data
result_pred <- data.table(zip_code = c(14602L, 14603L, 14604L, 14605L, 14606L, 14607L, 14608L, 14609L, 
                                       14610L, 14611L, 14612L, 14613L, 14614L, 14615L, 14616L, 14617L, 
                                       14618L, 14619L, 14620L, 14621L, 14622L, 14623L, 14624L, 14625L, 
                                       14626L, 14627L, 14692L))
result_metrics <- data.table()
result_matrix <- list()

# Start Parallel
cl <- makeCluster(3)
registerDoParallel(cl)

####### Pipline Process 
for (i in 1:10) {
  # Loading each data file
  data_file <- data_dir[i]
  data_add <- paste(data_folder, data_file, sep = '')
  data <- fread(data_add, drop = 'V1')
  data$count <- (data$count > 0)*1
  
  # Rearrange Data to Model Format
  data <- data[,c(1, 73, 2:72)]
  # data <- data[,c(1, 142, 2:141)]
  roc_data <- data[zip %in% roc_zip]
  data <- data[!(zip %in% roc_zip)]
  
  # model_data <- data[,2:ncol(joint_data)]
  model_data <- na.exclude(data)
  
  
  ################### Data Reampling ###################
  nLoc <- nrow(model_data[count > 0])
  Loc <- model_data[count > 0]
  nonLoc <- model_data[count == 0]
  Loc$count <- factor(Loc$count)
  nonLoc$count <- factor(nonLoc$count)
  
  ZipRange = 40
  training <- WillSample(Loc, nonLoc, ZipRange = ZipRange, nTime = nTime, poissonModifier = poissonModifier)
  training <- na.exclude(training)
  print(paste('Finished Job Sequence - Data Resampling for ', data_file))
  
  ################### Feature Selection ###################
  # feat_Cols <- names(FindFeature(training))
  # nFeature <- length(feat_Cols)
  # new_training <- data.frame(training)
  # new_training <- new_training[,c('count', feat_Cols)]
  # print(paste('Finished Job Sequence - Feature Selection for ', data_file))
  # 
  ################### Model Tuning ###################
  nFeature <- length(colnames(training))
  model <- AutoTune(training, mtry_range = ceiling(nFeature/10):ceiling(nFeature/6))
  
  print(paste('Finished Job Sequence - Model Tuning for ', data_file))
  
  ## Getting Predictions
  pred <- predict(model, roc_data, type = "prob")
  model_result <- data.table(zip_code = as.numeric(roc_data$zip), round(pred[,'X1'],4))
  colnames(model_result)[2] <- data_file
  
  ## Getting Metrics
  report_metrics <- data.table(getTrainPerf(model))
  report_metrics$data <- data_file
  report_metrics <- report_metrics[,c(4,1:3)]
  
  ## Getting Confusion Matrix
  confuse_table <- confusionMatrix(model)$table
  
  result_pred <- merge(result_pred, model_result, by = 'zip_code')
  result_metrics <- rbind(result_metrics, report_metrics)
  result_matrix[[i]] <- confuse_table
  
  print(paste('Finished ', data_file))
}

# Shut down Parallel
stopCluster(cl)

find_NA <- function(df) {
  new_df <- df[rowSums(is.na(df)) > 0,]
  return(new_df)
}
