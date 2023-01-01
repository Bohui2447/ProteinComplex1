
### R version 3.5.1
## devtools::install_github("paodan/funcTools")
library(funcTools)
library(keras) # keras version 2.2.5.0
library(dplyr)
library(ggplot2)
library(purrr)
library(ROCR)

### loading data
# First Download the data and save them in the current working directory
load("./TrainingDataset.Rdata")
load("./ValidationDataset.Rdata") 

#### Hyperparameter searching (this will take a lot of time)
# In this study, instead of using the for loop, we submitted slurm batch jobs 
# to the HPC of Utrecht University to dramatically reduce the training time.
set.seed(1234)
seeds = sample(1:1000000, 40, replace = F)
for(seed in seeds){
  set.seed(seed) # 186720
  N =  50
  params = list(dropOut1 = round(runif(N, 0, 0.5), digits = 3),
                dropOut2 = round(runif(N, 0, 0.5), digits = 3),
                dropOut3 = round(runif(N, 0, 0.5), digits = 3),
                units1 = sample(50*(4:15), size = N, replace = T),
                units2 = sample(10*(5:20), size = N, replace = T),
                units3 = sample(5*(2:10), size = N, replace = T))
  
  for(ni in 1:N){
    print(ni)
    {
      model = keras_model_sequential()  %>%
        layer_dense(units = params$units1[ni], input_shape = dim(trainNorm)[2], activation = "relu") %>% # 512
        layer_dropout(params$dropOut1[ni]) %>%
        layer_dense(units = params$units2[ni], activation = "relu") %>%
        layer_dropout(params$dropOut2[ni]) %>%
        layer_dense(units = params$units3[ni], activation = "relu") %>% 
        layer_dropout(params$dropOut3[ni]) %>%
        layer_dense(units = 1, activation = "sigmoid") %>% 
        compile(
          optimizer = optimizer_rmsprop(lr = 0.002),
          loss = "binary_crossentropy",
          metrics = list('accuracy')
        )
      model  %>%  summary(model)
      
      # training the model
      set.seed(123)
      history = model %>% fit(
        trainNorm, as.numeric(trainlables),
        epochs = 50,
        verbose = 1,
        batch_size = 2048,
        validation_data = list(valNorm, as.numeric(vallabels))
      )
    }
  }
}


### The paramaters for the best model used in our paper
params = list(dropOut1 = 0.438,
              dropOut2 = 0.214,
              dropOut3 = 0.037,
              units1 = 350,
              units2 = 140,
              units3 = 25)

t0 = system.time({
  model = keras_model_sequential()  %>%
    layer_dense(units = params$units1, input_shape = dim(trainNorm)[2], activation = "relu") %>% 
    layer_dropout(params$dropOut1) %>%
    layer_dense(units = params$units2, activation = "relu") %>%
    layer_dropout(params$dropOut2) %>%
    layer_dense(units = params$units3, activation = "relu") %>% 
    layer_dropout(params$dropOut3) %>%
    layer_dense(units = 1, activation = "sigmoid") %>% 
    compile(
      optimizer = optimizer_rmsprop(lr = 0.002),
      loss = "binary_crossentropy",
      metrics = list('accuracy')
    )
  model  %>%  summary(model)
  
  # training the model
  history = model %>% fit(
    trainNorm, as.numeric(trainlables),
    epochs = 50,
    verbose = 0,
    batch_size = 2048,
    validation_data = list(valNorm, as.numeric(vallabels))
  )
})


### Prediction and evaluation
# Loading the test dataset
load('./TestDataset.Rdata')

# Loading the best model
model = load_model_hdf5('TheBestModel.h5')


# model prediction
t1 = system.time({
  pre_prob = model %>% predict(testNorm_combine) %>% as.vector()
  pred <- prediction(predictions = pre_prob, labels = as.factor(test_combineLables))
})
cat("Time of training model:", t0[1], "\nTime of prediction:", t1[1])


# Function to calculate F1 score
F1 = function(precision, recall){
  2*precision*recall/(precision+recall)
}

# Function to plot precision-recall curve
plot.prc = function(prc.perf, main = NULL, xlim = c(0, 1), ylim = c(0,1), 
                    plotF1 = TRUE, nudge_x = 0.02){
  data_prc = data.frame(prc.perf@x.values, prc.perf@y.values)
  colnames(data_prc) = c(prc.perf@x.name, prc.perf@y.name)
  if(is.na(data_prc[1,1])){
    data_prc[1,1] = 0
  } 
  if(is.na(data_prc[1,2])){
    data_prc[1,2] = 1
  } 
  data_prc$F1 = F1(data_prc[,1], data_prc[,2])
  
  id = which.max(data_prc$F1)[1]
  data_F1 = data_prc[id,]
  g = ggplot(data_prc, aes_string(x = prc.perf@x.name, y = prc.perf@y.name))+
    geom_line()+
    geom_point(data = data_F1)+
    coord_cartesian(xlim = xlim, ylim = ylim)+
    ggtitle(main)+
    theme_Publication()
  if(plotF1) g = g +
    geom_text(aes(label = paste0("F1=", round(F1, 3))), data = data_F1, 
              hjust = "outward", nudge_x = nudge_x, color = "brown")
  return(g)
}

# precision, recall, and F1 score
prec_rec = performance(pred, "prec", "rec")
prec_rec_plot = plot.prc(prec_rec)
prec_rec_plot
F1_value = max(prec_rec_plot$data$F1)
F1_value
