library(terra) # raster manipulation
library(sits)  # classify eo cube time series
library(sf)    # vector manipulation
library(raster)
library(rgdal)
library(dplyr)

#devtools::install_github("rolfsimoes/sits@dev2")
devtools::install_github("oldlipe/sits@dev2")
devtools::install_github("e-sensing/sits@dev")
devtools::install_github("oldlipe/sits@feature/blocks_optimal")

#shp_dir <- "/data/RO_Teste"
shp_file <- "/data/Rondonia_231068/shp/BR_UF_2019.shp"
sf_shp <- st_read(shp_file)
plot(sf_shp)
rondonia <- subset(sf_shp, SIGLA_UF == "RO")
plot(rondonia)
#sf_C62L04 <- sf_shp[1,]

######Criar Cubo####
#Create data cube -62,553 -11,631
cube_L5_231068 <- sits_cube(
  source = "MPC",
  collection = "LANDSAT-C2-L2",
  roi = rondonia,
  #roi = c(
 # xmin = -62.2010,
 # ymin = -10.7784,
 # xmax = -64.075,
 # ymax = -10.984,
 # crs = 4326
 #),
  start_date = "1988-01-01",
  end_date = "1988-12-31",
  bands = c("BLUE", "GREEN","RED","NIR08","SWIR16","SWIR22","CLOUD")
)

#Alterar o nome do sensor
cube_L5_231068_reg <- dplyr::mutate(cube_L5_231068_reg, sensor = "TM-ETM")

saveRDS(cube_L5_231068_copy, "/data/cube_L5_231068_copy.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####Regulariza o cubo#####
# regularize a data cube
cube_L5_231068_reg <- sits_regularize(
  cube = cube_L5_231068,
  period = "P16D",
  res = 30,
  output_dir	= "/data/Rondonia_231068/reg/1989",
  multicores = 60,
)

saveRDS(cube_L5_231068_reg, "/data/RO_Teste/RO_231068/samples/cube_L5_231068_reg.rds")

#plot image
plot(cube_L5_231068_reg, red = "SWIR16", green = "NIR08", blue = "BLUE",
     date = c("1988-09-25"))
plot(cube_L5_231068_reg, red = "SWIR22", green = "NIR08", blue = "BLUE",
     date = c("1985-01-27"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####Cria índices no cubo regularizado#####
cube_L5_231068_reg <- sits_apply(
  data = cube_L5_231068_reg,
  NDVI = (NIR08 - RED)/(NIR08 + RED),
  memsize = 240,
  multicores = 60,
  output_dir = "/data/Rondonia_231068/reg/1989"
)

cube_L5_231068_reg <- sits_apply(
  data = cube_L5_231068_reg,
  EVI = 2.5 * ((NIR08 - RED) / (RED + 6 * GREEN - 7.5 * BLUE + 1)),
  memsize = 240,
  multicores = 60,
  output_dir = "/data/Rondonia_231068/reg/1989"
)


cube_L5_231068_reg <- sits_apply(
  data = cube_L5_231068_reg,
  NBR = (NIR08 - SWIR16) / (NIR08 + SWIR16),
  memsize = 240,
  multicores = 60,
  output_dir = "/data/Rondonia_231068/reg/1989"
)

cube_L5_231068_reg <- sits_apply(
  data = cube_L5_231068_reg,
  NDWI = (NIR08 - GREEN) / (NIR08 + GREEN),
  memsize = 240,
  multicores = 60,
  output_dir = "/data/Rondonia_231068/reg/1988"
)

#Amostras para o MLME
endmembers_spectra <- #LANDSAT5
  tibble::tibble(
    type = c("solo", "sombra", "veg"),
    RED = c(10048, 8730, 8125),
    BLUE = c(8854, 8148, 7803),
    NIR08 = c(14361, 8713, 16142),
    GREEN = c(9632, 8966, 8444),
    SWIR16 = c(16090, 7706, 11307),
    SWIR22 = c(12178, 7721, 8579)
  )
#cria o mlme no cubo
mixture_cube <- sits_mixture_model(
  data = cube_L5_231068_reg,
  endmembers = endmembers_spectra,
  rmse_band = FALSE,
  memsize = 120,
  multicores = 60,
  output_dir = "/data/Rondonia_231068/reg/1986",
  progress = TRUE
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####RECUPERA qualquer cubo criado nos praâmetros do SITS#####
#Recuperando o cubo
cube_L5_231068_reg <- sits_cube(
  source = "MPC",
  collection = "LANDSAT-C2-L2",
  data_dir = "/data/Rondonia_231068/reg/1988",
  parse_info = c("X1", "tile","band","date"),
  multicores = 10
)

#seleciona somente alguns tiles de todo o cubo
cube_L5_231068_reg_select <- sits_select(
  data = cube_L5_231068_reg,
  tiles = c("231067","231068","232067","232066")
)

#Datas das imagens disponíveis
sits_timeline(cube_L5_231068_reg)
#bandas disponíveis no cubo
sits_bands(cube_L5_231068_reg)

 #####Lê as anostras e mostra sua resposta ao longo do cubo####
#get data
samples_231068 <- sits_get_data(
  cube = cube_L5_231068_reg,
  samples = "/data/RO_Estado/tiles/samples/samples_1988.csv",
  #samples = "/data/Rondonia_231068/samples/1989/samples_231068_1989.csv",
  multicores = 40,
)

#samples_231068 <- saveRDS(cube_L5_231068_reg, "/data/RO_Teste/RO_231068/samples/samples_231068_1985-08-01_getData.rds")

#Mapa de SOM
som <- sits_som_map(
  data = samples_231068,
  grid_xdim = 5,
  grid_ydim = 5,
  alpha = 1,
  rlen = 100,
  distance = "euclidean",
  som_radius = 2,
  mode = "online"
)
plot(som)

#"Apara" as amostras por meio do SOM
new_samples <- sits_som_clean_samples(
  som_map = som,
  prior_threshold = 0.6,
  posterior_threshold = 0.6,
  keep = c("clean", "analyze")
)
# print the new sample distribution
sits_labels_summary(new_samples)

# Mapa de SOM com amostras "aparadas"
new_cluster <- sits_som_map(
  data = new_samples,
  grid_xdim = 4,
  grid_ydim = 4,
  alpha = 1.0,
  distance = "euclidean",
)
plot(new_cluster)

#Gráfico de confusão
patterns <- sits_patterns(samples_231068)
plot(patterns)

#plota as amostras
plot(sits_select(samples_231068, "NDVI"))
#sits_view(samples_231068)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####Processo de classificação####
# Train the samples with Random Forests model.
rfor_model <- sits_train(
  samples = samples_231068,
  ml_method = sits_rfor(num_trees = 100)
)

# plot the most important variables of the model
plot(rfor_model)

#Train the samples with TempCNN model.
tempcnn_model <- sits_train(
  samples = samples_231068,
  ml_method = sits_tempcnn(
    optimizer            = torchopt::optim_adamw,
    cnn_layers           = c(128, 128, 128),
    cnn_kernels          = c(7, 7, 7),
    cnn_dropout_rates    = c(0.2, 0.2, 0.2),
    epochs               = 100,
    validation_split     = 0.2,
    verbose = FALSE
  )
)
# plot the most important variables of the model
plot(tempcnn_model)

#Probabilidade de cada classe
samples_231068_probs <- sits_classify(
  data = cube_L5_231068_reg,
  ml_model = rfor_model,
  multicores = 120,
  memsize = 60,
  output_dir = "/data/Rondonia_231068/results/1988/classify/v4_mateus",
  version = 'v4_rfor_Mateus'
)
# plot the probability cube for class vs
plot(samples_231068_probs)

#v4_rfor_mateus: umida, cr_mineracao, cr_queimada, cr_solo, cr_veg, floresta,
#hidrografia, urbana

#Filtro smooth da propabilidade
samples_231068_bayes <- sits_smooth(
  cube = samples_231068_probs,
  window_size = 7,
  neigh_fraction = 0.5,
  smoothness = c(20, 20, 20, 20, 20, 20, 0, 20),
  #smoothness = 20,
  multicores = 64,
  memsize = 120,
  output_dir = "/data/Rondonia_231068/results/1988/bayes/v4_mateus/",
  version = 'v4_rfor_Mateus'
)
plot(samples_231068_bayes)

#Classificação em si
samples_231068_maps <- sits_label_classification(
  cube = samples_231068_bayes,
  multicores = 64,
  memsize = 120,
  output_dir = "/data/Rondonia_231068/results/1988/maps/v4_mateus",
  version = 'v4_rfor_Mateus'
)
plot(samples_231068_maps, title = "RO Classification Map")

#Seleciona as imagens de um cubo de dados
mosaic_rgb <- sits_select(
  cube_L5_231068_reg_select,
  bands = c('SWIR16', 'NIR08','BLUE'),
  start_date = '1988-09-25',
  end_date = '1988-09-25'
)
#datas de um cubo de imagens
sits_timeline(mosaic_rgb)

#Mosaico de imagens. Pode ser tanto para imagens de um cubo de dados/n
#quanto para mosaico de classificações (maps)./n
#O parametro roi serve para recortar a área de interesse
mosaic_cube <-sits_mosaic(
  cube= samples_231068_maps,
  #cube = mosaic_rgb,
  crs = 4326,
  #roi = rondonia,
  #roi = c(
    #xmin = -62,553,
    #ymin = -11,630,
    #xmax = -60,900,
    #ymax = -11,631,
    #crs = '4326'
  #),
  multicores = 50,
  output_dir = "/data/Rondonia_231068/results/1989/mosaic",
  version = "v1_rfor_email",
  progress = TRUE
)
# show the location of the classification file
samples_231068_maps$file_info[[1]]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####Processo de máscara do PRODES#####
# Open mask map
prodes2021 <- sits_cube(
  source = "USGS",
  collection = "LANDSAT-C2L2-SR",
  data_dir = "/data/PRODES",
  parse_info = c("X1","date","tile","band","version"),
  #bands = "class",
  labels = c(as.character(1:100), "NoForest")
)
prodes2021 <- dplyr::mutate(prodes2021, sensor = "ETM")

# Open classification map
ro_class <- sits_cube(
  source = "MPC",
  collection = "LANDSAT-C2-L2",
  data_dir = "/data/Rondonia_231068/results/1988/mosaic/v4_mateus",
  version = "v4",
  parse_info = c("X1", "X2", "tile", "start_date", "end_date","band", "version"),
  bands = "class",
  labels = c("umida","cr_mineracao","cr_queimada", "cr_solo", "cr_veg",
             "floresta", "hidrografia","area_urbana")
)
ro_class <- dplyr::mutate(ro_class, sensor = "ETM")
plot(ro_class, title = "RO Classification Map")

# Reclassify cube
ro_mask <- sits_reclassify(
  #cube = cropped_class
  cube = ro_class,
  mask = prodes2021,
  rules = list(
    "NonForest" = mask == "NoForest"
  ),
  memsize = 10,
  multicores = 1,
  output_dir = "/data/Rondonia_231068/results/1988/mosaic/v4_mateus/masked",
  version = "masked"
)

plot(ro_mask, palette = "Spectral")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####Precisão das amostras#####
val_rfor <- sits_kfold_validate(
  samples = samples_231068,
  folds = 5,
  ml_method = sits_rfor(),
  multicores = 5
)
# print the validation statistics
sits_accuracy_summary(val_rfor)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####Entropia para calcular os pontos mais incertos#####
bayes_uncert_entropy <- sits_uncertainty(
  cube = samples_231068_probs,
  type = "entropy",
  #window_size = 5,
  multicores = 60,
  memsize = 200,
  output_dir = "/data/Rondonia_231068/results/1988/uncert",
  version = "v0_rfor"
)

plot(bayes_uncert_entropy)

#Cria shp com os pontos mais incertos
new_samples_uncert <- sits_uncertainty_sampling(
  uncert_cube = bayes_uncert_entropy,
  n = 200,
  min_uncert = 0.7,
  sampling_window = 10
)

sits_view(new_samples_uncert)
sits_to_csv(new_samples_uncert, "/data/Rondonia_231068/results/1988/uncert/incertos_1988.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####Medição da Precisão#####
results <- list()
# Give a name to the results of the random forest model (see above)
val_rfor$name <- "rfor"
# store the rfor results in a list
results[[length(results) + 1]] <- val_rfor

# Temporal CNN
val_tcnn <- sits_kfold_validate(
  samples = samples_231068,
  ml_method = sits_tempcnn(
    optimizer = torchopt::optim_adamw,
    opt_hparams = list(lr = 0.001)
  ),
  folds = 5,
  multicores = 10
)

# Give a name to the result
val_tcnn$name <- "TempCNN"
# store the results in a list
results[[length(results) + 1]] <- val_tcnn

# Save to an XLS file
#xlsx_file <- "data/RO_Teste/RO_231068/results/class_prodes/1985-1986/model_comparison.xlsx"
sits_to_xlsx(results, file = "/data/Rondonia_231068/samples/1986/rfor.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####Melhorando as amostras#####
clusters <- sits_cluster_dendro(
  samples = samples_231068,
  bands = c("NDVI", "EVI"),
  dist_method = "dtw_basic",
  linkage = "ward.D2"
)

sits_cluster_frequency(clusters)

clusters_new <- dplyr::filter(clusters, cluster != 2)
clean <- sits_cluster_clean(clusters_new)
sits_cluster_frequency(clean)

#####Validação da classificação
samples_231068_maps <- sits_cube(
  source = "MPC",
  collection = "LANDSAT-C2-L2",
  data_dir = "/data/Rondonia_231068/results/1988/maps/v4_mateus",
  parse_info = c("X1", "X2", "tile","start_date","end_date","band","version","X3","X4"),
  bands = "class",
  labels = c("Area_Umida","CR_Mineracao","CR_Queimada","CR_Solo","CR_Veg","Floresta","Hidrografia","Urbano"),
  version = "v4",
  multicores = 10
)
# get ground truth points
valid_csv <- "/data/Rondonia_231068/samples/1988/samples_1988_valid_mateus.csv"
# calculate accuracy according to Olofsson's method
area_acc <- sits_accuracy(data = samples_231068_maps,
                          validation = valid_csv)
# print the area estimated accuracy
area_acc
#write.table(area_acc, "/data/Rondonia_231068/samples/1988/accuracy_classification_1988_mateus.csv",sep = "\t", row.names = FALSE)



# retrieve the metadata for the classified cube
# the files are stored as Dropbox links
cerrado_classif_rds <- system.file("extdata/Cerrado/cerrado_classif_dropbox.rds",
                                   package = "sitsdata"
)
# read the cube metadata
cerrado_classif <- readRDS(cerrado_classif_rds)
# plot one tile of the classification
plot(cerrado_classif, tile = "044048")

valid_csv <- system.file("extdata/csv/cerrado_lc8_validation.csv",
                         package = "sitsdata"
)

area_acc <- sits_accuracy(cerrado_classif,
                          validation_csv = valid_csv
)
