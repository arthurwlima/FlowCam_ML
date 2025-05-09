################################################################################
####################### Arthur Lima 2025.01.13 #################################
# Conjunto de rotinas pra processar os resultados da FlowCam                   #
# 1 - Calcular a contagem de particulas por classe e atributos medios de cada  # 
# classe dentro de cada amostra                                                #
# 2 - Concatenar as contagens de cada amostra em uma planilha unica por projeto#
# 3 - Classificar particulas em novas amostras a partir dos filtros de bibliote# 
#cas especificas 
################################################################################
library(MASS)

require(ggplot2)
require(dplyr)
require(patchwork)
require(tidyr)
library(e1071)
library(stringi)

setwd("/home/arthurw/Documents/FitoPK/FlowCam")

################################################################################
# 1 - Calcular a contagem de particulas por classe e atributos medios de cada  # 
################################################################################
{
rm(list=ls())
listArq <- list.files("./", pattern="*_data.csv", full.names = F, recursive=T) 
cont=0
print(paste("Total de arquivos: ", length(listArq)))
for(i in listArq){
  cont= cont+1
  print(paste(cont, as.character(i)))
  d <- read.table(i, header=T, sep=',')
  names(d)
  dd <- d %>% 
    group_by(Class) %>% 
    dplyr::summarise(across(c("Circularity", "Circularity..Hu.", "Diameter..ESD.",
                            "Average.Blue", "Average.Green", "Average.Red",
                            "Length", "Roughness", "Transparency", "Volume..ESD.", "Width"), 
                          ~ mean(., na.rm = TRUE)))
  dd[,2:length(names(dd))] <- round(dd[,2:length(names(dd))], 2)
  dd$contagem <- table(d$Class)
  
  sample_id0 <- strsplit(i, "/")[[1]][2]
  sample_id <- strsplit(sample_id0, "_")[[1]][1]
  dd$sample_id <- sample_id

  write.table(dd, paste(sample_id, "_resumo.csv", sep=''), row.names=F, sep=',')
}}

################################################################################
# 2 - Concatenar os resumos de cada amostra em uma planilha unica por categoria#
################################################################################
{
rm(list=ls())
listArq <- list.files("./", pattern="*_resumo.csv", full.names = F, recursive=T) 

my_data <- list()
for (i in listArq) {
  my_data[[i]] <- read.table(file = i, header=T, sep=',')
}


dd <- Reduce(function(...) merge(..., all= TRUE), my_data)

names(dd)
dd.1 <- select(dd, c("Class", "sample_id", "contagem"))

dd.1$Class = stri_trans_general(str = dd.1$Class, id = "Latin-ASCII")

dd_wide <- pivot_wider(dd.1, names_from = "Class", values_from = "contagem", values_fill = 0)

names(dd_wide)
write.table(dd_wide, "Contagem_Clases_Amostra.csv", row.names=F, sep=',')
}
####################### 20250213 #######################################
# Arthur Lima
# Script para classificar amostras de fitoplancton coletadas pela flowcam
# 1 - Organizar as bibliotecas de microalgas de interesse
# 2 (?) - Criar uma regra de classificacao para identificar as microalgas (DCA)
# 3 - Classificar novas particulas com algoritmo de KNN
########################################################################
# Filtro de chattonela no arquivo, valores separados por |
# /home/arthurw/Documents/FitoPK/FlowCam/Chattonella subsalsa 40 dias- CCMR0025-x20.lst
# Primeiras 59 linhas com definicoes do dicionario de dados
# Cabecalho insluido no arquivo
# /home/arthurw/Documents/FitoPK/FlowCam/Chattonella subsalsa 40 dias- CCMR0025-x20.header.lst
########################################################################
# 2025-02-28
# Filtro de chattonela com valores bem diferentes dos observados nas amostras
# Criar um 'novo filtro' a partir dos atributos de cada grupo em amostras classificadas
########################################################################
# 2025-05-09
# Classificacao por SVM funcionando bem para a mesma amostra, mas sem funcionar para uma segunda amostra
# testar o treinamento do modelo com uma combinacao maior de amostras

d.samp <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20240221_data.csv", header=T, sep=',')
d.samp2 <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20230906_data.csv", header=T, sep=',')
d.samp$Class = stri_trans_general(str = d.samp$Class, id = "Latin-ASCII")
d.samp2$Class = stri_trans_general(str = d.samp2$Class, id = "Latin-ASCII")



d.header <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/Chattonella subsalsa 40 dias- CCMR0025-x20.header.2.lst", header=T, sep='|')

d.samp$chaetoceros <- 0
d.samp$chaetoceros [d.samp$Class =='Chaetoceros' ] <- 1

d.samp2$chaetoceros <- 0
d.samp2$chaetoceros [d.samp2$Class =='Chaetoceros' ] <- 1

d.header$chaetoceros <- ''

intersect(names(d.header), names(d.samp))
d.samp <- select(d.samp, intersect(names(d.header), names(d.samp)))
d.samp2 <- select(d.samp2, intersect(names(d.header), names(d.samp2)))

#dd <- merge(d.chat, d.samp, all=T)
#names(dd)



# Reduz o tamanho da amostra
#  n = 10000
#  amostra <- sample(1:nrow(dd), size = n)
#  dd2 <- dd[amostra, ] 
#
  
fit = svm(factor(d.samp$chaetoceros) ~ ., data = d.samp[,2:23], scale = T, kernel = "radial")
y_hat = predict(fit, newdata=d.samp2[,2:23])
table(y_hat)
table(d.samp2$chaetoceros)
tab <- table(y_hat, d.samp2$chaetoceros)
caret::confusionMatrix(tab)
"Confusion Matrix and Statistics


y_hat     0     1
0 42676     0
1     0   579

Accuracy : 1          
95% CI : (0.9999, 1)
No Information Rate : 0.9866     
P-Value [Acc > NIR] : < 2.2e-16 "

################################################################################
# euglena
d.samp <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20240221_data.csv", header=T, sep=',')
d.samp2 <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20230906_data.csv", header=T, sep=',')
d.samp$Class = stri_trans_general(str = d.samp$Class, id = "Latin-ASCII")
d.samp2$Class = stri_trans_general(str = d.samp2$Class, id = "Latin-ASCII")

d.samp$euglena <- 0
d.samp$euglena [d.samp$Class =='Euglena' ] <- 1
table(d.samp$euglena)
d.samp2$euglena <- 0
d.samp2$euglena [d.samp2$Class =='Euglena' ] <- 1
table(d.samp2$euglena)

d.header$euglena <- ''

intersect(names(d.header), names(d.samp))
d.samp <- select(d.samp, intersect(names(d.header), names(d.samp)))
d.samp2 <- select(d.samp2, intersect(names(d.header), names(d.samp2)))

#dd <- merge(d.chat, d.samp, all=T)
#names(dd)



# Reduz o tamanho da amostra
#  n = 10000
#  amostra <- sample(1:nrow(dd), size = n)
#  dd2 <- dd[amostra, ] 
#

fit = svm(factor(d.samp$euglena) ~ ., data = d.samp[,2:23], scale = T, kernel = "radial", cost=100000)
y_hat = predict(fit, newdata=d.samp2[,2:23])
table(y_hat)
table(d.samp2$euglena)
tab <- table(y_hat, d.samp2$euglena)
caret::confusionMatrix(tab)


################################################################################
# Ciano e fito filamentoso
d.samp <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20240221_data.csv", header=T, sep=',')
d.samp2 <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20230906_data.csv", header=T, sep=',')
d.samp$Class = stri_trans_general(str = d.samp$Class, id = "Latin-ASCII")
d.samp2$Class = stri_trans_general(str = d.samp2$Class, id = "Latin-ASCII")

d.samp$euglena <- 0
d.samp$euglena [d.samp$Class =='Cianobact�rias e fito filamentoso' ] <- 1
table(d.samp$euglena)
d.samp2$euglena <- 0
d.samp2$euglena [d.samp2$Class =='Cianobact�rias e fito filamentoso' ] <- 1
table(d.samp2$euglena)

d.header$euglena <- ''

intersect(names(d.header), names(d.samp))
d.samp <- select(d.samp, intersect(names(d.header), names(d.samp)))
d.samp2 <- select(d.samp2, intersect(names(d.header), names(d.samp2)))


fit = svm(factor(d.samp2$euglena) ~ ., data = d.samp2[,2:23], scale = T, kernel = "radial", cost = 100000)
y_hat = predict(fit, newdata=d.samp[,2:23])
table(y_hat)
table(d.samp2$euglena)
tab <- table(y_hat, d.samp$euglena)
caret::confusionMatrix(tab)

################################################################################
# Cylindrotheca
d.samp <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20240221_data.csv", header=T, sep=',')
d.samp2 <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20230906_data.csv", header=T, sep=',')
d.samp$Class = stri_trans_general(str = d.samp$Class, id = "Latin-ASCII")
d.samp2$Class = stri_trans_general(str = d.samp2$Class, id = "Latin-ASCII")

table(d.samp2$Class)
d.samp$euglena <- 0
d.samp$euglena [d.samp$Class =='Cylindroteca' ] <- 1
table(d.samp$euglena)
d.samp2$euglena <- 0
d.samp2$euglena [d.samp2$Class =='Cylindrotheca' ] <- 1
table(d.samp2$euglena)

d.header$euglena <- ''

intersect(names(d.header), names(d.samp2))
d.samp <- select(d.samp, intersect(names(d.header), names(d.samp)))
d.samp2 <- select(d.samp2, intersect(names(d.header), names(d.samp2)))
d.samp_comp <- rbind(d.samp, d.samp2)

fit = svm(factor(d.samp_comp$euglena) ~ ., data = d.samp_comp[,2:23], scale = T, kernel = "radial", cost = 100000)
y_hat = predict(fit, newdata=d.samp2[,2:23])
table(y_hat)
table(d.samp2$euglena)
tab <- table(y_hat, d.samp2$euglena)
caret::confusionMatrix(tab)
############################################################################
###
d.samp <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20240221_data.csv", header=T, sep=',')
d.samp2 <- read.table("/home/arthurw/Documents/FitoPK/FlowCam/RawData/PIER HANGAR 20231001_data.csv", header=T, sep=',')
d.samp$Class = stri_trans_general(str = d.samp$Class, id = "Latin-ASCII")
d.samp2$Class = stri_trans_general(str = d.samp2$Class, id = "Latin-ASCII")

d.samp_comp <- rbind(d.samp, d.samp2)
d.header$Class <- ''
intersect(names(d.header), names(d.samp_comp))
d.samp_comp <- dplyr::select(d.samp_comp, intersect(names(d.header), names(d.samp_comp)))

# Reduz o tamanho da amostra
n = 10000
amostra <- sample(1:nrow(d.samp_comp), size = n)
dd2 <- d.samp_comp[amostra, ] 
dd2[,2:23] <- scale(dd2[,2:23], center=T)
?scale
flow.lda <- lda(Class ~ ., data= dd2[,2:24])

svg("flow.lda.svg", 10, 10)
plot(flow.lda)
dev.off()


summary(flow.lda)
ld1 <- flow.lda$scaling[,1]
ld1[order(abs(ld1), decreasing=T)]

ld2 <- flow.lda$scaling[,2]
ld2[order(abs(ld2), decreasing=T)]

ggplot(dd2, aes(x=Intensity, y=Diameter..ESD., col=Class))+
  geom_point(alpha=0.7)
