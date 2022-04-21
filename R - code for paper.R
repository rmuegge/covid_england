library(reshape2)
library(rgdal)
library(leaflet)
library(dplyr)
library(xtable)
library(ggplot2)
library(cluster)
library(factoextra)
library(gridExtra)
library(scales) 
library(RColorBrewer)
library(htmltools)
library(tidyverse)
library(coda)
library(ggmcmc)
library(ggpubr)
library(fossil)
library(spdep)

#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load in the datafile
covid <- read.csv("coviddata.csv")

###################
#Fitting the model#
###################
#Note: the file coviddata.csv already contains the estimated risks we obtained from fitting the model below
#The code is included solely for reproducibility
#W is the neighbourhood matrix defined in Section 2.2. (Methods)
W <- readRDS("W_nbhd_matrix.rds")

#Note: fitting the model takes several hours (~8hrs) to run
fit1 <- ST.CARar(deaths~offset(log(expecteddeaths)),family="poisson",
                 data=covid.ordered,W=W,AR=1,MALA=TRUE,burnin=200000,n.sample=2200000,thin=1000,verbose=TRUE)

#Note: fitting the model takes several hours (~8hrs) to run
fit2 <- ST.CARar(deaths~offset(log(expecteddeaths)),family="poisson",
                 data=covid.ordered,W=W,AR=2,MALA=TRUE,burnin=200000,n.sample=2200000,thin=1000,verbose=TRUE)
covid$fit_ar2 <- fit2$`fitted.values`
covid$est_risk <- covid$fit_ar2/covid$expecteddeaths

#Results Section
##################
####lockdown 1####
##################
#weeks after lockdown 1 (starts 20_14)
#include 1 week prior as baseline level, and up to the 13th week after intro of lockdown
w_vec <- c(paste0("20_",13:26)) 
#actual weeks of lockdown:
l_int <- c(paste0("20_",14:20))

#weeks_since_lockdown is 0 for last week before lockdown, then 1,2,3,...
w_vec <- as.data.frame(w_vec)
w_vec$weeks_since_ld <- 0:(nrow(w_vec)-1)
colnames(w_vec) <- c("week","weeks_since_ld")

#obtain baseline risk (from last week before lockdown) by area
baseline_risk <- covid[which(covid$week==w_vec$week[1]),]
baseline_risk <- baseline_risk[,c("area","est_risk")]
colnames(baseline_risk) <- c("area","baseline_risk")

#merge weeks after lockdown with other data for that time
wal <- covid[which(covid$week %in% w_vec$week),]
wal <- merge(w_vec,wal)
wal$weeks_since_ld <- as.factor(wal$weeks_since_ld)

#merge data with baseline_risk
wal <- merge(wal,baseline_risk)
#scaled: estimated/baseline
wal$sc_risk <- wal$est_risk/wal$baseline_risk
wal <- wal[order(wal$area,wal$weeks_since_ld),]
#add abbreviations (used for plot legend)
abbr <- as.data.frame(cbind(c("East Midlands","East of England","London","North East",
                "North West","South East","South West","West Midlands",
                "Yorkshire and The Humber"),
              c("EM","EE","L","NE","NW","SE","SW","WM","YH")))
colnames(abbr) <- c("region","abbr")
wal <- merge(wal,abbr)
wal$abbr <- as.factor(wal$abbr)

#Figure 2(a)
bp_est_ld1 <- ggplot(data=wal,aes(weeks_since_ld,est_risk)) + geom_boxplot() +
  ylab("Estimated risk") + xlab("") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept=1,linetype="dashed",colour="red") + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + ggtitle("(a) first lockdown, 26/03/20-12/05/20") +
  geom_hline(yintercept=median(baseline_risk$baseline_risk),linetype="dashed",colour="blue") +
  geom_boxplot()

#week vs estimated risk by region
reg_wal <- wal 
#population-based weighted averages by region
reg_wal$weighted_risk <- reg_wal$est_risk*reg_wal$pop
reg_wal <- reg_wal[,c("abbr","weighted_risk","pop","week")]
reg_wal <- aggregate(reg_wal[,c("weighted_risk","pop")],by=list(reg_wal$abbr,reg_wal$week),FUN=sum)
colnames(reg_wal) <- c("abbr","week","weighted_risk","pop")
reg_wal$weighted_risk <- reg_wal$weighted_risk/reg_wal$pop
#add weeks_since_ld to these weighted averages
reg_wal <- merge(reg_wal,w_vec)
reg_wal$weeks_since_ld <- as.factor(reg_wal$weeks_since_ld)

#choose color palette
myColors <- brewer.pal(9,"Set1")
names(myColors) <- levels(reg_wal$abbr)

#Figure 3(a) left with legend
line_reg_ld1 <- ggplot() + geom_line(data=reg_wal,aes(weeks_since_ld,weighted_risk,group=abbr,colour=abbr)) + xlab("") + 
  ylab("Estimated risk") +
  labs(fill="Region", colour="Region") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept = 1, linetype="dashed",colour="red") +
  geom_line(data=reg_wal,aes(weeks_since_ld,weighted_risk,group=abbr,colour=abbr)) + ggtitle("(a) First lockdown, 26/03/20-12/05/20") +
  scale_colour_manual(name="Region",values=myColors,limits=c("EE","EM","L","NE","NW","SE","SW","WM","YH")) + theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.key.size = unit(3,"line")) + ylim(0,10) +
  guides(color = guide_legend(override.aes = list(size = 2)))

#extract legend to add to the combined plots later
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#extract legend for the regional risk plots
legend_reg_plot <- get_legend(line_reg_ld1)

#Figure 3(a) left without legend (as in paper)
line_reg_ld1 <- ggplot() + geom_line(data=reg_wal,aes(weeks_since_ld,weighted_risk,group=abbr,colour=abbr)) + xlab("") + 
  ylab("Estimated risk") +
  labs(fill="Region", colour="Region") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept = 1, linetype="dashed",colour="red") +
  geom_line(data=reg_wal,aes(weeks_since_ld,weighted_risk,group=abbr,colour=abbr)) + ggtitle("(a) First lockdown, 26/03/20-12/05/20") +
  scale_colour_manual(name="Region",values=myColors,limits=c("EE","EM","L","NE","NW","SE","SW","WM","YH")) + theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.key.size = unit(3,"line"),legend.position = "none") + ylim(0,10) +
  guides(color = guide_legend(override.aes = list(size = 2)))

#week vs scaled estimated risk by region
reg_wal_sc <- wal 
#compute population-based weighted averages 
reg_wal_sc$weighted_sc_risk <- reg_wal_sc$sc_risk*reg_wal_sc$pop
reg_wal_sc <- reg_wal_sc[,c("abbr","weighted_sc_risk","pop","week")]
reg_wal_sc <- aggregate(reg_wal_sc[,c("weighted_sc_risk","pop")],by=list(reg_wal_sc$abbr,reg_wal_sc$week),FUN=sum)
colnames(reg_wal_sc) <- c("abbr","week","weighted_sc_risk","pop")
reg_wal_sc$weighted_sc_risk <- reg_wal_sc$weighted_sc_risk/reg_wal_sc$pop
#add weeks_since_ld to these weighted scaled risks
reg_wal_sc <- merge(reg_wal_sc,w_vec)
reg_wal_sc$weeks_since_ld <- as.factor(reg_wal_sc$weeks_since_ld)

#Figure 3(a) right
line_reg_sc_ld1 <- ggplot() + geom_line(data=reg_wal_sc,aes(weeks_since_ld,weighted_sc_risk,group=abbr,colour=abbr)) + xlab("") + 
    ylab("Scaled estimated risk") +
    labs(fill="Region", colour="Region") +
    annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
    geom_hline(yintercept = 1, linetype="dashed",colour="blue") +
    geom_line(data=reg_wal_sc,aes(weeks_since_ld,weighted_sc_risk,group=abbr,colour=abbr)) + ggtitle("") +
    scale_colour_manual(name="Region",values=myColors,limits=c("EE","EM","L","NE","NW","SE","SW","WM","YH")) + theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.key.size = unit(3,"line"),legend.position = "none") + ylim(0,10) +
    guides(color = guide_legend(override.aes = list(size = 2)))

# Determine number of clusters
LA_risk <- wal[,c("area","week","est_risk")]
LA_risk <- LA_risk[order(LA_risk$area,LA_risk$week),]
LA_risk <- LA_risk[which(LA_risk$week %in% w_vec$week[2:nrow(w_vec)]),]
#from long to wide format
LA_risk <- reshape(LA_risk, idvar = "area", timevar = "week", direction = "wide")

#only risks for average silh. widths and wss
risks <- LA_risk[,2:ncol(LA_risk)]

#function to compute the average silhouette width
avg_sil <- function(k) {
  km.res <- kmeans(risks, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(risks))
  mean(ss[, 3])
}
# Compute and plot within-cluster sum of squares for k = 2 to k = 10
k.values <- 2:10
# extract avg silhouette for 2-10 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)
avg_sil_values <- as.data.frame(cbind(k.values,avg_sil_values))

#compute the within-cluster sum of squares
wcss <- (nrow(LA_risk)-1)*sum(apply(risks,2,var))
for(i in 2:10){
  wcss[i] <- sum(kmeans(risks,centers=i,nstart=25)$withinss)
} 
wcss <- as.data.frame(cbind(1:10,wcss))
colnames(wcss) <- c("no_clusters","wcss")

#create the average silhouette width plot
sil_ld1 <- ggplot() + geom_line(data=avg_sil_values,aes(x=k.values,y=avg_sil_values,group=1)) + geom_point(data=avg_sil_values,aes(x=k.values,y=avg_sil_values,group=1),colour="#003865") +
  scale_x_continuous(breaks=2:10) + theme(panel.grid.minor.x = element_blank()) + xlab("") + ylab("Average Silhouette Width") + ggtitle("(a) Lockdown 1") + ylim(0,0.6) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),plot.title = element_text(size=16))

#create the within-cluster sum of squares plot
wcss_ld1 <- ggplot() + geom_line(data=wcss,aes(x=no_clusters,y=wcss,group=1)) + geom_point(data=wcss,aes(x=no_clusters,y=wcss,group=1),colour="#003865") +
  scale_x_continuous(breaks=1:10) + theme(panel.grid.minor.x = element_blank()) + xlab("") + ylab("Within-Cluster Sum of Squares") + ggtitle("") + ylim(0,6005) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),plot.title = element_text(size=16))

#plot the average silhouette width and wcss (see the supplementary materials)
grid.arrange(sil_ld1,wcss_ld1,nrow=1)

#output above suggests 2 clusters
#apply k-means algorithm with k=2 clusters from 25 starting points
clust <- kmeans(risks,2,nstart=25)
LA_risk$cluster <- clust$cluster
covid_risk <- LA_risk[,c("area","cluster")]
covid_risk <- merge(covid_risk,wal)
#repeat the 4 lines above if cluster 1 is not the smaller cluster
table(covid_risk$cluster)
#rename the cluster groups so they indicate the number of LADs in that cluster
covid_risk$cluster[which(covid_risk$cluster==1)] <- c(paste0("1 (",table(covid_risk$cluster)[1]/nrow(w_vec),")")) 
covid_risk$cluster[which(covid_risk$cluster==2)] <- c(paste0("2 (",table(covid_risk$cluster)[2]/nrow(w_vec),")")) 
covid_risk$cluster <- as.factor(covid_risk$cluster)

#computed the weighted risk by cluster, according to population sizes
covid_risk$weighted_risk <- covid_risk$est_risk*covid_risk$pop
covid_avg_risk <- aggregate(covid_risk[,c("weighted_risk","pop")],by=list(covid_risk$cluster,covid_risk$week),FUN=sum)
colnames(covid_avg_risk) <- c("cluster","week","sum_risk","sum_pop")
covid_avg_risk$weighted_risk <- covid_avg_risk$sum_risk/covid_avg_risk$sum_pop 
covid_avg_risk <- merge(covid_avg_risk,w_vec)

covid_avg_risk$weeks_since_ld <- as.factor(covid_avg_risk$weeks_since_ld)
covid_risk$weeks_since_ld <- as.factor(covid_risk$weeks_since_ld)

#Figure 4(a)
#Note: colors might be switched
clust_line_ld1 <- ggplot() + geom_line(data=covid_avg_risk,aes(weeks_since_ld,weighted_risk,group=cluster,colour=cluster)) + xlab("") + 
  ylab("Average estimated risk") + ylim(0,9.5) + ggtitle("(a) First lockdown, 26/03/20-12/05/20") +
  labs(fill="Cluster", colour="Cluster") + geom_hline(yintercept = 1, linetype="dashed",colour="red") + 
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_line(data=covid_avg_risk,aes(weeks_since_ld,weighted_risk,group=cluster,colour=cluster)) +
  scale_colour_manual(name="Cluster (# LAD)",values=c("#e41a1c","#377eb8")) + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16),legend.title = element_text(size = 16),
        legend.text = element_text(size=16),plot.title = element_text(size=16))

high_clust <- covid_avg_risk[which(covid_avg_risk$cluster==paste0("1 (",table(covid_risk$cluster)[1]/nrow(w_vec),")")),]
max(high_clust$weighted_risk)
low_clust <- covid_avg_risk[which(covid_avg_risk$cluster==paste0("2 (",table(covid_risk$cluster)[2]/nrow(w_vec),")")),]
max(low_clust$weighted_risk)

#extract only LADs with their assigned clusters
LA_clusters1 <- unique(LA_risk[,c("area","cluster")])

#get numerical summaries for IMD
imd <- read.csv("IMD_index.csv")
#Buckinghamshire is not included in this dataset and needs to be aggregated from the areas listed below
buck <- mean(imd$IMD...Average.rank[which(imd$Local.Authority.District.code..2019. %in%c("E07000004","E07000005","E07000006","E07000007"))])
imd <- imd[,c(1,3)]
buck <- as.data.frame(cbind("E06000060",buck))
colnames(buck) <- colnames(imd)
imd <- rbind(imd,buck)
colnames(imd) <- c("area","avg_rank")
imd$avg_rank <- as.numeric(imd$avg_rank)
#order the LADs according to their average IMD rank
imd <- imd[order(imd$avg_rank),]
#assign the new rank for only those LADs in our data (islands were removed)
imd$imd_rank <- c(nrow(imd):1)

#add imd ranks to the cluster assignment data
LA_clusters <- merge(LA_clusters1,imd)
#select only those LADs in the higher risk cluster
LA_high <- LA_clusters[which(LA_clusters$cluster==1),]
LA_high <- merge(LA_high,covid[,c("area","region")])
LA_high <- unique(LA_high)
#get average IMD rank
mean(LA_high$imd_rank)

LA_low <- LA_clusters[which(LA_clusters$cluster==2),]
LA_low <- merge(LA_low,covid[,c("area","region")])
LA_low <- unique(LA_low)
#get average IMD rank
mean(LA_low$imd_rank)

#show which regions the LADs in the higher risk cluster belong to
table(LA_high$region)

##################
####lockdown 2####
##################
#weeks after lockdown 2 (starts 20_46)
#include 1 week prior as baselevel, and 8 weeks after intro of lockdown
w_vec <- c(paste0("20_",45:53))
#actual weeks of lockdown:
l_int <- c(paste0("20_",46:49))

#weeks_since_lockdown is 0 for last week before lockdown, then 1,2,3,...
w_vec <- as.data.frame(w_vec)
w_vec$weeks_since_ld <- 0:(nrow(w_vec)-1)
colnames(w_vec) <- c("week","weeks_since_ld")

#obtain baseline risk (from last week before lockdown) by area
baseline_risk <- covid[which(covid$week==w_vec$week[1]),]
baseline_risk <- baseline_risk[,c("area","est_risk")]
colnames(baseline_risk) <- c("area","baseline_risk")

#merge weeks after lockdown with other data for that time
wal <- covid[which(covid$week %in% w_vec$week),]
wal <- merge(w_vec,wal)
wal$weeks_since_ld <- as.factor(wal$weeks_since_ld)

#merge data with baseline_risk
wal <- merge(wal,baseline_risk)
#scaled: estimated/baseline
wal$sc_risk <- wal$est_risk/wal$baseline_risk
wal <- wal[order(wal$area,wal$weeks_since_ld),]
#add abbreviations for plot legend
abbr <- as.data.frame(cbind(c("East Midlands","East of England","London","North East",
                              "North West","South East","South West","West Midlands",
                              "Yorkshire and The Humber"),
                            c("EM","EE","L","NE","NW","SE","SW","WM","YH")))
colnames(abbr) <- c("region","abbr")
wal <- merge(wal,abbr)
wal$abbr <- as.factor(wal$abbr)

#Figure 2(b)
bp_est_ld2 <- ggplot(data=wal,aes(weeks_since_ld,est_risk)) + geom_boxplot() +
  ylab("Estimated risk") + xlab("") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept=1,linetype="dashed",colour="red") + ylim(0,15) + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + ggtitle("(b) second lockdown, 05/11/20-02/12/20") +
  geom_hline(yintercept=median(baseline_risk$baseline_risk),linetype="dashed",colour="blue") +
  geom_boxplot()

#week vs estimated risk by region
reg_wal <- wal 
#compute population-based weighted average risk
reg_wal$weighted_risk <- reg_wal$est_risk*reg_wal$pop
reg_wal <- reg_wal[,c("abbr","weighted_risk","pop","week")]
reg_wal <- aggregate(reg_wal[,c("weighted_risk","pop")],by=list(reg_wal$abbr,reg_wal$week),FUN=sum)
colnames(reg_wal) <- c("abbr","week","weighted_risk","pop")
reg_wal$weighted_risk <- reg_wal$weighted_risk/reg_wal$pop
#add weeks_since_ld
reg_wal <- merge(reg_wal,w_vec)
reg_wal$weeks_since_ld <- as.factor(reg_wal$weeks_since_ld)

#Figure 3(b) left plot
line_reg_ld2 <- ggplot() + geom_line(data=reg_wal,aes(weeks_since_ld,weighted_risk,group=abbr,colour=abbr)) + xlab("") + 
  ylab("Estimated risk") +
  labs(fill="Region", colour="Region") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept = 1, linetype="dashed",colour="red") +
  geom_line(data=reg_wal,aes(weeks_since_ld,weighted_risk,group=abbr,colour=abbr)) + ggtitle("(b) Second lockdown, 05/11/20-02/12/20") +
  scale_colour_manual(name="Region",values=myColors,limits=c("EE","EM","L","NE","NW","SE","SW","WM","YH")) + theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.key.size = unit(3,"line"),legend.position = "none") + ylim(0,10.5) +
  guides(color = guide_legend(override.aes = list(size = 2)))

#week vs scaled estimated risk by region
reg_wal_sc <- wal 
reg_wal_sc$weighted_sc_risk <- reg_wal_sc$sc_risk*reg_wal_sc$pop
#compute weighted scaled estimated risk
reg_wal_sc <- reg_wal_sc[,c("abbr","weighted_sc_risk","pop","week")]
reg_wal_sc <- aggregate(reg_wal_sc[,c("weighted_sc_risk","pop")],by=list(reg_wal_sc$abbr,reg_wal_sc$week),FUN=sum)
colnames(reg_wal_sc) <- c("abbr","week","weighted_sc_risk","pop")
reg_wal_sc$weighted_sc_risk <- reg_wal_sc$weighted_sc_risk/reg_wal_sc$pop
#add weeks_since_ld
reg_wal_sc <- merge(reg_wal_sc,w_vec)
reg_wal_sc$weeks_since_ld <- as.factor(reg_wal_sc$weeks_since_ld)

#Figure 3(b) right plot
line_reg_sc_ld2 <- ggplot() + geom_line(data=reg_wal_sc,aes(weeks_since_ld,weighted_sc_risk,group=abbr,colour=abbr)) + xlab("") + 
  ylab("Scaled estimated risk") +
  labs(fill="Region", colour="Region") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept = 1, linetype="dashed",colour="blue") +
  geom_line(data=reg_wal_sc,aes(weeks_since_ld,weighted_sc_risk,group=abbr,colour=abbr)) + ggtitle("") +
  scale_colour_manual(name="Region",values=myColors,limits=c("EE","EM","L","NE","NW","SE","SW","WM","YH")) + theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.key.size = unit(3,"line"),legend.position = "none") + ylim(0,10) +
  guides(color = guide_legend(override.aes = list(size = 2)))

##################
####lockdown 3####
##################
#weeks after lockdown 3 (starts 21_01)
#include 1 week prior as baselevel, and 13 weeks after intro of lockdown
w_vec <- c("20_53",paste0("21_0",1:9),paste0("21_",10:13))
#actual weeks of lockdown:
l_int <- c(paste0("21_0",1:9),paste0("21_",10:12))

#weeks_since_lockdown is 0 for last week before lockdown, then 1,2,3,...
w_vec <- as.data.frame(w_vec)
w_vec$weeks_since_ld <- 0:(nrow(w_vec)-1)
colnames(w_vec) <- c("week","weeks_since_ld")

#obtain baseline risk (from last week before lockdown) by area
baseline_risk <- covid[which(covid$week==w_vec$week[1]),]
baseline_risk <- baseline_risk[,c("area","est_risk")]
colnames(baseline_risk) <- c("area","baseline_risk")

#merge weeks after lockdown with other data for that time
wal <- covid[which(covid$week %in% w_vec$week),]
wal <- merge(w_vec,wal)
wal$weeks_since_ld <- as.factor(wal$weeks_since_ld)

#merge data with baseline_risk
wal <- merge(wal,baseline_risk)
#scaled: estimated/baseline
wal$sc_risk <- wal$est_risk/wal$baseline_risk
wal <- wal[order(wal$area,wal$weeks_since_ld),]
#add abbreviations for legend in plot
abbr <- as.data.frame(cbind(c("East Midlands","East of England","London","North East",
                              "North West","South East","South West","West Midlands",
                              "Yorkshire and The Humber"),
                            c("EM","EE","L","NE","NW","SE","SW","WM","YH")))
colnames(abbr) <- c("region","abbr")
wal <- merge(wal,abbr)
wal$abbr <- as.factor(wal$abbr)

#Figure 2(c)
bp_est_ld3 <- ggplot(data=wal,aes(weeks_since_ld,est_risk)) + geom_boxplot() +
  ylab("Estimated risk") + xlab("Weeks after lockdown") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept=1,linetype="dashed",colour="red") + ylim(0,15) + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + ggtitle("(c) third lockdown, 05/01/21-28/03/21") +
  geom_hline(yintercept=median(baseline_risk$baseline_risk),linetype="dashed",colour="blue") +
  geom_boxplot()

#Figure 2
#grid.arrange(bp_est_ld1,bp_ld1,bp_est_ld2,bp_ld2,bp_est_ld3,bp_ld3,nrow=3,layout_matrix = cbind(c(1,3,5), c(2,4,6),c(2,4,6), c(2,4,6)))
grid.arrange(bp_est_ld1,bp_est_ld2,bp_est_ld3,nrow=3)

#week vs estimated risk by region
reg_wal <- wal 
#compute population-based weighted averages
reg_wal$weighted_risk <- reg_wal$est_risk*reg_wal$pop
reg_wal <- reg_wal[,c("abbr","weighted_risk","pop","week")]
reg_wal <- aggregate(reg_wal[,c("weighted_risk","pop")],by=list(reg_wal$abbr,reg_wal$week),FUN=sum)
colnames(reg_wal) <- c("abbr","week","weighted_risk","pop")
reg_wal$weighted_risk <- reg_wal$weighted_risk/reg_wal$pop
#add weeks_since_ld
reg_wal <- merge(reg_wal,w_vec)
reg_wal$weeks_since_ld <- as.factor(reg_wal$weeks_since_ld)

#Figure 3(c) left plot
line_reg_ld3 <- ggplot() + geom_line(data=reg_wal,aes(weeks_since_ld,weighted_risk,group=abbr,colour=abbr)) + xlab("Weeks after lockdown") + 
  ylab("Estimated risk") +
  labs(fill="Region", colour="Region") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept = 1, linetype="dashed",colour="red") +
  geom_line(data=reg_wal,aes(weeks_since_ld,weighted_risk,group=abbr,colour=abbr)) + ggtitle("(c) Third lockdown, 05/01/21-28/03/21") +
  scale_colour_manual(name="Region",values=myColors,limits=c("EE","EM","L","NE","NW","SE","SW","WM","YH")) + theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.key.size = unit(3,"line"),legend.position = "none") + ylim(0,10.5) +
  guides(color = guide_legend(override.aes = list(size = 2)))

#week vs scaled estimated risk by region
reg_wal_sc <- wal 
#compute population-based weighted scaled risks
reg_wal_sc$weighted_sc_risk <- reg_wal_sc$sc_risk*reg_wal_sc$pop
reg_wal_sc <- reg_wal_sc[,c("abbr","weighted_sc_risk","pop","week")]
reg_wal_sc <- aggregate(reg_wal_sc[,c("weighted_sc_risk","pop")],by=list(reg_wal_sc$abbr,reg_wal_sc$week),FUN=sum)
colnames(reg_wal_sc) <- c("abbr","week","weighted_sc_risk","pop")
reg_wal_sc$weighted_sc_risk <- reg_wal_sc$weighted_sc_risk/reg_wal_sc$pop
#add weeks_since_ld
reg_wal_sc <- merge(reg_wal_sc,w_vec)
reg_wal_sc$weeks_since_ld <- as.factor(reg_wal_sc$weeks_since_ld)

#Figure 3(c) right plot
line_reg_sc_ld3 <- ggplot() + geom_line(data=reg_wal_sc,aes(weeks_since_ld,weighted_sc_risk,group=abbr,colour=abbr)) + xlab("Weeks after lockdown") + 
  ylab("Scaled estimated risk") +
  labs(fill="Region", colour="Region") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_hline(yintercept = 1, linetype="dashed",colour="blue") +
  geom_line(data=reg_wal_sc,aes(weeks_since_ld,weighted_sc_risk,group=abbr,colour=abbr)) + ggtitle("") +
  scale_colour_manual(name="Region",values=myColors,limits=c("EE","EM","L","NE","NW","SE","SW","WM","YH")) + theme(axis.text=element_text(size=14),axis.title=element_text(size=16),legend.key.size = unit(3,"line"),legend.position = "none") + ylim(0,10) +
  guides(color = guide_legend(override.aes = list(size = 2)))

#Figure 3
grid.arrange(line_reg_ld1,line_reg_sc_ld1,line_reg_ld2,line_reg_sc_ld2,line_reg_ld3,line_reg_sc_ld3,legend_reg_plot,nrow=3,layout_matrix = rbind(c(1,2,7),c(3,4,7),c(5,6,7)),
             widths = c(5,5,1))

# Determine number of clusters
LA_risk <- wal[,c("area","week","est_risk")]
LA_risk <- LA_risk[which(LA_risk$week %in% w_vec$week[2:nrow(w_vec)]),]
LA_risk <- LA_risk[with(LA_risk,order(area,week)),]
LA_risk <- reshape(LA_risk, idvar = "area", timevar = "week", direction = "wide")

#only risks for average silh. widths and wss
risks <- LA_risk[,2:ncol(LA_risk)]

avg_sil_values <- map_dbl(k.values, avg_sil)
avg_sil_values <- as.data.frame(cbind(k.values,avg_sil_values))

wcss <- (nrow(LA_risk)-1)*sum(apply(risks,2,var))
for(i in 2:10){
  wcss[i] <- sum(kmeans(risks,centers=i,nstart = 25)$withinss)
} 
wcss <- as.data.frame(cbind(1:10,wcss))
colnames(wcss) <- c("no_clusters","wcss")

sil_ld3 <- ggplot() + geom_line(data=avg_sil_values,aes(x=k.values,y=avg_sil_values,group=1)) + geom_point(data=avg_sil_values,aes(x=k.values,y=avg_sil_values,group=1),colour="#003865") +
  scale_x_continuous(breaks=2:10) + theme(panel.grid.minor.x = element_blank()) + xlab("Number of clusters") + ylab("Average Silhouette Width") + ggtitle("(b) Lockdown 3") + ylim(0,0.6) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),plot.title = element_text(size=16))

wcss_ld3 <- ggplot() + geom_line(data=wcss,aes(x=no_clusters,y=wcss,group=1)) + geom_point(data=wcss,aes(x=no_clusters,y=wcss,group=1),colour="#003865") +
  scale_x_continuous(breaks=1:10) + theme(panel.grid.minor.x = element_blank()) + xlab("Number of clusters") + ylab("Within-Cluster Sum of Squares") + ggtitle("") + ylim(0,6005) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),plot.title = element_text(size=16))

grid.arrange(sil_ld1,wcss_ld1,sil_ld3,wcss_ld3,nrow=2)

#output above suggests 2 clusters
clust <- kmeans(risks,2,nstart=25)
LA_risk$cluster <- clust$cluster
covid_risk <- LA_risk[,c("area","cluster")]
covid_risk <- merge(covid_risk,wal)
#repeat the 4 lines above if cluster 1 is not the smaller cluster
table(covid_risk$cluster)
covid_risk$cluster[which(covid_risk$cluster==1)] <- c(paste0("1 (",table(covid_risk$cluster)[1]/nrow(w_vec),")")) 
covid_risk$cluster[which(covid_risk$cluster==2)] <- c(paste0("2 (",table(covid_risk$cluster)[2]/nrow(w_vec),")")) 

#computed the weighted risk by cluster, according to population sizes
covid_risk$weighted_risk <- covid_risk$est_risk*covid_risk$pop
covid_avg_risk <- aggregate(covid_risk[,c("weighted_risk","pop")],by=list(covid_risk$cluster,covid_risk$week),FUN=sum)
colnames(covid_avg_risk) <- c("cluster","week","sum_risk","sum_pop")
covid_avg_risk$weighted_est_risk <- covid_avg_risk$sum_risk/covid_avg_risk$sum_pop 
covid_avg_risk <- merge(covid_avg_risk,w_vec)

covid_avg_risk$weeks_since_ld <- as.factor(covid_avg_risk$weeks_since_ld)

#Figure 4(b)
#Note: colors might be switched
clust_line_ld3 <- ggplot() + geom_line(data=covid_avg_risk,aes(weeks_since_ld,weighted_est_risk,group=cluster,colour=cluster)) + xlab("Weeks after lockdown") + 
  ylab("Average estimated risk") + ylim(0,9.5) + ggtitle("(b) Third lockdown, 05/01/21-28/03/21") +
  labs(fill="Cluster", colour="Cluster") + geom_hline(yintercept = 1, linetype="dashed",colour="red") +
  annotate("rect",xmin=unique(wal$weeks_since_ld[which(wal$week==l_int[1])]),xmax=unique(wal$weeks_since_ld[which(wal$week==l_int[length(l_int)])]),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_line(data=covid_avg_risk,aes(weeks_since_ld,weighted_est_risk,group=cluster,colour=cluster)) +
  scale_colour_manual(name="Cluster (# LAD)",values=c("#e41a1c","#377eb8")) + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16),legend.title = element_text(size = 16),legend.text = element_text(size=16),plot.title = element_text(size=16))

grid.arrange(clust_line_ld1,clust_line_ld3,nrow=2)

#get peak risks for the two clusters
high_clust <- covid_avg_risk[which(covid_avg_risk$cluster==paste0("1 (",table(covid_risk$cluster)[1]/nrow(w_vec),")")),]
max(high_clust$weighted_est_risk)
low_clust <- covid_avg_risk[which(covid_avg_risk$cluster==paste0("2 (",table(covid_risk$cluster)[2]/nrow(w_vec),")")),]
max(low_clust$weighted_est_risk)

#plot clusters on a map
LA_clusters3 <- unique(LA_risk[,c("area","cluster")])

#add IMD data
LA_clusters <- merge(LA_clusters3,imd)
#look only at higher risk cluster
LA_high <- LA_clusters[which(LA_clusters$cluster==1),]
LA_high <- merge(LA_high,covid[,c("area","region")])
LA_high <- unique(LA_high)
#get average IMD rank
mean(LA_high$imd_rank)

LA_low <- LA_clusters[which(LA_clusters$cluster==2),]
LA_low <- merge(LA_low,covid[,c("area","region")])
LA_low <- unique(LA_low)
#get average IMD rank
mean(LA_low$imd_rank)

#show which regions the LADs in higher risk cluster belong to
table(LA_high$region)

######
#Contingency tables
cluster_table_13 <- data.frame("Cluster_Lockdown_1"=LA_clusters1$cluster,"Cluster_Lockdown_3"=LA_clusters3$cluster)
table(cluster_table_13)
#rand index between clusterings from lockdowns 1 and 3
rand.index(LA_clusters1$cluster,LA_clusters3$cluster)
