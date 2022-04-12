library(ggplot2)
library(spdep)

#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

covid <- read.csv("coviddata.csv")

#Data Section
#show regions on map (see supplementary material)
one_week <- covid[which(covid$week=="21_03"),]
Bdry <- readOGR(dsn="Local_Authority_Districts__May_2020__Boundaries_UK_BFC.shp")
Covid.Bdry <- merge(x=Bdry,y=one_week,by.x="LAD20CD",by.y="area",all.x=FALSE)
Covid.Bdry.ll <- spTransform(Covid.Bdry,CRS("+proj=longlat +datum=WGS84 +no_defs"))
variable <- Covid.Bdry.ll@data$region
colours <- colorFactor(palette=c("#377eb8","#e41a1c","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999"),domain=variable,reverse=FALSE)
leaflet(data=Covid.Bdry.ll) %>% addTiles() %>% addPolylines(opacity=1,weight=1,color="black") %>%
  addPolygons(fillColor=~colours(variable),color="",fillOpacity=0.7,
              weight=1,smoothFactor = 0,opacity=1.0) %>%
  addLegend(pal=colours,values=variable,opacity=1,title="Regions") %>%
  addScaleBar(position="bottomleft")

#compute SMR
covid$smr <- covid$deaths/covid$expecteddeaths

#extract week start dates for plot axis
week_starts <- gsub("-.*","",covid$week_span)
week_starts <- paste0(week_starts,gsub("_.*","",covid$week))
covid$starts <- week_starts
ticks <- unique(covid[,c("week","starts")])
ticks <- ticks[seq(from=1,to=67,by=4),]

#get average SMR by week
avg_SMR <- aggregate(covid$smr,by=list(covid$week),FUN=mean)
colnames(avg_SMR) <- c("week","avg_SMR")

#### plot SMR
#Figure 1(a)
ggplot() + geom_point(data=covid,aes(x=week,y=smr)) + 
  ylab("SMR") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(name = "Week starting on", breaks=ticks$week, labels=ticks$starts) +
  annotate("rect",xmin=c("20_14"),xmax=c("20_20"),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00")  +
  annotate("rect",xmin=c("20_46"),xmax=c("20_49"),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  annotate("rect",xmin=c("21_01"),xmax=c("21_12"),ymin=0,ymax=Inf,alpha=0.2,fill="#E69F00") +
  geom_line(data=avg_SMR,aes(x=week,y=avg_SMR,group=1),colour="black",size=1)

#create the function below to modify leaflet
##############
#this function allows to plot the legend in opposite order
addLegend_decreasing <- function (map, position = c("topright", "bottomright", "bottomleft", 
                                                    "topleft"), pal, values, na.label = "NA", bins = 7, colors, 
                                  opacity = 0.5, labels = NULL, labFormat = labelFormat(), 
                                  title = NULL, className = "info legend", layerId = NULL, 
                                  group = NULL, data = getMapData(map), decreasing = FALSE) {
  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors)) 
      stop("You must provide either 'pal' or 'colors' (not both)")
    if (missing(title) && inherits(values, "formula")) 
      title <- deparse(values[[2]])
    values <- evalFormula(values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && col2rgb(na.color, alpha = TRUE)[[4]] == 
        0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins)) 
      warning("'bins' is ignored because the palette type is not numeric")
    if (type == "numeric") {
      cuts <- if (length(bins) == 1) 
        pretty(values, bins)
      else bins	
      
      if (length(bins) > 2) 
        if (!all(abs(diff(bins, differences = 2)) <= 
                 sqrt(.Machine$double.eps))) 
          stop("The vector of breaks 'bins' must be equally spaced")
      n <- length(cuts)
      r <- range(values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1])/(r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE){
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      }else{
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")
      
    }
    else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n])/2
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }
      
    }
    else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- quantile(values, probs = p, na.rm = TRUE)
      mids <- quantile(values, probs = (p[-1] + p[-n])/2, 
                       na.rm = TRUE)
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    }
    else if (type == "factor") {
      v <- sort(unique(na.omit(values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE){
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      }else{
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    }
    else stop("Palette function not supported")
    if (!any(is.na(values))) 
      na.color <- NULL
  }
  else {
    if (length(colors) != length(labels)) 
      stop("'colors' and 'labels' must be of the same length")
  }
  legend <- list(colors = I(unname(colors)), labels = I(unname(labels)), 
                 na_color = na.color, na_label = na.label, opacity = opacity, 
                 position = position, type = type, title = title, extra = extra, 
                 layerId = layerId, className = className, group = group)
  invokeMethod(map, data, "addLegend", legend)
}

#Plot average SMR by LAD on a map
##############
#take average SMR by LAD
avg_SMR <- aggregate(x=covid[,c("smr")],by=list(covid$area),FUN=mean) 
colnames(avg_SMR) <- c("area","avg_SMR")
Covid.Bdry <- merge(x=Bdry,y=avg_SMR,by.x="LAD20CD",by.y="area",all.x=FALSE)
Covid.Bdry.ll <- spTransform(Covid.Bdry,CRS("+proj=longlat +datum=WGS84 +no_defs"))

variable <- Covid.Bdry.ll@data$avg_SMR
colours <- colorNumeric(palette="RdBu",domain=c(min(variable),1,max(variable)),reverse=TRUE)
#plot average SMR on map
#Figure 1(b)
leaflet(data=Covid.Bdry.ll) %>% addTiles() %>% 
  addPolygons(fillColor=~colours(variable),color="",fillOpacity=0.7,
              weight=1,smoothFactor = 0.5,opacity=1.0) %>%
  addLegend_decreasing(pal=colours,values=variable,opacity=1,title="SMR",decreasing=TRUE)

#Moran's I test
weeks <- unique(covid$week)
stats <- c(rep(NA,length(weeks)))
pvals <- c(rep(NA,length(weeks)))

week1 <- covid[which(covid$week==weeks[1]),]
Covid.Bdry <- merge(x=Bdry,y=week1,by.x="LAD20CD",by.y="area",all.x=FALSE)
Covid.Bdry.ll <- spTransform(Covid.Bdry,CRS("+proj=longlat +datum=WGS84 +no_defs"))
#create binary adjacency neighbourhood matrix
W.nb <- poly2nb(Covid.Bdry.ll,row.names=Covid.Bdry.ll@data$LAD20CD)
W.list <- nb2listw(W.nb,style="B")

#compute the p-value for Moran's I at each week
#this takes a few minutes
for(i in 1:length(weeks)){
  weeki <- covid[which(covid$week==weeks[i]),]
  Covid.Bdry <- merge(x=Bdry,y=weeki,by.x="LAD20CD",by.y="area",all.x=FALSE)
  Covid.Bdry.ll <- spTransform(Covid.Bdry,CRS("+proj=longlat +datum=WGS84 +no_defs"))
  mori <- moran.mc(x=Covid.Bdry.ll$smr,listw=W.list,nsim=10000)
  stats[i] <- as.numeric(mori$statistic)
  pvals[i] <- as.numeric(mori$p.value)
}

stats_new <- stats[!is.na(stats)]
c(min(stats_new),max(stats_new),mean(stats_new),median(stats_new))
pvals_new <- pvals[!is.na(stats)]
c(min(pvals_new),max(pvals_new),mean(pvals_new),median(pvals_new))
#number of p-values for which we reject the null hypothesis of no spatial autocorrelation
sum(pvals_new<0.00077)
#proportion of weeks with at least one death for which the null-hypothesis is rejected
sum(pvals_new<0.00077)/65

#temporal autocorrelation test
areas <- unique(covid$area)
pvals <- c(rep(NA,length(areas)))
for(i in 1:length(areas)){
  one.area <- covid[covid$area==areas[i],]
  lb <- Box.test(x=one.area$smr,lag=10,type="Ljung-Box")  
  pvals[i] <- lb$p.value
}
#number of p-values for which we reject the null hypothesis of no temporal autocorrelation
sum(pvals<0.00016)
#proportion of LADs for which the null-hypothesis is rejected
sum(pvals<0.00016)/312
