library(scales)

pal37 = c('208 210 85',
         '100 100 80',
         '74 73 19',
         '52 18 8',
         '45 12 11',
         '147 34 30',
         '223 58 37',
         '242 179 85',
         '250 224 80',
         '196 41 44',
         '128 46 92',
         '129 68 122',
         '132 74 152',
         '96 44 134',
         '72 29 107',
         '60 11 58',
         '94 160 186',
         '32 33 69',
         '80 113 141',
         '45 22 42',
         '28 37 79',
         '31 64 114',
         '26 64 53',
         '21 51 18',
         '17 45 14',
         '103 149 81',
         '150 206 94',
         '89 79 56',
         '46 52 36',
         '179 216 99',
         '14 34 12',
         '120 149 45',
         '223 234 77',
         '157 74 91',
         '104 91 209',
         '74 152 200',
         '40 99 175')
length(pal37)
pal37H = sapply(strsplit(pal37, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

#show_col(pal37H)



pal2 = c('186 137 65',
         '161 116 61',
         '156 125 84',
         '45 67 114',
         '60 81 121',
         '52 119 182',
         '17 48 96',
         '107 108 145',
         '82 102 145',
         '117 155 214',
         '209 48 99',
         '226 83 76',
         '239 155 80',
         '234 177 72',
         '214 49 106',
         '169 40 33',
         '191 90 39',
         '150 60 35',
         '165 51 77',
         '87 39 145',
         '80 41 133',
         '91 33 115',
         '140 61 153',
         '162 64 166',
         '138 203 237',
         '59 135 199',
         '34 77 137',
         '13 34 83',
         '148 33 20',
         '83 50 25',
         '227 110 46',
         '234 152 57',
         '241 194 75',
         '122 75 130')
length(pal2)

pal34H = sapply(strsplit(pal2, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

#show_col(pal34H)

#' Adds transparency to colour
#'
#' @param cols Vector of colours.
#' @param alphas Single value or vector of alphas
#' @param ... Passed to rgb
#' @return rgb colours with transparency set.
colAlpha = function(cols,alphas,...) {
  if(length(alphas)==1)
    alphas = rep(alphas,length(cols))
  tmp = col2rgb(cols)
  sapply(seq_len(ncol(tmp)),function(e) rgb(tmp[1,e],tmp[2,e],tmp[3,e],alphas[e]*255,maxColorValue=255,...))
}


col25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  #"black", 
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)



col22=c("#9cb169","#a75acb","#59b648","#636edd","#99b534","#d5439a","#51bf7f","#d73f52","#4cc8c6","#cd5d2a","#609dd8","#d0a63d","#5e63a9","#687327","#bd8dd9","#377b40","#9c4c89","#479e80","#df80ae","#a8793e","#a74a57","#e28674")


theme_classic_2 = theme_classic(base_size = 15) + theme(panel.border = element_rect(fill=F),axis.line = element_blank())



#' Save raster and non-raster versions of plots
#'
#' Saves a pdf and png version of everything.
#'
#' @param baseName Name of file to save.  File name without extension.
#' @param plotFun Function to generate the plot.  Must take two arguments, noFrame and noPlot which control if the plot part / frame part are drawn.
#' @param width Width in inches.
#' @param heights Height in inches.
#' @param res Resolution in pixels per inch for rasterised image.
#' @param rawData If provided, save the raw data used to make the plot.
#' @param row.names Passed to write.table
#' @param col.names Passed to write.table
#' @param ... Passed to plotFun
saveFig = function(baseName,plotFun,width=4,height=3,res=300,rawData=NULL,row.names=FALSE,col.names=TRUE,...){
  #Save the base versions
  pdf(paste0(baseName,'.pdf'),width=width,height=height,useDingbats=FALSE)
  plotFun(noFrame=FALSE,noPlot=FALSE,...)
  dev.off()
  png(paste0(baseName,'.png'),width=width,height=height,res=300,units='in')
  plotFun(...)
  dev.off()
  #And the various bits of frame in vector format
  pdf(paste0(baseName,'_frame.pdf'),width=width,height=height,useDingbats=FALSE)
  plotFun(noPlot=TRUE,noFrame=FALSE,...)
  dev.off()
  #And the various bits of the actual plot in raster format
  png(paste0(baseName,'_plot.png'),width=width,height=height,res=300,units='in')
  plotFun(noPlot=FALSE,noFrame=TRUE,...)
  dev.off()
  #Save raw data if provided
  if(!is.null(rawData)){
    write.table(rawData,paste0(baseName,'_rawData.tsv'),
                quote=FALSE,
                sep='\t',
                row.names=row.names,
                col.names=col.names)
  }
}



  
  
  
  
  
  
  
  
  
  
  

  
