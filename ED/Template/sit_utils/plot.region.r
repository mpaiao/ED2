#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#


#----- Paths. -----------------------------------------------------------------------------#
main    = "/n/moorcroftfs2/mlongo/EDBRAMS/final_ed/trop_southam/pve+sas" # Main path.
here    = file.path(main,"sit_utils")          # This path
srcdir  = c( "/Users/mlongo/Util/Rsc"          #    R will select the first one
           , "/n/home00/mlongo/Util/Rsc"       #
           )#end c                             #
outroot = here                                 # Directory for figures
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      List files to be read.  The current and the previous check, and the path with the   #
# job order.                                                                               #
#------------------------------------------------------------------------------------------#
joborder  = file.path(main,"joborder.txt" )
lastcheck = file.path(here,"lastcheck.txt")
mycheck   = file.path(here,"mycheck.txt"  )
#------------------------------------------------------------------------------------------#




#----- Plot options. ----------------------------------------------------------------------#
outform        = c("png")  # Formats for output file.  Supported formats are:
                           #   - "X11" - for printing on screen
                           #   - "eps" - for postscript printing
                           #   - "png" - for PNG printing
                           #   - "pdf" - for PDF printing
depth          = 96        # PNG resolution, in pixels per inch
paper          = "letter"  # Paper size, to define the plot shape
ptsz           = 14        # Font size.
lwidth         = 2.5       # Line width
inset          = 0.01      # inset between legend and edge of plot region.
fracexp        = 0.40      # Expand the y axis by this amount.
ncolours       = 20        # Number of colours to split the real variables
mtext.xoff     = -8.50     # Offset for the x label
mtext.yoff     = -1.00     # Offset for the y label
mtext.xadj     =  0.50     # Offset for the x label
mtext.yadj     =  0.65     # Offset for the y label
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#




#----- Loading some packages and scripts. -------------------------------------------------#
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#----- Set how many formats we must output. -----------------------------------------------#
outform = tolower(outform)
nout    = length (outform)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
dummy = dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Read job order.  We always use this file to build the array.                         #
#------------------------------------------------------------------------------------------#
cat0(" + Read ",basename(joborder),".")
names.jobs       = scan(file=joborder,skip=1,nlines=1,what="character",quiet=TRUE)
names.jobs       = gsub(pattern="_",replacement=".",x=tolower(names.jobs))
jobs             = read.table(file=joborder,skip=3,header=FALSE,comment.char=""
                             ,stringsAsFactors=FALSE)
names(jobs)      = names.jobs
njobs            = nrow(jobs)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Read the last and the current check.  For the current check, we normally skip the    #
# last line to avoid trouble, unless the file is complete.                                 #
#------------------------------------------------------------------------------------------#
names.check   = c("run","lon","lat","year","month","day","hhmm","runt"
                 ,"agb","bsa","lai","scb")
cat0(" + Read ",basename(lastcheck),".")
last          = read.table(file=lastcheck,skip=0,header=FALSE,comment.char=""
                          ,col.names=names.check,stringsAsFactors=FALSE)
cat0(" + Read ",basename(mycheck),".")
ncurr         = length(readLines(mycheck))
if (ncurr == njobs){
   curr = read.table(file=mycheck,skip=0,header=FALSE,comment.char=""
                    ,col.names=names.check,stringsAsFactors=FALSE)
}else if (ncurr > 0){
   curr = read.table(file=mycheck,skip=0,nrows=ncurr-1,header=FALSE,comment.char=""
                    ,col.names=names.check,stringsAsFactors=FALSE)
}else{
   curr = data.frame(rep(NA,times=length(names.check)),names=names.check)
   curr = curr[-1,]
}#end f
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Find the dimensions.                                                                 #
#------------------------------------------------------------------------------------------#
lon             = unique(sort(jobs$lon        ))
lat             = unique(sort(jobs$lat        ))
dlon            = median(diff(lon))
dlat            = median(diff(lat))
n.lon           = length(lon        )
n.lat           = length(lat        )
key.lon         = sprintf("%+06.2f",lon)
key.lat         = sprintf("%+06.2f",lat)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Initialise the variables.                                                            #
#------------------------------------------------------------------------------------------#
template   = array( data     = NA
                  , dim      = c(n.lon,n.lat)
                  , dimnames = list(key.lon,key.lat)
                  )#end array
datum      = list ( agb    = template
                  , lai    = template
                  , bsa    = template
                  , scb    = template
                  , status = template
                  , yearn  = template
                  )#end list
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Copy the information from last check to joborder.                                    #
#------------------------------------------------------------------------------------------#
jobs$status = rep("INITIAL",times=njobs)
jobs$yearn  = jobs$yeara
jobs$agb    = rep(NA,times=njobs)
jobs$bsa    = rep(NA,times=njobs)
jobs$lai    = rep(NA,times=njobs)
jobs$scb    = rep(NA,times=njobs)

il = match(last$run,jobs$run); l.sel = ! is.na(il)
ic = match(curr$run,jobs$run); c.sel = ! is.na(ic)
jobs$status[il[l.sel]] = last$runt[l.sel]  ; jobs$status[ic[c.sel]] = curr$runt[c.sel]
jobs$agb   [il[l.sel]] = last$agb [l.sel]  ; jobs$agb   [ic[c.sel]] = curr$agb [c.sel]
jobs$lai   [il[l.sel]] = last$lai [l.sel]  ; jobs$lai   [ic[c.sel]] = curr$lai [c.sel]
jobs$bsa   [il[l.sel]] = last$bsa [l.sel]  ; jobs$bsa   [ic[c.sel]] = curr$bsa [c.sel]
jobs$scb   [il[l.sel]] = last$scb [l.sel]  ; jobs$scb   [ic[c.sel]] = curr$scb [c.sel]
jobs$yearn [il[l.sel]] = last$year[l.sel]  ; jobs$yearn [ic[c.sel]] = curr$year[c.sel]





keep            = names(jobs) %in% c("lon","lat","iata","yeara","yearz","yearn","status"
                                     ,"agb","bsa","lai","scb")
jobs            = jobs[,keep]
weird           = is.finite(jobs$lai) & abs(jobs$lai) > 20
jobs$lai[weird] = NA
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Find the indices to map the data to the arrays.                                      #
#------------------------------------------------------------------------------------------#
i.lon               = match(jobs$lon        ,lon        )
i.lat               = match(jobs$lat        ,lat        )
index               = cbind(i.lon,i.lat)
datum$agb   [index] = jobs$agb
datum$lai   [index] = jobs$lai
datum$bsa   [index] = jobs$bsa
datum$scb   [index] = jobs$scb
datum$status[index] = jobs$status
datum$yearn [index] = jobs$yearn
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Run the matrices.                                                                    #
#------------------------------------------------------------------------------------------#
yr.range              = range(c(jobs$yeara,jobs$yearz))
yr.cut                = pretty(yr.range,n=10)
yr.keep               = yr.cut > yr.range[1] & yr.cut < yr.range[2]
yr.brks               = c(-Inf,c(yr.range[1],yr.cut[yr.keep],yr.range[2]))
n.cut                 = length(yr.cut)-1
yr.cut                = cut(datum$yearn,yr.brks)
yr.level              = levels(yr.cut)
n.level               = length(yr.level)
#------------------------------------------------------------------------------------------#



#------ Find the year colours. ------------------------------------------------------------#
datum$yr.idx          = match(yr.cut,yr.level) + 0 * datum$yearn
initial               = datum$status == "INITIAL"
crashed               = datum$status == "CRASHED"
hydfail               = datum$status == "HYDFAIL"
bad.met               = datum$status == "BAD_MET"
metmiss               = datum$status == "METMISS"
stopped               = datum$status == "STOPPED"
extinct               = datum$status == "EXTINCT"
ststate               = datum$status == "STSTATE"
the.end               = datum$status == "THE_END"
datum$yr.idx[initial] = 0
datum$yr.idx[the.end] = n.level + 1
datum$yr.idx[ststate] = n.level + 2
datum$yr.idx[extinct] = n.level + 3
datum$yr.idx[stopped] = n.level + 4
datum$yr.idx[metmiss] = n.level + 5
datum$yr.idx[bad.met] = n.level + 6
datum$yr.idx[crashed] = n.level + 7
datum$yr.idx[hydfail] = n.level + 8
yr.cscheme            = c("grey89",iatlas(n=n.level),"royalblue4","steelblue3","purple3"
                         ,"mediumpurple1","deepskyblue","firebrick4","hotpink","red3")
ybottom               = rep(0,times=n.level+8)
ytop                  = rep(1,times=n.level+8)
xleft                 = seq(from=-1,to=n.level+6)
xright                = seq(from= 0,to=n.level+7)
xat                   = seq(from=-1,to=n.level+6)+0.5
xbrks                 = seq(from=-1,to=n.level+7)+0.5
xlabel                = c("Initial",yr.brks[-1],"Finish","StState","Extinct"
                         ,"Stopped","MetMiss","Bad Met","Crashed","HydFail")
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Limits for longitude and latitude.                                                  #
#------------------------------------------------------------------------------------------#
limlon = c(min(lon)-0.5*dlon,max(lon)+0.5*dlon)
limlat = c(min(lat)-0.5*dlat,max(lat)+0.5*dlat)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Create the status map for all sites.                                                #
#------------------------------------------------------------------------------------------#
cat0(" Plot the current status.")

#------ Make plot annotation. -------------------------------------------------------------#
letitre      = paste("Polygon status")
#------------------------------------------------------------------------------------------#

for (o in sequence(nout)){
   #----- Open file. ----------------------------------------------------------------------#
   fichier = file.path(here,paste0("stt_region.",outform[o]))
   dummy   = open.plot( fichier = fichier
                      , outform = outform[o]
                      , size    = size
                      , depth   = depth
                      , ptsz    = ptsz
                      )#end open.plot
   #---------------------------------------------------------------------------------------#


   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig = par(no.readonly = TRUE)
   mar.orig = par.orig$mar
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Split the plotting window.                                                       #
   #---------------------------------------------------------------------------------------#
   layout( mat     = rbind(2,1)
         , heights = c(3,1)
         )#end layout
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     First, let's plot the legend.                                                     #
   #---------------------------------------------------------------------------------------#
   par(mar=c(3,2,2,2)+0.1)
   plot.new()
   plot.window(xlim=range(xleft,xright),ylim=range(ybottom,ytop),xaxs="i",yaxs="i")
   rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,col=yr.cscheme)
   box()
   axis(side=1,at=xat,srt=45,labels=FALSE)
   text(x=xat,y=par("usr")[3]-0.6,labels=xlabel,srt=30,adj=1,xpd=TRUE,cex=1.00)
   title(main="Status",ylab="",xlab="")
   #---------------------------------------------------------------------------------------#





   #----- Set the window. -----------------------------------------------------------------#
   par(mar = c(3.1,3.1,4.1,2.1))
   plot.new()
   plot.window(xlim=limlon,ylim=limlat,xaxs="i",yaxs="i")
   axis(side=1)
   axis(side=2)
   box()
   title(main=letitre,xlab="",ylab="")
   image(x=lon,y=lat,z=datum$yr.idx,col=yr.cscheme,breaks=xbrks,add=TRUE)
   southammap()
   #---------------------------------------------------------------------------------------#



   #----- Close the device. ---------------------------------------------------------------#
   dummy = close.plot(outform=outform[o])
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Create parameter space maps for all other variables.                                #
#------------------------------------------------------------------------------------------#
cat0(" Plot the current properties...","\n")
key.var  = c("lai","bsa","agb","scb")
desc.var = c( desc.unit( desc = "Leaf area index"     , unit = untab$m2lom2)
            , desc.unit( desc = "Basal area"          , unit = untab$cm2om2)
            , desc.unit( desc = "Above-ground biomass", unit = untab$kgcom2)
            , desc.unit( desc = "Soil carbon"         , unit = untab$kgcom2)
            )#end c
n.var    = length(key.var)
for (v in sequence(n.var)){
   cat0("   - ",desc.var[v],".")
   #----- Collapse realisations using the median. -----------------------------------------#
   this.var       = datum[[key.var[v]]]
   rien           = ! is.finite(this.var)
   this.var[rien] = NA
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Break the data into bins.                                                        #
   #---------------------------------------------------------------------------------------#
   if (all(is.na(this.var))){
     var.brks    = c(-1,0,1)
     n.brks      = length(var.brks)
     var.cut     = cut(as.numeric(this.var),breaks=var.brks)
   }else{
     var.brks    = pretty(this.var,n=ncolours)
     n.brks      = length(var.brks)
     var.cut     = cut(this.var,breaks=var.brks)
   }#end if
   var.lev     = levels(var.cut)
   var.idx     = match(var.cut,var.lev) + 0 * this.var
   var.cscheme = iatlas(n=n.brks-1)
   #---------------------------------------------------------------------------------------#


   #----- Make the edges. -----------------------------------------------------------------#
   xleft       = var.brks[-n.brks]
   xright      = var.brks[     -1]
   ybottom     = rep(0,times=n.brks)
   ytop        = rep(1,times=n.brks)
   xat         = var.brks
   #---------------------------------------------------------------------------------------#


   #------ Find the soil texture key and description. -------------------------------------#
   letitre      = paste(desc.var[v])
   #---------------------------------------------------------------------------------------#

   for (o in sequence(nout)){
      #----- Open file. -------------------------------------------------------------------#
      fichier = file.path(here,paste0(key.var[v],"_region.",outform[o]))
      dummy   = open.plot( fichier = fichier
                         , outform = outform[o]
                         , size    = size
                         , depth   = depth
                         , ptsz    = ptsz
                         )#end open.plot
      #------------------------------------------------------------------------------------#


      #----- Save the margins to avoid losing the data. -----------------------------------#
      par.orig = par(no.readonly = TRUE)
      mar.orig = par.orig$mar
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Split the plotting window.                                                    #
      #------------------------------------------------------------------------------------#
      layout( mat     = rbind(2,1)
            , heights = c(4,1)
            )#end layout
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     First, let's plot the legend.                                                  #
      #------------------------------------------------------------------------------------#
      par(mar=c(3,3,2,2)+0.1)
      plot.new()
      plot.window(xlim=range(xleft,xright),ylim=range(ybottom,ytop),xaxs="i",yaxs="i")
      rect(xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,col=var.cscheme)
      box()
      axis(side=1,at=xat)
      title(main=desc.var[v],xlab="",ylab="")
      #------------------------------------------------------------------------------------#



      #----- Set the window. --------------------------------------------------------------#
      par(mar = c(3.1,3.1,4.1,2.1))
      plot.new()
      plot.window(xlim=limlon,ylim=limlat,xaxs="i",yaxs="i")
      axis(side=1)
      axis(side=2)
      box()
      title(main=letitre,xlab="",ylab="")
      image(x=lon,y=lat,z=this.var,col=var.cscheme,breaks=var.brks,add=TRUE)
      southammap()
      #------------------------------------------------------------------------------------#



      #----- Close the device. ------------------------------------------------------------#
      dummy = close.plot(outform=outform[o])
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#
