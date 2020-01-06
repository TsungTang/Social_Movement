list.of.packages <- c("ggplot2", "Rcpp" , 'data.table' , 'pipeR' ,
                      'haven' ,'igraph','useful','stringr','factoextra','maps' ,
                      'cluster','lubridate','RColorBrewer','readxl','plotly','tnet')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)
library(stringr)
library(igraph)
library(pipeR)
library(lubridate)
library(RColorBrewer)
library(readxl)

# copy xlsx to working directory

  root.dir = str_extract(getwd() , ".+?(?=/R.analysis)" )
  lastdate.pos = which.max(as.numeric(str_extract(list.files(root.dir) , '[0-9]*')) )
  
  last.data.dir = paste0(root.dir , '/' , list.files(root.dir)[lastdate.pos]) 
  
  last.rawdata = str_subset(list.files(last.data.dir) ,'xlsx' )
  if(length(last.rawdata) > 1) {
    last.rawdata = last.rawdata[-1]
  }
  
  lastdate = last.rawdata %>>% str_extract('[0-9]*')
  
  file.copy(paste0(last.data.dir , '/' ,last.rawdata) , 
            paste0(getwd() ,'//data//' , lastdate ,'.xlsx' ) ,overwrite = T )
  
  dir = paste0('data//' , lastdate , '.xlsx')
  read_xlsx( dir , sheet = 1)  %>>% 
    setDT() %>>% 
    fwrite(paste0('data//',lastdate , 'pecktime.csv'))
  
  read_xlsx( dir , sheet = 2)  %>>% 
    setDT() %>>% 
    fwrite(paste0('data//' ,lastdate , 'namelist.csv'))
  



## Read the latest data
csv.list <- str_subset(list.files('data') ,'csv') 
latest.date <- max(str_extract(csv.list , '[0-9]+') %>>% as.numeric()) %>>% as.character()
latest.data <- str_subset(csv.list , latest.date)

for(i in seq_along(latest.data)) {
  
  obj.name <- switch(str_extract(latest.data[i], '[a-z]+') , 
                         'namelist' = 'name.list' , 
                         'pecktime' = 'peck.time' , 
                         'timebar'  = 'timebar')
  
  assign(obj.name , fread(paste0('data//',latest.data[i]) ,  encoding = 'UTF-8' ,colClasses = 'character' ))

  }


setnames(setDT(name.list) , c("event" , paste0('group' , 1: (length(name.list) -1) ) ) )


# remove trim
cols = names(name.list)
name.list <- name.list[,lapply(.SD, str_trim, side = "both"), .SDcols = cols]
setkey(name.list , event)

peck.time <- peck.time[,lapply(.SD, str_trim , side = 'both') ]
setnames(peck.time , c('event' ,'peckt' , 'regular' , 'traditional' , 'main.issue'  ,'time.scale' , 'china' , 'labor',
                       'gender','environment' ,'land'))
peck.time[,peckt := peckt %>>% as.Date( format='%Y.%m.%d')]



## include regular event or not ----------------------------------------------------------#
event.sel  = function(event = 'all'){
  events <<- event
  if(event == 'all') {
    print("do nothing")
  }else if(event == 'nonregular'){
    name.list <<- name.list[!name.list$event %>>% str_detect('[0-9]{2}')]
    peck.time <<- peck.time[!peck.time$event %>>% str_detect('[0-9]{2}')]
    cat('non_regular event only')
  }else if(event == 'regular') {
    name.list <<- name.list[name.list$event %>>% str_detect('[0-9]{2}')]
    peck.time <<- peck.time[peck.time$event %>>% str_detect('[0-9]{2}')]
    cat('regular only')
    
  }
}

event.sel(event = 'all')

  
#--------------------------------------------------------------------------------------###

##----#


netowrkdt.dur <- function(start.t = '20120101' , end.t = '20181231') {
  start.t <- ymd(start.t ) 
  end.t   <- ymd(end.t ) 
  
  # event list & group list
  event.dur.t <- peck.time[peckt <= end.t & peckt >= start.t]
  
  name.list.t <- name.list[event %in% event.dur.t$event]
  
  group.list.t <- unlist(name.list.t[,-1] ,use.names = F )
  group.list.t <- unique(group.list.t[str_length(group.list.t)>0])
  
  all.list <- c(event.dur.t$event , group.list.t)
  
  # edge list dt
  all.edge <- list()
  for(i in 2:ncol(name.list.t)){
    all.edge[[i]] <-  name.list.t[,c(1,i) , with = F]
  } ; edge.list <- rbindlist(all.edge)
  
  edge.list <- edge.list[str_length(group1)>0,]
  edge.list <- unique(edge.list[,c("event", "group1")])  # drop duplicate
  all.t <- graph_from_edgelist(as.matrix(edge.list) , directed = F)
  
  if(length(all.list[! all.list %in% V(all.t)$name]) > 0 ){
    all.t <- all.t +   all.list[! all.list %in% V(all.t)$name]
  }
  V(all.t)$type <- bipartite_mapping(all.t)$type 
  
  # drop duplicate
  all.t <- simplify(all.t)
  
  return(list('igraph' = all.t , 'all.list' = all.list , 'event.list' = event.dur.t$event , 'group.list' = group.list.t , 'edge.list' = edge.list))
} 

rawlist2012 <- netowrkdt.dur(start.t = '20120101' , '20121231')
rawlist2013 <- netowrkdt.dur(start.t = '20130101' , '20131231')
rawlist2014 <- netowrkdt.dur(start.t = '20140101' , '20141231')
rawlist2015 <- netowrkdt.dur(start.t = '20150101' , '20151231')
rawlist2016 <- netowrkdt.dur(start.t = '20160101' , '20161231')
rawlist2017 <- netowrkdt.dur(start.t = '20170101' , '20171231')
rawlist2018 <- netowrkdt.dur(start.t = '20180101' , '20181231')
rawlistt1v1 <- netowrkdt.dur(start.t = '20120101' , '20140318')
rawlistt2v1 <- netowrkdt.dur(start.t = '20140319' , '20160501')
rawlistt3v1 <- netowrkdt.dur(start.t = '20160502' , '20181231')

net.2012 <- rawlist2012$igraph
net.2013 <- rawlist2013$igraph
net.2014 <- rawlist2014$igraph
net.2015 <- rawlist2015$igraph
net.2016 <- rawlist2016$igraph
net.2017 <- rawlist2017$igraph
net.2018 <- rawlist2018$igraph
net.t1v1 = rawlistt1v1$igraph
net.t2v1 = rawlistt2v1$igraph
net.t3v1 = rawlistt3v1$igraph

rawlist2012.2018 <- netowrkdt.dur()
net.2012.2018 <- rawlist2012.2018$igraph

    
    
    

#-- calculate activity time--------------------------------------------------------------------####
setkey(rawlist2012.2018$edge.list , 'event')
raw.groupt <- rawlist2012.2018$edge.list[peck.time][order(peckt)][,.SD[c(1,.N)] ,by = 'group1'][!is.na(group1)][,.(group1 , peckt)]
raw.groupt[,dur := c('start','end')]
group.act.t <- dcast.data.table(raw.groupt , group1~dur ,value.var = 'peckt')

group.event.time <- merge( peck.time  ,group.act.t, by.x = 'event' , by.y = 'group1' ,all =T )
setnames(group.event.time , 'event' ,'name')
group.event.time[,comb:= ifelse(!is.na(peckt) , as.character(peckt) , paste0(start, ' to ' ,end))]
#=====================================================================================#


#---------------------------------network character  --------------------------------####
net.char <- function(input.g , type = 'all' , cohesive.block = FALSE){
  
  # individual char
  output.vertex <- switch(type , 
                          'event' = !V(input.g)$type ,
                          'group' = V(input.g)$type ,
                          'all'   =  rep( TRUE, V(input.g) %>>% length))
  
  ver.list <- V(input.g)$name[output.vertex]
  
  
  type.v <- ifelse(V(input.g)$type == F , 'event' , 'group')
  type.v <- type.v[output.vertex]
  
  act.time = group.event.time[match(ver.list , group.event.time$name)]$comb
  
  degree.v <-  degree(input.g)[output.vertex]
  betweenness.v <- round(betweenness(input.g,directed = F,normalized = T)[output.vertex]*100,4) # rescale 0-100
  constraint.v <- round((1-constraint(input.g)[output.vertex])*100,4) # reverse & rescale 0-100 愈高代表愈居於結構洞位置
  eigen.cen.v <- round(eigen_centrality(input.g,directed = F)$vector[output.vertex] * 100 ,4) # rescale 0-100
  
  bicomponent_list <- biconnected_components(input.g)
  largest.bicomp <- bicomponent_list$components[[which.max(sapply(bicomponent_list$components, length)) ]]$name
  largest.bicomp <-  ifelse(ver.list %in% largest.bicomp ,1,0)
  
  cutpoint <- ifelse(ver.list %in% igraph::articulation.points(input.g)$name ,1,0)
  
  # community
  edge.clust <- cluster_edge_betweenness (input.g ,directed = F  )
  edge.clust <- edge.clust$membership[output.vertex]
  szieof.edge.com = edge.clust %>>% table()
  szieof.edge.com = szieof.edge.com[match(edge.clust , names(szieof.edge.com))] %>>% as.numeric()
  
  walktrap.com <- walktrap.community(input.g)$membership[output.vertex]
  szieof.walktrap.com = walktrap.com %>>% table()
  szieof.walktrap.com = szieof.walktrap.com[match(walktrap.com , names(szieof.walktrap.com))] %>>% as.numeric()
  
  louvain.com <- cluster_louvain(input.g)$membership[output.vertex]
  szieof.louvain.com = louvain.com %>>% table()
  szieof.louvain.com = szieof.louvain.com[match(louvain.com , names(szieof.louvain.com))] %>>% as.numeric()
  
  # time 
  # temp = group.event.time[name %in% ver.list]
  # temp[match(name , ver.list)]
  timedt <- group.event.time[match(ver.list , name ,nomatch = 0)][,.(peckt, regular ,traditional, main.issue,time.scale,china,labor,gender,
                                                                  environment ,land)]
  timedt <- timedt[,lapply(.SD , as.character)]
  # timedt <- group.event.time[match(ver.list , name ,nomatch = 0)]$comb
  # timedt = timedt %>>% as.character()
  
  #embeddedness   low performance
  if(cohesive.block == TRUE){
    gBlocks <- igraph::cohesive.blocks(input.g)
    embeddedness <- igraph::max_cohesion(gBlocks)[match(gBlocks$labels , ver.list , nomatch = 0)]
    ind.char  <- as.data.table(cbind(ver.list ,type.v , act.time , degree.v ,betweenness.v , constraint.v ,
                                     eigen.cen.v ,largest.bicomp ,cutpoint ,embeddedness , edge.clust,
                                     szieof.edge.com , walktrap.com,szieof.walktrap.com, louvain.com , szieof.louvain.com,
                                     timedt))
  }else if(cohesive.block == FALSE) {
    ind.char  <- as.data.table(cbind(ver.list ,type.v , act.time , degree.v ,betweenness.v , constraint.v , 
                                     eigen.cen.v ,largest.bicomp ,cutpoint ,edge.clust, szieof.edge.com,
                                     walktrap.com ,szieof.walktrap.com , louvain.com , szieof.louvain.com,timedt))
  }
  num.var <- lapply( ind.char, function(x) all(str_detect(x , '^[0-9]+$')==TRUE) )
  num.var <- names(num.var[num.var==T & !is.na(num.var)])
  ind.char[, (num.var) := lapply(.SD ,function(x)  round(as.numeric(x),3) ) , .SDcols = num.var]
  return(ind.char)

}

# eventchar2012.dt <- net.char(net.2012 , cohesive.block = T)
# groupchar2012.dt <- net.char(net.2012 ,type = 'group' , cohesive.block = T)
allver2012.dt <- net.char(net.2012 ,type = 'all' , cohesive.block = T)
allver2013.dt <- net.char(net.2013 ,type = 'all' , cohesive.block = T)
allver2014.dt <- net.char(net.2014 ,type = 'all' , cohesive.block = T)
allver2015.dt <- net.char(net.2015 ,type = 'all' , cohesive.block = T)
allver2016.dt <- net.char(net.2016 ,type = 'all' , cohesive.block = T)
allver2017.dt <- net.char(net.2017 ,type = 'all' , cohesive.block = T)
allver2018.dt <- net.char(net.2018 ,type = 'all' , cohesive.block = T)
allvert1v1.dt <- net.char(net.t1v1 ,type = 'all' , cohesive.block = T)
allvert2v1.dt <- net.char(net.t2v1 ,type = 'all' , cohesive.block = T)
allvert3v1.dt <- net.char(net.t3v1 ,type = 'all' , cohesive.block = T)

allver2012.2018.dt <- net.char(net.2012.2018 ,type = 'all' , cohesive.block = T)



# largest k component ---------------------------------------------------------------####
# source('G://我的雲端硬碟//drive//180302.for.learning//function//k.component.R')

# # k connected matrix
# netdt.list = ls() %>>% str_subset('net\\.') %>>% str_subset('[0-9]')
# netdt.list = netdt.list[netdt.list  %>>% str_length() < 10]
# k.connected.matrix = list()
# 
# cat( ' start time : ', as.character(st.time) ,'\n' )
# ptm <- proc.time()
# for(i in seq_along(netdt.list)){
#   print(netdt.list[i])
#   k.connected.matrix[[i]] = k.connected(get(netdt.list[i]))
# }
# 
# cat(' end time : ' ,as.character (Sys.time() ) ,'\n' ,
#     'time used : ' ,proc.time() - ptm , '\n')
# 
# ### 
# 
# largesit.kcomp = function(input.g) {
#   
#   st.time = Sys.time()
#   cat( ' start time : ', as.character(st.time) ,'\n' )
#   
#   mat = k.connected.matrix[[which(str_detect(netdt.list , input.g ))]]
#   
#   all.ver = V(get(input.g))$name 
#   comp1 = membership.k.component(mat , 1)
#   cat(' progress : 20%' , '\n')
#   comp2 = membership.k.component(mat , 2)
#   cat(' progress : 40%'  , '\n')
#   comp3 = membership.k.component(mat , 3)
#   cat(' progress : 60%' , '\n' )
#   comp4 = membership.k.component(mat , 4)
#   cat(' progress : 80%'  , '\n')
#   comp5 = membership.k.component(mat , 5)
#   cat(' progress : 100%'  , '\n')
#   
#   temp  = c(all.ver , comp1 , comp2 , comp3 ,comp4 ,comp5)
#   cat(' end time : ' ,as.character (Sys.time() ) ,'\n' ,
#       'time used : ' ,Sys.time() - st.time , '\n')
#   
#   return((temp %>>% table) - 1)
# }
# 
# 
# lcomp.n2012 = largesit.kcomp('net.2012')
# lcomp.n2013 = largesit.kcomp('net.2013')
# lcomp.n2014 = largesit.kcomp('net.2014')
# lcomp.n2015 = largesit.kcomp('net.2015')
# lcomp.n2016 = largesit.kcomp('net.2016')
# lcomp.n2017 = largesit.kcomp('net.2017')
# lcomp.n2018 = largesit.kcomp('net.2018')
# lcomp.nt1v1 = largesit.kcomp('net.t1v1')
# lcomp.nt2v1 = largesit.kcomp('net.t2v1')
# lcomp.nt3v1 = largesit.kcomp('net.t3v1')
# 
# 
# 
# ## sort
tmp.list = ls() %>>% str_subset('lcomp') %>>% str_sub(8,11)
# 
# for(i in 1:length(ls() %>>% str_subset('lcomp')) ) {
#    
#   
#   temp = get((ls() %>>% str_subset('lcomp'))[i] )
#   temp2 = get((ls() %>>% str_subset('^.{10,15}$') %>>% str_subset('allver'))[i])
# 
#   assign(paste0('lcomp.n' , tmp.list[i]) , temp[match(temp2$ver.list ,names(temp) )])
#   
# }
# 
# 
# allver2012.dt = cbind(allver2012.dt , 'largest.k.comp' = as.numeric(lcomp.n2012))
# allver2013.dt = cbind(allver2013.dt , 'largest.k.comp' = as.numeric(lcomp.n2013))
# allver2014.dt = cbind(allver2014.dt , 'largest.k.comp' = as.numeric(lcomp.n2014))
# allver2015.dt = cbind(allver2015.dt , 'largest.k.comp' = as.numeric(lcomp.n2015))
# allver2016.dt = cbind(allver2016.dt , 'largest.k.comp' = as.numeric(lcomp.n2016))
# allver2017.dt = cbind(allver2017.dt , 'largest.k.comp' = as.numeric(lcomp.n2017))
# allver2018.dt = cbind(allver2018.dt , 'largest.k.comp' = as.numeric(lcomp.n2018))
# allvert1v1.dt = cbind(allvert1v1.dt , 'largest.k.comp' = as.numeric(lcomp.nt1v1))
# allvert2v1.dt = cbind(allvert2v1.dt , 'largest.k.comp' = as.numeric(lcomp.nt2v1))
# allvert3v1.dt = cbind(allvert3v1.dt , 'largest.k.comp' = as.numeric(lcomp.nt3v1))


#-----------------------------------------------save data----------------------------------------#
tmp.list2 = c("2012","2013" ,"2014"  ,"2015","2016","2017" , "2018","t1v1" ,'t2v1' , "t3v1","2012.2018")

for(i in seq_along(tmp.list2) ) {
  
  dir = paste0(getwd() , '//data//clean.dt//'  ,events ,'//')
           
  temp1 = get(paste0('allver',tmp.list2[i] , '.dt'))
  temp2 = get(paste0('rawlist',tmp.list2[i] ))
  
  cat( dir , '\n' , paste0(' allver',tmp.list2[i] , '.dt') , '\n',
       paste0('rawlist',tmp.list2[i]) , '\n\n\n')  

  saveRDS(temp1,  paste0(dir ,'allver',tmp.list2[i] , '.dt' , '.rds'))
  saveRDS(temp2,  paste0(dir ,'rawlist',tmp.list2[i]  , '.rds'))

}

# fwrite(allver2012.2018.dt , 'data//clean.dt//allver2012.2018.dt.csv')

# slide windows --------------------------------------------#
# time length 6 month; speed 1 month
year = c(rep(2012,12),rep(2013,12),rep(2014,12),rep(2015,12),rep(2016,12),
         rep(2017,7))
x = c(1:12,1:12,1:12,1:12,1:12,1:7)
month = ifelse(str_count(x)==1, paste0('0',x) , x)
inpute.date = paste0(as.character(year), month , '01')
for(i in seq_along(inpute.date)){
  # time period
  start.date = ymd(inpute.date[i]) 
  end.date =  start.date %m+% months(12)-1
  start.char = as.character(start.date) %>>% str_remove_all('-')
  end.char =  as.character(end.date) %>>% str_remove_all('-')
  
  #assign
  assign(paste0('rawlist' ,start.char,'.',end.char),netowrkdt.dur(start.t = start.char , end.t = end.char))
  assign(paste0('net.',start.char,'.',end.char) ,
         get(paste0('rawlist' ,start.char,'.',end.char))$igraph
  )
  assign(paste0('allver' , start.char ,'.' , end.char , '.dt') ,
         net.char(get(paste0('net.',start.char,'.',end.char)),type = 'all' , cohesive.block = T )
         )
  
  
  #save
  dir = paste0(getwd() , '//data//clean.dt//'  ,events ,'//')
  cat( dir , '\n' , paste0('allver' , start.char ,'.' , end.char , '.dt') , '\n',
       paste0('rawlist' ,start.char,'.',end.char) , '\n\n\n')  
  
  saveRDS(get(paste0('allver' , start.char ,'.' , end.char , '.dt')),  
          paste0( dir , 'allver' , start.char ,'.' , end.char , '.dt' , '.rds'))
  saveRDS(get(paste0('rawlist' ,start.char,'.',end.char)),  
          paste0( dir ,'rawlist' ,start.char,'.',end.char, '.rds'))
  
}

#------------------------------------------------------------


#----------------------------------------------------------------------------#

#------------- which kind of data you want to run --------------------------------#
# # 
# read.data = function(events = 'all'){
# 
#   dir = paste0(getwd() , '//data//clean.dt//'  ,events ,'//')
# 
#   dir.file.list = list.files(dir) %>>% str_subset("rds")
#   cat('data directory : ' , dir , '1\n\n')
# 
#   for(i in seq_along(dir.file.list)) {
# 
#     assign( (dir.file.list[i] %>>% str_remove('.rds') ),
#             readRDS(paste0(dir , dir.file.list[i] ) )  , envir =  .GlobalEnv)
#     cat('read successful ' , '\n' ,
#         dir.file.list[i] , ' as ' , dir.file.list[i] %>>% str_remove('.rds') ,
#         '\n')
#   }
# 
# 
# }
# 
# events = 'all'
# read.data('all')
###


# 
# degree.10 <- quantile(allver2012.2018.dt$degree.v, probs = c(.9) )
# 
# #--------------------------------------------- igraph plot ----------------------------####
# igrapg.plot <- function(net.in , char.in){
#   V(net.in)$color <- ifelse(V(net.in)$type, "lightblue", "salmon")
#   V(net.in)$shape <- ifelse(V(net.in)$type, "circle", "square")
#   V(net.in)$size      <- ifelse(V(net.in)$type,  0.3 , 1.5)
#   # V(net.2016)$bicom    <- allver2016.$largest.bicomp
#   # groups(biconnected.components(net.2016)$components[[145]])
#   # vertenet.in.label.cex.in      <- ifelse(char.in$degree.v >= quantile(char.in$degree.v, probs = c(.95) ), 2 , 1)
#   V(net.in)$label.cex <- ifelse(char.in$degree.v >= degree.10 , 1.5, .2)
#   plot(cohesive.blocks(net.in ) , net.in, col = char.in$edge.clust,  vertex.label.color = "black",
#        vertex.size = 1)
#   
# }
# # 2012~2018 network each year
# igrapg.plot(net.2012 , allver2012.dt)
# igrapg.plot(net.2013 , allver2013.dt)
# igrapg.plot(net.2014 , allver2014.dt)
# igrapg.plot(net.2015 , allver2015.dt)
# igrapg.plot(net.2016 , allver2016.dt)
# igrapg.plot(net.2017 , allver2017.dt)
# igrapg.plot(net.2018 , allver2018.dt)


#------------------------------------------------------------------------------------------

switch.embtocol <- Vectorize(function(x) {
  switch(x ,
        '1'= brewer.pal(n =11, name = "Spectral")[7],
        '2'= brewer.pal(n =11, name = "Spectral")[9],
        '3'= brewer.pal(n =11, name = "Spectral")[11]) 
  }
)
    
switch.comptocol <- Vectorize(function(x) {
  switch(x ,
         '0'= brewer.pal(n =11, name = "Spectral")[6],
         '1'= brewer.pal(n =11, name = "Spectral")[5],
         '2'= brewer.pal(n =11, name = "Spectral")[4],
         '3'= brewer.pal(n =11, name = "Spectral")[3],
         '4'= brewer.pal(n =11, name = "Spectral")[2],
         '5'= brewer.pal(n =11, name = "Spectral")[1]
         ) 
}
)            

switch.communitycol = Vectorize(function(x) {
  switch(x , 
         '1' =  brewer.pal(n =11, name = "Spectral")[1],
         '2' =  brewer.pal(n =11, name = 'Spectral')[2],
         '3' =  brewer.pal(n =11, name = 'Spectral')[3],
         '4' =  brewer.pal(n =11, name = 'Spectral')[4],
         '5' =  brewer.pal(n =11, name = 'Spectral')[5],
         '6' =  brewer.pal(n =11, name = 'Spectral')[6],
         '7' =  brewer.pal(n =11, name = 'Spectral')[7],
         '8' =  brewer.pal(n =11, name = 'Spectral')[8],
          brewer.pal(n =11, name = 'Spectral')[9]
  )
}
)

#----------------------------------- visnetwork ---------------------------------------####
library(visNetwork)
# temp <- allver2012.dt[match(rawlist2012$all.list,allver2012.dt$ver.list )]
# 
# temp[,embeddedness.col := switch.embtocol(as.character(embeddedness))]
# 
# nodes <- data.table(id=rawlist2012$all.list , label=rawlist2012$all.list,
#                     group = c(rep('event',length(rawlist2012$event.list))  ,rep('group',length(rawlist2012$group.list))),
#                     shadow = c(rep(TRUE,length(rawlist2012$event.list))  ,rep(FALSE,length(rawlist2012$group.list))),
#                     value = c(rep(.5,length(rawlist2012$event.list))  ,rep(.01,length(rawlist2012$group.list))),
#                     color = temp$embeddedness.col,
#                     title =  paste0(
#                       "<p>Name: ", temp$ver.list , "<br>",
#                       "Type: ", temp$type , '<br>' ,
#                       "Date: ", temp$comb , '<br>' ,
#                       "Degree: ", temp$degree.v , "<br>",
#                       "walktrap community: ", temp$walktrap.com , '<br>' ,
#                       "Betweenness: ", temp$betweenness.v , '<br>' ,
#                       "Constraint: ", temp$constraint.v ,'<br>',
#                       "Eigenvector: ", temp$eigen.cen.v ,'<br>',
#                       "Largest BiComponent: ", ifelse(temp$largest.bicomp==1 ,'Yes' ,'No') ,'<br>',
#                       "Embeddedness: ", temp$embeddedness ,'<br>'
#                     )
#                     )
# edges <- data.table(from = rawlist2012$edge.list[[1]] , to = rawlist2012$edge.list[[2]])
# 
# visNetwork(nodes, edges, width = "100%" , height = "800px") %>>% visPhysics(stabilization = FALSE) %>>%
#   visEdges(smooth = FALSE) %>>%
#   visNodes(scaling = list(label = list(enabled = T))) %>>% visEdges(smooth = FALSE) %>>%
#   visIgraphLayout()

plot.visnetwork <- function(input.rawlist = rawlist2012 , input.char = allver2012.dt) {
  temp <- input.char[match(input.rawlist$all.list,input.char$ver.list )]
  # temp[,embeddedness.col := switch.comptocol(as.character(embeddedness))]
  # temp[,larg.kcomp.col := switch.comptocol(as.character(largest.k.comp))]
  temp[,com.col := switch.communitycol(as.character(edge.clust))]
  
  
  nodes <- data.table(id=input.rawlist$all.list , label=input.rawlist$all.list, 
                      shadow = c(rep(TRUE,length(input.rawlist$event.list))  ,rep(FALSE,length(input.rawlist$group.list))),
                      value = c(rep(.5,length(input.rawlist$event.list))  ,rep(.01,length(input.rawlist$group.list))),
                      color = temp$com.col,
                      title =  paste0(
                        "<p>Name: ", temp$ver.list , "<br>",
                        "Type: ", temp$type , '<br>' ,
                        "Date: ", temp$act.time , '<br>' ,
                        "Degree: ", temp$degree.v , "<br>",
                        "Community: ", temp$edge.clust , '<br>' ,
                        "Size of Community: ", temp$szieof.edge.com , '<br>' ,
                        'Tag China: ' ,ifelse(temp$china == 1 ,'Yes' ,'No') , '<br>' ,
                        'Tag Labor: ' , ifelse(temp$labor == 1 ,'Yes' ,'No')  , '<br>' ,
                        'Tag Gender: ' , ifelse(temp$gender  == 1 ,'Yes' ,'No'), '<br>' ,
                        'Tag Environment: ' , ifelse(temp$environment == 1 ,'Yes' ,'No'), '<br>',
                        'Tag Land: ' , ifelse(temp$land == 1 ,'Yes' ,'No') , '<br>' , 
                        'Tag Traditional' , ifelse(temp$traditional == 1 ,'Yes' ,'No') , '<br>' ,
                        "Betweenness: ", temp$betweenness.v , '<br>' ,
                        "Constraint: ", temp$constraint.v ,'<br>',
                        "Eigenvector: ", temp$eigen.cen.v ,'<br>',
                        'Largest bicomponent: ' , ifelse(temp$largest.bicomp == 1 ,'Yes' ,'No') , '<br>' ,
                        # "Largest K Component: ", temp$largest.k.comp ,'<br>',
                        # "Embeddedness: ", temp$embeddedness ,'<br>',
                        "Main issue: " , temp$main.issue ,'<br>' ,
                        "Time Scale: " , temp$time.scale ,'<br>'  
                      )
  )
  edges <- data.table(from = input.rawlist$edge.list[[1]] , to = input.rawlist$edge.list[[2]])
  
  visNetwork(nodes, edges, width = "100%" , height = "800px") %>>% visPhysics(stabilization = FALSE) %>>%
    visEdges(smooth = FALSE) %>>% 
    visNodes(scaling = list(label = list(enabled = T))) %>>% visEdges(smooth = FALSE) %>>% 
    visIgraphLayout()
}

# export html
for(i in seq_along(tmp.list)){
  dir = paste0(getwd() , '//plot//'  ,events ,'//')
  
  dttime = tmp.list[i]
  
  temp1 = get((ls() %>>% str_subset(paste0('allver' , dttime ,'.dt'))))
  temp2 = get((ls() %>>% str_subset( paste0('rawlist' , dttime , '$' ) ) ) )
  
  visSave(plot.visnetwork(temp2 , temp1)  ,file = paste0(dir , 'net',dttime ,'.html') )
  cat(dir , '\n',dttime , '\n')
}

# windows data type export
date.list = str_subset(ls() , 'rawlist')[(str_subset(ls() , 'rawlist') %>>% str_count() >20)] %>>%
  str_sub(8,24)

for( i in seq_along(date.list)) {
  dir = paste0(getwd() , '//plot//'  ,events ,'//')
  dttime = date.list[i]
  
  temp1 = get(paste0('allver' , dttime ,'.dt'))
  temp2 = get(paste0('rawlist' , dttime ))

  visSave(plot.visnetwork(temp2 , temp1)  ,file = paste0(dir , 'net',dttime ,'.html') )
  cat(dir , '\n',dttime , '\n')
  
}


# visSave(plot.visnetwork(rawlist2012 , allver2012.dt)  ,file = "net2012.html")
# visSave(plot.visnetwork(rawlist2013 , allver2013.dt)  ,file = "net2013.html")
# visSave(plot.visnetwork(rawlist2014 , allver2014.dt)  ,file = "net2014.html")
# visSave(plot.visnetwork(rawlist2015 , allver2015.dt)  ,file = "net2015.html")
# visSave(plot.visnetwork(rawlist2016 , allver2016.dt)  ,file = "net2016.html")
# visSave(plot.visnetwork(rawlist2017 , allver2017.dt)  ,file = "net2017.html")
# visSave(plot.visnetwork(rawlist2018 , allver2018.dt)  ,file = "net2018.html")
