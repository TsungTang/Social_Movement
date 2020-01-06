
library(data.table)
library(stringr)
library(igraph)
library(pipeR)
library(lubridate)
library(RColorBrewer)
library(readxl)
#------------- which kind of data you want to run --------------------------------#

read.data = function(events = 'all'){
  
  dir = paste0(getwd() , '//data//clean.dt//'  ,events ,'//')
  
  dir.file.list = list.files(dir) %>>% str_subset("rds")
  cat('data directory : ' , dir , '\n\n')
  
  for(i in seq_along(dir.file.list)) {
    
    
    assign( (dir.file.list[i] %>>% str_remove('.rds') ), 
            readRDS(paste0(dir , dir.file.list[i] ) )  , envir =  .GlobalEnv)
    cat('read successful ' , '\n' , 
        dir.file.list[i] , ' as ' , dir.file.list[i] %>>% str_remove('.rds') ,
        '\n')
  }
  
  
}


read.data('all')


allver2012.2018.dt[type.v=='event' , act.time.date:=as.Date(act.time) ]

rawlist.name = str_subset(ls()  , 'rawlist')[str_subset(ls()  , 'rawlist') %>>% str_count() > 20]
allver.name = str_subset(ls()  , 'allver')[str_subset(ls()  , 'allver') %>>% str_count() > 20]

rawlist.name2 = str_subset(ls()  , 'rawlistt')[str_subset(ls()  , 'rawlistt') %>>% str_count() > 10]
allver.name2 = str_subset(ls()  , 'allvert')[str_subset(ls()  , 'allvert') %>>% str_count() > 10]

alllist2 = list()
for(g in seq_along(rawlist.name2)) {
  # per group
  cat('time: ' , g , '\n')
  
  neibor.list = list()
  for( i in get(rawlist.name2[g])$group.list){
    neib.event = neighborhood(get(rawlist.name2[g])$igraph , nodes =i )
    neib.event = names(neib.event[[1]])
    temp = get(allver.name2[g])[ver.list %in% neib.event & type.v == 'event']
    temp[,lapply(.SD , as.numeric) , .SDcols = c('degree.v', 'regular','traditional' , 
                                                 'china', 'labor', 'gender', 'environment', 'land')]
    temp = temp[,lapply(.SD , as.numeric) , .SDcols = c('degree.v' ,'edge.clust', 'szieof.edge.com', 'walktrap.com', 
                                                        'szieof.walktrap.com' , 'regular','traditional' , 
                                                        'china', 'labor', 'gender', 'environment', 'land')]
    temp = temp[,.('ver.list' = i , 'n.edg.com' = length(unique(edge.clust)) , 'mean.edg.comsize'=  mean(szieof.edge.com), 'n.walk.com' = length(unique(walktrap.com)) ,
                   'mean.walk.comsize'  =  mean(szieof.walktrap.com)  , 'mean.event.degree' = mean(degree.v) ,
                   'n.reg'  = sum(regular) , 'n.trad' = sum(traditional), 'n.china'= sum(china) , 'n.labor' =sum(labor) ,
                   'n.enviro' = sum(environment) , 'n.land' = sum(land) , time = g ) ]
    
    neibor.list[[i]] = temp
  }
  neibor.list = rbindlist(neibor.list)
  
  alllist2[[g]] = merge(get(allver.name2[g])[type.v=='group'] , neibor.list ,by= 'ver.list')
}

alllist = rbindlist(alllist)
alllist[,ver.list.t := paste0(ver.list ,'.t' ,time)]
alllist[,time.scale:=ifelse(time<=16 , 'Time 1' ,
                            ifelse(time < 55 , 'Time 2' , 'Time 3' ) )]

alllist[,embeddedness:= as.ordered(embeddedness)]
alllist[, nested := 1 ][embeddedness==1 & largest.bicomp == 1, nested := 2 ][
  embeddedness>=2 & largest.bicomp == 1, nested := 3 ]

alllist[,nested:=as.ordered(nested)]

# model
library(MASS)
fit1 = polr(nested ~ degree.v + n.edg.com + mean.edg.comsize + n.land + 
              n.labor + n.china + n.enviro  + n.reg + factor(time.scale) , data = alllist, Hess=TRUE)
# summary(fit1)
(ctable <- round(coef(summary(fit1)) , 4)  )
# ## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
# ## combined table
(ctable <- cbind(ctable, "p value" = round(p,4) ))
# 
fit1 = polr(nested ~ degree.v + n.edg.com  + n.land + n.labor + n.china + n.enviro  + n.reg  , 
            data = alllist[time.scale=='Time 3'], Hess=TRUE ,  method = "logistic")
# summary(fit1)
ctable <- coef(summary(fit1))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
## combined table
(ctable <- cbind(ctable, "p value" = round(p,4) ))

#-------------------------

alllist2 = rbindlist(alllist2)
alllist2[,ver.list.t := paste0(ver.list ,'.t' ,time)]

###
g = decompose.graph(net.2012  , min.vertices = 1)
g


outlist = list()
for(g in seq_along(rawlist.name2)) {
  cat('Time: ' , g , '---------------------------' , '\n\n\n\n')
  bc = cohesive_blocks(get(rawlist.name2[g])$igraph)
  allposib.col = cbind(get(rawlist.name2[g])$group.list,get(rawlist.name2[g])$group.list)
  test = expand.grid(get(rawlist.name2[g])$group.list , get(rawlist.name2[g])$group.list , stringsAsFactors = F)
  testg = graph_from_edgelist(as.matrix(test) , directed = F)
  testg = simplify(testg , remove.loops = T)
  test = get.adjedgelist(testg , mode = 'all')
  test  =get.edgelist(testg  )
  test = as.data.table(test)
  
  
  allblock = list()
  for(i in 1:(bc %>>% length) ) {
    allblock[[i]] =  names(bc$blocks[[i]])
  }
  
  cohesion.sim = vector()
  for( i in 1:nrow(test)) {
    cohesion.vec = bc$cohesion[sapply(allblock , function(x) all(test[i,] %>>% unlist %in% x) )]
    cohesion.sim[i] = max(cohesion.vec)
    if( i %% 100 ==0){
      cat((i/nrow(test))*100  , '%', '\n')
    }
  }
  cat('End Time: ' , g , '---------------------------' , '\n\n\n\n')
  outlist[[g]] = cbind(test , cohesion.sim)
  
}

forcomb.dt = alllist2[time==1]
setnames(forcomb.dt , names(forcomb.dt) , paste0('V1.', names(forcomb.dt)))

test=merge( outlist[[1]], forcomb.dt ,
       by.x = 'V1' , by.y = 'V1.ver.list' , all.x = T)
forcomb.dt = alllist2[time==1]
setnames(forcomb.dt , names(forcomb.dt) , paste0('V2.', names(forcomb.dt)))
test=merge( test, forcomb.dt ,
            by.x = 'V2' , by.y = 'V2.ver.list' , all.x = T)
  
test[,cohesion.sim:= as.ordered(cohesion.sim)]

# model ------------------ #
fit1 = polr(cohesion.sim ~ V2.n.china + V1.n.china ,
            data = test, Hess=TRUE)

(ctable <- coef(summary(fit1)))
# ## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
# ## combined table
(ctable <- cbind(ctable, "p value" = round(p,4) ))

fit1 = lm(cohesion.sim ~ V1.degree.v + V2.degree.v +V1.n.edg.com + V2.n.edg.com + V1.mean.edg.comsize  + 
            V2.mean.edg.comsize + V2.n.enviro + V1.n.enviro + V2.n.enviro * V1.n.enviro + V2.n.china + 
            V1.n.china + V2.n.china * V1.n.china + V1.n.labor + V2.n.labor + V1.n.labor*V2.n.labor , data = test)

summary(fit1)



#

names(bc$blocks[2][[1]])
bc$cohesion
names(bc$blocks[2][[1]]) %in% allver20130301.20140228.dt[largest.bicomp==1]$ver.list 
allver20130301.20140228.dt[largest.bicomp==1]$ver.list %in% names(bc$blocks[2][[1]]) 


tmp = allver2012.2018.dt[type.v == 'event' & gender == 1]$ver.list

tmpg = delete_vertices(net.20120101.20121231, V(net.20120101.20121231)$name [V(net.20120101.20121231)$name %in% tmp])

#remove isolate vertex
tmpg = delete_vertices(tmpg , V(tmpg)$name[degree(tmpg)==0])
tmpg.c = cohesive.blocks(tmpg)
tmpg.c$cohesion
tmpg.c$blocks

bipartite.projection(rawlist20120101.20121231$igraph )


allver20120101.20121231.dt
V(rawlist20120101.20121231$igraph)$name

bc = cohesive_blocks(rawlist20130501.20140430$igraph)
bc
bc$cohesion
names(bc$blocks[[7]])
 
bc2 = cohesive_blocks(rawlist20130601.20140531$igraph)
bc2[9]
bc2$cohesion
names(bc2$blocks[[9]])

allver20130501.20140430.dt[embeddedness>1 & type.v=='group']$ver.list
allver20130601.20140531.dt[embeddedness>1 & type.v=='group']$ver.list %in% allver20130501.20140430.dt$ver.list

# match.group 
t.glist = allver20121101.20131031.dt$ver.list
tp1.glist = allver20121201.20131130.dt$ver.list

t.dt = allver20121101.20131031.dt[ver.list%in%tp1.glist]
tp1.dt = allver20121201.20131130.dt[match( t.glist , tp1.glist, nomatch = 0)]

tp1.dt[(tp1.dt$embeddedness - t.dt$embeddedness) != 0 ]

allver20130501.20140430.dt[match( allver20130601.20140531.dt[type.v=='group']$ver.list,allver20130501.20140430.dt[type.v=='group']$ver.list  )]

allver20130601.20140531.dt[type.v=='group' & ver.list %in% allver20130501.20140430.dt$ver.list ,]


allver20130601.20140531.dt[embeddedness>1 & type.v=='group']
allver20141001.20150930.dt[embeddedness>1 & type.v=='group']
allver20150901.20160831.dt[embeddedness>1 & type.v=='group']
allver20160501.20170430.dt[embeddedness>1 & type.v=='group']
