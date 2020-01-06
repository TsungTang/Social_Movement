library(data.table)
library(stringr)
library(igraph)
library(pipeR)
library(lubridate)
library(RColorBrewer)
library(readxl)

group.past.dt = list()
# 每一時間點觀察組織過去參加的事件狀況
for(z in seq_along(rawlist.name)){  # 從第一個時間點到最後時間點
  posi1 = z
  tmp.grouplist = get(rawlist.name[[z]])$group.list
  tmpgraph = get(rawlist.name[posi1])$igraph
  
  # community
  tmpgraph2 = bipartite.projection(tmpgraph)$proj2
  c2.g = cluster_fast_greedy(tmpgraph2, weights =E(tmpgraph2)$weight )
  c2.g2 = multilevel.community(tmpgraph2, weights =E(tmpgraph2)$weight )
  c2.g3 = walktrap.community(tmpgraph2, weights =E(tmpgraph2)$weight)
  c2.g4 = leading.eigenvector.community(tmpgraph2, weights =E(tmpgraph2)$weight)
  
  tmp.dt = data.table('ver.list' = c2.g$names, 'fg.com' = c2.g$membership ,'ml.com'= c2.g2$membership ,'wt.com' =c2.g3$membership ,
                      'eig.com'=c2.g4$membership )
  tmp.dt[,fg.com.N:=.N , by ='fg.com']
  tmp.dt[,ml.com.N:=.N , by ='ml.com']
  tmp.dt[,wt.com.N:=.N , by ='wt.com']
  tmp.dt[,eig.com.N:=.N , by ='eig.com']
  ##
  
  # past event##
  group.past = list()
  for(ix in seq_along(tmp.grouplist)){
    test.g = make_ego_graph(tmpgraph , nodes = tmp.grouplist[ix])
    
    group.past[[ix]] = get(allver.name[posi1])[ver.list %in% V(test.g[[1]])$name & type.v=='event'][  #ego 過去參加的事件類型
      ,list(tmp.grouplist[ix]  ,posi1 ,max(degree.v) ,min(degree.v) , round(mean(degree.v),3) ,.N,sum(as.numeric(regular)) ,
            sum(as.numeric(traditional)) ,sum(as.numeric(china)),sum(as.numeric(labor)),
            sum(as.numeric(gender)) ,sum(as.numeric(environment)),sum(as.numeric(land))  )]
  }
  group.past = rbindlist(group.past)
  #---#
  
  tmp.dt = cbind(tmp.dt , group.past[,-1] )
  group.past.dt[[z]] = tmp.dt
  cat(rawlist.name[z] , '\n')
}

group.past.dt = rbindlist(group.past.dt)

col.vec = c('time' ,'max.event.size' , 'min.event.size' , 'mean.event.size','pati.n','pati.regular' ,'pati.trad' ,'pati.china' ,'pati.labor','pati.gender','pati.envi','pati.land')
setnames(group.past.dt ,names(group.past.dt)[10:21], col.vec)

# fwrite(group.past.dt , 'data/clean.dt/group.past.dt.csv')
#  1-16  ;  17-42 ; 43-67
five.att = c('pati.china' , 'pati.labor' ,'pati.gender' ,'pati.envi' ,'pati.land')
comm.dt = list()
for (z in five.att) {
  count.list = list()
  for (i in 1:67) {
    test = group.past.dt[time==i & get(z)>1 , fg.com]
    count.table = sort(table(test) , decreasing = T)  
    count.table  = as.vector(count.table)
    
    if(length(count.table) > 3 ){
      count.table = c(count.table[1:3] ,
                      sum(count.table[4:length(count.table)]))
    }
    
    if(length(count.table) <= 3 ){
      add.0=4 - length(count.table)
      count.table = c(count.table  , rep(0 , add.0 ))
    }
    count.list[[i]] = count.table
  }
  count.list = do.call(rbind , count.list)
  comm.dt[[z]] = count.list
}

each.att.com = list()
for (i in five.att) {
  tmp.mat1 =apply(comm.dt[[i]][1:16,] , 2 , sum) 
  tmp.mat2 =apply(comm.dt[[i]][17:42,] , 2 , sum) 
  tmp.mat3 =apply(comm.dt[[i]][43:67,] , 2 , sum) 
  tmp.mat =apply(comm.dt[[i]] , 2 , sum) 
  
  each.att.com[[i]] = cbind(rep(i ,4) , c(1,2,3,'ALL'),
                            rbind(paste0(tmp.mat1 , '(',round(tmp.mat1  /sum(tmp.mat1),2) ,')' )  ,
                                  paste0(tmp.mat2 , '(',round(tmp.mat2  /sum(tmp.mat2),2) ,')' ) ,
                                  paste0(tmp.mat3 , '(',round(tmp.mat3  /sum(tmp.mat3),2) ,')' )  ,
                                  paste0(tmp.mat , '(',round(tmp.mat  /sum(tmp.mat),2) ,')' ) )
                            )

}
each.att.com = do.call(rbind , each.att.com )
write.csv(each.att.com ,file = 'plot/暫定採用圖表/each.att.com.csv')

# round(tmp.mat/sum(tmp.mat),3)

round(comm.dt$pati.china[1:16,]/rowSums(comm.dt$pati.china[1:16,]),3)

tmp.mat = apply( comm.dt$pati.china, 2 , sum) 
tmp.mat/sum(tmp.mat)

library(plotly)

y <- c('giraffes', 'orangutans', 'monkeys')
SF_Zoo <- c(20, 14, 23)
LA_Zoo <- c(12, 18, 29)
data <- data.frame(y, SF_Zoo, LA_Zoo)

p <- plot_ly(data, x = ~SF_Zoo, y = ~y, type = 'bar', orientation = 'h', name = 'SF Zoo',
             marker = list(color = 'rgba(246, 78, 139, 0.6)',
                           line = list(color = 'rgba(246, 78, 139, 1.0)',
                                       width = 3))) %>%
  add_trace(x = ~LA_Zoo, name = 'LA Zoo',
            marker = list(color = 'rgba(58, 71, 80, 0.6)',
                          line = list(color = 'rgba(58, 71, 80, 1.0)',
                                      width = 3))) %>%
  layout(barmode = 'stack',
         xaxis = list(title = ""),
         yaxis = list(title =""))
