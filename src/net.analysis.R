## 合作關係

library(plotly)
library(tnet)

# explore

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


  # creat random graph with fixed group degree
make.random.graph = function(in.el , fixed = 'fixed.one.edge.group' , size = 30) {  # first column: event 

  in.el = as.matrix(in.el)
  event.list = unique(in.el[,1])
  group.list = unique(in.el[,2])
  allcomb  =expand.grid(event.list,group.list )
  
  g= graph_from_edgelist(in.el)
  
  if(fixed == 'fixed.one.edge.group'){

     # note: All group degree must be greater than or equal to 1
     # numeric vertex name
     event.list.num = 1:length(event.list)
     group.list.num = max(event.list.num)+1:length(group.list)
     
     # length of multiple edge groups & name of multiple edge groups
     multiple.edge.group.length = length(names(table(in.el[,2]))[table(in.el[,2])!=1])
     multipe.edge.group.num = max(event.list.num)+1:multiple.edge.group.length
     # all combination of multiedge group
     multi.allcomb = expand.grid(event.list.num , multipe.edge.group.num)
     random.el = list()
     for( i in 1: size) {
       
       # ONE ARC FOR RANDOM EVENT
    one.arc.el  = cbind(sample( event.list.num ,length(group.list.num) , replace = T ), group.list.num)
    one.arc.el  = as.data.frame(one.arc.el);names(one.arc.el) = c('Var1' ,'Var2')
    
    
  # other edge (for multiple edge groups)

  # all combination of multiedge group minus existed edge(one.arc.el)
  multi.allcomb.2 = dplyr::anti_join(multi.allcomb , one.arc.el , by = c("Var1", "Var2"))
  p.mul = (length(E(g)) - nrow(one.arc.el)) / nrow(multi.allcomb.2) # probability
  sample.row.mul = vector()
  for(q in 1:nrow(multi.allcomb.2)) {
    sample.row.mul[q] = sample(c(T,F) , 1 , replace = T  , prob =c(p.mul,1-p.mul)  )
  }
  multi.allcomb.2 = multi.allcomb.2[sample.row.mul,]
 
  # combine two type edge
  tmp.g  = graph_from_edgelist( as.matrix(rbind(multi.allcomb.2 , one.arc.el)) , directed = F)
  V(tmp.g)$type <- V(tmp.g) %in% event.list.num 
  random.el[[i]] = get.edgelist(tmp.g )
  
     }
     
   #-----------#
     
  }else if(fixed == 'full.random') {
    
    random.netlist = list()
    random.el = list()
    for( i in 1:size){
      random.netlist[[i]] =  sample_bipartite(length(event.list) ,length(group.list), type = 'gnp', p = (E(g) %>>% length) / nrow(allcomb))
      random.el[[i]] = get.edgelist(random.netlist[[i]])
    }
    
    tmp = V(random.netlist[[1]])[V(random.netlist[[1]])$type]
    if(all(random.el[[1]][,2] %in% tmp) == F  ) {
      stop('something wrong in event group column')
    } 
  
  }
  return(random.el)
}
  
  
analysis1 = function(in.el , period = '' ,fixed = 'fixed.one.edge.group', size = 30 , seed = 12221){
  set.seed(seed)
    in.el=as.matrix(in.el)
   
    # random graph with fixed degree
    random.el = make.random.graph(in.el , fixed = fixed , size = size)
    n.event = unique(in.el[,1]) %>>% length()
    n.group = unique(in.el[,2]) %>>% length()
    n.edge.O = nrow(in.el)
    
    # density
    den.O = edge_density(graph_from_edgelist(in.el , F))
    
    tmp = sapply(random.el , function(x) edge_density(graph_from_edgelist(x , F) ) )
    den.R.d = mean(tmp) - sd(tmp) * 2
    den.R.u = mean(tmp) + sd(tmp) * 2
    den.R = mean(tmp)
    den.R.sd = sd(tmp)
    
    # CC
    CC.O = clustering_tm(net =in.el)
    tmp = sapply(random.el , clustering_tm )
    CC.R.d = mean(tmp) - sd(tmp) * 2
    CC.R.u = mean(tmp) + sd(tmp) * 2
    CC.R = mean(tmp)
    CC.R.sd = sd(tmp)
    
    # averge path
    avg.path.O = average.path.length(graph_from_edgelist(in.el ,F) )
    tmp = sapply(random.el , function(x) average.path.length(graph_from_edgelist(x , F) ) )
    avg.path.R.d = mean(tmp) - sd(tmp) * 2
    avg.path.R.u = mean(tmp) + sd(tmp) * 2
    avg.path.R = mean(tmp)
    avg.path.R.sd = sd(tmp)
    
    
    # cohesion
    # k.con.O  = k.connected(graph_from_edgelist(in.el , F))
    # k.con.T  = k.connected(graph_from_edgelist(random.el.T , F))
    # k.con.F  = k.connected(graph_from_edgelist(random.el.F , F))
    # tmp = sapply(random.el , function(x) k.connected(graph_from_edgelist(x , F) ) )
    
    #bi
    bicomponent_list <- biconnected_components(graph_from_edgelist(in.el , F))
    n.bi.o = sapply(bicomponent_list$components, length) %>>% max
    
    tmp = sapply(random.el , function(x) {
      tmp.bicomponent_list <- biconnected_components(graph_from_edgelist(x , F))
      tmp.n.bi = sapply(tmp.bicomponent_list$components, length) %>>% max
      return(tmp.n.bi)
    } )
    bi.R.d = mean(tmp) - sd(tmp) * 2
    bi.R.u = mean(tmp) + sd(tmp) * 2
    bi.R = mean(tmp)
    bi.R.sd = sd(tmp)
    # biratio
    bi.obs.ran.ratio = n.bi.o/bi.R
    
    # modularity
    raw.modularity = modularity(cluster_edge_betweenness(graph_from_edgelist(in.el , F)))
    scale.modularity =  raw.modularity/length(V(graph_from_edgelist(in.el , F)))
    g = graph_from_edgelist(in.el , F)
    V(g)$type <- V(g)$name %in% in.el[,2]
    g = bipartite.projection(g)$proj2
    raw.modularity2 = modularity(cluster_fast_greedy(g))
    scale.modularity2 =  raw.modularity2/length(V(g))
    
    
    c('period'=period, 'n.event'=n.event,'n.group'=n.group,'n.edge.O' =n.edge.O,'den.O'= den.O,
       'den.R.d'= den.R.d ,'den.R.u' = den.R.u, 'den.R'=den.R,'den.R.sd'=den.R.sd,'avg.path.O'=avg.path.O, 
      'avg.path.R.d'=avg.path.R.d,'avg.path.R.u'=avg.path.R.u,'avg.path.R'=avg.path.R,'avg.path.R.sd'=avg.path.R.sd,
      'CC.O'=CC.O, 'CC.R.d'=CC.R.d , 'CC.R.u'=CC.R.u,'CC.R'=CC.R,'CC.R.sd'=CC.R.sd,'n.bi.O'=n.bi.o,
      'bi.R.d'=bi.R.d,'bi.R.u'=bi.R.u,'bi.R'=bi.R,'bi.R.sd'=bi.R.sd,'bi.obs.ran.ratio'=bi.obs.ran.ratio,
      'raw.modularity' = raw.modularity , 'scale.modularity' = scale.modularity ,'raw.modularity2'=raw.modularity2,
      'scale.modularity2' = scale.modularity2
      )
}


  
ana.2012 = analysis1(rawlist2012[[5]] , period='2012' , fixed ='full.random' )
ana.2013 = analysis1(rawlist2013[[5]] , period='2013' , fixed ='full.random')
ana.2014 = analysis1(rawlist2014[[5]] , period='2014' , fixed ='full.random')
ana.2015 = analysis1(rawlist2015[[5]] , period='2015' , fixed ='full.random')
ana.2016 = analysis1(rawlist2016[[5]] , period='2016' , fixed ='full.random')
ana.2017 = analysis1(rawlist2017[[5]] , period='2017' , fixed ='full.random')
ana.2018 = analysis1(rawlist2018[[5]] , period='2018' , fixed ='full.random')
ana.byyear = rbind(ana.2012,ana.2013,ana.2014,ana.2015,ana.2016,ana.2017,ana.2018)
ana.byyear=as.data.table(ana.byyear)
ana.byyear = ana.byyear[,lapply(.SD ,function(x) round(as.numeric(x),3)  )]

#--- windows
rawlist.name = str_subset(ls()  , 'rawlist')[str_subset(ls()  , 'rawlist') %>>% str_count() > 20]
ana.list = list()
for( i in seq_along(rawlist.name)) {
  date.char = rawlist.name %>>% str_sub(8) %>>% {.[i]}
  remove.event = get(paste0('allver' , date.char ,'.dt'))[traditional==0]$ver.list
  in.edge.list = get(rawlist.name[i])[[5]][!event %in% remove.event]

  ana.list[[i]] =
    analysis1(in.edge.list ,  period = str_sub(rawlist.name[i],8,str_count(rawlist.name[i])) ,
              fixed = 'fixed.one.edge.group')
  cat('data ', i , '\n')
}
ana.byyear.w = do.call(rbind ,ana.list)
ana.byyear.w = as.data.table(ana.byyear.w)
ana.byyear.w = ana.byyear.w[,names(ana.byyear.w)[-1] := lapply(.SD ,function(x) round(as.numeric(x),5)  ) , .SDcols = names(ana.byyear.w)[-1]]
ana.byyear.w[,scale.modularity:=scale.modularity*100]
ana.byyear.w[,time:=1:.N]
ana.byyear.w[,period:= str_replace(period , '[.]' , '-')]

# top 30% events & regular events
top.30.event = allver2012.2018.dt[degree.v>19&type.v=='event']$ver.list
regular.event = allver2012.2018.dt[ regular == 1 & type.v=='event']$ver.list
china.event = allver2012.2018.dt[ china == 1 & type.v=='event']$ver.list
labor.event = allver2012.2018.dt[ labor == 1 & type.v=='event']$ver.list
gender.event = allver2012.2018.dt[ gender == 1 & type.v=='event']$ver.list
environment.event = allver2012.2018.dt[ environment == 1 & type.v=='event']$ver.list
land.event = allver2012.2018.dt[ land == 1 & type.v=='event']$ver.list


date.list = str_subset(ls() , 'rawlist')[(str_subset(ls() , 'rawlist') %>>% str_count() >20)] %>>%
  str_sub(8,24)

# add event label
top.30.event.list = vector()
regular.event.list = vector()
china.event.list = vector()
labor.event.list = vector()
gender.event.list = vector()
environment.event.list = vector()
land.event.list = vector()

for(i in seq_along(date.list)) {
  top.30.event.list[i] = table(get(paste0('rawlist' ,date.list[i] ))$event.list %in% top.30.event)[2]
  regular.event.list[i] = table(get(paste0('rawlist' ,date.list[i] ))$event.list %in% regular.event)[2]
  china.event.list[i] = prop.table(table(get(paste0('rawlist' ,date.list[i] ))$event.list %in% china.event))[2]
  labor.event.list[i] = prop.table(table(get(paste0('rawlist' ,date.list[i] ))$event.list %in% labor.event))[2]
  gender.event.list[i] = prop.table(table(get(paste0('rawlist' ,date.list[i] ))$event.list %in% gender.event))[2]
  environment.event.list[i] = prop.table(table(get(paste0('rawlist' ,date.list[i] ))$event.list %in% environment.event))[2]
  land.event.list[i] = prop.table(table(get(paste0('rawlist' ,date.list[i] ))$event.list %in% land.event))[2]
}


ana.byyear.w[, c('top.30' ,'regular' ,'china' ,'labor','gender','environment','land'):=list(top.30.event.list , regular.event.list , china.event.list ,labor.event.list,gender.event.list,
                                                                                            environment.event.list,land.event.list)]

# add degree by event label
china.deg.list = vector()
labor.deg.list = vector()
gender.deg.list = vector()
environment.deg.list = vector()
land.deg.list = vector()
deg.list = vector()

for(i in seq_along(date.list)){
  deg.list[i] = sum(get(paste0( 'allver', date.list[i] , '.dt'))[type.v=='event' , degree.v ])
  china.deg.list[i] = sum(get(paste0( 'allver', date.list[i] , '.dt'))[type.v=='event' & china ==1, degree.v ])/deg.list[i]
  labor.deg.list[i] = sum(get(paste0( 'allver', date.list[i] , '.dt'))[type.v=='event' & labor ==1, degree.v ])/deg.list[i]
  gender.deg.list[i] = sum(get(paste0( 'allver', date.list[i] , '.dt'))[type.v=='event' & gender ==1, degree.v ])/deg.list[i]
  environment.deg.list[i] = sum(get(paste0( 'allver', date.list[i] , '.dt'))[type.v=='event' & environment ==1, degree.v ])/deg.list[i]
  land.deg.list[i] = sum(get(paste0( 'allver', date.list[i] , '.dt'))[type.v=='event' & land ==1, degree.v ])/deg.list[i]
}

ana.byyear.w[, c('china.d' ,'labor.d','gender.d','environment.d','land.d' ,'deg'):=list(china.deg.list , labor.deg.list , gender.deg.list ,environment.deg.list,land.deg.list,
                                                                                        deg.list)]


# saveRDS(ana.byyear.w , 'data/clean.dt/ana.byyear.w.rds')
# ana.byyear.w[, (n.group + n.event)/(mean(n.group)+mean(n.event))]
ana.byyear.w[,.(period, n.bi.O/n.bi.R*((n.group + n.event)/(mean(n.group)+mean(n.event))) )]
ana.byyear.w[,.(CC.O,CC.R,(CC.O)/(CC.R))]

#------#
# legend setting
set.legend <- function(x =.5, y = .5){
  list(
    font = list(
      family = "Garamond",
      size = 16,
      color = "#000"),
    bgcolor = "#E2E2E2",
    bordercolor = "#FFFFFF",
    borderwidth = 4 , x =x , y = y)
}

# v line
vline <- function(x = 0, color = "red" , dash = 'line') {
  list(
    type = 'line',
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color , dash =dash)
  )
}
# add text
addtext = function(x = 0 , y = 0 , text = '', ax = 95 , ay = -60) {
  list(
    x = x,
    y = y,
    text = paste0('<b>' , text , '<b>'),
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 1,
    ax = ax,
    ay = ay,
    font = list(color = '#000000', size = 28)
  )
}

# font
f <- list(
  size = 25
)

#  y axis 2
ay = function(title = "" , range = NA) {
  list(
    tickfont = list(color = "red" ,size= 50),
    overlaying = "y",
    side = "right",
    title = title,
    titlefont  = f ,
    tickwidth  =4 ,
    showline = TRUE,
    zeroline = TRUE ,
    range = range
  )
}

ayo = function(title = '', range = NA){
  list(
    title  =title ,
    showline = TRUE,
    zeroline = TRUE,
    titlefont = f ,
    tickwidth  =4,
    tickfont  = list(size= 50),
    range = range
    
  )
} 

ax <- list(
  title = "Timeline",
  showline = F,
  zeroline = TRUE,
  showticklabels = FALSE,
  showgrid = FALSE,
  titlefont  = f
)


# descriptive
# p1 <- plot_ly(ana.byyear, x = ~period, y = ~n.event , name = 'Number of Events', type = 'scatter', mode = 'lines+markers') %>>%
#   add_trace(y = ~n.bi.O, name = 'Bicomponent size', mode = 'lines+markers' ) %>>%
# layout(title = "Descriptive",
#        xaxis = list(title = ""),
#        yaxis = list (title = "Count"))
# 
# p2 <- plot_ly(ana.byyear, x = ~period, y = ~n.group , name = 'Number of Groups', type = 'scatter', mode = 'lines+markers') %>>%
#   add_trace(y = ~n.edge.O, name = 'Number of Arcs', mode = 'lines+markers' ) %>>%
#   layout(title = "Descriptive",
#          xaxis = list(title = ""),
#          yaxis = list (title = "Count"))


# htmlwidgets::saveWidget(as.widget(subplot(p1,p2 ,nrows = 2)), "descriptive.html")

#window 
p2 <- plot_ly(ana.byyear.w, x = ~period, y = ~n.event , name = 'Number of Events', type = 'scatter', mode = 'lines+markers',line= list(width =9 ,dash = 'dot')) %>>%
  add_trace(y = ~n.group, name = 'Number of Groups', mode = 'lines+markers',line= list(width = 3.5)) %>>%
  add_trace(y = ~top.30, name = 'Top 30%', mode = 'lines+markers',line= list(width = 3.5)) %>>%
  add_trace(y = ~regular, name = 'Regular', mode = 'lines+markers' ,line= list(width = 3.5)) %>>%
  add_trace(y = ~china, name = 'China', mode = 'lines' ,yaxis = "y2",line= list(width =9 , dash ='solid')) %>>%
  add_trace(y = ~labor, name = 'Labor', mode = 'lines' ,yaxis = "y2",line= list(width = 9 , dash ='solid')) %>>%
  add_trace(y = ~gender, name = 'Gender', mode = 'lines' ,yaxis = "y2",line= list(width = 9 , dash ='solid')) %>>%
  add_trace(y = ~environment, name = 'Environment', mode = 'lines',yaxis = "y2" ,line= list(width = 9, dash ='solid')) %>>%
  add_trace(y = ~land, name = 'Land', mode = 'lines',yaxis = "y2" ,line= list(width = 9 , dash ='solid')) %>>%
  layout(title = "", yaxis2 = ay("" , range = c(0,.4)) ,
         xaxis = ax,
         yaxis = ayo("" , range = c(16,32)),
         margin = list(l = 100, r = 40, b = 30, t = 20),
         legend = set.legend(1.2,1),
         shapes = list(vline(15 ,dash = 'dash') , vline(41 , dash='dot') ) ,
         annotations = list(addtext(x = 15 , y =30 , paste0('2014-03-18' , '<br>' , '太陽花學運爆發')) ,
                            addtext(x = 41 , y =30 , paste0('2016-05-20' , '<br>' , '民進黨政府上台'))))
htmlwidgets::saveWidget(as_widget(p2), "event.descrip.w.html")

p1 <- plot_ly(ana.byyear.w, x = ~period, y = ~ deg , name = 'Median Degree of Events', type = 'scatter', mode = 'lines+markers',line= list(width =9 ,dash = 'dot') ) %>>%
  add_trace(y = ~china.d, name = 'China', mode = 'lines' ,line= list(width = 9 , dash ='solid')  ,yaxis = "y2") %>>%
  add_trace(y = ~labor.d, name = 'Labor', mode = 'lines' ,line= list(width = 9 , dash ='solid'),yaxis = "y2") %>>%
  add_trace(y = ~gender.d, name = 'Gender', mode = 'lines',line= list(width = 9 , dash ='solid'),yaxis = "y2") %>>%
  add_trace(y = ~environment.d, name = 'Environment', mode = 'lines' ,line= list(width = 9 , dash ='solid'),yaxis = "y2") %>>%
  add_trace(y = ~land.d, name = 'Land', mode = 'lines' ,line= list(width = 9 , dash ='solid'),yaxis = "y2") %>>%
  layout(title = "",yaxis2 = ay("" ,range = c(0,.45) ) ,
         xaxis = ax,
         yaxis = ayo("" ,range = c(200,550) ),
         margin = list(l =110, r = 40, b = 30, t = 20),
         legend = set.legend(1.2,1),
         shapes = list(vline(15 ,dash = 'dash') , vline(41 , dash='dot') ) ,
         annotations = list(addtext(x = 15 , y =500 , paste0('2014-03-18' , '<br>' , '太陽花學運爆發')) ,
                            addtext(x = 41 , y =500 , paste0('2016-05-20' , '<br>' , '民進黨政府上台'))))



###

#The default order will be alphabetized unless specified as below:


# p3 <- plot_ly(ana.byyear, x = ~period, y = ~n.bi.O, name = 'Size of Bicomponent', type = 'scatter', mode = 'lines+markers',
#              line = list(color = 'rgb(205, 12, 24)', width = 4)) %>%
#   # add_trace(y = ~tri.ration.F, name = 'Tricomponent ratio ', line = list(color = 'rgb(22, 96, 167)', width = 4)) %>%
#   add_trace(y = ~bi.R.u, name = 'Random Graph 95% CI UPPER LIMIT', line = list(color = 'rgb(205, 12, 24)', width = 4, dash = 'dash')) %>%
#   add_trace(y = ~bi.R.d, name = 'Random Graph 95% CI LOWER LIMIT', line = list(color = 'rgb(205, 12, 24)', width = 4, dash = 'dash')) %>%
#   layout(title = "Size of Bicomponent",
#          xaxis = list(title = "Years"),
#          yaxis = list (title = "Ratio"),legend = l)


# window
p3 <- plot_ly(ana.byyear.w, x = ~period, y = ~((n.bi.O)/(bi.R)), name = 'Bicomponent% in Obs. / Bicomponent% in Random', type = 'scatter', mode = 'lines+markers',
              line = list(color = 'rgb(205, 12, 24)', width = 4)) %>>% 
  add_trace(y = ~scale.modularity, name = 'Modularity / ln(Network Size)', line = list(color = 'rgb(205, 12, 24)', width = 1, dash = 'dash' ), mode = 'lines' , yaxis = "y2") %>>%
  layout(title = "Size of Bicomponent",
         xaxis = list(title = ""),
         yaxis = list (title = "Odds Ratio"),margin = list(l = 35, r = 30, b = 130, t = 25),
         legend = l ,  yaxis2 = ay() , 
         shapes = list(vline(15 ,dash = 'dash') , vline(41 , dash='dot') ) ,
         annotations = list(addtext(x = 15 , y =1 , paste0('2014-03-18' , '<br>' , '太陽花學運爆發')) ,
                            addtext(x = 41 , y =1 , paste0('2016-05-20' , '<br>' , '民進黨政府上台'))))
htmlwidgets::saveWidget(as_widget(p3), "bicomponent.modularity.w.html")

p3.1 <- plot_ly(ana.byyear.w, x = ~period, y = ~n.bi.O, name = 'Bicomponent Size', type = 'scatter', mode = 'lines',
              line = list(color = 'rgb(205, 12, 24)', width = 4)) %>>% 
  # add_trace(y = ~scale.modularity, name = 'Modularity / ln(Network Size)', line = list(color = 'rgb(205, 12, 24)', width = 4, dash = 'dash' ), mode = 'lines' , yaxis = "y2") %>%
  add_trace(y = ~(CC.O/CC.R) , name = 'CC Observation/CC Random' ,  line = list(color = 'rgb(205, 12, 24)', width = 4, dash = 'dot' , connectgaps = TRUE) , mode = 'lines' , yaxis = "y2") %>>%
  layout(title = "",
         xaxis = ax,
         yaxis = ayo(''),margin = list(l = 35, r = 30, b = 30, t = 25), 
         legend = set.legend(1.2,1) , yaxis2 = ay() ,
         shapes = list(vline(15 ,dash = 'dash') , vline(41 , dash='dot') ) ,
         annotations = list(addtext(x = 15 , y =120 , paste0('2014-03-18' , '<br>' , '太陽花學運爆發')) ,
         addtext(x = 41 , y =120 , paste0('2016-05-20' , '<br>' , '民進黨政府上台'))
         ))
htmlwidgets::saveWidget(as_widget(p3.1), "bicomponent.modularity.w.2.html")

p4 = plot_ly(ana.byyear.w, x = ~period, y = ~CC.O/CC.R, name = 'Clustering Coefficient', type = 'scatter', mode = 'lines+markers',
                line = list(color = 'rgb(205, 12, 24)', width = 4)) %>>% 
  add_trace(y = ~scale.modularity, name = 'Modularity / ln(Network Size)', line = list(color = 'rgb(205, 12, 24)', width = 1, dash = 'dash' ), mode = 'lines' , yaxis = "y2") %>%
  layout(title = "",
         xaxis = list(title = ""),
         yaxis = list (title = "CC"),margin = list(l = 35, r = 30, b = 130, t = 25),
         legend = set.legend(1.05,1) , yaxis2 = ay() ,
         shapes = list(vline(15 ,dash = 'dash') , vline(41 , dash='dot') ) ,
         annotations = list(addtext(x = 15 , y =.6 , paste0('2014-03-18' , '<br>' , '太陽花學運爆發')) ,
                            addtext(x = 41 , y =.6 , paste0('2016-05-20' , '<br>' , '民進黨政府上台'))
         ))
htmlwidgets::saveWidget(as_widget(p4), 'CC.modularity.html')

  
p4.1 = plot_ly(ana.byyear.w, x = ~period, y = ~(1-CC.O)/(1-CC.R), name = 'Clustering Coefficient', type = 'scatter', mode = 'lines+markers',
             line = list(color = 'rgb(205, 12, 24)', width = 4)) %>>% 
  add_trace(y = ~avg.path.O, name = 'Averge Path Length', line = list(color = 'rgb(205, 12, 24)', width = 1, dash = 'dash' ), mode = 'lines' , yaxis = "y2") %>%
  layout(title = "",
         xaxis = list(title = ""),
         yaxis = list (title = "CC"),margin = list(l = 35, r = 30, b = 130, t = 25), legend = l , yaxis2 = ay() ,
         shapes = list(vline(15 ,dash = 'dash') , vline(41 , dash='dot') ) ,
         annotations = list(addtext(x = 15 , y =.6 , paste0('2014-03-18' , '<br>' , '太陽花學運爆發')) ,
                            addtext(x = 41 , y =.6 , paste0('2016-05-20' , '<br>' , '民進黨政府上台'))
         ))
htmlwidgets::saveWidget(as_widget(p4), 'CC.modularity.html')

 # p2.2 <- plot_ly(ana.byyear, x = ~period, y = ~bi.ration.F, name = 'Bicomponent ratio ', type = 'scatter', mode = 'lines+markers',
#               line = list(color = 'rgb(205, 12, 24)', width = 4)) %>%
#   add_trace(y = ~tri.ration.F, name = 'Tricomponent ratio ', line = list(color = 'rgb(22, 96, 167)', width = 4)) %>%
#   layout(title = "Ratio of Bicomponent and Tricomponent",
#          xaxis = list(title = "Years"),
#          yaxis = list (title = "Ratio"))



# table 3 period
ana.t1v1 = analysis1(rawlistt1v1[[5]] ,  period = 't1v1' ,  fixed = 'fixed.one.edge.group')
ana.t2v1= analysis1(rawlistt2v1[[5]] , period='t2v1' ,  fixed = 'fixed.one.edge.group')
ana.t3v1= analysis1(rawlistt3v1[[5]] , period='t3v1' ,   fixed = 'fixed.one.edge.group')
ana.byperiod = rbind(ana.t1v1,ana.t2v1,ana.t3v1)
ana.byperiod=as.data.table(ana.byperiod)
ana.byperiod[,names(ana.byperiod)[-1] := lapply(.SD ,function(x) round(as.numeric(x),3)) , .SDcols = names(ana.byperiod)[-1]]
ana.byperiod[,CC.ratio:=round(CC.O/CC.R,3)]
ana.byperiod[,avg.path.ratio:=round( avg.path.O/avg.path.R,3) ]

column1 = c('' , 'Events' , 'Groups' , 'Small-World Parameters' ,' Clustering Coefficient' , '  (Random expected)' , ' Ratio of Observed to Random' ,
            ' Average Path Length', '  (Random expected)' ,' Ratio of Observed to Random' , 'Size of Largest Bicomponent' , ' Observed' , 
            '  Ratio of Bicomponet' , ' Random Group Assignment', '  Ratio of Bicomponet' ,' Ratio of Observed to Random')

column2 = c('2012-20140318' , ana.t1v1[2] , ana.t1v1[3] , '' )
column2_4 = matrix(nrow = 16 , ncol = 3)
column2_4[1,] = c('2012-20140318' , '20140319-20160501' , '20160502-20180631')
column2_4[2,] = ana.byperiod$n.event
column2_4[3,] = ana.byperiod$n.group
column2_4[4,] = ""
column2_4[5,] = ana.byperiod$CC.O
column2_4[6,] =  paste0(ana.byperiod$CC.R , '(' , ana.byperiod$CC.R.sd,')')
column2_4[7,] = ana.byperiod$CC.ratio
column2_4[8,] = ana.byperiod$avg.path.O
column2_4[9,] = paste0(ana.byperiod$avg.path.R , '(' , ana.byperiod$avg.path.R.sd,')')
column2_4[10,] = ana.byperiod$avg.path.ratio
column2_4[11,] = ""
column2_4[12,] = ana.byperiod$n.bi.O
column2_4[13,] = round(ana.byperiod$n.bi.O/ (ana.byperiod$n.event + ana.byperiod$n.group),3)
column2_4[14,] = paste0(ana.byperiod$bi.R , '(' ,ana.byperiod$bi.R.sd  , ")")
column2_4[15,] = round(ana.byperiod$bi.R/ (ana.byperiod$n.event + ana.byperiod$n.group),3)

column2_4[16,] = round(ana.byperiod$n.bi.O/ana.byperiod$bi.R,3)

out.table1 = cbind(column1, column2_4)

write.csv(out.table1 , 'table1.csv' ,row.names = F )


# des table2
allvert1v1.dt[,degree.v , by = 'type.v' ]

###--------------------- degree distribution in three period #------------------####
deg.t1v1 = degree(net.t1v1)
deg.t2v1 = degree(net.t2v1)
deg.t3v1 = degree(net.t3v1)
deg.2012.2018 = degree(net.2012.2018)

deg.t1v1.g = deg.t1v1[names(deg.t1v1) %in% rawlistt1v1$group.list]
deg.t2v1.g = deg.t2v1[names(deg.t2v1) %in% rawlistt2v1$group.list]
deg.t3v1.g = deg.t3v1[names(deg.t3v1) %in% rawlistt3v1$group.list]
deg.2012.2018.g = deg.2012.2018[names(deg.2012.2018) %in% rawlist2012.2018$group.list]

deg.t1v1.e = deg.t1v1[names(deg.t1v1) %in% rawlistt1v1$event.list]
deg.t2v1.e = deg.t2v1[names(deg.t2v1) %in% rawlistt2v1$event.list]
deg.t3v1.e = deg.t3v1[names(deg.t3v1) %in% rawlistt3v1$event.list]
deg.2012.2018.e = deg.2012.2018[names(deg.2012.2018) %in% rawlist2012.2018$event.list]

bicom.size = function(net) {
  bicomponent_list <- biconnected_components(net)
  n.bi.o = sapply(bicomponent_list$components, length) %>>% max
  return(n.bi.o)
}

# bicomponent sensitivity analysis--------------------------------------------------------------#
#find largest bicomponent ----------------------------------#
detect.largest.vertex = function(net = net.2012 ,rawlist = rawlist2012 , intype = 'event'){
  
  deg.net= degree(net)
  if(intype == 'event') {
    deg.net = deg.net[names(deg.net) %in% rawlist$event.list]  
  }else if(intype == 'group') {
    deg.net = deg.net[names(deg.net) %in% rawlist$group.list]  
  }
  max.deg = max(deg.net)
  larg.ver = names(which.max(deg.net))
  
  cat(larg.ver , ' degree = ' , max.deg , '\n')
  
  net.tmp = delete_vertices(net, which(V(net)$name %in% larg.ver))
  return( list('g' = net.tmp , 'degree'= max.deg , 'vertex' = larg.ver))
}

bicomponent.sensitivity = function(net , rawlist , intype){
  sizeofbicom = bicom.size(net)
  netlist.tmp = detect.largest.vertex(net ,rawlist = rawlist , intype= intype)
  
  tmp.list = data.frame(matrix(ncol = 5))
  dtname = c('Iteration' ,'Vertext' ,'Vertex.degree' , 'Netsize' , 'Bicom.size')
  names(tmp.list) = dtname
  
  tmp.list[1,] = c(1 , netlist.tmp$vertex , netlist.tmp$degree , length(V(netlist.tmp$g)) , sizeofbicom)
 
  it = 2
  if(intype == 'event'){
      i = length(V(net)[!V(net)$type])-1
  }else if(intype == 'group'){
      i = length(V(net)[V(net)$type])-1
  }
  while(i > 1 ){
    netlist.tmp = detect.largest.vertex(netlist.tmp$g ,rawlist = rawlist , intype)
    sizeofbicom = bicom.size(netlist.tmp$g)
    
      tmp.list[it,] = c(it , netlist.tmp$vertex , netlist.tmp$degree , length(V(netlist.tmp$g)) , sizeofbicom)
    it = it +1
    i = i-1
  }
  return(tmp.list)
}

tmpdt.list = list()
netname = list('t1v1' , "t2v1" , 't3v1' ,'2012.2018')
for (i in netname) {
  cat('\n' , i , '\n\n')
  tmpdt.e = bicomponent.sensitivity(get(paste0('net.',i)), get(paste0('rawlist' , i)) , 'event')
  tmpdt.g = bicomponent.sensitivity(get(paste0('net.',i)), get(paste0('rawlist' , i)) , 'group')
  
  setDT(tmpdt.e)[,type := 'Event'] ; setDT(tmpdt.g)[,type := 'Group'] 
  
  tmpdt.list[[i]] = rbind(tmpdt.e ,tmpdt.g)
  
}
  # legend setting

  l.set = function(x = .2,y = .95 , showlegend = FALSE){
    l <- list(
      font = list(
        family = "Garamond",
        size = 30,
        color = "#000"),
      bgcolor = "#E2E2E2",
      bordercolor = "#FFFFFF",
      borderwidth = 4 , x =x , y = y,
      showlegend = showlegend)
  }
plot.list = list()
  for (i in netname) {
    if(i == '2012.2018'){
     l  = l.set(showlegend = T)
    }else{
     l  = l.set()
    }
    
  plot.list[[i]]  = plot_ly(data = tmpdt.list[[i]], x = ~(V(get(paste0('net.' ,i)) ) %>>% length-as.numeric(Iteration)), y = ~as.numeric(Bicom.size),
            mode = 'lines+markers', symbol = ~type , type = 'scatter' ,
            text= ~paste('Vertex: ', Vertext , '<br>','Vertex.degree: '  ,Vertex.degree)) %>>%
      layout( xaxis = list(title = 'Vertices remain' ,range = list(650,1050) ),
              yaxis = list(title = 'Size of Bicomponent'  ,range = list(0,250)) , 
              legend = l)
  }


  p1 =  plot_ly(data = tmpdt.list[[4]], x = ~(V(net.2012.2018) %>>% length-as.numeric(Iteration)), y = ~as.numeric(Bicom.size),
          mode = 'lines+markers', symbol = ~type , type = 'scatter' , marker = list(size = 14) ,
         text= ~paste('Vertex: ', Vertext , '<br>','Vertex.degree: '  ,Vertex.degree)) %>>%
  layout( 
    xaxis = list(
    title = 'Vertices remain',
    showline = T,
    zeroline = TRUE,
    titlefont  = list(size=30),
    titlefont  = list(size=30),
    tickfont = list(size=30),range = list(650,1050) 
  ) ,
    yaxis =  list(
            title = 'Size of Bicomponent',
            tickwidth  =4 ,
            showline = TRUE,
            zeroline = TRUE,
            titlefont  = list(size=25),
            tickfont = list(size=25) , type ='log'
          ) , 
          legend = l.set(showlegend = T) )
htmlwidgets::saveWidget(as.widget(p1), "Sensitivity.bicomp.html")



###----#


p2 <- plot_ly(alpha = 0.6) %>%
  add_histogram(x = ~deg.2012.2018.g , name = 'Groups' ) %>%
  add_histogram(x = ~deg.2012.2018.e , name = 'Events') %>%
  layout(barmode = "overlay" ,
          xaxis = list(domain = c(0.6, 0.95)),
         yaxis = list(domain = c(0.6, 0.95)))

ay <- list(
  tickfont = list(color = "red"),
  overlaying = "y",
  side = "right",
  title = "second y axis"
)
p <-  plot_ly(data = tmp.dt[Iteration<300,], x = ~300-as.numeric(Iteration), y = ~as.numeric(Bicom.size),
          mode = 'lines+markers', symbol = ~type ,
         text= ~paste('Vertex: ', Vertext , '<br>','Vertex.degree: '  ,Vertex.degree) , yaxis = "y2") %>%
  add_histogram(x = ~deg.2012.2018.g , name = 'Groups' ) %>%
  add_histogram(x = ~deg.2012.2018.e , name = 'Events') %>%
  layout(
    title = "Double Y Axis", yaxis2 = ay,
    xaxis = list(title="x")
  )

# plot_ly(x = 0:120, y = deg.2012.2018.e.den, mode = "lines", fill = "tozeroy", yaxis = "y2", name = "Density") %>% 
#   add_trace(x = deg.2012.2018.e, y = deg.2012.2018.e.den, mode = "lines", fill = "tozeroy", yaxis = "y2", name = "Density") %>% 
#   layout(yaxis2 = list(overlaying = "y", side = "right"))

## -- degree plot by catgorical # group only ----------------------------------####



degobj = str_subset(ls() , 'deg.t....g|deg.2012.2018.g')
for(i in 1:4){
  assign(degobj[i] ,ifelse(get(degobj[i]) >3 , '>=4' , get(degobj[i])) )
}


f <- list(
  family = "Courier New, monospace",
  size = 18,
  color = "#7f7f7f"
)
x <- list(
  title = "Degree",
  titlefont = f,
  categoryorder = 'array', categoryarray = c('1', '2', '3','>=4')
)
y <- list(
  title = "Count",
  titlefont = f
)
p <- plot_ly(alpha = 0.5) %>>%
  add_histogram(x = ~deg.t1v1.g , name = "T1") %>>%
  add_histogram(x = ~deg.t2v1.g , name = "T2") %>>%
  add_histogram(x = ~deg.t3v1.g, name = "T3") %>>%
  layout(xaxis = x, yaxis = y, barmode = 'group')



#---------------------------------------------------------------------------------####



### 事件分類附表
event.des = allver2012.2018.dt[type.v=='event', .(ver.list , peckt ,china , labor, gender , environment , land ,regular) ]
event.des = cbind(event.des[,1:2],
      event.des[,lapply(.SD , function(x) ifelse(x==1,'Y' ,'') ) , .SDcols = c( 'regular','china','labor','gender','environment','land')]
)
event.des = event.des[order(peckt %>>% ymd())]
setnames(event.des ,c('事件名稱' ,'高峰時間' ,'例行性事件' , '中國' , '勞工' ,'性別' ,'環境' ,'土地'))
write.csv(event.des , 'plot/暫定採用圖表/事件描述附表.csv' , row.names = FALSE  )
###









####

plot_ly(ana.byyear, x = ~period, y = ~CC.O , name = 'CC.O ', type = 'scatter', mode = 'lines+markers' , yaxis = "y1") %>>%
  add_trace( y = ~avg.path.O, name = 'avg.path.O', mode = 'lines+markers' ,yaxis = "y2") %>>%



p <- plot_ly(ana.byyear, x = ~period, y = ~den.O, name = 'den.O', type = 'scatter', mode = 'lines') %>%
  add_trace(y = ~CC.O, name = 'CC.O', mode = 'lines+markers') %>%
  add_trace(y = ~bi.ration.F, name = 'bi.ration.F', mode = 'lines+markers')

#--------------------------- k component--------------------------------#
clustering_tm(edgelist)
[1] 0.6558297
> average.path.length(rawlist2015[[1]])
[1] 3.959552
> membership.k.component(k.connected(rawlist2015[[1]]) , 2)

library(ggplot2)
library(plotly)

table(allver2012.dt$largest.k.comp , allver2012.dt$type.v) 
table(allver2013.dt$largest.k.comp , allver2013.dt$type.v) 
table(allver2014.dt$largest.k.comp , allver2014.dt$type.v) 
table(allver2015.dt$largest.k.comp , allver2015.dt$type.v) 
table(allver2016.dt$largest.k.comp , allver2016.dt$type.v) 
table(allver2017.dt$largest.k.comp , allver2017.dt$type.v) 
table(allver2018.dt$largest.k.comp , allver2018.dt$type.v) 

mean(allver2012.dt$largest.k.comp)/12
mean(allver2014.dt$largest.k.comp)/12
mean(allver2017.dt$largest.k.comp)/12
mean(allver2018.dt$largest.k.comp)/6

allver2012.dt[,year := 2012]
allver2013.dt[,year := 2013]
allver2014.dt[,year := 2014]
allver2015.dt[,year := 2015]
allver2016.dt[,year := 2016]
allver2017.dt[,year := 2017]
allver2018.dt[,year := 2018]


allvar = rbindlist(list(allver2012.dt , allver2013.dt , allver2014.dt , 
               allver2015.dt , allver2016.dt , allver2017.dt ,allver2018.dt))
# allvar[,year:=as.character(year)]

dcast.data.table(allvar , ver.list ~ largest.k.comp  , value.var = 'year' )

with(allvar , 
     table(year , largest.k.comp))
test1 = allvar[type.v == 'event',mean(largest.k.comp) , by = 'year']
test2 = allvar[type.v == 'group',mean(largest.k.comp) , by = 'year']
test3 = allvar[ ,mean(largest.k.comp) , by = c('year')]


p <- plot_ly(test, x = ~year, y = ~V1 , mode = 'line')

p <- plot_ly(data, x = ~x, y = ~trace_0, name = 'trace 0', type = 'scatter', mode = 'lines') %>%
  add_trace(y = ~trace_1, name = 'trace 1', mode = 'lines+markers') %>%
  add_trace(y = ~trace_2, name = 'trace 2', mode = 'markers')



