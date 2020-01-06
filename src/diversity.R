


# modeling diversity of each edge and event diversity


# analysis level  # group or arc?
library(data.table)
library(stringr)
library(igraph)
library(pipeR)
library(lubridate)
library(RColorBrewer)
library(readxl)

# type of edge:
# 1. consolidation  # same communication
#    1. new consolidation
#    2. repeat consolidation
# 2. bridge         # cross comunication 
#    1. new bridge
#    2. repeat bridge
# 3. extention    
# 4. encounter  # two vertex are new


#-------------------fun----------------------------------#
# Determine if range of vector is FP 0.
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}
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
# -------------------------------------------------#


# new event
allver2012.2018.dt[type.v=='event' , act.time.date:=as.Date(act.time) ]
begin.event.list  = allver20120101.20121231.dt[type.v=='event',ver.list]
other.event.list = rawlist2012.2018$event.list[!rawlist2012.2018$event.list %in% begin.event.list ]

rawlist.name = str_subset(ls()  , 'rawlist')[str_subset(ls()  , 'rawlist') %>>% str_count() > 20]
allver.name = str_subset(ls()  , 'allver')[str_subset(ls()  , 'allver') %>>% str_count() > 20]

# data wrangling per "other.event.list"
edge.list = list()
event.table = list()
for(i in seq_along(other.event.list)){
  input.event.name = other.event.list[i]
  evend.infor = allver2012.2018.dt[ver.list==input.event.name , ] 
  
  evend.infor$act.time %>>% str_remove_all('-')
  
  # find time position
  num.list = rawlist.name %>>% str_sub(-8,-1) %>>% as.numeric()
  num.list = num.list - as.numeric(evend.infor$act.time %>>% str_remove_all('-') ) 
  posi = which(abs(0-num.list)==min(abs(0-num.list)))
  
  g =get(rawlist.name[posi])$igraph
  
  # neighbor of new event
  ego.g = make_ego_graph(g, nodes = input.event.name)
  ego.edge.list = data.table(bipartite.projection(ego.g[[1]])$proj2 %>>% as_edgelist() )
  ego.edge.list[,time:= posi]
  standard1 = apply(ego.edge.list , 1 , 
              function(x) {
                # dose vertex exist before?
            match.posi  =  which(as.character(x) %in% V(get(rawlist.name[posi-1])$igraph)$name )
            out = length(match.posi)
            if(length(match.posi)==0){
              out= 0
            }
            return(out)
              } 
  )
  
  # 0 = ecounter  無-無
  # 1 = extension 有-無
  # 2             有-有
  

  ### 合併團體過去參加事件attribute
  # 組織合作的community(one mode) & 過去參加事件類型
  tmp.dt = group.past.dt[time == posi-1]
  
  ego.edge.list[,exist.before:= standard1]
  ego.edge.list[exist.before==0,connected.type:='ecounter'][
    exist.before==1,connected.type:='extension']
  
  ego.edge.list = rbind(ego.edge.list[,c(1,2,3,4)] , ego.edge.list[,.('V1'=V2,'V2' =V1, time ,exist.before )])
  
  ## test---
  testdt1= tmp.dt[match(ego.edge.list[exist.before==2 ,]$V1, tmp.dt$ver.list)]  # ego
  testdt2=tmp.dt[match(ego.edge.list[exist.before==2 ,]$V2, tmp.dt$ver.list)]   # alter
  setnames(testdt1 , paste0('E.',names(testdt1)));setnames(testdt2 , paste0('A.',names(testdt2) ))
  combtestdt = cbind(testdt1  ,testdt2  )
  
  combtestdt[,c('fg.com.cross','ml.com.cross' , 'wt.com.cross' , 'eig.com.cross'):=list(!(E.fg.com == A.fg.com ) ,!(E.ml.com == A.ml.com ) , 
                                                                      !(E.wt.com == A.wt.com ) , !(E.eig.com == A.eig.com ))]
  
  ego.edge.list[exist.before==2 , names(combtestdt):= as.list(combtestdt)]
  
  # if(any(standard1==2)){
  #   standard2 = 
  #     apply(ego.edge.list[exist.before == 2,.(V1,V2)] ,1 ,function(x)  {
  #       tmp.dt[ver.list %in% as.character(x)]$community %>>% 
  #         zero_range()
  #     }  
  #     )  # FALSE = bridge ; TRUE = consolidation
  #   
  #   ego.edge.list[exist.before==2 , cross.com:= !standard2]
  #   ego.edge.list[exist.before==2 , connected.type:=
  #                   ifelse(cross.com== T, 'bridge','consolidation')]
  # }
  if(nrow(ego.edge.list)==0){
    ego.edge.list=  rbind(ego.edge.list , NA , fill =T)
    ego.edge.list[,x:=NULL]
  }
  edge.list[[i]] = ego.edge.list[,names(evend.infor):= as.list(evend.infor)]
  # event.table[[i]] =  ego.edge.list$connected.type %>>% table()
 cat(i,'\n') 
}
 
edge.dt = rbindlist(edge.list , fill = T)

# merge with ana.byyear.w

edge.dt = merge(edge.dt , ana.byyear.w[,.(n.event,n.group, time)] , by = 'time' , all.x = T)

cross.community.type = function(dt = edge.dt  , community.type='fg.com'){
  edge.dt[,cross.com:=NULL]
  edge.dt[,cross.com:=get(paste0(community.type ,'.cross' ))]
  edge.dt[exist.before==2 , connected.type:=
                               ifelse(cross.com== T, 'bridge','consolidation')]
}
cross.community.type(community.type='fg.com')
edge.dt[connected.type == 'consolidation', consolidation.count := .N ,by = 'ver.list'][is.na(consolidation.count),consolidation.count:=0]
edge.dt[connected.type == 'bridge', bridge.count:=.N ,by = 'ver.list'][is.na(bridge.count),bridge.count:=0]
edge.dt[connected.type == 'ecounter', ecounter.count:=.N ,by = 'ver.list'][is.na(ecounter.count),ecounter.count:=0]
edge.dt[connected.type == 'extension', extension.count:=.N ,by = 'ver.list'][is.na(extension.count),extension.count:=0]
# 
# edge.dt[ ,connected.type2:=connected.type]
# edge.dt[connected.type %in% c('extension','ecounter') ,connected.type2:='other']
# edge.dt[ ,connected.type3:=ifelse(connected.type %in% c('bridge') , 1 , 0)]

#dv ref group
edge.dt[,connected.type:= relevel(factor(connected.type), ref = "ecounter")]
#main issue ==1 > other
# edge.dt[main.issue %in% names(edge.dt$main.issue %>>% table)[(edge.dt$main.issue %>>% table) < 250],
#         main.issue:= '其它']
#issue ref group = gender
# edge.dt[,main.issue:= relevel(factor(main.issue), ref = "性別")]

# time level <2014-03-18 | < 2016-05-20 | >= 2016-05-20
edge.dt[,timelevel:= ifelse(act.time.date < '2014-03-18' , 
                            'Time 1' , ifelse( act.time.date < '2016-05-20',
                                               'Time 2' , 'Time 3') 
                            )  ]



fit = list()
for(i in 1:3){
  intime = paste0('Time ' , i)
  fit[[i]] = glm(factor(cross.com) ~ 
              factor(china) + factor(labor) + factor(gender) + factor(environment) +
              factor(land)  + factor(regular) + degree.v + E.fg.com.N +  E.pati.n + n.event + n.group , data = edge.dt[timelevel == intime] ,  family = "binomial" )
  
}
fit[[4]] =  glm(factor(cross.com) ~ 
                  factor(china) + factor(labor) + factor(gender) + factor(environment) +
                  factor(land)  + factor(regular) + degree.v + E.fg.com.N +  E.pati.n + n.event + n.group +  relevel(factor(timelevel), ref = "Time 3") ,
                data = edge.dt ,  family = "binomial" )

library(stargazer)
stargazer(fit[[1]] , fit[[2]] ,fit[[3]] ,fit[[4]] , title = 
            "Logistic Regression of Between Community Tie" ,dep.var.labels = 'Between Community Tie' , style = 'asr' ,type="html", out="diversity.model.3.html",out.header = FALSE,
          star.char = c("", "a"), star.cutoffs = c(1, 0.01), notes.append = FALSE, notes = "<sup>a</sup> Value is significant at the p < .01 ",
          column.labels   = c("~2014-03-18 ", "~2016-05-20 " , ' 2016-05-21~' , 'All'),
          covariate.labels = c( 'China' , 'Labor' , 'Gender' , 'Environment' ,'Land' , 'Regular Event' ,'N of SMGs' , 'Size of Community' , 'Number of events involved in the past'  ,
                                'Number of past events' ,'Number of past groups' , 'Before March 18 Movement' , 'Before Third Party Alternation'))

# interaction
#E.pati.china E.pati.labor E.pati.gender E.pati.envi E.pati.land

# all interaction
ego.var= c('E.pati.china', 'E.pati.labor', 'E.pati.gender', 'E.pati.envi', 'E.pati.land')
edge.dt[,(ego.var) := lapply(.SD , function(x)ifelse(is.na(x) , NA , ifelse(x>=1 , 1, 0 ) ) ) , .SDcols = ego.var] # recode to 1;0

event.var = c('factor(china)' , 'factor(labor)' , 'factor(gender)' , 'factor(environment) ', 
    'factor(land)'  )
exp.df = expand.grid(paste0('factor(',ego.var,')') , event.var)
exp.df[,2] = paste0(exp.df[,1] , '*' , exp.df[,2])

#based formula
tmp.f = 'factor(cross.com) ~ 
            factor(china) + factor(labor) + factor(gender) + factor(environment) + 
            factor(land)  + factor(regular) + degree.v + E.fg.com.N +  E.pati.n + n.event + n.group '


get.coef= list()
out.coef = list()
for(i in 1:3){
  intime = paste0('Time ' , i)
  for (x in 1:nrow(exp.df)) {
    tmp.f2 = paste0(tmp.f ,'+', exp.df[x,1],'+',exp.df[x,2])
    fit2 = glm( as.formula(tmp.f2), data = edge.dt[timelevel == intime] ,  family = "binomial" )
    
    get.coef[[x]] = coef(summary(fit2))[nrow(coef(summary(fit2))),c(1,2,4)]
    get.coef[[x]]  = round(get.coef[[x]],3)
  }
  out.coef[[i]] = do.call( rbind, get.coef)

}

get.coef= list()
for (x in 1:nrow(exp.df)) {
  tmp.f2 = paste0(tmp.f ,'+  relevel(factor(timelevel), ref = "Time 3")' ,'+', exp.df[x,1],'+',exp.df[x,2])
  fit2 = glm( as.formula(tmp.f2), data = edge.dt,  family = "binomial" )
  
  # 有多少這樣的組織
  var1 = exp.df[x,1] %>>% str_extract( "\\([^()]+\\)") %>>%  # group
                          { str_sub(.,1+1,str_length(.)-1)}
  var2 = exp.df[x,2] %>>% str_extract_all( "\\([^()]+\\)") %>>%{.[[1]][2]}%>>% # event
                          { str_sub(.,1+1,str_length(.)-1)}
  nub.event = fit2$data[eval(parse(text = var1))==1 & eval(parse(text = var2)) == 1 ,E.ver.list] %>>% 
    length()
  
  
  get.coef[[x]] = c(coef(summary(fit2))[nrow(coef(summary(fit2))),c(1,2,4)] ,'nub.event' = nub.event)
  get.coef[[x]]  = round(get.coef[[x]],3)
}
out.coef[[4]] = do.call( rbind, get.coef)


library(reshape2)
library(ggplot2)

# heatm = list()
# for (i in 1:4) {
i=4
  coe.mat=matrix(out.coef[[i]][,1], nrow = 5)
  p.mat=matrix(out.coef[[i]][,3] , nrow = 5)
  n.mat=matrix(out.coef[[i]][,4] , nrow = 5)
  
  
  rownames(coe.mat) = c('china' , 'labor' , 'gender' ,'environment' ,'land')
  colnames(coe.mat) = c('china' , 'labor' , 'gender' ,'environment' ,'land')
  coe.mat[p.mat>=.01] = NA
  
  melted_coe.mat <- melt(coe.mat, na.rm = FALSE)
  melted_coe.mat$paste.txt = paste0(ifelse(is.na(melted_coe.mat$value) ,'',melted_coe.mat$value),'\n','N = ',melt(n.mat)[,3])
  names(melted_coe.mat)[1:2] = c('組織曾參與之事件屬性' ,'新事件屬性' )
 heatm[[i]]  = ggplot(data = melted_coe.mat, aes(`組織曾參與之事件屬性`, `新事件屬性`, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-3,3), space = "Lab", 
                         name="Coefficients")+theme_minimal()+
    geom_text(aes(`組織曾參與之事件屬性`, `新事件屬性`, label = paste.txt), color = "black", size = 4)
  
# }
library("gridExtra")grid.arrange(p1, p2, nrow = 1)

grid.arrange(heatm[[1]], heatm[[2]], heatm[[3]] , heatm[[4]])
library(gplots)
heatmap.2(as.matrix(melted_coe.mat))
#--------------------------------------------------------------------------------------------------###


edge.dt[time==2,]

#-----------------#


library(nnet)
function(informula = "factor(connected.type)~degree.v + szieof.edge.com  + 
           factor(china) + factor(labor) + factor(gender) + factor(environment) +
           factor(land)  + factor(regular)" , time = c(1,2,3) , family = 'multinomial' ) {
  
  if(all(time %in% c(1,2,3))) {
    informula = paste0(informula , '+ factor(timelevel)')
  }
  
  fit = multinom(as.formula(informula) , data = edge.dt[timelevel %in% paste0('Time ' , time)])
  z <- summary(test)$coefficients/summary(test)$standard.errors
  
  # 2-tailed z test
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  ifelse(p < .001 , 'sig' , "N") 
  
  exp(coef(test))
  
}


test = glm(factor(cross.com) ~ degree.v + szieof.edge.com + 
      factor(china) + factor(labor) + factor(gender) + factor(environment) +
      factor(land)  + factor(regular) +   relevel(factor(timelevel), ref = "Time 3"), data = edge.dt ,  family = "binomial" )



test <- multinom(factor(connected.type2) ~ degree.v + szieof.edge.com+ factor(timelevel) + 
                   factor(china) + factor(labor) + factor(gender) + factor(environment) +
                   factor(land)  + factor(regular), data = edge.dt)
z <- summary(test)$coefficients/summary(test)$standard.errors
z
# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
ifelse(p < .001 , 'sig' , "N") 
exp(coef(test))

test1 <- multinom(factor(connected.type) ~ degree.v + szieof.edge.com + 
                   factor(china) + factor(labor) + factor(gender) + factor(environment) +
                   factor(land)  + factor(regular), data = edge.dt[timelevel == 'Time 3'])
z1 <- summary(test1)$coefficients/summary(test1)$standard.errors
z1
# 2-tailed z test1
p1 <- (1 - pnorm(abs(z1), 0, 1)) * 2
p1
exp(coef(test1))



g= rawlist2012.2018$igraph
g2 = bipartite.projection(g)$proj2

c1 = cluster_fast_greedy(g2)
modularity(c1)
c1$membership %>>% table()
c2 = cluster_leading_eigen(g2)
modularity(c2)
c2$membership %>>% table()
c3 = cluster_edge_betweenness(g2)
modularity(c3)
c3$membership %>>% table()
c4=cluster_walktrap(g2)
modularity(c4)
c4$membership %>>% table()

