library(spatstat)
library(rstan)
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())
source("utils.r")

load("WashingtonDC-agls.rdata")

library(docopt)

'Usage:
model-mismatch.r [--start start] [--end end] [--bwSpace bwSpace] [--bwTime bwTime] [--nonseparable] [--duplicates] [--model model]

Options:
--start start [default: 2010-01-01]
--end end [default: 2010-01-31]
--bwSpace bw-space [default: 1.609]
--nonseparable
--duplicates
--bwTime bw-time [default: 14]
--model model [default: hawkes-model.stan]
]' -> doc

opts <- docopt(doc)

start_date = as.Date(opts$start)
end_date = as.Date(opts$end) 
bw_space = as.numeric(opts$bwSpace)
bw_time = as.numeric(opts$bwTime)
separable = !opts$nonseparable
remove.duplicates = !opts$duplicates

fname = sprintf("data/hawkes-model-%s-%s-%.02f-%.02f-%d-%d.rdata", start_date,end_date,bw_space,bw_time,separable,remove.duplicates)
load(fname)

llmask = colMeans(out$ll) > quantile(colMeans(out$ll),1-mean(out$a))
llshapes = rep(22,length(llmask))
llshapes[llmask == T] = 20
llalphas = rep(.25,length(llmask))
llalphas[llmask == T] = .6

library(ggplot2)
df = data.frame(x=data$Space[,1],y=data$Space[,2],ll=factor(llmask,c("Background","Excitatory")),
                year=cut(data$Time,c(-1,365,365*2,365*3+2), labels=c("2010","2011","2012")))
g = ggplot(data=df,aes(x,y))
g = g + geom_point(aes(shape=llshapes,alpha=llalphas),size=2) 
g = g + facet_wrap(~ year)
g = g + scale_shape_identity() 
g = g + scale_alpha_identity() 
g = g+  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line.x=element_blank(),
    axis.line.y=element_blank(),
    axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      plot.background=element_blank())
g
ggsave("shots.pdf",g,width=12,height=6)

