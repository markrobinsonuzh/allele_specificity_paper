
#Make diagram to explain rationale of weight
library(ggplot2)


#We show 4 cases

#### Score plots ####
MM <- c(4,8,0,0)
UU <- c(4,0,0,4)
UM <- c(0,0,4,4)
MU <- c(0,0,4,0)
eps <- 1

calc_weight <- function(MM, UU, beta=0.5, a=.2) {
  s1 <- beta+MM
  s2 <- beta+UU
  stats::pbeta(.5+a, shape1=s1, shape2=s2)-stats::pbeta(.5-a,
                                                        shape1=s1, shape2=s2)
}
w <- calc_weight(MM,UU)
withw <- abs(log10(((MM+eps)*(UU+eps) ) / ( (MU+eps)*(UM+eps) ) ) * w)
without <- abs(log10(((MM+eps)*(UU+eps) ) / ( (MU+eps)*(UM+eps) ) )  )
x <- data.frame(ASMtuple = c(without, withw), 
                Weight =c(rep("No",4),rep("Yes",4)),
                Case = rep(1:4,2),
                Tuple = 1.5)

p2a <- ggplot(x[x$Weight == "No",]) +
  geom_point(aes(Tuple,ASMtuple), size = 6, color="red") +
  geom_segment(aes(x=Tuple, xend=Tuple, y=0, yend=ASMtuple),
               size = 1, color="red") +
  ylab("|log-odds|") +
  facet_grid(~Case) +
  theme_bw() +
  scale_y_continuous(limits=c(-0.1,1.6)) +
  theme(text = element_text(size = 10), 
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_blank()) 

p2b <- ggplot(x[x$Weight == "Yes",]) +
  geom_point(aes(Tuple,ASMtuple), color = "blue", size = 6) +
  geom_segment(aes(x=Tuple, xend=Tuple, y=0, yend=ASMtuple),
               size = 1, color="blue") +
  ylab("|log-odds*wi|") + 
  facet_grid(~Case) +
  theme_bw() +
  scale_y_continuous(limits=c(-0.1,1.6)) +
  theme(text = element_text(size = 10), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_blank()) 


#### plot PDF for each weight ####

quant <- seq(0,1, 0.05)
cdfs <- c(dbeta(quant, 4.5,4.5),
       dbeta(quant, 8.5,0.5),
       dbeta(quant, 0.5,0.5),
       dbeta(quant, 0.5,4.5))

x <- data.frame(q = rep(quant,4),
                PDF = cdfs,
                Case = rep(1:4,each = length(quant)))

p3 <- ggplot(x, aes(q,PDF)) +
  geom_area(aes(x = ifelse(q>0.30 & q< 0.70 , q, NA)), fill = "pink") +
  geom_line(size = 1.5) + 
  facet_grid(~Case) + 
  theme_bw() +
  theme(text = element_text(size = 10), 
        strip.text.x = element_blank()) +
  #xlab(expression("0.5"*epsilon)) + 
  ylab("PDF") 
  #xlab("")

#### Reads ####
xpoint <- data.frame(CpG=c(1:2),
                     value=
                       c(
                       rep(c(21,21,19,19,19,21,21,19),2),
                       rep(c(21,21,19,19,19,21,21,19),2),
                       rep(c(19,19,19,19,21,19,21,21),2),
                       rep(c(19,19,19,19,21,19,21,21),2)),
                     read=rep(1:8,each = 8),
                     case=rep(1:4, each = 2),
                     stringsAsFactors = FALSE)
xread <- data.frame(read = 1:8,
                    start = 1,
                    end = 2)


p1 <- ggplot() +
  scale_shape_identity() +
  theme_bw() +
  theme(text = element_text(size = 10), 
        axis.text=element_blank(),
        #axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_segment(data=xread,aes(x=start-0.5,y=read,xend=end+0.5,yend=read), color= "gray") +
  geom_point(data=xpoint, aes_(x=~CpG, y=~read, shape=~value), fill="white",
             size=5) +
  scale_x_continuous(limits = c(-1,4)) +
  scale_y_continuous(limits=c(0,9)) +
  facet_grid(~case) +
  ylab("Read") +
  xlab("Tuple")

#### grid-plot ####

p4 <- cowplot::plot_grid(p1, p3, p2a,p2b, ncol=1, nrow = 4, 
                         rel_heights = c(2, 1.5,1,1),
                         labels = c("A","B","C",""),
                         align = "v",
                         axis = "lr")

p4

ggplot2::ggsave("curvesNscatters/weight_diagram.png",p4,
                width = 7, height = 8)
