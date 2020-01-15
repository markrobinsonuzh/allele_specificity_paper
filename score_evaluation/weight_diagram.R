
#Make diagram to explain rationale of weight
library(ggplot2)


#We show 4 cases
#bottom
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
p2 <- ggplot(x) +
  geom_point(aes(Tuple,ASMtuple,color = Weight), size = 6) +
  scale_color_manual(values = c("red","blue")) +
  facet_grid(~Case) +
  theme_bw() +
  scale_y_continuous(limits=c(-0.1,1.6)) +
  theme(text = element_text(size = 10), 
        axis.text.x=element_blank(),
        #axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_blank()) 

#plot CDF (and PDF?) for each weight

quant <- seq(0.3,0.7, 0.1)
cdfs <- c(pbeta(quant, 4.5,4.5),
       pbeta(quant, 8.5,0.5),
       pbeta(quant, 0.5,0.5),
       pbeta(quant, 0.5,4.5))

x <- data.frame(q = rep(quant,4),
                CDF = cdfs,
                Case = rep(1:4,each = 5))

p3 <- ggplot(x) +
  geom_line(aes(q,CDF), size = 2) +
  facet_grid(~Case) + 
  theme_bw() +
  theme(text = element_text(size = 10), 
        strip.text.x = element_blank()) +
  xlab(expression("0.5-+"*epsilon))

#top
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


p1 <- ggplot() +
  scale_shape_identity() +
  theme_bw() +
  theme(text = element_text(size = 10), 
        axis.text=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_point(data=xpoint, aes_(x=~CpG, y=~read, shape=~value), fill="white",
             size=5) +
  scale_x_continuous(limits = c(-1,4)) +
  scale_y_continuous(limits=c(0,9)) +
  facet_grid(~case) +
  ylab("Read")

p4 <- cowplot::plot_grid(p1, p2, p3, ncol=1, nrow = 3, 
                         labels = c("A","B","C"),
                         align = "v",
                         axis = "lr")
p4

ggplot2::ggsave("curvesNscatters/weight_diagram.png",p4,
                width = 6, height = 6.5)
