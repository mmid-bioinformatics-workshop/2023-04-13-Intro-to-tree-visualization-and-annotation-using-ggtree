## Introduction to ggTree
## April 13th, 2023
## by Taylor Davedow
##********************************************************

##********************************************************
## Exercise 1 Answer ----
## HINT: use as.factor() around a continuous variable 
## to read as a discrete scale
##********************************************************


exercise1 <- gg_flip +
  geom_tiplab(aes(label = serovar),
              offset = 0.0001, 
              size = 5)+
  geom_tippoint(aes(color = as.factor(st)), 
                shape = 18,
                size = 4, 
                alpha = 0.5)+
  scale_color_manual(values = c("pink", "blue", "green", "orange", "purple"),
                     name = "Sequence Type")+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom",
        legend.spacing = unit(0, "cm"), 
        panel.border = element_blank(),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

##********************************************************
## END
##******************************************************** 