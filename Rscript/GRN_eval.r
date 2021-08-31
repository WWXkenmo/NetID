#################################################
##collect EPR and plot
blood_GENIE3_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/blood/500HVG_TFs_raw/GENIE3/NN_res.csv")
blood_GRNBOOST2_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/blood/500HVG_TFs_raw/GRNBOOST2/NN_res.csv")
blood_PIDC_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/blood/500HVG_TFs_raw/PIDC/NN_res.csv")

blood_GENIE3_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/blood/500HVG_TFs_netID/GENIE3/NN_res.csv")
blood_GRNBOOST2_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/blood/500HVG_TFs_netID/GRNBOOST2/NN_res.csv")
blood_PIDC_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/blood/500HVG_TFs_netID/PIDC/NN_res.csv")


EPR <- data.frame(EPR = log2(c(blood_GENIE3_raw$EPR, blood_GRNBOOST2_raw$EPR, blood_PIDC_raw$EPR,
                          blood_GENIE3_netID$EPR, blood_GRNBOOST2_netID$EPR, blood_PIDC_netID$EPR)+1),
                  Methods = c(rep("raw",18),rep("netID",18)), Network = c(blood_GENIE3_raw$X, blood_GRNBOOST2_raw$X, blood_PIDC_raw$X,
                          blood_GENIE3_netID$X, blood_GRNBOOST2_netID$X, blood_PIDC_netID$X),
                  GRN = c(rep(c(rep("GENIE3",6),rep("GRNBOOST2",6),rep("PIDC",6)),2)))

p1 <- ggplot(data = log2(EPR+1), aes(x=Methods, y=EPR, col=Network, shape=Network)) +
  geom_point(size=2) +
  theme_bw(base_size = 14) +
  facet_wrap(~GRN, ncol=4) +
  labs(y = "EPR", 
       x = "Methods") +
  scale_color_manual(name="Network", values=c("#46b4d3", "#34393c", "#db494b", "#eba432", "#86aa6d", "#ac3d1f")) +
  scale_shape_manual(name="Network", values=c(0, 1, 2, 3, 4, 8)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))




ESC_GENIE3_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/ESC/500HVG_TFs_raw/GENIE3/NN_res.csv")
ESC_GRNBOOST2_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/ESC/500HVG_TFs_raw/GRNBOOST2/NN_res.csv")
ESC_PIDC_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/ESC/500HVG_TFs_raw/PIDC/NN_res.csv")

ESC_GENIE3_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/ESC/500HVG_TFs_netID/GENIE3/NN_res.csv")
ESC_GRNBOOST2_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/ESC/500HVG_TFs_netID/GRNBOOST2/NN_res.csv")
ESC_PIDC_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/ESC/500HVG_TFs_netID/PIDC/NN_res.csv")

EPR <- data.frame(EPR = log2(c(ESC_GENIE3_raw$EPR, ESC_GRNBOOST2_raw$EPR, ESC_PIDC_raw$EPR,
                          ESC_GENIE3_netID$EPR, ESC_GRNBOOST2_netID$EPR, ESC_PIDC_netID$EPR)+1),
                  Methods = c(rep("raw",18),rep("netID",18)), Network = c(ESC_GENIE3_raw$X, ESC_GRNBOOST2_raw$X, ESC_PIDC_raw$X,
                          ESC_GENIE3_netID$X, ESC_GRNBOOST2_netID$X, ESC_PIDC_netID$X),
                  GRN = c(rep(c(rep("GENIE3",6),rep("GRNBOOST2",6),rep("PIDC",6)),2)))

p2 <- ggplot(data = EPR, aes(x=Methods, y=EPR, col=Network, shape=Network)) +
  geom_point(size=2) +
  theme_bw(base_size = 14) +
  facet_wrap(~GRN, ncol=4) +
  labs(y = "EPR", 
       x = "Methods") +
  scale_color_manual(name="Network", values=c("#46b4d3", "#34393c", "#db494b", "#eba432", "#86aa6d", "#ac3d1f")) +
  scale_shape_manual(name="Network", values=c(0, 1, 2, 3, 4, 8)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))



HSC_GENIE3_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/HSC/500HVG_TFs_raw/GENIE3/NN_res.csv")
HSC_GRNBOOST2_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/HSC/500HVG_TFs_raw/GRNBOOST2/NN_res.csv")
HSC_PIDC_raw <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/HSC/500HVG_TFs_raw/PIDC/NN_res.csv")

HSC_GENIE3_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/HSC/500HVG_TFs_netID/GENIE3/NN_res.csv")
HSC_GRNBOOST2_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/HSC/500HVG_TFs_netID/GRNBOOST2/NN_res.csv")
HSC_PIDC_netID <- read.csv("/mnt/data1/weixu/NetID/Beeline-master/outputs/HSC/500HVG_TFs_netID/PIDC/NN_res.csv")


EPR <- data.frame(EPR = log2(c(HSC_GENIE3_raw$EPR, HSC_GRNBOOST2_raw$EPR, HSC_PIDC_raw$EPR,
                          HSC_GENIE3_netID$EPR, HSC_GRNBOOST2_netID$EPR, HSC_PIDC_netID$EPR)+1),
                  Methods = c(rep("raw",18),rep("netID",18)), Network = c(HSC_GENIE3_raw$X, HSC_GRNBOOST2_raw$X, HSC_PIDC_raw$X,
                          HSC_GENIE3_netID$X, HSC_GRNBOOST2_netID$X, HSC_PIDC_netID$X),
                  GRN = c(rep(c(rep("GENIE3",6),rep("GRNBOOST2",6),rep("PIDC",6)),2)))

p3 <- ggplot(data = EPR, aes(x=Methods, y=EPR, col=Network, shape=Network)) +
  geom_point(size=2) +
  theme_bw(base_size = 14) +
  facet_wrap(~GRN, ncol=4) +
  labs(y = "EPR", 
       x = "Methods") +
  scale_color_manual(name="Network", values=c("#46b4d3", "#34393c", "#db494b", "#eba432", "#86aa6d", "#ac3d1f")) +
  scale_shape_manual(name="Network", values=c(0, 1, 2, 3, 4, 8)) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10))
