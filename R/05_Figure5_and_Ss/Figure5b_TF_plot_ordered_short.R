####run script for 5a
df_sorted_gyn_combined=df_sorted_gyn
df_sorted_gyn_combined$speices=rep('gyn',nrow(df_sorted_gyn_combined))

df_sorted_bid_combined=df_sorted_bid
df_sorted_bid_combined$speices=rep('bid',nrow(df_sorted_bid_combined))

df_sorted_combined=rbind(df_sorted_gyn_combined,df_sorted_bid_combined)

gyn_TF_plot_order=levels(df_sorted_gyn$TF)
bid_TF_plot_order=levels(df_sorted_bid$TF)
#matching_TFs=intersect(gyn_TF_plot_order,bid_TF_plot_order)
bid_only_TFs= bid_TF_plot_order[!bid_TF_plot_order %in% gyn_TF_plot_order ]
#gyn_only_TFs= gyn_TF_plot_order[!gyn_TF_plot_order %in% bid_TF_plot_order ]

new_order_bid=c(bid_only_TFs,gyn_TF_plot_order)

df_sorted_combined$TF=factor(df_sorted_combined$TF,levels = new_order_bid)
df_sorted_combined$speices=factor(df_sorted_combined$speices,levels = c('gyn','bid'))

ggplot(df_sorted_combined,aes(x=Freq,y=TF,fill=color))+
  geom_bar(stat='identity',position = position_dodge(width = 0.9), width = 0.7)+
  scale_fill_manual(values=c("#238443", "#1d91c0"))+facet_grid(.~speices)+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y =  element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


ggplot(df_sorted_combined,aes(x=Freq,y=TF,fill=color))+
  geom_bar(stat='identity',position = position_dodge(width = 0.9), width = 0.7)+
  scale_fill_manual(values=c("#238443", "#1d91c0"))+facet_grid(color~speices)+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y =  element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

##############
library(dplyr)
df_sorted_combined_filter <- df_sorted_combined %>%
    group_by(TF, color) %>%
    # Conditionally apply filters based on the 'color' column
    filter((color == "bs" & !all(Freq_bs == 0)) | 
             (color == "meso" & !all(Freq_meso == 0)) | 
             (color != "bs" & color != "meso")) %>%
    ungroup()  # Ungroup after filtering
  
df_sorted_combined_filter$color=factor(df_sorted_combined_filter$color,levels=c( "bs", "meso" ))
ggplot(df_sorted_combined_filter,aes(x=Freq,y=TF,fill=color))+
    geom_bar(stat='identity',position = position_dodge(width = 0.9), width = 0.7)+
    scale_fill_manual(values=rev(c("#238443", "#1d91c0")))+facet_grid(color~speices)+
    theme_bw() +
    theme(axis.line = element_line(colour = "grey"),
          axis.text.x = element_text(size = 8, face = "bold"),
          axis.text.y = element_text(size = 8, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y =  element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  

  
##########
df_sorted_combined_bs=df_sorted_combined[ which(df_sorted_combined$color=='bs'),]
df_sorted_combined_bs_filtered <- df_sorted_combined_bs %>%
  group_by(TF) %>%
  filter(!all(Freq_bs == 0)) %>%
  ungroup()

not_kept_bs_1=as.character(df_sorted_combined_bs_filtered[which(df_sorted_combined_bs_filtered$Freq_bs==0),]$TF)
not_kept_bs_2=as.character(names(table(df_sorted_combined_bs_filtered$TF)[which(table(df_sorted_combined_bs_filtered$TF)==1)]))
not_kept_bs=c(not_kept_bs_1,not_kept_bs_2)
df_sorted_combined_bs_filtered=df_sorted_combined_bs_filtered[which(! df_sorted_combined_bs_filtered$TF %in% not_kept_bs),]


p1=ggplot(df_sorted_combined_bs_filtered,aes(x=Freq,y=TF,fill=color))+
  geom_bar(stat='identity',position = position_dodge(width = 0.9), width = 0.7)+
  scale_fill_manual(values=c( "#1d91c0"))+facet_grid(.~speices)+xlim(0, 8)+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y =  element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#bdbdbd"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


##########
df_sorted_combined_meso=df_sorted_combined[ which(df_sorted_combined$color=='meso'),]
df_sorted_combined_meso_filtered <- df_sorted_combined_meso %>%
  group_by(TF) %>%
  filter(!all(Freq_meso == 0)) %>%
  ungroup()%>%
  arrange(
    desc(speices == "gyn"),  # Sort 'gyn' first
    desc(Freq_meso),         # Sort Freq_meso in descending order
    desc(speices == "bid")   # Then sort 'bid'
  )

df_sorted_combined_meso_filtered$TF=factor(df_sorted_combined_meso_filtered$TF,levels=rev(as.character(unique(df_sorted_combined_meso_filtered$TF))))

not_kept_meso_1=as.character(df_sorted_combined_meso_filtered[which(df_sorted_combined_meso_filtered$Freq_meso==0),]$TF)
not_kept_meso_2=as.character(names(table(df_sorted_combined_meso_filtered$TF)[which(table(df_sorted_combined_meso_filtered$TF)==1)]))
not_kept_meso=c(not_kept_meso_1,not_kept_meso_2)
df_sorted_combined_meso_filtered=df_sorted_combined_meso_filtered[which(! df_sorted_combined_meso_filtered$TF %in% not_kept_meso),]

p2=ggplot(df_sorted_combined_meso_filtered,aes(x=Freq,y=TF,fill=color))+
  geom_bar(stat='identity',position = position_dodge(width = 0.9), width = 0.7)+
  scale_fill_manual(values=c("#238443"))+facet_grid(.~speices)+xlim(0, 8)+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y =  element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#bdbdbd"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


p2/p1
