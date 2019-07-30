library(ggplot2)

data <- read.csv(file = "/Users/messyada/Desktop/OSU/committee_courses/NRT/research/viral_bac_counts/bac_vir_stats.csv", header = TRUE)
View(data)

ggplot(data = data, aes(x=Reef.Type, fill=factor(Microbe), y=Counts)) + geom_boxplot() + scale_fill_manual(values=c("#ffffff", "#000000"), name = "Microbe") + facet_wrap(~Location) + xlab("Reef Type") + ylab("Counts") + ggtitle("Viral and Bacterial Counts Across Reef Types and Island Location") + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12), plot.title = element_text(size=14,hjust = 0.5))

mypath <- file.path("/Users/messyada/Desktop/",paste("bac_vir_boxplot", ".png", sep = "")) 
png(filename=mypath, width=13, height=15, units="in", res=300, pointsize=20)

ggplot(data = data, aes(x=Reef.Type, fill=factor(Microbe), y=Counts)) + geom_boxplot() + scale_fill_manual(values=c("#ffffff", "#000000"), name = "Microbe") + facet_wrap(~Location) + xlab("Reef Type") + ylab("Counts") + ggtitle("Viral and Bacterial Counts Across Reef Types and Island Location") + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12), plot.title = element_text(size=14,hjust = 0.5))

#MAKE SURE TO SPLIT UP BY SIDE OF ISLAND

