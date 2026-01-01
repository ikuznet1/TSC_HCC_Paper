library(ggplot2)

data <- read.csv('~/Documents/TSC_Paper/240111_plasmaAFPELISA_analyzed.csv')

p <- ggplot(data, aes(x=Genotype, y=AFP)) + geom_boxplot()



#######################################
#############  FIGURE 3C  #############
#######################################

data <- read.csv('~/Documents/TSC_Paper/vv011_Key.csv')
data <- data[2:69,]

key <- data.frame(mouse=data$'Mouse..',genotype=data$'END')


data <- read.csv('~/Documents/TSC_Paper/CPC17100/CPC19000_Arany-Li_ER_04-23-2024.csv')
data <- data[1:84,]

a <- data$'X.3'
b <- data$'X.4'
c <- data$'X.5'
a <- !unlist(lapply(a,'str_detect','n/a'))
b <- !unlist(lapply(b,'str_detect','n/a'))
c <- unlist(lapply(c,'str_detect','tumor'))

onco <- data.frame(mouse = as.numeric(gsub("([0-9]+).*$","\\1",data[3:84,1])), image = data[3:84,1],adenoma = a[3:84],carcinoma = b[3:84],sarcoma = c[3:84])

worst_onco <- list()
for (i in unique(onco$mouse)){
	idx = match(onco$mouse,i)
	idx = !is.na(idx)

	worst_onco[i] <- TRUE %in% onco$carcinoma[idx]
	if (!worst_onco[[i]]){worst_onco[i] <- TRUE %in% onco$adenoma[idx]}

}

worst_onco <- unlist(worst_onco)

oncframe<-data.frame(mouse <- unique(onco$mouse),genotype <- key$genotype[match(unique(onco$mouse),key$mouse)],onco = worst_onco)
colnames(oncframe) <- c('mouse','genotype','onco')

onccounts<-table(oncframe$genotype,oncframe$onco)
#752 is a TSC KO mouse with a sarcoma but no HCC
#753 is a TSC KO mouse with a sarcoma but no HCC
oncperc <- onccounts / rowSums(onccounts)

toplot <- data.frame(rownames(oncperc),oncperc[,2])
colnames(toplot) <- c('geno','perc')
toplot$geno <- factor(toplot$geno,levels=c("WT","TSC KO","TFE3 KO","DKO "))

pdf('~/Documents/TSC_Paper/Oncobygeno_bar.pdf',width=5,height=5)
ggplot(toplot, aes(x=geno, y=perc)) +  geom_bar(stat="identity",width=0.6) + theme_classic() + xlab("Genotype") + ylab("Frequency") 
dev.off()

