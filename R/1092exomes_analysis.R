kinship0_1000 <- '/home/raony/mendel/projects/1000genomes/integrated_call_sets/1092exomes/1092exomes_sureselect.allinfo.kin0'

data <- read.table(kinship0_1000, head=T)
head(data)



data$pop1 <- gsub('.+_', '', data$ID1)
data$pop2 <- gsub('.+_', '', data$ID2)

data$poppop <- paste(data$pop1, '_',data$pop2, sep='')
str(data)
list_pop <- split(data,data$poppop)

mean_p <- lapply(list_pop,function(df.f)mean(as.numeric(df.f$Kinship)))
head(mean_p)

order_meanp <- data.frame(pop=names(mean_p), mean=unlist(mean_p))
head(order_meanp)
order_meanp <- order_meanp[order(order_meanp$mean),]

list_pop <- list_pop[order_meanp$pop]
names(list_pop)

boxplot(lapply(list_pop,function(df.f)df.f$Kinship), las=2, cex.axis=0.5)
abline(h=0, lty=2)
abline(v=1:300, lty=2, col='gray')
?abline

boxplot(split(as.numeric(data$V8),data$poppop), las=2, cex.axis=0.5)
boxplot(split(as.numeric(data$V8),data$poppop), las=2, cex.axis=0.5)

?boxplot
