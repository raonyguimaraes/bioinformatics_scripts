### Load the package or install if not present
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

### Set the display a 2 by 2 grid
par(mfrow=c(2,2))

### Show all the colour schemes available
display.brewer.all()

### Generate random data matrix
rand.data <- replicate(8,rnorm(100,100,sd=1.5))

### Draw a box plot, with each box coloured by the 'Set3' palette
boxplot(rand.data,col=brewer.pal(8,"Set3"))

### Draw plot of counts coloured by the 'Set3' pallatte
br.range <- seq(min(rand.data),max(rand.data),length.out=10)
results <- sapply(1:ncol(rand.data),function(x) hist(rand.data[,x],plot=F,br=br.range)$counts)
plot(x=br.range,ylim=range(results),type="n",ylab="Counts")
cols <- brewer.pal(8,"Set3")
lapply(1:ncol(results),function(x) lines(results[,x],col=cols[x],lwd=3))

### Draw a pie chart
table.data <- table(round(rand.data))
cols <- colorRampPalette(brewer.pal(8,"Dark2"))(length(table.data))
pie(table.data,col=cols)