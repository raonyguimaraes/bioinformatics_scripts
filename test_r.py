#test R in python

import scipy
import numpy as np
#Description: VCF summary
from optparse import OptionParser
import os
import rpy2.robjects as ro

from rpy2.robjects import r

import rpy2.robjects.numpy2ri


#R functions
plot = ro.r.plot
summary = ro.r.summary
table = ro.r.table
rnorm = ro.r.rnorm
dataf = ro.DataFrame({})





array_test = np.array(['1','2','3','4','5','6','7','3','2','1','4','2'])

array_size = array_test.size
#convert row to column
array_test = array_test.reshape(array_size,1)
#print array_test
tlist = ro.StrVector(array_test)
#print table(tlist)
keyWordArgs = {'row.names':ro.StrVector(("seed"))}
x = ro.r['as.data.frame'](table(tlist))
ro.r['print'](x)

#print table to plot

ro.r.library("plotrix")

ro.r('par(mar=c(5,4,4,18))')
r.assign('x',x)
ro.r('addtable2plot(x, cex=0.8, bty="o",display.rownames=F,hlines=T,vlines=T, title="Variantes com Phylop score")')



raw_input("Press ENTER to exit")
