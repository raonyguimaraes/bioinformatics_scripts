"""
This script will predict the gender of your samples by looking at X-Chromosome 
data and requires that a marker map be applied to column name headers of the spreadsheet.

It is required that you only have the X-Chromosome data active in your spreadsheet before you 
run the script. To activate only the X-Chromosome, inactivate all chromosomes by going to 
Select > Column > Inactivate All Columns; then go to Select > Activate By Chromosomes, and 
uncheck all chromosomes except X.

Author: Dr. Bo Peng, MD Anderson Cancer Center, Carol Hall and Greta Linse, Golden Helix, Inc.
Last updated: 2-19-2009
"""
def hetero(row):
    'Calculate heterozygosity of this row'
    het = 0
    missing = 0
    for item in row:
        if item == None or len(item)!=3:
            missing += 1
        elif item[0] != item[2]:
            het += 1
    if len(row) == missing:
        return (0,0,missing,1)
    else:
        return (het, het*1.0/(len(row)-missing),missing, missing*1.0/len(row))

def val(v):
    if v < 0.1:
        return 1
    else:
        return 2
        
def geno_gender_check():
    ss = ghi.getCurrentObject()
    nR = ss.numRows()
    nC = ss.numCols()
    activeRows = []
    activeCols = []
    numRows = 0
    numCols = 0
    
    count0 = 0
    progress0 = ghi.progressDialog("Scanning Spreadsheet", nR+nC)
    
    if not ss.hasMarkerMap():
        ghi.message("The spreadsheet must have a marker map applied.")
        return
        
    for r in range(1,nR+1):
        if ss.getRowState(r)==1:
            numRows += 1
            activeRows.append(r)
        if progress0.wasCanceled()==True:
            return
        count0 += 1
        progress0.setProgress(count0)
        
    if len(activeRows) == 0:
        ghi.message("You need to have at least one active row.")
        return
    chrName = ss.getMarkerMapChromosome()
    
    for c in range(1,nC+1):
        if ss.getColType(c)== ghi.const.TypeGenotypic and ss.getColState(c)!= ghi.const.StateInactive:
            if chrName[c-1] in ['X','x',23,'23']:
                numCols += 1
                activeCols.append(c)
        if progress0.wasCanceled():
            return
        count0 += 1
        progress0.setProgress(count0)
        
    if len(activeCols) == 0:
        ghi.message("You need to have some active columns from the X chromosome.")
        return
        
    progress0.finish()

    count1 = 0
    progress1 = ghi.progressDialog("Calculating Heterozygosity and Missing Rates", numRows)
    
    newSS = ghi.dataSetBuilder("Heter_Sex", numRows)
    newSS.addRowLabels('WGS Number', [ss.cell(rr,0) for rr in activeRows])
    
    progress1.show()
    
    het = [0 for i in range(numRows)]
    missing = [0 for i in range(numRows)]
    het_rate = [0.0 for i in range(numRows)]
    missing_rate = [0.0 for i in range(numRows)]
    Sex_cat = [0 for i in range(numRows)]
    Sex_bin = [0 for i in range(numRows)]
    
    if progress1.wasCanceled()==True:
        return
    
    index = 0
    
    for r in activeRows:
        row = []
        for c in activeCols:
            row.append(ss.cell(r,c))
        h, hr, m, mr = hetero(row)
        het[index] = h
        het_rate[index] = hr
        missing[index] = m
        missing_rate[index] = mr
        index += 1
        
        if progress1.wasCanceled()==True:
            return
        count1 += 1
        progress1.setProgress(count1)
    
    progress1.finish()
    
    count2 = 0
    progress2 = ghi.progressDialog("Imputing Gender", len(het_rate))
    
    Sex = [val(x) for x in het_rate]
    
    if progress2.wasCanceled()==True:
        return
        
    for i in range(len(Sex)):
        if Sex[i] == 1:
            Sex_cat[i] = 'M'
            Sex_bin[i] = 0
        elif Sex[i] == 2:
            Sex_cat[i] = 'F'
            Sex_bin[i] = 1
        if progress2.wasCanceled()==True:
            return
        count2 += 1
        progress2.setProgress(count2)
        
    progress2.finish()
    
    newSS.addIntColumn('Heterozygosity', het)
    newSS.addRealColumn('Heterozygosity Rate', het_rate)
    newSS.addIntColumn('Missing', missing)
    newSS.addRealColumn('Missing Rate', missing_rate)
    newSS.addCategoricalColumn('Sex', Sex_cat)
    newSS.addBoolColumn('Sex', Sex_bin)
        
    newSS.finish(ss.getID())
    
try:
    geno_gender_check()
except:
    ghi.error()