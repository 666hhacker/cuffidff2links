import re
import cookielib, urllib2, urllib
import poster
import sys
import os
from sets import Set
import csv

def EnrichrLink(genesStr,clusterInfo=''):
    #post a gene list to enrichr server and get the link.
    cj = cookielib.CookieJar()
    opener = poster.streaminghttp.register_openers()
    opener.add_handler(urllib2.HTTPCookieProcessor(cookielib.CookieJar()))

    params = {'list':genesStr,'description':clusterInfo}
    datagen, headers = poster.encode.multipart_encode(params)
    url = "http://amp.pharm.mssm.edu/Enrichr/enrich"
    request = urllib2.Request(url, datagen,headers)
    urllib2.urlopen(request)


    x = urllib2.urlopen("http://amp.pharm.mssm.edu/Enrichr/share")
    responseStr = x.read()
    splitPhrases = responseStr.split('"')
    linkID = splitPhrases[3]
    shareUrlHead = "http://amp.pharm.mssm.edu/Enrichr/enrich?dataset="
    enrichrLink = shareUrlHead + linkID
    return enrichrLink

def sig2pcainput(groupPath,genePlusId,fileCategory):
    sigTest_id = [item[1] for item in genePlusId]
    with open(groupPath) as cf:
        cder = csv.reader(cf,delimiter='\t') #cder: abbr. for csv reader
        header = cder.next()
        trackingIdIdx = header.index('tracking_id')
        conditionIdx = header.index('condition')
        FPKMIdx = header.index('FPKM')
        replicateIdx = header.index('replicate')
        group = [((row[trackingIdIdx],row[conditionIdx],row[replicateIdx]),row[FPKMIdx])
			for row in cder if row[trackingIdIdx] in sigTest_id]
    group = dict(group)

    label = [(item[1],item[2]) for item in group.keys()]
    label = list(Set(label))
    label.sort()

    header = '\t'.join([perLabel[0] for perLabel in label])
    header = 'Genes\t' + header

    wStr = ''
    for perGene in genePlusId:
        line = '\t'.join([group[(perGene[1],perLabel[0],perLabel[1])]
         if (perGene[1],perLabel[0],perLabel[1])in group.keys()
         else '' for perLabel in label])
        wStr = wStr + perGene[0] + '\t' + line + '\n'
    wStr = header + '\n' + wStr
    dirName = os.path.dirname(groupPath)
    with open(os.path.join(dirName,'EnrichrLinks', fileCategory + '_pcainput.txt'),'w') as f:
        f.write(wStr)

        
def readperdiff(path):
    logStr = '';
    print '\n'
    logStr = logStr + '\n'
    #ignore empty file
    if os.path.getsize(path) < 1000:
        print 'Empty file:'+path+'\n'
        return 'Empty file:'+path+'\n'

    print 'Processing:'+path
    logStr = logStr + 'Processing:'+path + '\n'
    
    baseName = os.path.basename(path)
    fileName = os.path.splitext(baseName)[0]
    dirName = os.path.dirname(path)
    if dirName == '':
	dirName = os.getcwd()
    writeDir = os.path.join(dirName,'EnrichrLinks')
    if not os.path.exists(writeDir):
	os.makedirs(writeDir)

    f = open(path)
    content = f.read()
    f.close()
    lines = re.split('\n',content)
    linesCount = len(lines)

    #remove vacant lines at the tail of the file if any
    while len(lines[linesCount-1])<100:
	lines.pop()
	linesCount = linesCount-1
	if linesCount == 0:
	    break

    
    header = lines[0]
    dat = [lines[i] for i in range(1,linesCount)]

    header = re.split('\t',header)
    dat = [re.split('\t',datum) for datum in dat]

    #Filter out significant rows
    sigMetaIdx = header.index('significant')
    sigDat = [ row for row in dat if row[sigMetaIdx] == 'yes']
    geneMetaIdx = header.index('gene')

    #Generate Matlab PCA input for this .diff file
	#test if its group tracking file exist first
    underscoreIdx = fileName.rfind('_')
    if underscoreIdx == -1:
        perLogStr = 'invalid .diff file. name should be xxx_exp.diff format'
        print perLogStr
        logStr = logStr + perLogStr
    else:
        fileCategory = fileName[:underscoreIdx]
        if(fileCategory[len(fileCategory)-1] == 's'):
            groupTrackFileName = fileCategory+'.read_group_tracking'
        else:
            groupTrackFileName = fileCategory+'s.read_group_tracking'
        groupTrackPath = os.path.join(dirName,groupTrackFileName)
        if(os.path.exists(groupTrackPath)):
           test_idMetaIdx = header.index('test_id')
           sigTest_id = [(row[geneMetaIdx],row[test_idMetaIdx]) for row in sigDat]
           sigTest_idUnique = list(Set(sigTest_id))
           sigTest_idUnique.sort()
           sig2pcainput(groupTrackPath,sigTest_idUnique,fileCategory)
        else:
           perLogStr = ('read_group_tracking file dose not exist for "' + path + 
                    '". FPKM values are unavailable and no PCA input file will be generated')
           print perLogStr
           logStr = logStr + perLogStr
       

    #seperate significant rows into different comparison groups
    sample1MetaIdx = header.index('sample_1')
    sample2MetaIdx = header.index('sample_2')
    comparisons = []
    comparisonsDat = []
    for row in sigDat:
	rowComparison = [row[sample1MetaIdx],row[sample2MetaIdx]]
	if rowComparison not in comparisons:
	    comparisons.append(rowComparison)
	    comparisonsDat.append([])
	comparisonsIdx = comparisons.index(rowComparison)
	comparisonsDat[comparisonsIdx].append(row)

    #seperate rows in each comparison goup into up and down subgroup
    try:
	log2FoldMetaIdx = header.index('log2(fold_change)')
    except ValueError:
	return logStr +'"' + path + '"' + ' is not a valid file: log2(fold_change) field does not find\n'
    comparisonsUpDat = [[row for row in perComparisonDat
			 if float(row[log2FoldMetaIdx])>0]
			for perComparisonDat in comparisonsDat]
    comparisonsDownDat = [[row for row in perComparisonDat
			   if float(row[log2FoldMetaIdx])<0]
			  for perComparisonDat in comparisonsDat]

    #extract genes in each subgroup
    comparisonsUpGenes = [ [row[geneMetaIdx] for row in perComparisonUpDat]
			   for perComparisonUpDat in comparisonsUpDat]
    comparisonsDownGenes = [[row[geneMetaIdx] for row in perComparisonDownDat]
			    for perComparisonDownDat in comparisonsDownDat]


    #enrichr post and write the response links into txt file
    comparisonsCount = len(comparisons)
    wStr = ''
    for i in range(comparisonsCount):
	upInfo = comparisons[i][0]+','+comparisons[i][1]+'\tUp Genes'
	if len(comparisonsUpGenes[i])==0:
	    upLink = ''
	else:
	    for j in range(len(comparisonsUpGenes[i])):
		perGene = comparisonsUpGenes[i][j]
		match = re.search(r'[\w\-@./]+',perGene)
		if match.group(0) is not perGene:
		    print '"' + perGene + '"' + ' is not a valid gene symbol and converted to ' + '"' + match.group(0) + '"'
		    logStr = logStr + '"' + perGene + '"' + ' is not a valid gene symbol and converted to ' + '"' + match.group(0) + '"\n'

		    comparisonsUpGenes[i][j] = match.group(0)
		    
	    upLink = EnrichrLink('\n'.join(comparisonsUpGenes[i]), upInfo)
	
	downInfo = comparisons[i][0]+','+comparisons[i][1]+'\tDown Genes'
	if len(comparisonsDownGenes[i]) == 0:
	    downLink = ''
	else:
	    for j in range(len(comparisonsDownGenes[i])):
		perGene = comparisonsDownGenes[i][j]
	       
		match = re.search(r'[\w\-@./]+',perGene)
		if match.group(0) is not perGene:
		    print '"' + perGene + '"' + ' is not a valid gene symbol and converted to ' + '"' + match.group(0)+ '"'
		    logStr = logStr + '"' + perGene + '"' + ' is not a valid gene symbol and converted to ' + '"' + match.group(0)+ '"\n' 

		    comparisonsDownGenes[i][j] = match.group(0)
		    
	    downLink = EnrichrLink('\n'.join(comparisonsDownGenes[i]), downInfo)
	
	wStr = wStr + upInfo+'\t'+ upLink +'\n'
	wStr = wStr + downInfo+ '\t' + downLink+'\n'
		
    f = open(os.path.join(writeDir,fileName+'_enrichrLinks.txt'),'w')
    f.write(wStr)
    f.close()


    #below is printing up and down genes of each comparison to a single txt file
    rowMax = 0
    for perComparisonUpGenes in comparisonsUpGenes:
	if rowMax < len(perComparisonUpGenes):
	    rowMax = len(perComparisonUpGenes)
	    
    for perComparisonDownGenes in comparisonsDownGenes:
	if rowMax < len(perComparisonDownGenes):
	    rowMax = len(perComparisonDownGenes)

    colMax = 2*len(comparisons)
    wStr = ''
    for i in range(colMax):
	if i%2 == 0:
	    wStr = wStr + comparisons[i/2][0]+','+ comparisons[i/2][1]+'\t'
	else:
	    wStr = wStr + '\t'
    wStr.rstrip('\t')
    wStr =wStr + '\n'

    secondRow = '\t'.join(['Up Genes','Down Genes']*(colMax/2))
    wStr = wStr + secondRow + '\n'

    for i in range(rowMax):
	for j in range(len(comparisons)):
	    try :
		wStr = wStr + comparisonsUpGenes[j][i] + '\t'
	    except IndexError:
		wStr = wStr + '\t'
	    try :
		wStr = wStr + comparisonsDownGenes[j][i] + '\t'
	    except IndexError:
		wStr = wStr + '\t'
	wStr.rstrip('\t')
	wStr =wStr + '\n'

    f = open(os.path.join(writeDir,fileName+'_updown.txt'),'w')
    f.write(wStr)
    f.close()
    return logStr
    
def readdiff(path,regPattern=r'.+\.diff'):
    
    if os.path.isfile(path):
	 strpart = readperdiff(path)
	 dirName = os.path.dirname(path);
	 with open(os.path.join(dirName,r'Enrichr_log.txt'),'w') as f:
	    f.write(strpart)
	 print '\nlog file is at' + '"'+ dirName + '"' + '\n'
    else:
	logStr = '';
	for directory, dirnames,filenames in os.walk(path):
	    strparts = [readperdiff(os.path.join(directory,perFile)) for perFile in filenames
	       if re.match(regPattern,perFile) != None]
	    logStr = logStr + ' '.join(strparts)
	with open(os.path.join(path,r'Enrichr_log.txt'),'w') as f:
	    f.write(logStr)
	print '\nlog file is at' + '"'+ path + '"' + '\n'

def mergeFiles(pathx,pathxs,identifier):
    folderName = os.path.basename(pathx)
    resDir = os.path.join(pathx,'mergedFiles')
    if(not os.path.exists(resDir)):
        os.mkdir(resDir)
    wFile = os.path.join(resDir,folderName+'_'+identifier)
    if(os.path.exists(wFile)):
        os.remove(wFile)
    f = open(wFile,'a')
    logStr = ''
    for perPath in pathxs:
        perLogStr = perPath + '\n'
        print perLogStr
        logStr = logStr + perLogStr
        with open(perPath) as rf:
                content = rf.read()
        while(content[len(content)-1] != '\n'):
                content.pop()
        f.write(content)
    f.close()
    return logStr
        
def merge(pathx,regPattern=r'.+\.diff'):
    lonelyDiff = []
    coupleDiff = []
    global logStr
    logStr = ''
    def logPrintInfo(perLogStr):
        print perLogStr
        global logStr
        logStr = logStr + perLogStr
    for directory,dirnames,filenames in os.walk(pathx):
        for perFile in filenames:
            if(directory!='mergedFiles' and re.match(regPattern,perFile) != None):
                fileNameSplits = perFile.split('_')
                fileCategory = fileNameSplits[len(fileNameSplits)-2]
                if(fileCategory[len(fileCategory)-1] == 's'):
                    trackFile = fileCategory+'.read_group_tracking'
                else:
                    trackFile = fileCategory+'s.read_group_tracking'
                if(trackFile in filenames):
                    coupleDiff.append((directory,perFile,trackFile))
                else:
                    lonelyDiff.append((directory,perFile))
    if(not lonelyDiff):
        logPrintInfo('All .diff files have their read_group_tracking files\n')
        diffFiles = [os.path.join(item[0],item[1]) for item in coupleDiff]
        logPrintInfo('Now merging all the .diff files...\n')
        logStr = logStr + mergeFiles(pathx,diffFiles,'Merged_with_FPKM_exp.diff')
        trackingFiles = [os.path.join(item[0],item[2]) for item in coupleDiff]
        logPrintInfo('Now merging read_group_tracking files...\n')
        logStr = logStr + mergeFiles(pathx,trackingFiles,'Merged_with_FPKMs.read_group_tracking')
    elif(not coupleDiff):
        logPrintInfo('None of the .diff files has their read_group_tracking files\n')
        totalDiff = [os.path.join(item[0],item[1]) for item in lonelyDiff]
        logPrintInfo('Now merging all the .diff files...\n')
        logStr = logStr + mergeFiles(pathx,totalDiff,'Merged_Total_exp.diff')
    else:
        logPrintInfo('Not all .diff files have their read_group_tracking files\n')
        totalDiff = ([os.path.join(item[0],item[1]) for item in coupleDiff] +
                            [os.path.join(item[0],item[1]) for item in lonelyDiff])
        logPrintInfo('Now merging all the .diff files...\n')
        logStr = logStr + mergeFiles(pathx,totalDiff,'Merged_Total_exp.diff')
        diffFiles = [os.path.join(item[0],item[1]) for item in coupleDiff]
        logPrintInfo('Now merging only the .diff files that have read_group_tracking files...\n')
        logStr = logStr + mergeFiles(pathx,diffFiles,'Merged_with_FPKM_exp.diff')
        trackingFiles = [os.path.join(item[0],item[2]) for item in coupleDiff]
        logPrintInfo('Now merging read_group_tracking files...\n')
        logStr = logStr + mergeFiles(pathx,trackingFiles,'Merged_with_FPKMs.read_group_tracking')
    logPrintInfo('finish\n')
    with open(os.path.join(pathx,'Merge_log.txt'),'w') as lf:
        lf.write(logStr)
