# -*- coding: utf-8 -*-
'''
Created on 2012-7-23

@author: CHENGCHAO
'''

### import necessary libraries
#import re
import os
import sys
import logging
import time
import datetime
import commands

### checking out the python version
if sys.version_info < (2, 6) or sys.version_info >= (3, 0):
    print ("Python Version error: must use phthon 2.6 or greater (not supporting python 3 yet)")
    sys.exit(-1)
#else:
#    print "Passed version test"
#    sys.exit(0)


############################ parameter variables
### ___ required values
gff = ''     # gff file
samFile = ''     # sam file
bamFile = ''     # bam file
allsj = ''    # junction file(.bed)
faFile = ''      # fasta file
outDir = ''    # output directory


for paramIndex in range(1, len(sys.argv)):  # going through all parameters
    if (sys.argv[paramIndex] == '-gff'):  # gff file
        paramIndex += 1  # increase index
        gff = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-allsj'):
        paramIndex += 1  # increase index
        allsj = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-sam'):
        paramIndex += 1  # increase index
        samFile = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-bam'):
        paramIndex += 1  # increase index
        bamFile = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-fa'):
        paramIndex += 1  # increase index
        faFile = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-outDir'):  # outDir
        paramIndex += 1  # increase index
        outDir = sys.argv[paramIndex]


### checking out the required arguments
if (outDir == '' or gff == '' or samFile == ''or bamFile == '' or faFile == '' or allsj == ''):  # at least one required param is missing
    print ('Not enough arguments!!')
    print ('Usage (with fastq files):\n\tpython pipeline.py -gff gffFile -allsj sjFile -sam samFile -bam bamFile -fa faFile -outDir outDir')
    sys.exit()


def listToString(x):
    rVal = ''
    for a in x:
        rVal += a + ' '
    return rVal

os.system('mkdir ' + outDir)
oFile = open(outDir + '/commands.txt', 'a')  # file that will contain list of commands excuted here

### setting up the logging format
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outDir + '/log.runCASC.' + str(datetime.datetime.now()),
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv))
startTime = time.time()

scriptPath = os.path.abspath(os.path.dirname(__file__))  # absolute script path
binPath = scriptPath + '/bin'  # absolute bin path
outPath = os.path.abspath(outDir)  # absolute output path

tempPath = outPath + '/temp'
expressionPath = outPath + '/Expression'
sjPath = outPath + '/SJ'
asPath = outPath + '/AS'
os.system('mkdir ' + tempPath)
os.system('mkdir ' + expressionPath)
os.system('mkdir ' + sjPath)
os.system('mkdir ' + asPath)

########## functions here... ############


def doPrepare():  # do some preparation work
    '''
    do some preparation work
    '''
    logging.debug("make Annofile: 将gff文件的feature一列用数字代替后的文件")
    cmd = 'perl ' + binPath + '/make_gff_at.pl -i ' + gff + ' -o ' + tempPath + '/annofile'
    oFile.write('###### make Annofile #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("makeAnnofile is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in makeAnnofile %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get chrLengthFile: 从fasta文件中提取chrLength")
    cmd = 'perl ' + binPath + '/getChrLength.pl ' + faFile + ' ' + tempPath + '/chrLengthFile'
    oFile.write('###### get chrLengthFile #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get chrLengthFile is done with status %s" % status)
    if (int(status) != 0): # it did not go well
        logging.debug("error in get chrLengthFile %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get known sj from gffFile: 从gff文件中提取已知的splice junction")
    cmd = 'perl ' + binPath + '/gff2sj.pl -gff ' + gff + ' -o ' + tempPath + '/gffsj'
    oFile.write('###### get gff sj #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get gff sj is done with status %s" % status)
    if (int(status) != 0): # it did not go well
        logging.debug("error in get gff sj %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

######## end of doPrepare() ##############


def doSj():  # count SJ expression,get total SJ, get novel SJ...
    logging.debug("get new sj: 获取所有新的sj")
    cmd = 'perl ' + binPath + '/get_known_new_sj.pl -gffsj ' + tempPath + '/gffsj' + \
        ' -sj ' + allsj + ' -on ' + sjPath + '/novel_sj' + ' -ok ' + sjPath + '/known_sj'
    cmd += " && more " + allsj + " | wc -l | more > " + tempPath + "/allsjNum"
    cmd += " && more " + sjPath + "/novel_sj" + " | wc -l | more > " + tempPath + "/novelsjNum"
    cmd += " && more " + sjPath + "/known_sj" + " | wc -l | more > " + tempPath + "/knownsjNum"
    oFile.write('###### get new sj #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get new sj is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get new sj %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)
    allsjNumFile = file(tempPath + '/allsjNum')
    allsjNum = int(allsjNumFile.readline())-1
    novelsjNumFile = file(tempPath + '/novelsjNum')
    novelsjNum = int(novelsjNumFile.readline())
    knownsjNumFile = file(tempPath + '/knownsjNum')
    knownsjNum = int(knownsjNumFile.readline())
    f = file(sjPath + '/sjNum', 'w')
    f.write("allsjNum\t" + str(allsjNum) + "\n")
    f.write("novelsjNum\t" + str(novelsjNum) + "\n")
    f.write("knownsjNum\t" + str(knownsjNum) + "\n")
    allsjNumFile.close()
    novelsjNumFile.close()
    knownsjNumFile.close()
    f.close()
######## end of doSj() ##############


def detectAS():  # detectAS without Intro...
    logging.debug("get_raw_novel_AS_events")
    cmd = 'perl ' + binPath + '/ABlas.pl -gff ' + gff
    cmd += ' -allsj ' + allsj + ' -sj ' + sjPath + '/novel_sj' + ' -o ' + asPath + '/novel_as'
    oFile.write('###### get_raw_AS_events #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_raw_AS_events is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_raw_AS_events %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("detectIntroR")
    # cmd = 'perl5.14 ' + binPath + '/run_ABlas_IntroR_unstrand_huangweibing.pl -gff ' + gff + ' -bam ' + bamFile + \
    # ' -chrlen ' + tempPath + '/chrLengthFile' + ' -o ' + asPath + '/raw_introR_events' + \
    # ' -anno ' + tempPath + '/annofile -od ' + asPath + ' -fa ' + faFile
    cmd = 'perl ' + binPath + '/ABlas_IntroR.pl -gff ' + gff + ' -sam ' + samFile + \
    ' -chrlen ' + tempPath + '/chrLengthFile' + ' -o ' + asPath + '/raw_introR_events' + \
    ' -anno ' + tempPath + '/annofile -od ' + expressionPath
    oFile.write('###### do EXP and detectIntroR #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("do EXP and detectIntroR is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in do EXP and detectIntroR %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get_uniq_AS")
    cmd = "cd " + asPath + " && sh " + binPath + "/novel_as_uniq.sh " + asPath
    cmd += " && sh " + binPath + "/novel_as_sj.sh " + binPath + "/ABlas_quantify.pl " + allsj
    oFile.write('###### get_uniq_AS #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_uniq_AS is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_uniq_AS %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get_available_AS")
    cmd = "cd " + asPath + ' && perl5.14 ' + binPath + '/intronR_quantify.pl -gff novel_as_IR_uniq_sj -fa ' + faFile
    cmd += ' -bam ' + bamFile + ' -o ' + 'novel_as_IR_uniq_sj2' + ' && perl5.14 ' + binPath + '/intronR_quantify.pl -gff raw_introR_events_uniq_sj_N -fa ' + faFile
    cmd += ' -bam ' + bamFile + ' -o ' + 'raw_introR_events_uniq_sj2'
    cmd += " && sh " + binPath + "/filter.sh"
    oFile.write('###### get_available_AS #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_available_AS is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_available_AS %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    # logging.debug("stat all")
    # cmd = "cd " + asPath + "&&cat availabe_as_events_novel availabe_introR_events_novel > AS_events_novel&&cat availabe_as_events_known availabe_introR_events_known > AS_events_known&&cat availabe_as_events availabe_introR_events > AS_events&&less -S AS_events_known | awk '{print $5}' | awk -F '|' '{print $1}'|sort|uniq -c > AS_events_known_result&&less -S AS_events_novel | awk '{print $5}' | awk -F '|' '{print $1}'|sort|uniq -c > AS_events_novel_result&&less -S AS_events | awk '{print $5}' | awk -F '|' '{print $1}'|sort|uniq -c > AS_events_result"
    # oFile.write('###### stat all #####\n' + cmd + '\n#\n')
    # oFile.flush()
    # status, output = commands.getstatusoutput(cmd)
    # logging.debug("stat all is done with status %s" % status)
    # if (int(status) != 0):  # it did not go well
    #     logging.debug("error in stat all %s" % status)
    #     logging.debug("error detail: %s" % output)
    #     raise Exception()
    # logging.debug(output)
######## end of detectAS() ##############

######## end of functions ##############


################## actual process ##############

####
#### 1. do some preparation work
####
logging.debug("1. do some preparation work")
logging.debug("start preparation work..")
try:
    doPrepare()
    # pass
except:
    logging.debug("There is an exception in preparation")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    sys.exit(-1)
logging.debug("done preparation")


####
#### 2. doSJ
####
logging.debug("2. count SJ expression,get total SJ, get novel SJ...")
logging.debug("start counting..")
try:
    doSj()
    # pass
except:
    logging.debug("There is an exception in counting")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    sys.exit(-1)
logging.debug("done doSj..")


####
#### 3. detectAS_without_intro
####
logging.debug("3. detectAS")
logging.debug("start detectAS..")
try:
    detectAS()
   # pass
except:
    logging.debug("There is an exception in detectAS")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    sys.exit(-1)
logging.debug("done detectAS..")