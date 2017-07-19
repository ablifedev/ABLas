# -*- coding: utf-8 -*-
'''
detect novel as events
Created on 2012-7-23
@author: CHENGCHAO
'''

# import necessary libraries
#import re
import os
import sys
import logging
import time
import datetime
import commands

# checking out the python version
if sys.version_info < (2, 6) or sys.version_info >= (3, 0):
    print("Python Version error: must use phthon 2.6 or greater (not supporting python 3 yet)")
    sys.exit(-1)
# else:
#    print "Passed version test"
#    sys.exit(0)


# parameter variables
# ___ required values
gff = ''     # gff file
fa = ''     # fa file
mod = ''     # gene mod file
bamFile = ''     # sam file
allsj = ''    # junction file(.bed)
knownsj = ''    # junction file(.bed)
novelsj = ''    # junction file(.bed)
chrLengthFile = ''      # chr Length file
known_nir = ''     # known_nir
known_ir = ''     # known_ir
# readLength = ''     # readLength for iden_intronR
outDir = ''    # output directory
sp_type = 'animal'    # output directory


for paramIndex in range(1, len(sys.argv)):  # going through all parameters
    if (sys.argv[paramIndex] == '-gff'):  # gff file
        paramIndex += 1  # increase index
        gff = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-fa'):  # outDir
        paramIndex += 1  # increase index
        fa = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-mod'):  # outDir
        paramIndex += 1  # increase index
        mod = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-allsj'):  # outDir
        paramIndex += 1  # increase index
        allsj = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-bamFile'):  # outDir
        paramIndex += 1  # increase index
        bamFile = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-chrlen'):  # outDir
        paramIndex += 1  # increase index
        chrLengthFile = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-known_ir'):  # outDir
        paramIndex += 1  # increase index
        known_ir = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-known_nir'):  # outDir
        paramIndex += 1  # increase index
        known_nir = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-outDir'):  # outDir
        paramIndex += 1  # increase index
        outDir = sys.argv[paramIndex]
    elif (sys.argv[paramIndex] == '-spType'):  # outDir
        paramIndex += 1  # increase index
        sp_type = sys.argv[paramIndex]


# checking out the required arguments
# at least one required param is missing
if (outDir == '' or gff == '' or bamFile == '' or chrLengthFile == ''):
    print('Not enough arguments!!')
    # print ('Usage (with fastq files):\n\tpython run_CASC.py -s1 reads1_1[,reads1_2] \
    # -s2 reads2_1[,reads2_2] -gtf gtfFile -bi bowtieIndexBase -o outDir -t readType -len readLength [options]*')
    print ('Example\n\tpython run_CASC.py -gff ~/data/assource/hg18_knownGene_Final_sort_cuttran_addintron.gff -mod ~/data/assource/hg_gene_model_genesym -allsj \
~/data/CUGBP1/Splicemap/barcode6_CUG_BP1_0.fq/junction_nUM_color.bed -bamFile ~/data/CUGBP1/Splicemap/barcode6_CUG_BP1_0.fq/good_hits.sam -chrlen \
~/data/assource/chr_len -outDir ~/data/testout&\n')
    sys.exit()


def listToString(x):
    rVal = ''
    for a in x:
        rVal += a + ' '
    return rVal

os.system('mkdir -p ' + outDir)
# file that will contain list of commands excuted here
oFile = open(outDir + '/commands.txt', 'a')

# setting up the logging format
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outDir + '/log.runCASC.' +
                    str(datetime.datetime.now()),
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv))
startTime = time.time()

scriptPath = os.path.abspath(os.path.dirname(__file__))  # absolute script path
binPath = scriptPath + '/bin/'  # absolute bin path
outPath = os.path.abspath(outDir)  # absolute output path

tempPath = outPath + '/temp'
expressionPath = outPath + '/Expression'
sjPath = outPath + '/SJ'
asPath = outPath + '/AS'
os.system('mkdir -p ' + tempPath)
os.system('mkdir -p ' + expressionPath)
os.system('mkdir -p ' + sjPath)
os.system('mkdir -p ' + asPath)

########## functions here... ############


def doPre():  # do preparing work
    '''
    do preparing work
    '''
    # perl ~/work/AS/CASC_cotton/bin/Anno/addintro.pl -g
    # Graimondii_221_gene_exons.gff -o Graimondii_221_gene_exons_addintron.gff
    global gff
    logging.debug("add intron to gff: 在gff中增加intron信息")
    cmd = 'perl ' + binPath + '/Anno/addintro.pl -g ' + \
        gff + ' -o ' + tempPath + '/gff.addintron'
    oFile.write('###### add intron to gff #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("add intron to gff is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in add intron to gff %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)
    gff = tempPath + '/gff.addintron'  # 后面使用gff都是加了intron的gff

    # perl ~/work/AS/CASC_cotton/bin/Anno/make_gff_at.pl -i
    # ../Graimondii_221_gene_exons_addintron.gff -o
    # Graimondii_221_gene_exons_addintron.gff.anno
    logging.debug("make Annofile: 将gff文件的feature一列用数字代替后的文件")
    cmd = 'perl ' + binPath + '/Anno/make_gff_at.pl -i ' + \
        gff + ' -o ' + tempPath + '/annofile'
    oFile.write('###### make Annofile #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("makeAnnofile is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in makeAnnofile %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

######## end of doDist() ##############


def doSj():  # count SJ expression,get total SJ, get novel SJ...
    global knownsj
    global novelsj
    # perl ~/work/AS/CASC_cotton/bin/Anno/genemodel2sj.pl -gff
    # Graimondii_221_gene_exons.gff -mod Graimondii_221_gene_model -o
    # model2sj.bed
    logging.debug("genemodel to bed:提取genemodel中的junctions")
    cmd = 'perl ' + binPath + '/Anno/genemodel2sj.pl -gff ' + \
        gff + ' -mod ' + mod + ' -o ' + tempPath + '/annoSj.bed'
    oFile.write('###### genemodel to bed #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("genemodel to bed is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in genemodel to bed %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    # perl /users/chengc/work/AS/CASC_cotton/bin/Anno/gff2bed.pl -i
    # Graimondii_221_gene_exons.gff -o gff2sj.bed
    logging.debug("gff to bed:提取gff中的junctions")
    cmd = 'perl ' + binPath + '/Anno/gff2bed.pl -i ' + \
        gff + ' -o ' + tempPath + '/gffSj.bed'
    oFile.write('###### gff to bed #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("gff to bed is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in gff to bed %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    # perl ~/work/AS/CASC_cotton/bin/SJ/get_known_novel_sj.pl -asj
    # ../../../gff2sj.bed -allsj
    # /data4/project4/huangjiao/cotton/cotton_mRNA/cotton_mRNA_clean_0108/cotton-root-12h_all.fq_0108/junctions.bed
    # -o novel_junctions.bed -o2 known_junctions.bed
    logging.debug("get new sj: 获取所有新的sj")
    cmd = 'perl ' + binPath + '/SJ/get_known_novel_sj.pl -asj ' + tempPath + '/gffSj.bed' + \
        ' -allsj ' + allsj + ' -o ' + sjPath + '/new_sj' + ' -o2 ' + sjPath + '/known_sj'
    cmd += " ; cat " + allsj + " | sed '1d' | wc -l  > " + tempPath + "/allsjNum"
    cmd += " ; cat " + sjPath + "/new_sj" + \
        " | wc -l  > " + tempPath + "/novelsjNum"
    cmd += " ; cat " + sjPath + "/known_sj" + \
        " | wc -l  > " + tempPath + "/knownsjNum"
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
    allsjNum = int(allsjNumFile.readline())
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
    knownsj = sjPath + '/known_sj'
    novelsj = sjPath + '/new_sj'
######## end of doSj() ##############


def detectAS():  # detectAS without Intro...
    # 对已知的事件直接赋值，进行筛选：
    # perl ~/work/AS/CASC_cotton/bin/AS/nas/as_quantify.pl -as ../../../known_AS_NIR.txt -sj all_junctions.bed -o root_known_AS_NIR.txt.quantify &
    # perl ~/work/AS/CASC_cotton/bin/AS/nas/as_filter.pl -as
    # root_known_AS_NIR.txt.quantify -num_alt_min 2 -num_model_min 2
    # -num_sum_alt_model 10 -o root_known_AS_NIR.txt.quantify.filter
    logging.debug("known_as_quantify")
    cmd = 'perl ' + binPath + '/AS/nas/as_quantify.pl -as ' + known_nir
    cmd += ' -sj ' + allsj + ' -o ' + asPath + '/known_AS_NIR.txt.quantify ;'
    cmd += 'perl ' + binPath + '/AS/nas/as_filter.pl -as ' + asPath + \
        '/known_AS_NIR.txt.quantify -num_alt_min 2 -num_model_min 2 -num_sum_alt_model 10'
    cmd += ' -o ' + asPath + '/known_AS_NIR.txt.quantify.filter'
    oFile.write('###### known_as_quantify #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("known_as_quantify is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in known_as_quantify %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    # novel的事件用novel sj来找
    # perl ~/work/AS/CASC_cotton/bin/AS/nas/detectAS.pl -gff ../../../Graimondii_221_gene_exons.gff -mod ../../../Graimondii_221_gene_model -sj novel_junctions.bed -o root_novel_AS_NIR.txt
    # perl ~/work/AS/CASC_cotton/bin/AS/nas/as_quantify.pl -as root_novel_AS_NIR.txt -sj all_junctions.bed -o root_novel_AS_NIR.txt.quantify &
    # perl ~/work/AS/CASC_cotton/bin/AS/nas/as_filter.pl -as
    # root_novel_AS_NIR.txt.quantify -num_alt_min 2 -num_alt_ratio 0.15 -o
    # root_novel_AS_NIR.txt.quantify.filter
    logging.debug("get_novel_AS_events")
    cmd = 'perl ' + binPath + '/AS/nas/detectAS.pl -gff ' + gff + ' -mod '
    cmd += mod + ' -sj ' + novelsj + ' -o ' + asPath + '/novel_AS_NIR_all.txt ;'
    cmd += """less -S """ + asPath + """/novel_AS_NIR_all.txt | perl -ne 'chomp;@line=split(/\\t/);$tmp="$line[1];$line[2];$line[4];$line[5]";@pos=split(/;/,$tmp);@pos2=sort {$a<=>$b} @pos;print $_,"\\t",join(";",@pos2),"\\n";' > """ + asPath + """/novel_AS_NIR_all_pos.txt ;less -S """ + asPath + """/novel_AS_NIR_all_pos.txt | perl -ne 'chomp;@line=split(/\\t/);next if($nir{$line[-1]}==1);$nir{$line[-1]}=1;pop @line;print join("\\t",@line),"\\n";' > """ + asPath + """/novel_AS_NIR.txt;"""
    cmd += 'perl ' + binPath + '/AS/nas/as_quantify.pl -as ' + \
        asPath + '/novel_AS_NIR.txt'
    cmd += ' -sj ' + allsj + ' -o ' + asPath + '/novel_AS_NIR.txt.quantify ;'
    cmd += 'perl ' + binPath + '/AS/nas/as_filter.pl -as ' + asPath + \
        '/novel_AS_NIR.txt.quantify -num_alt_min 2 -num_alt_ratio 0.15'
    cmd += ' -o ' + asPath + '/novel_AS_NIR.txt.quantify.filter'
    oFile.write('###### get_novel_AS_events #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_novel_AS_events is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_novel_AS_events %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get_all_available_AS")
    cmd = "cat " + asPath + '/known_AS_NIR.txt.quantify.filter ' + asPath + \
        '/novel_AS_NIR.txt.quantify.filter > ' + asPath + '/available_AS_NIR_events ;'
    cmd += "cat " + asPath + '/available_AS_NIR_events |' + \
        "awk '{print $7}' |sort|uniq -c" + " > " + \
        asPath + "/available_AS_NIR_events_result ;"
    cmd += "cat " + asPath + '/known_AS_NIR.txt.quantify.filter |' + \
        "awk '{print $7}' |sort|uniq -c" + " > " + asPath + \
        "/available_known_AS_NIR_events_result ;"
    cmd += "cat " + asPath + '/novel_AS_NIR.txt.quantify.filter |' + \
        "awk '{print $7}' |sort|uniq -c" + " > " + \
        asPath + "/available_novel_AS_NIR_events_result "
    oFile.write('###### get_available_AS #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_available_AS is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_available_AS %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get_NIR_sj")
    cmd = 'perl ' + binPath + '/AS/nas/as_get_sj.pl -as ' + \
        asPath + '/available_AS_NIR_events'
    cmd += ' -sj ' + allsj + ' -o ' + asPath + '/available_AS_NIR_events_sj'
    oFile.write('###### get_NIR_sj #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_NIR_sj is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_NIR_sj %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

######## end of detectAS() ##############


def detectIntroR():  # detectIntroR..
    logging.debug("detectIntroR")
    cmd = 'perl5.16 ' + binPath + '/EXP/detect_exon_intron_exp_base.pl -bam ' + bamFile + ' -fa ' + fa + ' -chrlen ' + chrLengthFile + \
        ' -anno ' + tempPath + '/annofile -od ' + expressionPath + ' & \n'
    if sp_type == "plant":
        cmd += 'perl5.16 ' + binPath + '/AS/nas/detectIR_strand_plants.pl -gff ' + gff + ' -mod ' + mod + ' -fa ' + \
            fa + ' -bam ' + bamFile + ' -chrlen ' + chrLengthFile + \
            ' -o ' + asPath + '/available_AS_IR_events_all && '
        cmd += """less -S """ + asPath + """/available_AS_IR_events_all | perl -ne 'chomp;@line=split(/\\t/);$tmp="$line[0]:$line[1]:$line[2]:$line[3]";next if($ir{$tmp}==1);$ir{$tmp}=1;print $_,"\\n";' > """ + asPath + """/available_AS_IR_events & \n """
    elif sp_type == "animal":
        cmd += 'perl5.16 ' + binPath + '/AS/nas/detectIR_strand_animals.pl -gff ' + gff + ' -mod ' + mod + ' -fa ' + fa + ' -as ' + asPath + \
            '/available_AS_NIR_events_sj' + ' -bam ' + bamFile + ' -chrlen ' + \
            chrLengthFile + ' -o ' + asPath + '/available_AS_IR_events_all && '
        cmd += """less -S """ + asPath + """/available_AS_IR_events_all | perl -ne 'chomp;@line=split(/\\t/);$tmp="$line[0]:$line[1]:$line[2]:$line[3]";next if($ir{$tmp}==1);$ir{$tmp}=1;print $_,"\\n";' > """ + asPath + """/available_AS_IR_events & \n """
    cmd += 'wait \n'
    with open(tempPath + "/ir.sh", 'w') as o:
        o.writelines(cmd)
    oFile.write('###### do EXP and detectIntroR #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput('sh ' + tempPath + '/ir.sh')
    logging.debug("do EXP and detectIntroR is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in do EXP and detectIntroR %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get_known_novel_IntroR")
    cmd = 'perl ' + binPath + '/AS/analyze/get_known_novel_as.pl -i ' + \
        asPath + '/available_AS_IR_events '
    # cmd += '-gas /data2/genome/human/hg19/hg19_ASevents -ok ' + asPath + '/available_as_events_known -on ' + asPath + '/available_as_events_novel'
    cmd += '-gas ' + known_ir + ' -ok ' + asPath + \
        '/available_known_AS_IR_events -on ' + asPath + '/available_novel_AS_IR_events'
    oFile.write('###### get_known_novel_IntroR #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_known_novel_IntroR is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_known_novel_IntroR %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("stat_IntroR")
    cmd = "cat " + asPath + '/available_AS_IR_events |' + \
        "awk '{print $5}' |sort|uniq -c" + " > " + \
        asPath + "/available_AS_IR_events_result ;"
    cmd += "cat " + asPath + '/available_known_AS_IR_events |' + \
        "awk '{print $5}' |sort|uniq -c" + " > " + \
        asPath + "/available_known_AS_IR_events_result ;"
    cmd += "cat " + asPath + '/available_novel_AS_IR_events |' + \
        "awk '{print $5}' |sort|uniq -c" + " > " + \
        asPath + "/available_novel_AS_IR_events_result "
    oFile.write('###### stat_IntroR #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_available_IntroR is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_available_IntroR %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get_IR_sj")
    cmd = 'perl ' + binPath + '/AS/nas/as_get_sj_IR.pl -as ' + \
        asPath + '/available_AS_IR_events'
    cmd += ' -sj ' + allsj + ' -o ' + asPath + '/available_AS_IR_events_sj'
    oFile.write('###### get_IR_sj #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_IR_sj is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_IR_sj %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("get_AS_sj")
    cmd = "cd " + asPath + \
        " && cat available_AS_NIR_events_sj available_AS_IR_events_sj > AS_events_sj"
    oFile.write('###### get_AS_sj #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("get_AS_sj is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in get_AS_sj %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)

    logging.debug("stat all")
    cmd = "cd " + asPath + "&&cat available_AS_NIR_events_result available_AS_IR_events_result > AS_events_result&&cat available_known_AS_NIR_events_result available_known_AS_IR_events_result > AS_events_known_result&&cat available_novel_AS_NIR_events_result available_novel_AS_IR_events_result > AS_events_novel_result"
    oFile.write('###### stat all #####\n' + cmd + '\n#\n')
    oFile.flush()
    status, output = commands.getstatusoutput(cmd)
    logging.debug("stat all is done with status %s" % status)
    if (int(status) != 0):  # it did not go well
        logging.debug("error in stat all %s" % status)
        logging.debug("error detail: %s" % output)
        raise Exception()
    logging.debug(output)
######## end of detectIntroR ##############

######## end of functions ##############


################## actual process ##############

####
# 1. count Exon/Intro/mRNA/5UTR/3UTR/intergenic expression
####
logging.debug("1. preparing work")
logging.debug("start..")
try:
    doPre()
    # pass
except:
    logging.debug("There is an exception in counting")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    sys.exit(-1)
logging.debug("done doPre..")


####
# 2. doSJ
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
# 3. detectAS_without_intro
####
logging.debug("3. detectAS without Intro...")
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


####
# 4. detectIntroR
####
logging.debug("4. detectIntroR")
logging.debug("start detectIntroR..")
try:
    detectIntroR()
    # pass
except:
    logging.debug("There is an exception in detectIntroR")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    sys.exit(-1)
logging.debug("done detectIntroR..")
