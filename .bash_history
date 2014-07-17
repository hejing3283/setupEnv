ltt
less NRP2.KeyReg.o3334630
less NRP2.KeyReg.e3334630
ll
ll *.e*
less PDIK1L.KeyReg.e3334524
less ZCCHC24.KeyReg.e3334584
ll *.e*|wc
grep error *.e*
cd ../../temp-tgene_big.txt_3/log/
ll *.e*|wc
ll *.e*
less APPL1.KeyReg.e3336163
less APPL1.KeyReg.o3336163 
less TIA1.KeyReg.e3336267
less TIA1.KeyReg.o3336267 
..
ll
less APPL1_candidateRegs.txt
ll
less SLC9A6_candidateRegs.txt
less CCNT2_candidateRegs.txt
awk 'NR==4{print FILENAME, $0}' *txt
awk 'FNR==4{print FILENAME, $0}' *txt
awk 'FNR==4&&$2<=0.05{print FILENAME, $0}' *txt
cd ../temp-tgene_big.txt_2
awk 'FNR==4&&$2<=0.05{print FILENAME, $0}' *txt
cd ../temp-tgene_big.txt_1
awk 'FNR==4&&$2<=0.05{print FILENAME, $0}' *txt
awk 'FNR==5&&$2<=0.05{print FILENAME, $0}' *txt
cd --
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/temp-tgene_big.txt_3
awk 'FNR==5&&$2<=0.05{print FILENAME, $0}' *txt
..
ll temp-tgene_big.txt_*/*txt|wc -l
wc -l tgene_big.txt
vi ../RUN.sh 
cd temp-tgene_big.txt_1/log/
ll *e*
cat GATAD2B.KeyReg.o3330926
grep *END* |wc
grep *END* *o* |wc
grep "#---END----"  *o* |wc
ll *o* | wc
..
grep "#---END----" temp-tgene_*/log/*o*| head
grep "#---END----" temp-tgene_*/log/*o*| awk -F"/|:|." '{print $4}' |head
grep "#---END----" temp-tgene_*/log/*o*| head| awk -F"/|:|." '{print $4}' 
grep "#---END----" temp-tgene_*/log/*o*| head| awk -F"/|:|\." '{print $4}' 
grep "#---END----" temp-tgene_*/log/*o*| head| awk -F"/|:|" '{split($3,".",a); print a[1]}' 
grep "#---END----" temp-tgene_*/log/*o*| head| awk -F"/|:" '{split($3,".",a); print a[1]}' 
grep "#---END----" temp-tgene_*/log/*o*| head| awk -F"/|:" '{split($3,a,"."); print a[1]}' 
grep "#---END----" temp-tgene_*/log/*o*| awk -F"/|:" '{split($3,a,"."); print a[1]}' > jobs.done.05162014
wc -l jobs.done.05162014
wc -l ../tgene_*txt
cat  ../tgene_*txt > ../tgene_all.txt
grep -v -w -f jobs.done.05162014 ../tgene_all.txt 
grep -v -w -f jobs.done.05162014 ../tgene_all.txt |wc
grep -v -w -f jobs.done.05162014 ../tgene_all.txt > jobs.fail.05162014
grep -w -f jobs.fail.05162014 tgene_big.txt
grep -w -f jobs.fail.05162014 tgene_big.txt|wc
fg
vi ../RUN.sh 
vi /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg_v4.py
fg
vi /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-2_regKeyRegulators_v6.R
fg
ll
fg
../RUN.sh 
qstcnt
ll
qstcnt
qst
cd ~/DATA/
ll
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/dnaseBS/
ll
less wgEncodeRegDnaseClusteredV2.bed
wc -l wgEncodeRegDnaseClusteredV2.bed
qstcnt
qst
qstcnt
qst
qstcnt
qst
qstcnt
qst
cd $crnsc 
cd model
ll step3-4_greedyOptCorr.py 
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/
ll
grep expfile RUN.sh 
less RUN.sh 
grep TRAPPC8 /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix|lss
vi tgene_noexpdata.txt
cd runApr30/
ll
ltt
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/dnaseBS/
ll
head ucsc_human_hg19_regulation_uniformSNaselHS_MCF7.bed
less ucsc_human_hg19_regulation_uniformSNaselHS_MCF7.bed
cd ../networks/
ll
rdf /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
cd ../../other/
ll
qrsh -l mem=8g,time=8::
exit
screen -ls
exit
ssh login3
screen -ls
screen
tcga
ll
pwd
exit
ssh tcga
tcga
ll
ll *11A*
ll *11A*bam|wc
ll *11A*bam
exit
ssh tcga
tcga
ll
less input_gtBatch_v2.txt_part3
ll input_g*
ll input_gtBatch_v2.txt_part*
ll input_gtBatch_v2.txt_part?+
ll input_gtBatch_v2.txt_part?
ll input_gtBatch_v2.txt_part??
cd ~/DATA/
cd projMisc/varReg/
ll
cd data/other/
ll
pwd
ll
chmod g+r *
ll
..
chmod g+r other/
ll
..
ll
chmod g+r data/
ll
..
ll
chmod g+r varReg/
tcga
ll *bam
cd $fdata
ll
cd 02042014/
ll
cd wgs/
ll
cd done/
ll
ll input_wgsCall_part*
crn 
ll
cd result/
ll
cd som
ll
..
ll
cd data/
ll
cd wgs/
ll
..
ll
cd somaticMutation/
ll
..
ll
cd $fres
ll
cd 02022014/
ll
cd wgsVars/
ll
cd rawVars/
ll
ll *vcf
ll *csv.vcf
ll *csv.vcf|wc
ll *genome*csv.vcf|wc
ls -1 *genome*csv.vcf|cut -d'.' -f1|less
#ls -1 *genome*csv.vcf|cut -d'.' -f1 > ~/DATA/projMisc/varReg/data/
mkdir ~/DATA/projMisc/varReg/data/062014/wgsVar
ls -1 *genome*csv.vcf|cut -d'.' -f1 > ~/DATA/projMisc/varReg/data/062014/wgsVar/sample.list_called_Feb2014
head ~/DATA/projMisc/varReg/data/062014/wgsVar/sample.list_called_Feb2014
cd ~/DATA/projMisc/varReg/data/
ll
c 062014/
cd other/
ll
cut -d- -f1,2,3,4 ../062014/wgsVar/sample.list_called_Feb2014 
cut -d- -f1,2,3 ../062014/wgsVar/sample.list_called_Feb2014 > sample.list_called_Feb2014.pid 
grep -f sample.list_called_Feb2014.pid tcga_brca_all_primarySolidTumor_06072014.tsv
grep -f sample.list_called_Feb2014.pid tcga_brca_all_primarySolidTumor_06072014.tsv|grep "-10A" |less
grep -f sample.list_called_Feb2014.pid tcga_brca_all_primarySolidTumor_06072014.tsv|grep "10A" |less
ll
wc -l tcga_brca_all_blood_normal\)06072014.tsv 
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal |grep "10A" |less
ll
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |grep "10A" |less
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |grep "10A" |wc
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |grep "10A" |awk -F"\t" '{print $2}'
wc -l ~/SCRATCH/projFocus/ceRNA/data/wgs/brca_wgs_normalbam_summary_02232014.tsv 
wc -l ~/SCRATCH/projFocus/ceRNA/data/wgs/brca_wgs_bam_summary_02042014.tsv_NormalSample
ll
echo tcga
echo $tcga
grep tcga ~/.bashrc*
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*01*bam
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*10A*bam
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*10*bam
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*bam
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*11A*bam
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*11A*bam |awk -F"/" 
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*11A*bam |awk -F"/" '{print NF}' 
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*11A*bam |awk -F"/" '{print $NF}' 
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*11A*bam |awk -F"/" '{print $NF}' > sample.normal.list_downloaded_Feb2014.pid 
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*11A*bam |awk -F"/" '{gsub(".bam","",$NF);print $NF}' > sample.normal.list_downloaded_Feb2014.pid 
cat sample.normal.list_downloaded_Feb2014.pid 
grep sample.normal.list_downloaded_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv 
grep -f sample.normal.list_downloaded_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv 
grep TCGA-BH-A0DG-11A-43D-A12Q-09 tcga_brca_all_blood_normal_06072014.tsv
grep TCGA-BH-A0DG-11A tcga_brca_all_blood_normal_06072014.tsv
les tcga_brca_all_blood_normal_06072014.tsv
less tcga_brca_all_blood_normal_06072014.tsv
ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*10A*bam |awk -F"/" '{gsub(".bam","",$NF);print $NF}' > sample.normal.list_downloaded_Feb2014.pid 
less tcga_brca_all_blood_normal_06072014.tsv 
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv 
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |less
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |wc
wc -l sample.list_called_Feb2014.pid
colcount tcga_brca_all_blood_normal_06072014.tsv
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |awk '{print $2, $17}' |less
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |awk -F"\t" '{print $2, $17}' |less
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |awk -F"\t" '{print $2"\t"$17}' |less
#ls -1 /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*10A*bam |awk -F"/" '{gsub(".bam","",$NF);print $NF}' > sample.normal.list_downloaded_Feb2014.pid 
less /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part1
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |awk -F"\t" '{print $2"\t"$17}' |less
grep -f sample.list_called_Feb2014.pid tcga_brca_all_blood_normal_06072014.tsv |awk -F"\t" '{print $2"\t"$17}' >input_gtBatch_v2.tcga.brca.BN4calledSample.txt
ll
wc -l input_gtBatch_v2.tcga.brca.BN4calledSample.txt
catwhich splitByN 
#splitByN 
pwd
tcga
ll
ll */
ll
ll normal
md tcgaNB4called
cd tcgaNB4called/
splitByN /ifs/home/c2b2/ac_lab/jh3283/DATA/projMisc/varReg/data/other/input_gtBatch_v2.tcga.brca.BN4calledSample.txt 10 
ll
mv /ifs/home/c2b2/ac_lab/jh3283/DATA/projMisc/varReg/data/other/input_gtBatch_v2.tcga.brca.BN4calledSample.txt_* .
ll
df -u
df .
df G -B .
df -h .
..
ll .sh
ll *.sh
ll README 
cat README 
less ~/scripts/TCGA/gtBatch_v2.sh
#vi ~/scripts/TCGA/gtBatch_v2.sh
less ~/scripts/TCGA/gtBatch_v2.sh
ll /ifs/scratch/c2b2/TCGA/soft/GeneTorrent/mykey.pem
ll /ifs/scratch/c2b2/TCGA/soft/GeneTorrent/cghub.key 
vi tgene_noexpdata.txt
vi ~/scripts/TCGA/gtBatch_v2.sh
exit
amlx
ll
cd callVars/result_anno_filter_final/
ll
cd somFinal/
ll
cd finals/
ll
wc -l *Tu*
ll *Tu*|wc
echo "215/16"|bc
wc -l *Re*
screen 
ssh login2
screen -ls
screen -r 
ssh login2
screen 
ssh login2
screen -ls
screen -r
ssh login2
ssh login2
screen -ls
screen -r 1161
screen -ls
screen -r 28200.pts-26.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
ssh login2
screen -ls
screen -r 28200.pts-26.login2
ssh login2
screen -s
screen -ls
screen -r 28200.pts-26.login2
ssh login2
screen -r
screen -r 1161.pts-3.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
screen -ls
screen -r 28200.pts-26.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
ll /ifs/scratch/c2b2/ac_lab/rs3412/no1/brca
cat SNVB6-A0RE.e5818148
cd SNVB6-A0RE.e5818148
cd /ifs/scratch/c2b2/ac_lab/rs3412/no1/brca
cat SNVB6-A0RE.e5818148
qstat -u rs3412
qstat -u rs3412|ls
qstat -u rs3412|less
ll *e|wc
ll *e*|wc
ll
ll *.o* |wc -l
tcga
cd tcgaNB4called/
ll
df -h .
ssh login2
screen -ls
screen -r 1161
ssh login2
screen -ls
screen -r 28200.pts-26.login2
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/
ll
rdf pc_ceRNET.xlsx
chmod g+r pc_ceRNET.xlsx
;;
ll
rdf pc_ceRNET.xlsx
pwd
tcga
ll
cd tcgaNB4called/
ll
ll part_1/
ll part_2/
cd part_2/
rmempdir 
ll
..
ll
tcga
cd tcgaNB4called/
~/scripts/TCGA/gtBatch_v2.sh input_gtBatch_v2.tcga.brca.BN4calledSample.txt_1 
ll
wc -l logs.download
wc -l logs.downloaded 
mkdir input_gtBatch_v2.tcga.brca.BN4calledSample.txt_1 
mkdir part_1
mv *bam* part_1/
ll
rmempdir *
ll
rmempdir 
ll
rm *gto
ll
mkdir part_2
cd part_
cd part_2
#~/scripts/TCGA/gtBatch_v2.sh ../input_gtBatch_v2.tcga.brca.BN4calledSample.txt_2 
df -h .
~/scripts/TCGA/gtBatch_v2.sh ../input_gtBatch_v2.tcga.brca.BN4calledSample.txt_2 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~a~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~
~~~~
~~~~~
~~~~
ssh login2
screen -ls
screen -r 28200.pts-26.login2
labv
ssh login2
screen -ls
screen -r 1161.pts-3.login2
scree n-ls
screen -ls
screen -r 1161.pts-3.login2
screen -ls
screen -r 28200.pts-26.login2
ssh login2
ssh login2
screen -ls
screen -Dr 1161.pts-3.login2
ssh login2
screen -ls
 screen -r 28200.pts-26.login2
 screen -Dr 28200.pts-26.login2
ssh login2
screen -ls
screen -Dr 1161.pts-3.login2
screen -ls
screen -Dr 1161.pts-3.login2
screen -ls
screen -Dr 28200.pts-26.login2
screen -ls
ssh login3
ssh login2
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/
ll
pwd
head miBS_cupidPredict.hg19
pwd
ll
head ll
ll
rm ~\$JMCB_mirBinding.xlsx 
ll
less MC2014_human_mirBS.bed.txt
ll
less MC2014_human_mirBS.bed.csv
ll
mv MC2014_human_mirBS.bed.txt
rm MC2014_human_mirBS.bed.txt
rm MC2014_human_mirBS.bed.xls
mv JMCB_mirBinding_2.xls MC2014_mirBS.xls
ll
rm JMCB_mirBinding.xlsx
ll
vi README.txt
fg
vi README.txt 
ll
less miBS_cupidPredict.hg19
ll
less MC2014_human_mirBS.bed.csv
fg
vi README.txt 
vi /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/runs/runDataMiBS.sh
ll /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
ech $crnsc
echo $crnsc
fg
pwd
fg
ll
qrsh -l mem=8g,time=8::
screen -ls
screen -r 1161.pts-3.login2
screen -ls
screen -Dr 28200.pts-26.login2
ssh login2
screen -ls
screen -Dr 1161.pts-3.login2
ssh login2
screen -ls
screen -Dr 28200.pts-26.login2
tcga
cd tcgaNB4called/
ll
mkdir part_4
cd part_4
df -h .
~/scripts/TCGA/gtBatch_v2.sh ../input_gtBatch_v2.tcga.brca.BN4calledSample.txt_4
tcga
cd tcgaNB4called/
ll
ll *bam|wc
ll *bam|wc -l
df -h .
#~/scripts/TCGA/gtBatch_v2.sh input_gtBatch_v2.tcga.brca.BN4calledSample.txt_3 
ll
mkdir part_3
cd part_3
~/scripts/TCGA/gtBatch_v2.sh ../input_gtBatch_v2.tcga.brca.BN4calledSample.txt_3 
exit
ssh login2
screen -ls
screen -r 1161.pts-3.login2
ssh login2
screen -ls
screen -r 28200.pts-26.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
ssh login2
screen -ls
screen -r 1161.pts-3.login2
cd $fres
ll
cd 05012014/
ll
cd fucFilt/
ll
cd ../sigMut/
ll
cd runMay5/
ll
..
ll
cat RUN.sh 
ll /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/
less keyRegSummary_allRunMay_05212014_0.01
cd runMay5/
less keyRegSummary_allRunMay_05212014_0.01
ll
less keyRegSummary_allRunMay_05212014_0.01
..
ll
cd --
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/
ll
d tcgal2som/
ll
cd tcgal2som/
ll
grep genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-20-2014.matrix *sh
grep genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-20-2014.matrix ../som/*sh
ll test
cd ../../02022014/
ll
cd som/
ll
cd --
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/tcgal2som
grep genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-20-2014.matrix /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/runs/*sh
fg
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jun2014/
ll
chmod g+wr result_step2N3_numbersNtests_06132014.xlsx
chmod o+wr result_step2N3_numbersNtests_06132014.xlsx
ll
rm ~\$result_step2N3_numbersNtests_06132014.xlsx 
ll /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/
ll /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/
ll
ls -1
..
cd runMay5/
ll
wc -l optCorr.result_flex_max_1000.tsv.significant.summary.netfile.txt
wc -l optCorr.result_flex_max_1000.tsv.significant.summary.netfile_full.2col_May-21-2014.txt
less optCorr.result_flex_max_1000.tsv.significant.summary.netfile.txt
ll
rdf optCorr.result_flex_max_1000.tsv.significant.summary.netfile.txt
less keyRegSummary_allRunMay_05212014_0.01
rdf keyRegSummary_allRunMay_05212014_0.01
fg
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/
ll
cd runMay5/
ll
cat optCorr.result_flex_max_1000.tsv.significant.summary.netfile.txt
cd ../sigTest/
ll
less ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all
wc -l ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all
rdf ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all
cd ../../candiReg/
ll
./RUN.sh 
ll
fg
ll
cd ../sigMut/cancerGene/
ll
wc -l cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all
rdf cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all
fg
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/encode/
ll
l
ll
cd dnaseBS/
ll
cd T470
ll
cd T47D/
ll
cd ../HMEC/
ll
..
cd chipSeqMethy/
ll
ll */*
cd MCF7/
ll
pwd
ll
pwd
fg
qrsh -l mem=8g,time=8::
cd --
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/
ll
cd cg
ll
cd test/
ll
vi RUN.sh 
ll
rm fg
fg
rdf /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/brca_tcga_som.utr3p.matrix.sorted.uniq
head /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/tcgal2som/brca_tcga_som.utr3p.matrix.sorted.uniq
less /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/tcgal2som/brca_tcga_som.utr3p.matrix.sorted.uniq
fg
pwd
fg
chmod +x RUN.sh 
./RUN.sh 
fg
./RUN.sh 
fg
./RUN.sh 
ll
less header_mut_matrix.allsamples
fg
./RUN.sh 
ll
less cg_ceRNA_driver_greedy_1k.TFBindSite.06242014.mutmat
fg
./RUN.sh 
ll
rm .RUN.sh.swp
ll
vless cg_ceRNA_driver_greedy_1k.TFBindSite.06242014.mutmat
less cg_ceRNA_driver_greedy_1k.TFBindSite.06242014.mutmat
ll
qrsh -l mem=8g,time=8::
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/utrBS_to_genmoicsCoord/cupid/
ll
w
m
mem
free -m
qrsh -l mem=8g,time=8::
screen -ls
screen -r Dr 1161.pts-3.login2
screen -Dr 1161.pts-3.login2
ll
cd ~/DATA/database/UCSC/
ll
less refseq_human_hg19_allGene_promoter2kb_upstrem.fasta
qrsh -l mem=8g,time=8::
qrsh -l mem=8g,time=8::
ll
cd $crnsc 
ll
ll *search*
ll
ll *py
grep awk *py
grep slow *py
cp extractSNPmat.py grepf2f.py
vi grepf2f.py 
fg
vi grepf2f.py 
mv grepf2f.py ~/bin/grepf2f
grep promoter *py
grep promoter  processData/*py
ll extractSNPmat.py 
less extractSNPmat.py
grep tss *py
grep 2000 *py
grep 2000 */*py
less processData/getMutation4EachRegion.py
python processData/getMutation4EachRegion.py
fg
cd $crnsc 
grep promoter */*py
grep promoter2k runs/*sh
grep -C 5 promoter2k runs/runDataTcgasom.sh
ll $crnsc/processData/selectTargetRegionRow.py
rdf runs/runDataTcgasom.sh
cd model/
vi step4-5_sigMut_sigTest.r 
fg
cd /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/
ll
ll /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_lasso_Jun-13-2014.uniq.utr3p 
ll /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k
ll /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p 
ll /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p 
ll /ifs/scratch/c2b2/ac_lab/rs3412/no1/net/annovar/cohort/cohort.snv.res.somatic.nosegup.norep
less /ifs/scratch/c2b2/ac_lab/rs3412/no1/net/annovar/cohort/cohort.snv.res.somatic.nosegup.norep
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/
ll
cd step5-1_bsSite/
ll
ll cg_ceRNA_driver_greedy_*
mv  cg_ceRNA_driver_greedy_* cg/
ll
md bs_with_addedBS
cd bs_with_addedBS/
ll
less /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/brca_somlevel2_byPoint.matrix
rdf /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/brca_somlevel2_byPoint.matrix
ln -s /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix .
ll
less /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_greedy_Jun-13-2014.uniq
ln -s /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_greedy_Jun-13-2014.uniq .
ll
ln -s /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_lasso_Jun-13-2014.uniq .
ll
less brca_somlevel2_byPoint.matrix
ll
vi RUN>sh
ps
kill -kill 28873
ll
fg
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/
ll
head miBS_*
ll
rm miBS_cupidPredict.hg19
cd ../../../database/
cd UCSC/
ll
cd --
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/bs_with_addedBS
wc -l brca_somlevel2_byPoint.matrix 
less ceRNA_driver_lasso_mut.mat
less ceRNA_driver_lasso_mut.mat.uniq.bed 
wc -l ceRNA_driver_lasso_mut.mat.uniq.bed 
which bedtools
cd $fdata
ll
cd encode/
ll
cd dnaseBS/
ll
mkdir mcf7
mv mcf7 MCF7
ll *Mcf*
mv  *Mcf* MCF7/
ll
md HMEC 
md T47D
ll MCF7/*|wc
cd MCF7/
ll
cd ../../h3k27/
ll
mkdir MCF7
mv *Mcf* MCF7/
..
ll
mv h3k27/ chipSeqMethy/
cd chipSeqMethy/
ll
md HMEC
cd ../dnaseBS/
ll
cd MCF7/
ll
qrsh -l mem=8g,time=8::
tcga
cd tcgaNB4called/
ll
rdf part_1/
pwd
cd $fres 
ll
qrsh -l mem=8g,time=8::
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/
ll
cd cancerGene/
ll
wc -l cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all
cd ../runMay5/
ll
cd ../cancerGene/
ll
less cg_ceRNA_driver_greedy_Jun-13-2014.uniq
wc -l cg_ceRNA_driver_greedy_Jun-13-2014.uniq
cd ../../fucFilt/
ll
cd step5-1_bsSite/
ll
mv cg_rca_somlevel2_byPoint.matrix.optCorr_flex_max_1000.sig.netfilemutRegTFBS.txt cg_brca_somlevel2_byPoint.matrix.optCorr_flex_max_1000.sig.netfilemutRegTFBS.txt
wc -l cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014.genelist
wc -l *genelist
cat cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist
cat cg_ceRNA_driver_greedy_utr3.mirBindSite.mutgene.point_06172014.genelist
grep -wf cg_ceRNA_driver_greedy_utr3.mirBindSite.mutgene.point_06172014.genelist cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist 
less cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014.genelist
ll
less cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist.bed
cut -f4 cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist.bed|sort|uniq|wc
cut -f1,2 cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist.bed|sort|uniq|wc
less ../RUN.sh 
cd ../step5-2_dnaseSite/
ll
cd mutInDnas/
ll
wc -l cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv
cat cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv
cut -f1 cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv|sort|uniq |wc -l
wc -l cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv
ll
wc -l brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv
cut -1 brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv|sort|uniq|wc -l
cut -f1 brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv|sort|uniq|wc -l
ll
cat cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv
cd ../../step5-1_bsSite/
ll
ll *mut.point*
wc -l *mut.point*
head cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014
less cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014
cut -f1 cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014|sort|uniq|wc
cut -f1 cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014|sort|uniq
cut -f1,2 cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014|sort|uniq
cut -f2 cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014|sort|uniq
cut -f2 cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014|sort|uniq|wc
cut -f1,2 cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014|sort|uniq|sort -k2
ll
cut -f1,2 cg_ceRNA_driver_greedy_utr3.mirBindSite.point.mutgene_06172014|sort|uniq|sort -k2
cut -f1,2 cg_ceRNA_driver_greedy_utr3.mirBindSite.point.mutgene_06172014|sort|uniq|wc
less cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014
cut -f1 cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014|sort|uniq|wc -l
cut -f2 cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014|sort|uniq|wc -l
cut -f6 cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014|sort|uniq|wc -l
cut -f6 cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014|sort|uniq -c |sort -k 1nr
cut -f1,6 cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014|sort|uniq|cut -f 2|sort|uniq|sort -k 1nr
cut -f1,6 cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014|sort|uniq|cut -f 2|sort|uniq -c|sort -k 1nr
cut -f1,6 cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014|sort|uniq|cut -f 2|sort|uniq -c|sort -k 1nr
ll
grep ATL2 /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/brca_somaticMutation_all.annotation_Jun172014
grep -w ATL2 /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/brca_somaticMutation_all.annotation_Jun172014
grep cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014
grep 38527429 cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014
ll
cd ../step5-2_dnaseSite/
ll
cd mutInDnas/
ll
cat cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv 
cut -f1 WDR20
SETD1B
ARL8B
ATP1B1
BPTF
CBFB
ll
cut -f1 cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv|sort|uniq
cut -f1 cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv|sortxh|uniq
cut -f1 cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv|sortxh|uniq|wc -l 
cut -f1 cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv|sortxh|uniq -c |sort -k 1nr 
..
ll
vi RUN.sh 
fg
qrsh -l mem=8g,time=8::
qst 
qstats -u rs3412
qstat -u rs3412
qstat -u rs3412|less
ll /ifs/scratch/c2b2/ac_lab/rs3412/no1/brca
cd /ifs/scratch/c2b2/ac_lab/rs3412/no1/brca
ll
cd compareBlood/
ll
cd snv/
ll
less TCGA-AR-A0TX-01A-11D-A19H-09
ll
..
ll ../vcfs/
ll 
/ifs/scratch/c2b2/TCGA/data/lincrna/rnaseq/brca/fastq
ll 
/ifs/scratch/c2b2/TCGA/data/lincrna/rnaseq/brca/fastq
ll /ifs/scratch/c2b2/TCGA/data/lincrna/rnaseq/brca/fastq/*|wc
qrsh -l mem=8g,time=8::
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/utrBS_to_genmoicsCoord/cupid/
ll
..
ll
..
ll
cd /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/
ll
ltt
rm parseFasta.py
ll *bed*
ll *Bed*
ll *combine*
ll ../*combine*
mv ../combineBS_bed.py .
rdf combineBS_bed.py 
tcga
cd tcgaNB4called/
ll
ll part_3/*bam
ll part_3/*bam |wc
ls -1 part_3/*bam |wc
ls -1 part_3/*bam
ll part_3/*bam
ll part_4/*bam
ll
cd part_2/
cd ../part_3/
ll
cd --
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/utrBS_to_genmoicsCoord/cupid
ll
less CupidPred.scored_5.bed
ps
ll
rm CupidPred.scored_.bed.uniq
ll
qrsh -l mem=20g,time=2:: -n now
qrsh -l mem=20g,time=2:: -now n
ll
vi ../../RUN.sh 
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/
cd fucFilt/step5-1_bsSite/
ll
cd bs_with_addedBS/
ll
ls -1
wc -l ceRNA_driver_greedy_mut.mat ceRNA_driver_lasso_mut.mat
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS
cd utrBS_to_genmoicsCoord/cupid/
ll
wc -l CupidPred.scored_all.bed.uniq 
less CupidPred.scored_all.bed.uniq 
vi /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/mergeOverlapRecord_sorted.py
fg
..
ll
rdf miBS_clash_parclip_mc2014_cupid.uniq.bed
screen -ls
ps
exit
screen -ls
exit
screen -ls
ssh login2
ssh login5
screen 
screen -ls
screen 
ssh login2
tcga
cd tcgaNB4called/
ll
cd part_3
ll
rm fdb7b8a7-7c0e-4738-8766-d58163b3d29b*
rm -rf fdb7b8a7-7c0e-4738-8766-d58163b3d29b*
ll
cp ../input_gtBatch_v2.tcga.brca.BN4calledSample.txt_3 .
sed -i 1d input_gtBatch_v2.tcga.brca.BN4calledSample.txt_3 
head input_gtBatch_v2.tcga.brca.BN4calledSample.txt_3 
~/scripts/TCGA/gtBatch_v2.sh input_gtBatch_v2.tcga.brca.BN4calledSample.txt_3 
