$freec2bed = "/home/jimmy200340/opt/FREEC-11.6/scripts/freec2bed.pl";
my @sample = ("Cry095", "SC04", "V19");
my $path = "/home/jimmy200340/work/WES/0329/CNV_FREEC";
chdir "$path";
mkdir "gistic";
mkdir "gistic/gistic_output";
system ("cp /home/jimmy200340/opt/FREEC-11.6/freec2gistic.r gistic/");
system ("cp /home/jimmy200340/opt/FREEC-11.6/prepare.sh gistic/");
system ("cp /home/jimmy200340/opt/FREEC-11.6/ref_hg19.genome gistic/");

foreach $sample(@sample)
{
	chdir "$path/$sample";
	system ("perl $freec2bed -f $path/$sample/$sample.mark.sort.recal.bam_ratio.txt > $path/$sample/$sample.freec_segments.bed");
	system ("mv $path/$sample/$sample.freec_segments.bed $path/gistic");
	system ("cp $path/$sample/$sample.mark.sort.recal.bam_sample.cpn $path/gistic");
}

system ("cd $path/gistic");
$dir = "$path/gistic";
system ("Rscript $path/gistic/freec2gistic.r --work-dir \"$dir\"");
system ("chmod 775 $path/gistic/prepare.sh");
system ("cd $path/gistic && sh ./prepare.sh");
mkdir "/home/jimmy200340/opt/GISTIC_2_0_23/output";
system ("cp $path/gistic/samples.seg /home/jimmy200340/opt/GISTIC_2_0_23/output");
chdir "/home/jimmy200340/opt/GISTIC_2_0_23";
system ("./gistic2 -b output -seg output/samples.seg -refgene refgenefiles/hg19.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme && cp -r ./output/* $path/gistic/gistic_output");