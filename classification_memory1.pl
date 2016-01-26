#!/usr/bin/perl
print "Starting classify  $ARGV[0]\n";
if ( ! $ARGV[0] eq "")
{
chdir($ARGV[0]);
}


#get file metadata
%metadata=();
open (DAT, "settings.ini") or return "no settings.ini";
while (<DAT>){$line = $_; chomp($line); (@temp)=split(' = ',$line); $settings{@temp[0]}=@temp[1];} close (DAT);

###########################################################################
# USER DEFINED #
#$trainfile="forest_loss.pix"; #Training file name
#$outname="forest_loss"; #Classification output file name
#$export="N"; #Export to ERDAS IMG format: "Y" or "N"
$maxtrees=7;#$settings{treethreads}; #Number of bagged trees
$sampling=20; #Sampling rate for one tree
$resultfile="out.vrt";#$settings{resultfile};

$maskfile=$settings{maskfile};
$maskclass=$settings{maskclass};
###########################################################################
# SYSTEM #
$threads=$settings{threads};#10 #Number of parallel processes - Tiles
$treethreads=$settings{treethreads}; #Number of parallel processes - Trees
#$easi=settings{easi};#"C:\\PCI\\Geomatica2013\\exe\\easi.exe"; #Path to EASI
$cpp=$settings{cpp};#"C:\\MinGW\\bin\\x86_64-w64-mingw32-g++.exe"; #Path to C++ compiler
$memsize=18000000000; # $settings{memsize};   Stack memory allocation
$metrics=$settings{metricspath};#"C:\\krylov\\test\\tiles"; # settings{metrics} memsizetiles path
$proj=$settings{projpath};#tiles path;#"C:\\krylov\\test"; # settings{proj}tiles path
$treeJ="C:\\NextGIS_QGIS\\UMD\\tree.exe";
$PF=$settings{"ignore pf flags"};
@PF=split(",", $PF);
#print $PF;
#print @PF;
#die;
$subs=10;


###########################################################################
# REGION #
$region=$settings{region};#"samerica";
$ulxgrid=$settings{ulxgrid};#-4560000;
$ulygrid=$settings{ulygrid};#2400000;
$prolong=-60;
$tileside=$settings{tileside};#2000;
$tilebuffer=$settings{tilebuffer};#2;
$pixelsize=$settings{pixelsize};#30;

@tilelist=readpipe("dir $metrics /A:D /O:N /B"); foreach (@tilelist) {chomp;}
#&subset_tile_list;
###########################################################################
# Setup #
$|=1; use threads;
(@time)=localtime(); @time[5]=1900+@time[5]; @time[4]=1+@time[4]; $theTime = "@time[2]\:@time[1]\:@time[0] @time[4]/@time[3]/@time[5]";
open (LOG, ">log.txt"); print LOG "start time: $theTime\n";  print "start time: $theTime\n";
$treefolder="trees\_".$outname;
if (!-e "$cpp" or !-d "$metrics"){die ("ERROR: Check setup and data files\n");}
if (-e "$outname.out") {system("del $outname.out"); system("del $outname.i*");}
use Cwd; $workfolder=cwd; $workfolder =~ s/\//\\/g;
$cmetrics=$metrics; $cmetrics =~ s/\\/\\\\/g;
&get_metric_list;
#get extent
if (@tilelist==1){
$tile=@tilelist[0];

	$x1=substr $tile,0,3; $sx=substr $tile,3,1; if ($sx eq "W"){$centerx=($x1+0.5)*-1;} else {$centerx=$x1+0.5;}
	$y1=substr $tile,5,2; $sy=substr $tile,7,1; if ($sy eq "S"){$centery=($y1+0.5)*-1;} else {$centery=$y1+0.5;}
	$ulx=$centerx-0.5;
	$uly=$centery+0.5;
	$lrx=$centerx+0.5;
	$lry=$centery-0.5;
	$side=$tileside+2*$tilebuffer;


$hmin=$ulx; 
$vmin=$lry;
$hmax=$ulx; 
$vmax=$lry
#$hmin=substr $tile,0,3; 
#$vmin=substr $tile,4,3;
#$hmax=substr $tile,0,3; 
#$vmax=substr $tile,4,3;
}
else {
my @harray=(); my @varray=();
for $tile (@tilelist) {


	$x1=substr $tile,0,3; $sx=substr $tile,3,1; if ($sx eq "W"){$centerx=($x1+0.5)*-1;} else {$centerx=$x1+0.5;}
	$y1=substr $tile,5,2; $sy=substr $tile,7,1; if ($sy eq "S"){$centery=($y1+0.5)*-1;} else {$centery=$y1+0.5;}
	$ulx=$centerx-0.5;
	$uly=$centery+0.5;
	$lrx=$centerx+0.5;
	$lry=$centery-0.5;
	$side=$tileside+2*$tilebuffer;
    
	#$h=substr $tile,0,3;
	push(@harray, $ulx); #$v=substr $tile,4,3;
	push(@varray, $lry);

#$h=substr $tile,0,3; push(@harray, $h); $v=substr $tile,4,3; push(@varray, $v);

}
@hsort = sort { $a <=> $b } @harray; $hmin = shift(@hsort); $hmax = pop(@hsort);
@vsort = sort { $a <=> $b } @varray; $vmin = shift(@vsort); $vmax = pop(@vsort);
}
#$totalxsize=($hmax-$hmin+1)*2000; $totalysize=($vmax-$vmin+1)*2000;
$totalxsize=($hmax-$hmin+$tileside*$pixelsize)/$pixelsize; $totalysize=($vmax-$vmin+$tileside*$pixelsize)/$pixelsize;
$side=$tileside+2*$tilebuffer; $xsize=$side; $ysize=$side;
$smysize=$tileside/$subs+2*$tilebuffer;
$ystep=$tileside/$subs;
###########################################################################
# MAIN #
&clean;
&to_raster;
&export_training;
#die;
&trees;
#die;
&treeanalysis;
&class;
&add_vrt;
#&aggregate;
(@time)=localtime(); @time[5]=1900+@time[5]; @time[4]=1+@time[4]; $theTime = "@time[2]\:@time[1]\:@time[0] @time[4]/@time[3]/@time[5]";
print LOG "end time: $theTime\n";  print "end time: $theTime\n"; close (LOG);
print "\t-----------Process Complete-----------\n";
###########################################################################

############################################################################################
#Clean-up
sub clean {
readpipe("del sample*.txt 2>nul");
readpipe("del table*.txt 2>nul");
readpipe("del temp.txt 2>nul");
readpipe("del out.txt 2>nul");
readpipe("del *.train 2>nul");
readpipe("del *.sample* 2>nul");
readpipe("del script.eas 2>nul"); readpipe("del PRM.PRM 2>nul");
readpipe("del class.cpp 2>nul"); readpipe("del class.exe 2>nul");

readpipe ("del ???_???.bil");
if (-e "???_???.hdf") {readpipe ("del ???_???.hdf");}
if (-e "???_???.prj") {readpipe ("del ???_???.prj");}
if (-e "???_???.out") {system ("del ???_???.out");}
readpipe("del trainmask 2>nul");
}

sub to_raster {
#print "Query the number of tiles. \n";
for my $tile (@tilelist)
{
$h=substr $tile,0,3; # push(@harray, $h);
$v=substr $tile,4,3; #push(@varray, $v);
#print "$hmin, $vmin, $hmax, $vmax\n";
print "Rasterize background samples for $tile tile.\n";
$x1=substr $tile,0,3; $sx=substr $tile,3,1; if ($sx eq "W"){$centerx=($x1+0.5)*-1;} else {$centerx=$x1+0.5;}
$y1=substr $tile,5,2; $sy=substr $tile,7,1; if ($sy eq "S"){$centery=($y1+0.5)*-1;} else {$centery=$y1+0.5;}
$ulx=$centerx-0.5;

$uly=$centery+0.5;
$lrx=$centerx+0.5;
$lry=$centery-0.5;
#$ulx=$ulxgrid+$h*$tileside*$pixelsize;
#$uly=$ulygrid-$v*$tileside*$pixelsize;
#$lrx=$ulx+$tileside*$pixelsize;
#$lry=$uly-$tileside*$pixelsize;

print "gdal_rasterize -burn 1 -of EHdr -ot Byte -te $ulx $lry $lrx $uly -ts $tileside $tileside -l background $proj\\background.shp $proj\\$tile.bil  \n";
system("gdal_rasterize -burn 1 -of EHdr -ot Byte -te $ulx $lry $lrx $uly -ts $tileside $tileside -l background $proj\\background.shp $proj\\$tile.bil");
print "Rasterize target samples for $tile tile.\n";
print "gdal_rasterize -burn 2 -l target $proj\\target.shp $proj\\$tile.bil  \n";
system("gdal_rasterize -burn 2 -l target $proj\\target.shp $proj\\$tile.bil");
}
#die;
}

sub export_training {
#print @tilelist;
#C++ make samples
$pi_total=@trainiles; $pi_counter=1;
print "exporting training\n";
$procnumber=0; $threadinc=0; 
print "\tcompiling...\n";
open (OUT, ">class.cpp");
print OUT"#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/stat.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <string>
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <ogr_spatialref.h>


using namespace std;

int main(int argc, char* argv[])
{

GDALAllRegister();  //
GDALDataset  *INGDAL; //
GDALRasterBand  *INBAND; //


//arguments
if (argc != 2){cout << \"wrong argument\" <<endl; exit (1);}
string tilename=argv[1];
string filename;
//counters
int i, k, x, y;
int xloc, yloc;
//randomizer
int rvalue;
srand ( time(NULL) );
//size
int xsize=$side;
int ysize=$side;
int smysize=$smysize;
int bigloop, l1, l2;

long nonull;
int tempval, tempsum, avvalue;
ifstream TEMPIN;
ofstream Y[$maxtrees];
//int avemain[$num16u];
//int avediff[$numdiff];
//int avemain5[$num16u];
//int avediff5[$numdiff];
uint8_t tmp;
//source data
uint8_t training[ysize][xsize]; memset(training, 0, sizeof(training[0][0]) * xsize * ysize);
uint8_t smtraining[smysize][xsize]; memset(training, 0, sizeof(training[0][0]) * xsize * smysize);
uint8_t mask_t[ysize][xsize];
unsigned short pf[smysize][xsize];
";
foreach (@mainlist){print OUT "uint8_t $_\[smysize][xsize];\n";}
foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "unsigned short $out\[smysize][xsize];\n";}
foreach (@list16s) {print OUT "short $_\[smysize][xsize];\n";}
foreach (@list16u) {print OUT "unsigned short $_\[smysize][xsize];\n";}
foreach (@list32r) {print OUT "float $_\[smysize][xsize];\n";}



for($inc=1;$inc<=$maxtrees;$inc++){$kinc=$inc-1;
print OUT"filename=tilename+\"\.sample$inc\";
Y[$kinc].open(filename.c_str()); if(!Y[$kinc]) {cout << \"error making sample$inc\.txt file\" <<endl; exit (1);}
";}

print OUT"
filename=tilename+\"\.bil\";
TEMPIN.open(filename.c_str(), ios::binary); if(!TEMPIN) {cout << \"error reading bil $tile\" <<endl; exit (1);}
for(y=$tilebuffer; y<ysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++) {training[y][x]=(int) TEMPIN.get();}} TEMPIN.close();
";





if (! $maskfile eq ""){
$maskfile =~ s/\\|\//\\\\/g;
$maskfile=$maskfile . "\\\\";
print OUT"
filename=\"$maskfile\"+tilename+\"\.out\";
TEMPIN.open(filename.c_str(), ios::binary); if(!TEMPIN) {cout << \"error reading mask $tile\" <<endl; exit (1);}
for(y=$tilebuffer; y<ysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++) {mask_t[y][x]=(int) TEMPIN.get();}} TEMPIN.close();
for(y=$tilebuffer; y<ysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++) { if (mask_t[y][x]!=$maskclass)  {training[y][x]=0;}}}

";
}
print OUT" for(bigloop=0; bigloop<$subs; bigloop++) {
l1=bigloop*$ystep; l2=bigloop*$ystep+smysize;
for(y=l1; y<l2; y++) {for(x=0; x<xsize; x++){ 
smtraining[y-l1][x]=training[y][x];
}}

";


print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\pf.tif\";

INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, pf, xsize, smysize, GDT_UInt16, 0, 0); GDALClose(INGDAL);

for(y=$tilebuffer; y<smysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++){ 
"; 
foreach (@PF){
print OUT"if (pf[x][y]==$_) {smtraining[y][x]=0;}";
}
print OUT"
}}
";

foreach (@mainlist){
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\$_\.tif\";


INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, $_, xsize, smysize, GDT_Byte, 0, 0); GDALClose(INGDAL);



";}
foreach (@list16s){
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\$_\.tif\"; 


INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, $_, xsize, smysize, GDT_Int16, 0, 0); GDALClose(INGDAL);
";}


foreach (@list16u){
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\$_\.tif\";

INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, $_, xsize, smysize, GDT_UInt16, 0, 0); GDALClose(INGDAL);
";}


foreach (@list32r){
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\$_\.tif\";

INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, $_, xsize, smysize, GDT_Float32, 0, 0); GDALClose(INGDAL);
";}






foreach (@diff){($out,$in1,$in2)=split('\,',$_);
print OUT"for(y=0; y<smysize; y++) {for(x=0; x<xsize; x++) {
tempval=$in1\[y][x]-$in2\[y][x]+30000; if (tempval<0){tempval=0;} if (tempval>65535){tempval=65535;} $out\[y][x]=tempval;}}";}

print OUT"for(y=0; y<smysize; y++) {for(x=0; x<xsize; x++) {
if (smtraining[y][x]>0) {
nonull++;
for(k=0; k<$maxtrees; k++){
rvalue=rand() % 100 + 1;
if (rvalue<=$sampling){
";
=for_comment
$inc=0; foreach (@list16u){print OUT "
tempval=0; tempsum=0;
for (yloc=-1; yloc<=1; yloc++) {for(xloc=-1; xloc<=1; xloc++) {
if ((x+xloc)>0 and (y+yloc)>0 and (x+xloc)<xsize and (y+yloc)<ysize){++tempval; tempsum=tempsum+$_\[x+xloc][y+yloc];}}}
avemain[$inc]=static_cast<int>(((double)tempsum/tempval)+0.5);
"; ++$inc;}
#1
$inc=0; foreach (@list16u){print OUT "
tempval=0; tempsum=0;
for (yloc=-2; yloc<=2; yloc++) {for(xloc=-2; xloc<=2; xloc++) {
if ((x+xloc)>0 and (y+yloc)>0 and (x+xloc)<xsize and (y+yloc)<ysize){++tempval; tempsum=tempsum+$_\[x+xloc][y+yloc];}}}
avemain[$inc]=static_cast<int>(((double)tempsum/tempval)+0.5);
"; ++$inc;}
#1
=cut

=for_comment

$inc=0; foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "
tempval=0; tempsum=0;
for (yloc=-1; yloc<=1; yloc++) {for(xloc=-1; xloc<=1; xloc++) {
if ((x+xloc)>0 and (y+yloc)>0 and (x+xloc)<xsize and (y+yloc)<ysize){++tempval; tempsum=tempsum+$out\[x+xloc][y+yloc];}}}
avediff[$inc]=static_cast<int>(((double)tempsum/tempval)+0.5);
"; ++$inc;}
#2


$inc=0; foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "
tempval=0; tempsum=0;
for (yloc=-2; yloc<=2; yloc++) {for(xloc=-2; xloc<=2; xloc++) {
if ((x+xloc)>0 and (y+yloc)>0 and (x+xloc)<xsize and (y+yloc)<ysize){++tempval; tempsum=tempsum+$out\[x+xloc][y+yloc];}}}
avediff[$inc]=static_cast<int>(((double)tempsum/tempval)+0.5);
"; ++$inc;}
#2
=cut
print OUT"Y[k]<<(int)smtraining[y][x]";
$inc=0; foreach (@mainlist){print OUT "<<\",\"<<(int)$_\[y][x]";}
$inc=0; foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "<<\",\"<<(int)$out\[y][x]";#<<\",\"<<avediff[$inc]<<\",\"<<avediff5[$inc]"; 
++$inc;}
foreach (@list16s) {print OUT "<<\",\"<<$_\[y][x]";}
$inc=0; foreach (@list16u) {print OUT "<<\",\"<<$_\[y][x]";#<<\",\"<<avemain[$inc]<<\",\"<<avemain5[$inc]"; 
++$inc; }
foreach (@list32r) {print OUT "<<\",\"<<$_\[y][x]";}
print OUT"
<<endl;
}}}}}

//bigloop
}

";
print OUT"
for(k=0; k<$maxtrees\; k++){Y[k].close();}

//cout << endl ;

//cout << nonull ;

//cout << endl ;
return 0;
}
";
close (OUT);
system("$cpp -I\"C:/NextGIS_QGIS/GDAL/include\" -L\"C:/NextGIS_QGIS/GDAL/lib\" class.cpp -o class.exe -Wl,--stack,$memsize -static -lgdal");
if (!-e "class.exe") {die "ERROR: C++ make samples (compiler failure)";}
#die;
print "\tprocessing...\n";
@trainiles=@tilelist;
foreach $tile (@trainiles) {
#$pi_prc=int($pi_counter/$pi_total*100); print "  $pi_prc";
push @ClassThreads, threads->create(\&export_thread, $tile); ++$procnumber;
if ($procnumber==$threads or $threadinc==$#trainiles){foreach $thread (@ClassThreads)  {$thread->join();} @ClassThreads=(); $procnumber=0;}
++$threadinc; ++$pi_counter;} 
sub export_thread {system("class.exe $tile");}
#die;
system("del class.cpp"); system("del class.exe");
print "\n";

#aggregate to samples
print "aggregating samples\n";
for($inc=1;$inc<=$maxtrees;$inc++){
$pi_prc=int($inc/$maxtrees*100); print "$pi_prc ...";
@samples=();
for $tile (@trainiles) {
if (-e "$tile.sample$inc"){ open (DAT1, "$tile.sample$inc");
while (<DAT1>){push (@samples,$_);}
close (DAT1);
system("del $tile.sample$inc");}
}
open (OUT, ">sample$inc\.txt");
foreach (@samples) {print OUT "$_";}
close (OUT);
}
for($inc=1;$inc<=$maxtrees;$inc++){$size = (-s "sample$inc\.txt"); if ($size==0){die "ERROR: sample aggregation failure";}}
#system("del *.train");
#system("del trainmask");
@samples=();
print " done";
print "\n";
}

############################################################################################
# Build tree
sub trees {
print "building trees\n";
if (-d "$treefolder") {system "del $treefolder\\* /q";} else {mkdir "$treefolder";}

$procnumber=0; $inc=1; 
while ($inc<=$maxtrees) {
$pi_prc=int($inc/$maxtrees*100); print "$pi_prc ...";
push @ClassThreads, threads->create(\&TR, $inc);
++$procnumber;
if ($procnumber==$treethreads or $inc==$maxtrees){
foreach $thread (@ClassThreads)  { $thread->join(); }
@ClassThreads = (); $procnumber=0; }
++$inc; }
sub TR {
open (PARAM, ">param$inc.txt");
print PARAM "datafile = sample$inc\.txt\n";
print PARAM "header = false\n";
print PARAM "catcols = 1\n";
print PARAM "mincut = 1\n";
print PARAM "minsize = 2\n";
print PARAM "mindev = 0.001\n";
print PARAM "split = deviance\n";
close (PARAM);
@mylog=readpipe ("$treeJ param$inc.txt > $treefolder/tree$inc\.txt");
}
$checktrees=1;
for($inc=1;$inc<=$maxtrees;$inc++){if (!-e "$treefolder/tree$inc\.txt"){$checktrees=0;}}
if ($checktrees==0){
die ("ERROR: could not process trees\n");
}
system("del sample*");
system("del param*");
print " done";
print "\n";
}

############################################################################################
# Classification
sub class {
print "classification...\n";
print "\tcompiling...\n";
open (OUT, ">class.cpp");
print OUT"#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/stat.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <stdint.h>
#include <string>
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <ogr_spatialref.h>
using namespace std;

int main(int argc, char* argv[])
{

GDALAllRegister();  //
GDALDataset  *INGDAL; //
GDALRasterBand  *INBAND; //


//arguments
if (argc != 2){cout << \"wrong argument\" <<endl; exit (1);}
string tilename=argv[1];
string filename;
//counters
int x, y, xloc, yloc;
int xsize=$side;
int ysize=$side;
int tempval, tempsum;
//result
int result[$maxtrees];
int medposition=int (floor($maxtrees/2));
double array[$metcount]; array[0]=0; array[1]=0;
ofstream outfile;
//source data
ifstream TEMPIN;
int smysize=$smysize;
int bigloop, l1, l2;

uint8_t tmp;
//int avemain[$num16u];
//int avemain5[$num16u];
//int avediff[$numdiff];
//int avediff5[$numdiff];
uint8_t mask_t[ysize][xsize];
uint8_t smmask_t[smysize][xsize];memset(smmask_t, 0, sizeof(smmask_t[0][0])* smysize * xsize );
uint8_t output[ysize][xsize]; memset(output, 0, sizeof(output[0][0]) * ysize * xsize);
uint8_t smoutput[smysize][xsize]; memset(smoutput, 0, sizeof(smoutput[0][0])* smysize * xsize );
unsigned short pf[ysize][xsize];
";
foreach (@mainlist){print OUT "uint8_t $_\[smysize][xsize];\n";}
foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "unsigned short $out\[smysize][xsize];\n";}
foreach (@list16s) {print OUT "short $_\[smysize][xsize];\n";}
foreach (@list16u) {print OUT "unsigned short $_\[smysize][xsize];\n";}
foreach (@list32r) {print OUT "float $_\[smysize][xsize];\n";}






for($inc=1;$inc<=$maxtrees;$inc++){print OUT"int tree$inc(double inparam[$metcount])\;\n";}
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\pf.tif\";

INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, 0, xsize, ysize, pf, xsize, ysize, GDT_UInt16, 0, 0); GDALClose(INGDAL);

";


if (! $maskfile eq ""){
#$maskfile =~ s/\\|\//\\\\/g;
#$maskfile=$maskfile . "\\\\";
print OUT"
filename=\"$maskfile\"+tilename+\"\.out\";
TEMPIN.open(filename.c_str(), ios::binary); if(!TEMPIN) {cout << \"error reading $tile\" <<endl; exit (1);}
for(y=$tilebuffer; y<ysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++) {mask_t[y][x]=(int) TEMPIN.get();}} TEMPIN.close();
for(y=$tilebuffer; y<ysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++) {if (mask_t[y][x]==$maskclass) {mask_t[y][x]=1;} else {mask_t[y][x]=0;}}}
";
}

if ($maskfile eq ""){
$maskclass=1;
print OUT"
for(y=$tilebuffer; y<ysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++) {mask_t[y][x]=1;}}
";
}


print OUT "
for(y=$tilebuffer; y<ysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++) {";
foreach (@PF){
print OUT"if (pf[y][x]==$_) {mask_t[y][x]=0;}";
}
print OUT"
}}
";




print OUT" for(bigloop=0; bigloop<$subs; bigloop++) {
l1=bigloop*$ystep; l2=bigloop*$ystep+smysize;
for(y=l1; y<l2; y++) {for(x=0; x<xsize; x++){ 
smmask_t[y-l1][x]=mask_t[y][x];
}}

";



foreach (@mainlist){
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\$_\.tif\";


INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, $_, xsize, smysize, GDT_Byte, 0, 0); GDALClose(INGDAL);



";}
foreach (@list16s){
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\$_\.tif\"; 

INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, $_, xsize, smysize, GDT_Int16, 0, 0); GDALClose(INGDAL);
";}


foreach (@list16u){
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\$_\.tif\";

INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, $_, xsize, smysize, GDT_UInt16, 0, 0); GDALClose(INGDAL);
";}


foreach (@list32r){
print OUT"
filename=\"$cmetrics\\\\\"+tilename+\"\\\\$_\.tif\";

INGDAL = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly ); INBAND = INGDAL->GetRasterBand(1);
INBAND->RasterIO(GF_Read, 0, l1, xsize, smysize, $_, xsize, smysize, GDT_Float32, 0, 0); GDALClose(INGDAL);
";}




foreach (@diff){($out,$in1,$in2)=split('\,',$_);
print OUT"for(y=0; y<smysize; y++) {for(x=0; x<xsize; x++) {
tempval=$in1\[y][x]-$in2\[y][x]+30000; if (tempval<0){tempval=0;} if (tempval>65535){tempval=65535;} $out\[y][x]=tempval;}}";}

print OUT" cout<< 1;";
print OUT"for(y=0; y<smysize; y++) {for(x=0; x<xsize; x++) {
if (smmask_t[y][x]==1) {
";
=for_comment
$inc=0; foreach (@list16u){print OUT "
tempval=0; tempsum=0;
for (yloc=-1; yloc<=1; yloc++) {for(xloc=-1; xloc<=1; xloc++) {
if ((x+xloc)>0 and (y+yloc)>0 and (x+xloc)<xsize and (y+yloc)<ysize){++tempval; tempsum=tempsum+$_\[x+xloc][y+yloc];}}}
avemain[$inc]=static_cast<int>(((double)tempsum/tempval)+0.5);
"; ++$inc;}
#1
$inc=0; foreach (@list16u){print OUT "
tempval=0; tempsum=0;
for (yloc=-2; yloc<=2; yloc++) {for(xloc=-2; xloc<=2; xloc++) {
if ((x+xloc)>0 and (y+yloc)>0 and (x+xloc)<xsize and (y+yloc)<ysize){++tempval; tempsum=tempsum+$_\[x+xloc][y+yloc];}}}
avemain5[$inc]=static_cast<int>(((double)tempsum/tempval)+0.5);
"; ++$inc;}
#1
$inc=0; foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "
tempval=0; tempsum=0;
for (yloc=-1; yloc<=1; yloc++) {for(xloc=-1; xloc<=1; xloc++) {
if ((x+xloc)>0 and (y+yloc)>0 and (x+xloc)<xsize and (y+yloc)<ysize){++tempval; tempsum=tempsum+$out\[x+xloc][y+yloc];}}}
avediff[$inc]=static_cast<int>(((double)tempsum/tempval)+0.5);
"; ++$inc;}
#2
$inc=0; foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "
tempval=0; tempsum=0;
for (yloc=-2; yloc<=2; yloc++) {for(xloc=-2; xloc<=2; xloc++) {
if ((x+xloc)>0 and (y+yloc)>0 and (x+xloc)<xsize and (y+yloc)<ysize){++tempval; tempsum=tempsum+$out\[x+xloc][y+yloc];}}}
avediff5[$inc]=static_cast<int>(((double)tempsum/tempval)+0.5);
"; ++$inc;}
#2
=cut

$metnumber=2;
$inc=0; foreach (@mainlist){print OUT "array[$metnumber] = $_\[y][x]\;\n"; ++$metnumber; }
$inc=0; foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "array[$metnumber] = $out\[y][x]\;\n"; ++$metnumber; 
#print OUT "array[$metnumber] = avediff[$inc]\;\n"; ++$inc; ++$metnumber; print OUT "array[$metnumber] = avediff5[$inc]\;\n"; ++$metnumber;
}
foreach (@list16s) {print OUT "array[$metnumber] = $_\[y][x]\;\n"; ++$metnumber;}

$inc=0; foreach (@list16u) {print OUT "array[$metnumber] = $_\[y][x]\;\n"; ++$metnumber; 
#print OUT "array[$metnumber] = avemain[$inc]\;\n"; ++$inc; ++$metnumber; print OUT "array[$metnumber] = avemain5[$inc]\;\n"; ++$metnumber;
}
foreach (@list32r) {print OUT "array[$metnumber] = $_\[y][x]\;\n"; ++$metnumber;}



for($inc=1;$inc<=$maxtrees;$inc++){print OUT"result[$inc-1] = tree$inc(array)\;\n";}
print OUT"
std::sort(result, result + $maxtrees);

smoutput[y][x]=result[medposition];
}
else {smoutput[y][x]=0;}
}}

for(y=l1; y<l2; y++) {for(x=0; x<xsize; x++){ 
output[y][x]=smoutput[y-l1][x];
//output[y][x]=av1025_b3[y-l1][x]/256;
}}

//big_look
}







filename=tilename+\"\.out\";
outfile.open(filename.c_str(), ios::binary); if(!outfile) {cout << \"error making output file\" <<endl; exit (1);} 
for(y=$tilebuffer; y<ysize-$tilebuffer; y++) {for(x=$tilebuffer; x<xsize-$tilebuffer; x++) {
tmp=output[y][x];
outfile.put(tmp);}}
outfile.close();
return 0;
}
";
for($inc=1;$inc<=$maxtrees;$inc++){
print OUT"int tree$inc(double inparam[$metcount])
{
";

open (DAT1, "$treefolder\\tree$inc\.txt");
@tree=(<DAT1>);
close (DAT1);
$flag=0;
foreach (@tree) {
$line = $_; 
$fl=(substr $line, 0, 1);
if ($fl eq 1) {$flag=1; next;}
if  ($flag eq 1) {
#print $fl;
chomp($line);
#print $line;
 @line1=split(/\)/, $line);
# print @line1[1];
 @line2=split(' ', @line1[1]);
 if (@line1[0]%2==0) {print OUT ("if")} else {
			until (@line1[0]>$pr_node) {print OUT "}"; $pr_node=($pr_node-1)/2;}
			print OUT " else if";}
$pr_node=@line1[0];
$tempname=@line2[0];
$tempname =~ s/X//g;
$tempname=$tempname*1+1;
$tempname="inparam\[".$tempname."\]";
 print OUT" ($tempname @line2[1] @line2[2]) {";
 if (@line1[2] eq " *") {
 $proc=substr @line2[6], 1, 10 ;
 #print $proc;
 # $proc=1+int((1-$proc*1)+0.5);
 $proc=1+int((100-$proc*100)+0.5);
 print OUT "return($proc);} \n";
 }
else {print OUT "\n";} 
}}
until (7>$pr_node) {print OUT "}"; $pr_node=($pr_node-1)/2;}
print OUT "return 0\;\}\n";
}
close (OUT);

system("$cpp -I\"C:/NextGIS_QGIS/GDAL/include\" -L\"C:/NextGIS_QGIS/GDAL/lib\" class.cpp -o class.exe -Wl,--stack,$memsize -static -lgdal");
if (!-e "class.exe") {die "ERROR: C++ make samples (compiler failure)";}

print "\tprocessing...\n";
$procnumber=0; $threadinc=0; 
$pi_total=@tilelist; $pi_counter=1;
foreach $tile (@tilelist) {
$pi_prc=int((($pi_counter - 0.5)/$pi_total)*100); print "$pi_prc ...";
push @ClassThreads, threads->create(\&class_thread, $tile); ++$procnumber;
if ($procnumber==$threads or $threadinc==$#tilelist){foreach $thread (@ClassThreads)  {$thread->join();} @ClassThreads=(); $procnumber=0;}
++$threadinc; ++$pi_counter;} 
sub class_thread {system("class.exe $tile");}
#system("del class.cpp"); system("del class.exe");
foreach $tile (@tilelist) {if (!-e "$tile.out"){die "ERROR: classification failure on $tile";}}
print "done";
print "\n";
}

############################################################################################

sub add_vrt {
print "Creating result file (out.vrt)...\n";
open (DAT1, "$metrics\\@tilelist[0]\\metric8u.vrt");
foreach (<DAT1>) {
#print $_;
if ((substr $_, 0, 5) eq "<SRS>"){
chomp($_);
$projectstring=$_;
last;
}}
close (DAT1);





open (OUT, ">$resultfile");
$xsize=($hmax-$hmin+$tileside*$pixelsize)/$pixelsize;
$ysize=($vmax-$vmin+$tileside*$pixelsize)/$pixelsize;
$ulx=$hmin;
$uly=$vmax+$tileside*$pixelsize;
print OUT "<VRTDataset rasterXSize=\"$xsize\" rasterYSize=\"$ysize\">
<GeoTransform>$ulx, $pixelsize, 0, $uly, 0, -$pixelsize</GeoTransform>
$projectstring\n";
#$i=1; 
print OUT "<VRTRasterBand dataType=\"Byte\" band=\"1\">\n";
foreach $tile (@tilelist) {

$x1=substr $tile,0,3; $sx=substr $tile,3,1; if ($sx eq "W"){$centerx=($x1+0.5)*-1;} else {$centerx=$x1+0.5;}
$y1=substr $tile,5,2; $sy=substr $tile,7,1; if ($sy eq "S"){$centery=($y1+0.5)*-1;} else {$centery=$y1+0.5;}

$ulx=$centerx-0.5;
$uly=$centery+0.5;
$lrx=$centerx+0.5;
$lry=$centery-0.5;


$x=($ulx-$hmin)/$pixelsize; $y=($lry-$vmin)/$pixelsize;
#$h=substr $tile,0,3; $v=substr $tile,4,3;
#$x=($h-$hmin)*$tileside; $y=($v-$vmin)*$tileside;

print OUT"
    <SimpleSource>
      <SourceFilename relativeToVRT=\"1\">$tile.out</SourceFilename>
      <SourceBand>1</SourceBand>
      <SourceProperties RasterXSize=\"$tileside\" RasterYSize=\"$tileside\" DataType=\"Byte\" BlockXSize=\"$tileside\" BlockYSize=\"1\" />
      <SrcRect xOff=\"0\" yOff=\"0\" xSize=\"$tileside\" ySize=\"$tileside\" />
      <DstRect xOff=\"$x\" yOff=\"$y\" xSize=\"$tileside\" ySize=\"$tileside\" />
    </SimpleSource>
";}
print OUT "</VRTRasterBand>\n";
#$i++;
print OUT "</VRTDataset>\n";
close (OUT);
print "Building overview images... \n";
system("gdaladdo $resultfile 16 32 64 128 256 512 1024 2048");
}


###########################################################################
# Tree analysis
sub treeanalysis {
print "\ttree analysis...\n";
open (DAT, "metric_list.txt");
@metr=<DAT>; foreach (@metr) {chomp;} $temp = shift(@metr);
$number=$maxtrees;
$folder=$treefolder;

# Make node tables
for($treex=1;$treex<=$number;$treex++) { #loop for trees START
open (INP, "$folder/tree$treex.txt");
open (TEMP, ">temp.txt");
my @allinp=<INP>; close (INP);
$num=@allinp; $inc=0; $stop=0;
while ($stop < 1) {
$line = @allinp[$inc];
if ($line eq "	 * denotes terminal node\n") {$stop=1;}
++$inc; if ($inc>$num){die;}}

$line = @allinp[$inc];
($node,$metric,$temp2,$deviance,$temp3)=split('\ +',$line); $node =~ s/\)//;
print TEMP "$node,root,$deviance,0\n";
++$inc;

while ($inc < $num) { #create table START
$line = @allinp[$inc];
($temp1,$node,$metric,$sign,$rule,$temp2,$deviance,$temp3)=split('\ +',$line);
$node =~ s/\)//;
$find=index($line,"*");
$terminal=0; if ($find > 0) {$terminal=1;}
if ($terminal==1 and $deviance<0){$deviance=$deviance*-1;}
print TEMP "$node,$metric,$deviance,$terminal\n";
++$inc; } #create table END
close (TEMP);
open (TMP1, "temp.txt");
open (TMP2, ">table$treex.txt");
my @allinp=<TMP1>;
$num=@allinp;
$inc=0;

while ($inc < $num) { #tree loop START
$line = @allinp[$inc]; $line =~ s/\n//;
($node,$metric,$dev,$terminal)=split(',',$line);
if ($terminal==1){print TMP2 "$node,$metric,0,0,0\n";}
else { #if node not terminal START
$nch1=$node*2; $nch2=$node*2+1; $ch1dev=0; $ch2dev=0;
$inc1=$inc;
while ($inc1 < $num) { #search child START
$line1 = @allinp[$inc1];
$line1 =~ s/\n//;
($node1,$metric1,$dev1)=split(',',$line1);
if ($node1==$nch1) {$ch1dev=$dev1; $chmetric=$metric1;}
if ($node1==$nch2) {$ch2dev=$dev1;}
++$inc1; } #search child END
print TMP2 "$node,$chmetric,$dev,$ch1dev,$ch2dev\n";
} #if node not terminal END
++$inc; } #tree loop END
close (TMP2); } #loop for trees END
close (TMP1);

# Analyse node tables - create output table
open (TEMP, ">out.txt");
for($treex=1;$treex<=$number;$treex++) { #loop for tables START
open (INP, "table$treex.txt");
my @allinp=<INP>;
$line = @allinp[0];
$line =~ s/\n//;
($node,$metric1,$dev,$ch1dev,$ch2dev)=split(',',$line);
print TEMP "root\t $treex\t $dev\n";
close (INP);} #loop for tables END

$count=@metr; $incount=0; 
while ($incount < $count) { #loop for metrix START
$line=@metr[$incount];chomp($line);($name,$fullname)=split('\t',$line);
for($treex=1;$treex<=$number;$treex++) { #loop for tables START
$met=0;
open (INP, "table$treex.txt"); my @allinp=<INP>; $num=@allinp;
$inc=0;
while ($inc < $num) { #loop in table START
$line1 = @allinp[$inc];
$line1 =~ s/\n//;
($node,$metric1,$dev,$ch1dev,$ch2dev)=split(',',$line1);
$deviance=$dev-$ch1dev-$ch2dev;
if ($metric1 eq $name){$met=$met+$deviance;}
++$inc;} #loop in table END

print TEMP "$name\t $treex\t $met\n";
close (INP);} #loop for tables END
++$incount;} #loop for metrix END
close (TEMP);

# Analyse output table - create report
open (TEMP, "out.txt");
open (OUT, ">tree_report.txt");
my @data=<TEMP>;
print OUT "metric_ID\tmetric_name\tdeviance(decrease)\tpercent_decrease_of_root\n";
$count=@metr; $incount=0;
$inc1=0; $root=0;
while ($inc1 < $number) {
$line2 = @data[$inc1];
$line2 =~ s/\n//;
($name,$treex,$met)=split('\t',$line2);
$met =~ s/\ //;
$root=$root+$met;
++$inc1; }
$averoot=$root/$number;
print OUT "root\t\t$averoot\tNA\n";

while ($incount < $count) {
$instart=($incount+1)*$number;
$inend=$instart+$number;
$met1=0;
while ($instart < $inend) {
$line3 = @data[$instart];
$line3 =~ s/\n//;
($name,$treex,$met)=split('\t',$line3);
$met =~ s/\ //;
$met1=$met1+$met;
++$instart;}
$met1=$met1/$number;
$avemet=$met1/$averoot*100;
$line=@metr[$incount]; ($name1,$fullname)=split('\t',$line);
print OUT "$name\t$fullname\t$met1\t$avemet\n";
++$incount;}

close (OUT); close (TEMP); close (INP);
system("del table*.txt");
system("del temp.txt");
system("del out.txt");

#metric type analysis
open (DAT, "tree_report.txt");
@data=<DAT>;
close (DAT);
$header = shift(@data); $header = shift(@data); foreach (@data) {chomp;}

for $metric (@mainlist){
for $line (@data) {
($temp1,$name,$temp2,$value)=split('\t',$line);
$name =~ s/_3x3//;
if ($name eq $metric){
if (exists $hash{$metric}){$hash{$metric}=$hash{$metric}+$value;}
else {$hash{$metric}=$value;}
}}}

for $metric (@list16s){
for $line (@data) {
($temp1,$name,$temp2,$value)=split('\t',$line);
$name =~ s/_3x3//;
if ($name eq $metric){
if (exists $hash{$metric}){$hash{$metric}=$hash{$metric}+$value;}
else {$hash{$metric}=$value;}
}}}

for $metricline (@diff){
($metric,$in1,$in2)=split(',',$metricline);
for $line (@data) {
($temp1,$name,$temp2,$value)=split('\t',$line);
$name =~ s/_3x3//;
if ($name eq $metric){
$hash{$in1}=$hash{$in1}+$value/2;
$hash{$in2}=$hash{$in2}+$value/2;
}}}

open (OUT, ">metric_type_analysis.txt");
foreach $key (sort {$hash{$a} cmp $hash{$b}} keys %hash) {print OUT "$key\t$hash{$key}\n";}
close (OUT);
}

############################################################################################
# Tile list
sub subset_tile_list {
$list="
049_036
049_037
050_036
050_037
";
@tilelist=split('\n',$list); shift(@tilelist);
}

############################################################################################
# Metric list
sub get_metric_list {






$list="
";

if (-e "config/metric.txt") {
	open (DAT, "config/metric.txt") or return "no metrics list file";
	while (<DAT>){
		chomp($_);
		push (@mainlist,$_);
		}
	close (DAT);
	$nummain=@mainlist;
	}
 else {print ("Using default metrics list");
 @mainlist=split('\n',$list); shift(@mainlist); $nummain=@mainlist;}

 
 
$list="
diffav2550_b3,av2550_b3,av5075_b3
diffav2550_b3_NBR,av2550_b3_NBR,av5075_b3_NBR
diffav2550_b3_NDVI,av2550_b3_NDVI,av5075_b3_NDVI
diffav2550_b3_TERM,av2550_b3_TERM,av5075_b3_TERM
diffav2550_b4,av2550_b4,av5075_b4
diffav2550_b4_NBR,av2550_b4_NBR,av5075_b4_NBR
diffav2550_b4_NDVI,av2550_b4_NDVI,av5075_b4_NDVI
diffav2550_b4_TERM,av2550_b4_TERM,av5075_b4_TERM
diffav2550_b5,av2550_b5,av5075_b5
diffav2550_b5_NBR,av2550_b5_NBR,av5075_b5_NBR
diffav2550_b5_NDVI,av2550_b5_NDVI,av5075_b5_NDVI
diffav2550_b5_TERM,av2550_b5_TERM,av5075_b5_TERM
diffav2550_b7,av2550_b7,av5075_b7
diffav2550_b7_NBR,av2550_b7_NBR,av5075_b7_NBR
diffav2550_b7_NDVI,av2550_b7_NDVI,av5075_b7_NDVI
diffav2550_b7_TERM,av2550_b7_TERM,av5075_b7_TERM
diffav2550_NBR,av2550_NBR,av5075_NBR
diffav2550_NDVI,av2550_NDVI,av5075_NDVI
diffavmin10_b3,avmin10_b3,av90max_b3
diffavmin10_b3_NBR,avmin10_b3_NBR,av90max_b3_NBR
diffavmin10_b3_NDVI,avmin10_b3_NDVI,av90max_b3_NDVI
diffavmin10_b3_TERM,avmin10_b3_TERM,av90max_b3_TERM
diffavmin10_b4,avmin10_b4,av90max_b4
diffavmin10_b4_NBR,avmin10_b4_NBR,av90max_b4_NBR
diffavmin10_b4_NDVI,avmin10_b4_NDVI,av90max_b4_NDVI
diffavmin10_b4_TERM,avmin10_b4_TERM,av90max_b4_TERM
diffavmin10_b5,avmin10_b5,av90max_b5
diffavmin10_b5_NBR,avmin10_b5_NBR,av90max_b5_NBR
diffavmin10_b5_NDVI,avmin10_b5_NDVI,av90max_b5_NDVI
diffavmin10_b5_TERM,avmin10_b5_TERM,av90max_b5_TERM
diffavmin10_b7,avmin10_b7,av90max_b7
diffavmin10_b7_NBR,avmin10_b7_NBR,av90max_b7_NBR
diffavmin10_b7_NDVI,avmin10_b7_NDVI,av90max_b7_NDVI
diffavmin10_b7_TERM,avmin10_b7_TERM,av90max_b7_TERM
diffavmin10_NBR,avmin10_NBR,av90max_NBR
diffavmin10_NDVI,avmin10_NDVI,av90max_NDVI
difffirst_b3,first_b3,last_b3
difffirst_b4,first_b4,last_b4
difffirst_b5,first_b5,last_b5
difffirst_b7,first_b7,last_b7
difffirst_NBR,first_NBR,last_NBR
difffirst_NDVI,first_NDVI,last_NDVI
diffmeanfirst_b3,meanfirst_b3,meanlast_b3
diffmeanfirst_b4,meanfirst_b4,meanlast_b4
diffmeanfirst_b5,meanfirst_b5,meanlast_b5
diffmeanfirst_b7,meanfirst_b7,meanlast_b7
diffmeanfirst_NBR,meanfirst_NBR,meanlast_NBR
diffmeanfirst_NDVI,meanfirst_NDVI,meanlast_NDVI
diffmedfirst_b3,medfirst_b3,medlast_b3
diffmedfirst_b4,medfirst_b4,medlast_b4
diffmedfirst_b5,medfirst_b5,medlast_b5
diffmedfirst_b7,medfirst_b7,medlast_b7
diffmedfirst_NBR,medfirst_NBR,medlast_NBR
diffmedfirst_NDVI,medfirst_NDVI,medlast_NDVI
";

if (-e "config/calc.txt") {
	open (DAT, "config/calc.txt") or return "no calc metrics list file";
	while (<DAT>){
		chomp($_);
		push (@diff,$_);
		}
	close (DAT);
	$numdiff=@diff;
	}
 else {print ("Using default calc metrics list");
 @diff=split('\n',$list); shift(@diff);  $numdiff=@diff;}





$list="
dem
";

if (-e "config/metric16.txt") {
	open (DAT, "config/metric16.txt") or return "no 16-bit metrics list file";
	while (<DAT>){
		chomp($_);
		push (@list16s,$_);
		}
	close (DAT);
	$num16s=@list16s;
	}
 else {print ("Using default 16-bit metrics list");
@list16s=split('\n',$list); shift(@list16s); $num16s=@list16s;}

$list="
";
if (-e "config/metric16u.txt") {
	open (DAT, "config/metric16u.txt") or return "no 16-bit u metrics list file";
	while (<DAT>){
		chomp($_);
		push (@list16u,$_);
		}
	close (DAT);
	$num16u=@list16u;
	}
 else {print ("Using default 16-bit u metrics list");
@list16u=split('\n',$list); shift(@list16u); $num16u=@list16u;}
$list="
";

if (-e "config/metric32r.txt") {
	open (DAT, "config/metric32r.txt") or return "no 32-bit r metrics list file";
	while (<DAT>){
		chomp($_);
		push (@list32r,$_);
		}
	close (DAT);
	$num32r=@list32r;
	}
 else {print ("Using default 32-bit r metrics list");
@list32r=split('\n',$list); shift(@list32r); $num32r=@list32r;}





open (OUT, ">metric_list.txt");
print OUT"X1\ttraining\n"; $inc=2;
foreach (@mainlist){print OUT "X$inc\t$_\n"; ++$inc;}
foreach (@diff){($out,$in1,$in2)=split('\,',$_); print OUT "X$inc\t$out\n"; ++$inc;# print OUT "X$inc\t$out\_3x3\n"; ++$inc; print OUT "X$inc\t$out\_5x5\n"; ++$inc;
}
foreach (@list16s) {print OUT "X$inc\t$_\n"; ++$inc;}
foreach (@list16u) {print OUT "X$inc\t$_\n"; ++$inc; 
#print OUT "X$inc\t$_\_3x3\n"; ++$inc; print OUT "X$inc\t$_\_5x5\n"; ++$inc;
}
foreach (@list32r) {print OUT "X$inc\t$_\n"; ++$inc;}
close (OUT);
$metcount=$inc-1;
}
