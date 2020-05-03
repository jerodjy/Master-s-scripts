
#Ancestry specific population analysis

#This perl script takes RFMix pipeline's output such as the beagle, vit and fbk files in order to make an ancestry specific frequency matrix. 
#All of these files are very large and hard to manage with programs like R. In that way, this script reads all of this files and makes a frequency count of all SNPs of a specific ancestry
#for further ancestry specific analyses. The script generates a major allele frequency matrix per population, as well as minor allele frequency, total allele frequency and an ancestry 
#specific treemix input file. It from 20 to 80 minutes to run with ~600 individuals and ~250,000-700,000 SNPs. It also depends on the percentage of that specific ancestry.
#The script requires the next input files and parameters.

#Input files:
#	Beagle file
#	Vit file
#	Fbk file
#	Individual ID list with population data

#Population file should have the next columns:
#Individual_ID	Population
#It does not need to have the same order as the beagle and vit

#Input parameters:
#	Ancestry number to analyze
#	Fbk threshold
#	Output name files
#	Missingness flag

#Only SNPs with an fbk probability higher than the fbk threshold will be considered
#The default for the missingness flag is 0, if there are no counts in a SNP for at least one population, the script will ignore that SNP.
#If you want those SNPs to be considered (excluding those SNPs that did not have a single count in any population), write -miss 1

#Output files:
#	Major allele table
#	Minor allele table
#	Total count table
#	TreeMix input

#In all the output tables the rows represent the SNPs, while the columns represent populations. 
#This analysis counts the allele frequencies in all cites that have been clasified as a specific ancestry by RFMix
#Furthermore, this analysis is not individual specific, rather it takes into account all individuals and creates an individual average which represents their population.
#The output files will start with the --out string and will have the next endings: 
# "_maj-allele.txt"
# "_min-allele.txt"
# "_total-count.txt" 
# "_treemix.txt"

#Declare all variables
my (%vit_flag,%fbk_flag,%opts,%pops,%pops_pos)=();
my (@line,@ID_list,@pop_list,@beagle_inds,@line_snp_a,@line_snp_b,@line_snp_t,@treemix)=();
my ($exit,$ID,$i,$j,$row,$pop,$anc,$temp,$allele_a,$allele_b,$total_a,$total_b,$first_allele,$print_a,$print_b,$print_t,$monoallelic,$miss,$miss2,$inds_v,$inds_b,$inds_f,$snps_v,$snps_b,$snps_f)=0;

use strict;
use List::MoreUtils qw(uniq);
use List::Util qw(min max);
use Getopt::Long;
GetOptions(\%opts,'-beagle=s','-vit=s','-popinfo=s','-ancestry=i','-fbk=s','-fbkcut=f','-out=s','--miss=i');

#Default parameters

$opts{"fbkcut"} = $opts{fbkcut} || "0.9";
$opts{"miss"} = $opts{miss} || "0";
$opts{"out"} = $opts{out} || "ancestry_specific_count_output_anc$opts{ancestry}";

#Arguments read and used:
print "\n\tArguments listed \nBeagle: $opts{beagle}\nViterbi: $opts{vit}\nPopinfo: $opts{popinfo}\nAncestry to analyze: $opts{ancestry}\nFbk: $opts{fbk}\nFbk threshold: $opts{fbkcut}\nMissingness flag: $opts{miss}\nOutput: $opts{out}\n\n";

###################################################
## Validation of concordance between input files ##
###################################################

#Confirm files provided as arguments exist
if(`if [ ! -f $opts{beagle} ];then echo "1";fi`){
	print "Beagle file not found: $opts{beagle}\n";
	$exit=1;
}
if(`if [ ! -f $opts{vit} ];then echo "1";fi`){
	print "Vit file not found: $opts{vit}\n";
	$exit=1;
}
if(`if [ ! -f $opts{popinfo} ];then echo "1";fi`){
	print "Popinfo file not found: $opts{popinfo}\n";
	$exit=1;
}
if(`if [ ! -f $opts{fbk} ];then echo "1";fi`){
	print "Fbk file not found: $opts{fbk}\n";
	$exit=1;
}
if(!$opts{ancestry}){
	print "Ancestry number was not provided\n";
	$exit=1;
}
if($opts{miss}>1||$opts{miss}<0){
	print "Invalid flag option: $opts{miss}\n\t(Must be 0 or 1)\n";
	$exit=1;
}
if($opts{fbkcut}>1||$opts{fbkcut}<0){
	print "Invalid fbkcut value: $opts{fbkcut}\n\t(Must be a value between 0 and 1)\n";
	$exit=1;
}

#If something is missing the script will die
if($exit==1){
	die "\nAborted due to missing or invalid arguments\n";
}

#Comparison of total individuals
$inds_v=`wc -l $opts{vit}`;
$inds_v=$inds_v+0;	#I do this step to consider only the first number from the "wc -l" output
$inds_b=`awk '{print NF}' $opts{beagle} | head -n 1`;
$inds_b=$inds_b-2;

$temp=`awk '{print NF}' $opts{fbk} | head -n 1`;
$anc=`cut -f2 $opts{vit} | sort -rn | sed -n '1p'`;	#This command obtains the number of ancestries assuming that just by chance all ancestries will be found on the first ancestry column from the vit file. Thus obtaining the maximum number will be equivalent to getting the total number of ancestry references in RFMix.
$inds_f=`perl -E 'say $temp/$anc'`;

#If there is an input file with a different number of individuals, the script will die.
if($inds_v!=$inds_b&&$inds_b!=$inds_f){
	print "ERROR \nNumber of individuals is not the same in the vit, beagle nor fbk file: $inds_v, $inds_b and $inds_f\n\n";
	exit;
}elsif($inds_v!=$inds_b){
	print "ERROR \nNumber of individuals is not the same in the vit and beagle file: $inds_v and $inds_b \n\n";
	exit;
}elsif($inds_b!=$inds_f){
	print "ERROR \nNumber of individuals is not the same in the beagle and fbk file: $inds_b and $inds_f \n\n";
	exit;
}elsif($inds_v!=$inds_f){
	print "ERROR \nNumber of individuals is not the same in the vit and fbk file: $inds_v and $inds_f \n\n";
	exit;
}else{
	print "Total individuals are concordant in vit, beagle and fbk files: $inds_v individuals\n";
}

#Comparison of total of SNPs
$snps_v=`awk '{print NF}' $opts{vit} | head -n 1`;
$snps_v=$snps_v-1;
$snps_b=`wc -l $opts{beagle} `;
$snps_b=$snps_b-1;
$snps_f=`wc -l $opts{fbk} `;
$snps_f=$snps_f+0;

#If there is an input file with a different number of SNPs, the script will die.
if($snps_v!=$snps_b&&$snps_b!=$snps_f){
	print "ERROR \nNumber of genetic markers is not the same in the vit, beagle nor fbk file: $snps_v, $snps_b and $snps_f \n\n";
	exit;
}elsif($snps_v!=$snps_b){
	print "ERROR \nNumber of genetic markers is not the same in the vit and beagle file: $snps_v and $snps_b \n\n";
	exit;
}elsif($snps_b!=$snps_f){
	print "ERROR \nNumber of genetic markers is not the same in the beagle and fbk file: $snps_b and $snps_f \n\n";
	exit;
}elsif($snps_v!=$snps_f){
	print "ERROR \nNumber of genetic markers is not the same in the vit and fbk file: $snps_v and $snps_f \n\n";
	exit;
}else{
	print "Total genetic markers are concordant in vit, beagle and fbk files: $snps_v markers\n";
}

#############################################
## Analyzing SNP's ancestry and fbk values ##
#############################################

open(VIT,"$opts{vit}");

#First I will store all positions that have the specified ancestry in all individuals. This information will be stored in a hash, thus if the key exists, then that position should be taken into account for the allele counts.
#SNPs are shown in columns while individuals in rows. The first column provides the individual ID. There are no headers.
print "\nReading vit file...\n";
while(<VIT>){
	chomp($_);	
	@line=split("\t",$_);
	$ID=shift(@line);	#Get the individual ID from the first column.
	for($i=0;$i<$#line;$i++){
		if($line[$i]==$opts{ancestry}){
			$vit_flag{"$ID $i"}=1;		#All keys ("Individual_ID SNP") are stored in %vit_flag with the random value 1.
		}
	}
	
	push(@ID_list,$ID);
}

close(VIT);

@line=();

open(FBK,"$opts{fbk}");

#Then, once all positions of the specified ancestry are stored in a hash, they will be re-evaluated according to the fbkcut. In that way not all SNPs with the specified ancestry will be taken into account, just those that have a higher probability than the fbkcut threshold (the analyses usually use fbkcut = 0.9, thus only considering SNPs with an fbk-value > 0.9).
#There are (individuals x ancestries) columns, as this file provides the fbk value of each ancestry in every SNP. E.g. 10 individuals and 3 ancestries will result on 30 columns in the fbk file.
#All SNPs are shown in rows. There are no headers.
print "Reading fbk file...\n";

while(<FBK>){
	chomp($_);
	@line=split(" ",$_);
	for($i=0;$i<$#line;$i=$i+$anc){		#$i increases that way so that only the interest ancestry is taken into account. This avoids unnecessary comparisons.
		if(exists($vit_flag{"$ID_list[$i/$anc] $j"}) && $line[$i+$opts{ancestry}-1]>=$opts{fbkcut}){	#The first condition makes sure the SNP in a specific individual has the interest ancestry. The second one makes sure that ancestry is called with enough certainty: prob >= fbkcut.
			$fbk_flag{"$ID_list[$i/$anc] $j"}=1;	#All keys ("Individual_ID SNP") with the desired ancestry and fbk value are stored in %fbk_flag with the random value 1.
		}		
	}
	$j++;	#This variable will provide the row number/SNP number in the same way as the $i variable when reading the vit file.

}

close(FBK);

%vit_flag=();

#####################
## Reading popinfo ##
#####################

open(POP,"$opts{popinfo}");
print "Reading popinfo...\n";

#Obtaining popinfo data.
while(<POP>){
	chomp($_);
	@line=split("\t",$_);
	$pops{"$line[0]_A"}=$line[1];	#Individual ID (first column) is the key while the population label (second column) is the value.
	$pops{"$line[0]_B"}=$line[1];	#It assumes the haplotypes have exactly the same individual ID plus a "_A" and "_B" ending.
}

close(POP);

##########################################
## Reading beagle and generating output ##
##########################################

#Reading beagle file and generating output files.
open(BEAGLE,"$opts{beagle}");
print "Reading beagle and writing output files...\n";
open(MAJ,">$opts{out}_$opts{ancestry}_maj-allele.txt");
open(MIN,">$opts{out}_$opts{ancestry}_min-allele.txt");
open(TOT,">$opts{out}_$opts{ancestry}_total-count.txt");
open(TREEMIX,">$opts{out}_$opts{ancestry}_treemix.txt");

$row=-1;
#Obtaining population list without repeated elements.
@pop_list=uniq(values %pops);
$temp=join("\t",@pop_list);
#Printing the output files' headers.
print MAJ "$temp\n";
print MIN "$temp\n";
print TOT "$temp\n";
$temp=join(" ",@pop_list);
print TREEMIX "$temp\n";

while(<BEAGLE>){
	chomp($_);
	@line=split("\t",$_);

	if($row==-1){	#Obtaining list of indices per population according to beagle order. This condition is only true when reading the first row. 

		foreach $i (@pop_list){
			for($j=2;$j<$#line;$j++){
				if($pops{$line[$j]} eq $i){
					push(@{$pops_pos{$i}},$j);	#This array will contain the positions of all individuals belonging to every population. E.g. it will store the positions 2,3,4,5 and 6 as those columns correspond to individuals belonging to the first population.
				}
			}

		}

		#Stores all individual's ID in the beagle header into @beagle_inds. The first two elements are removed as they are not IDs. 
		@beagle_inds=@line;
		shift(@beagle_inds);
		shift(@beagle_inds);

		$temp=join("\t",@beagle_inds);
		
	}else{

		($first_allele,$total_a,$total_b,$print_a,$print_b,$print_t)=0;
		(@line_snp_a,@line_snp_b,@line_snp_t)=();
		foreach $i (@pop_list){	#Analyzing by individuals belonging to each population.
			($allele_a,$allele_b)=0;
			foreach $j (@{$pops_pos{$i}}){	#Analyzing per SNP
				$temp="$beagle_inds[$j] $row";
				if(exists($fbk_flag{"$beagle_inds[$j] $row"})){		#This evaluates if the specific individual and SNP were reliable according to the fbk threshold.
					if($first_allele!~m/[ACGT]/){			#If this is the first allele read in that specific allele, that allele will be considered as the allele A.
						$first_allele=$line[$j+2];
						$allele_a++;			#This counts the frequency of allele A on a specific population.
						$total_a++;			#This counts the frequency of allele A on all populations. It is not erased until changing row/SNP in the beagle file.
					}elsif($first_allele eq $line[$j+2]){		#Compares the new allele found with the first allele/allele A. 
						$allele_a++;	
						$total_a++;
					}else{
						$allele_b++;
						$total_b++;
					}
			
				}	
			}

			#Stores all population specific counts in an array. The length of these arrays will be the same as number of populations.
			push(@line_snp_a,$allele_a+0);
			push(@line_snp_b,$allele_b+0);
			push(@line_snp_t,$allele_a+$allele_b);
		
		}

		#Printing SNP counts

		@treemix=();

		if($first_allele=~m/[ACGT]/){	#If there was no allele information read in any individual at certain SNP, it is a case of complete missingness and will not be printed.
			$print_a=join("\t",@line_snp_a);
			$print_b=join("\t",@line_snp_b);
			$print_t=join("\t",@line_snp_t);

			if(0==min(@line_snp_t)){	#If there is at least one population without a single allele called it will be also discarded and will not be printed in the outputs.
				$miss2++;

				if($opts{miss}&&!(max(@line_snp_b)==0)){	#In case the flag miss is turned on, missing SNPs will be printed (except totally missing). (always excludes monoallelic SNPs)
					if($total_a>$total_b){	#The major allele will be the most common allele in the dataset.
						print MAJ "$print_a\n";
						print MIN "$print_b\n";
						print TOT "$print_t\n";
				
						for($i=0;$i<=$#line_snp_a;$i++){
							push(@treemix,"$line_snp_a[$i],$line_snp_b[$i]");
						}

						$temp=join(" ",@treemix);
						print TREEMIX "$temp\n";		
					}else{
						print MIN "$print_a\n";
						print MAJ "$print_b\n";
						print TOT "$print_t\n";

						for($i=0;$i<=$#line_snp_a;$i++){
							push(@treemix,"$line_snp_b[$i],$line_snp_a[$i]");
						}

						$temp=join(" ",@treemix);
						print TREEMIX "$temp\n";

					}
				}			


			}elsif(max(@line_snp_b)==0){	#If there is not even a single minor allele called at a SNP, it will be considered as monoallelic and will no be printed as it is not useful.
				$monoallelic++;
			}elsif($total_a>$total_b){	#The major allele will be the most common allele in the dataset.
				print MAJ "$print_a\n";
				print MIN "$print_b\n";
				print TOT "$print_t\n";
					
				for($i=0;$i<=$#line_snp_a;$i++){
					push(@treemix,"$line_snp_a[$i],$line_snp_b[$i]");
				}

				$temp=join(" ",@treemix);
				print TREEMIX "$temp\n";		
			}else{
				print MIN "$print_a\n";
				print MAJ "$print_b\n";
				print TOT "$print_t\n";

				for($i=0;$i<=$#line_snp_a;$i++){
					push(@treemix,"$line_snp_b[$i],$line_snp_a[$i]");
				}

				$temp=join(" ",@treemix);
				print TREEMIX "$temp\n";

			}

		}else{
			$miss++;
		}

	}
	
	$row++;		#This variable will provide the row number/SNP number in the same way as the $j variable when reading the fbk file.

}

close(BEAGLE);

(%vit_flag,%fbk_flag,%pops,%pops_pos)=();

if($miss>0){
	print "\n$miss complete missing data at a specific SNP were found and excluded\n";
}
if($miss2>0){
	if($opts{miss}){
		print "$miss2 missing data in at least one population at certain SNP (considered in the output files)\n";
	}else{
		print "$miss2 missing data in at least one population at certain SNP (excluded in the output files)\n";
	}
}
if($monoallelic>0){
	print "$monoallelic monoallelic SNPs (excluding monoallelic SNPs with missing data)\n";
}


$temp=`wc -l $opts{out}_$opts{ancestry}_total-count.txt | cut -f1 -d " "`;
chomp($temp);
$temp=$temp-1;
print "\n$temp SNPs remaining\n";

close(MAJ);
close(MIN);
close(TOT);
close(TREEMIX);
