
##perl script
##Using PGA to score risk genes and summary the resutls

system ("cd /Users/zhenwang/Documents/Project5_AD/4_PGAanlysis_new6/1_ALZ/1_stringent/EUR_PGA_v2.0.4 ");
system ("java -Xmx4g -cp src edu.yu.einstein.zdzlab.pga.GWASAnalyzer ALZ EUR --c 0.5");
system ("cd /Users/zhenwang/Documents/Project5_AD/4_PGAanlysis_new6/1_ALZ/2_unstringent/EUR_PGA_v2.0.4");
system ("java -Xmx4g -cp src edu.yu.einstein.zdzlab.pga.GWASAnalyzer ALZ EUR --c 0.5");


$score=21.4;
$score2=13.1;
$loci_blank=0;
$loci_all=0;
$loci_high=0;
$loci_non_strin_only=0;
$loci_strin=0;
$loci_length=0;
$loci_high_length=0;

open(AAA,"3_2_stringent_and_un_merge.infor.txt");
open(BBB,">3_3_stringent_and_un_merge.infor.highRisk.txt");
open(CCC,">>0_results_ALS.gene_region.summary.txt");

while($line=<AAA>){
    $n++;
    chomp $line;
    @arraynum1=split(/\t/,$line);
    $loci_all++;
    if(length($arraynum1[2])==0){$loci_blank++;}
    @pos=split(/:/,$arraynum1[0]);
    @sta_end=split(/\-/,$pos[1]);
    $length=$sta_end[1]-$sta_end[0]+1;
    $loci_length=$loci_length+$length;
    @gene=split(/\,/,$arraynum1[2]);
    @gene_non=split(/\,/,$arraynum1[5]);
    
    my %gene_high;
    my %gene_high_non;
    $gene_high_srt="";
    $gene_high_non_srt="";
    $gene_high_num=0;
    $gene_high_non_num=0;

    ###count the number of the high risk gene
    for $num(0..$#gene){
        $gene[$num]=~s/\s+//g;
        $gene[$num]=~s/\)//g;
        $gene_non[$num]=~s/\s+//g;
        $gene_non[$num]=~s/\)//g;
        @gene_one=split(/\(/,$gene[$num]);
        @gene_one_non=split(/\(/,$gene_non[$num]);
        $gene_one=~s/\s+//g;
        $gene_one_non=~s/\s+//g;
       if($gene_one[1]>=$score){ $gene_high_srt=$gene_high_srt.$gene[$num]."),"; $gene_high_num++;  }
       if($gene_one_non[1]>=$score2){ $gene_high_non_srt=$gene_high_non_srt.$gene_non[$num]."),"; $gene_high_non_num++;  }
       
    }
        
    if($gene_high_num>=1 or $gene_high_non_num>=1){
        $loci_high++;
        $loci_high_length=$loci_high_length+$length;
        if($gene_high_num==0){$gene_high_srt="NA"; $loci_non_strin_only++;}
        if($gene_high_non_num==0){$gene_high_non_srt="NA";}
        if($gene_high_num>=1){$loci_strin++;}
        print BBB "$arraynum1[0]\t$length\t$arraynum1[1]\t$gene_high_srt\t$gene_high_non_srt\n";
    }
}

print "All\tNotarget\tHighScore\tStringent\tnonStringentOnly\tLengthALL\tLengthHigh\n";
print "$loci_all\t$loci_blank\t$loci_high\t$loci_strin\t$loci_non_strin_only\t$loci_length\t$loci_high_length\n";
print CCC "All\tNotarget\tHighScore\tStringent\tnonStringentOnly\tLengthALL\tLengthHigh\n";
print CCC "$loci_all\t$loci_blank\t$loci_high\t$loci_strin\t$loci_non_strin_only\t$loci_length\t$loci_high_length\n";

close AAA;
close BBB;
close CCC;

