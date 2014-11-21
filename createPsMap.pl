#!/usr/bin/perl

# How to use:
# -----------
#
# perl createPsMap.pl -i Athaliana_peptide.fa -d phyloBlastDB_Drost_Gabel_Grosse_Quint.fa -p BLAST_Athaliana -r athaliana_blast_results -t 30 -a 64

use Bio::SeqIO;
use Getopt::Long;


my $refOrg = "Athaliana_peptide.fa";
my $db = "phyloBlastDB_Drost_Gabel_Grosse_Quint.fa";
my $save = "BLAST_Athaliana";
my $table = "athaliana_blast";
my $threshold = 30;
my $threads = 1;
my $evalue = 1e-5;

GetOptions (
"input=s" => \$refOrg,
"database=s" => \$db,
"prefix=s" => \$save,
"resultTable=s" => \$table,
"threshold=i" => \$threshold,
"a=i" => \$threads,
"evalue=s" => \$evalue
);

my @blastFiles = glob $db."*.phr";

if(scalar @blastFiles == 0){
    print "Creating BLAST database ...\n";
    system("formatdb -p T -i ".$db);
}

print "Loading ".$refOrg."...\n";

# Loading Genes
my $in = Bio::SeqIO->new(-format => "fasta",
-file => $refOrg);

# Loading and writing sequences into an Array
# makes it easier to restart at a specific sequences if
# the server crashes and the script has to be restarted

my @genome = ();

while(my $seq = $in->next_seq()){
    
    push(@genome,$seq);
}
print scalar(@genome)." sequences of ".$refOrg." will now be blasted against the database.\n\n";

#####################################################################################
#   Begin to BLAST, write the BLAST results into a MySQL-File for further studies   #
#   and create the phylostratigrap                                                  #
#   *only AS-sequences with a minimal length of $threshold were send to BLAST       #
#####################################################################################

my $geneNumber = scalar @genome -1;

for my $geneIndex (0..$geneNumber){
    
    my $gene = $genome[$geneIndex];
    
    my $orgName;
    
    if($gene->description =~ /\| \[(.*)\]\s\| \[(.*)\]/){
        $orgName = $1;
    }
    
    print "\nGene ->".$gene->display_id." --- ".$geneIndex." of ".$geneNumber."\n";
    
    if(length($gene->seq) <= $threshold){
        open(OUT,">>".$save."_short_queries(".$threshold.").fa") or die $!;
        print OUT $gene->display_id,"\n";
        print OUT $gene->seq,"\n";
        close(OUT);
        next;
    }
    
    $save = $save."_".$geneIndex;
    blastP($db,$gene,$save,$table,$orgName,$threads,$evalue);
}

####################################################
## blastP against local sequence database         ##
##                                                ##
## 1.Arg: BLAST database name                     ##
## 2.Arg: gene (Bio::Seq)                         ##
## 3.Arg: mysql file name                         ##
## 4.Arg: mysql table name                        ##
## 5.Arg: name of reference organism              ##
## 6.Arg: number of threads for BLAST             ##
## 7.Arg: BLAST e-value cutoff                    ##
####################################################
sub blastP{
    
    use Bio::Tools::Run::StandAloneBlast;
    use Bio::Tools::Run::StandAloneNCBIBlast;
    use Bio::SearchIO;
    use Bio::Search::HSP::HSPI;
    
    my $db = shift;
    my $gene = shift;
    my $save = shift;
    my $table = shift;
    my $refOrgName = shift;
    my $threads = shift;
    my $e_val = shift;
    
    my $psFile = $refOrgName;
    $psFile =~ s//"_"/g;
    
    my %psMap = ();
    
    my $refTreeLine = "";
    
    if($gene->description =~ /\| \[(.*)\]\s\| \[(.*)\]/){
        $refTreeLine = $2;
    }
    
    my $path = `pwd`;
    chomp($path);
    
    my $maxAlignments = 65200; # maximum number of alignments which can be handled with blastall
    
    @params = (-program  => 'blastp',
    -database =>  $path."/".$db,
    -expect => $e_val,
    -b => $maxAlignments,
    -a => $threads);
    
    $blast_obj = Bio::Tools::Run::StandAloneBlast->new(@params);
    $report_obj = $blast_obj->blastall($gene);
    
    #####################################
    while( my $result = $report_obj->next_result ) {
        
        my $previousAccession = "";
        
        while( my $hit = $result->next_hit ) {
            
            while( my $hsp = $hit->next_hsp ) {
                
                if(defined $hit and defined $hsp){
                    
                    
                    my $orgName = "";
                    my $hitTreeLine = "";
                    
                    if($hit->description =~ /\| \[(.*)\] \| \[(.*)\]/){
                        $orgName = $1;
                        $hitTreeLine = $2;
                        print $orgName," --> ",$hitTreeLine."\n";
                    }else{
                        next;
                    }
                    
                    if($orgName eq ""){ next; }
                    
                    if($hit->accession eq $previousAccession){
                        next;
                    }
                    
                    $previousAccession = $hit->accession;
                    
                    my @ps = ();
                    
                    if($refOrgName eq $orgName){
                        my @refOrgTree = split("; ",$refTreeLine);
                        # PS of reference organism is length of taxonomy (beginning at cellular organism up to
                        # the organisms taxid itselfs)
                        @ps = ($refOrgName, scalar(@refOrgTree)+2);
                        print join(" --> ",@ps),"\n";
                    }else{
                        @ps = get_PS($refTreeLine,$hitTreeLine);
                    }
                    
                    print $ps[0]."\t\t".$ps[1]."\t\t".$hit->accession."\t\t".$orgName."\n";
                    
                    
                    if(exists $psMap{$gene->display_id}){
                        if($psMap{$gene->display_id} > $ps[1]){
                            $psMap{$gene->display_id} = $ps[1];
                        }
                    }else{
                        $psMap{$gene->display_id} = $ps[1];
                    }
                    
                    open(OUT,">>".$path."/".$save.".sql") or die $!;
                    my $sql = 	"INSERT INTO ".$table." ".
                    "VALUES(".
                    "\'".$ps[0]."\',". # PS name
                    "\'".$ps[1]."\',". # PS number
                    "\'".$orgName."\',". # hit organism
                    "\'".$gene->display_id."\',". # query id
                    "\'".$hit->accession."\',". # hit id
                    "\'".$hsp->length('query')."\',".
                    "\'".$hsp->length('hit')."\',".
                    "\'".$hsp->percent_identity."\',".
                    "\'".$hsp->evalue."\',".
                    "\'".$hsp->significance."\',".
                    "\'".$hsp->score."\',".
                    "\'".$hsp->frac_identical('query')."\',".
                    "\'".$hsp->frac_identical('hit')."\');";
                    print OUT $sql."\n";
                    close(OUT);
                }
            }
        }
    }
    
    open(OUT,">".$path."/".$psFile."_psMap.csv") or die $!;
    print OUT "ID;PS";
    while(my($k,$v) = each %psMap){
        print OUT $k.";".$v."\n";
    }
    close(OUT);
}


sub get_PS{
    
    my $refTree = shift;
    my $hitTree = shift;
    
    my @refLineage = split("; ",$refTree);
    my @hitLineage = split("; ",$hitTree);
    
    my $MAXps = scalar @taxNodes;
    my @ps = ();
    
    
    # comparing taxonomy begins at the root (e.g. at node Archaea or Bacteria or Eukaryota)
    for(my $i=0, my $j=0;$i<scalar(@refLineage) && $j<scalar(@hitLineage);$i++,$j++){
        #print "i: ".$i."   ".$taxNodes[$i]."  <--> ".$hitLineage[$j]."\n";
        if($refLineage[$i] ne $hitLineage[$j]){
            if($i == 0 && $j == 0){
                $ps[0] = "Cellular organisms";
                $ps[1] = $i+1;
            }else{
                $ps[0] = $refLineage[$i-1];
                $ps[1] = $i+1;
            }
            last;
        }
        # should happen if two compared (different) organisms share the same taxonomy
        if($i == (scalar(@refLineage)-1) ){
            $ps[0] = $refLineage[$i];
            $ps[1] = $i+1;
            last;
        }
    }
    
    return @ps;
}
