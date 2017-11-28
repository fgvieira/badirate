#!/usr/bin/env perl
#
#    BadiRate - Estimating gene family turnover rates by likelihood-based methods
#      Copyright (C) 2012 Pablo Librado, Filipe Garrett Vieira, Julio Rozas and Universitat de Barcelona (UB)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see {http://www.gnu.org/licenses/}.


=head1 NAME
BadiRate.pl - DESCRIPTION 

=head1 SYNOPSIS
The basic command is:

    perl BadiRate.pl -tree NEWICK_FILE -fam FAMILY_SIZE_FILE [options] > output.bd
 
	-anc         = Reconstruct ancestral family sizes and the minimum number of gains/losses in each lineage   
	-bmodel      = Run local or free branch models 	   
	-sizefile    = Family Size File 	   
        -unobs       = Correct the likelihood for families absent in all extant species
	-h|help	     = Display this help	   
	-rmodel      = Family turnover rates to be estimated	   
	-out         = Set the output file	   
	-outlier     = Report families escaping the estimated stochastic process	   
	-ep          = Define the estimation procedure
	-print_ids   = Display nodes ids in Newick format	   
        -priorfile   = Prior File
        -n_max_int   = Modify the maximum number of family members in the internal phylogenetic nodes 
        -root_dist   = Estimation method for the root a priori distribution
	-seed	     = Seed of the pseudo-random number generator	   
        -start_val   = Starting values for the likelihood methods
	-treefile    = Phylogenetic tree in Newick format	      
	-version     = Report the BadiRate version	 


    -----OR------ 

    perl BadiRate.pl controlfile.bd > output.bd

=head1 DESCRIPTION

This script estimates the gain, birth, death and innovation gene (or DNA element) family rates. It implements three stochastic population models: 
1) BDI: Birth-Death-Innovation
2) LI: Lambda-Innovation, where Lambda is the birth-and-death parameter (equal birth and death rates are assumed)
3) GD: Gain-Death, where gain parameter accounts for all gains regardless they are originated by DNA duplications or represent an innovation
4) BD: Birth and Death model
5) L: Lambda, where Lambda is the birth-and-death parameter (equal birth and death rates are assumed)
In addition, for each family or subfamily, it also infers the most likely ancestral content and the (sub)families unlikely evolving under the estimated stochastic process.

=head1 AUTHOR - Pablo Librado Sanz

Email 

plibrado@ub.edu

=head1 CONTRIBUTORS

fgarret@ub.edu
jrozas@ub.edu

=cut

# Let the code begin...
use FindBin;                
use lib $FindBin::Bin."/lib";
use Inline (C =>'DATA');
use strict;
use warnings;
use Bio::TreeIO;
use Getopt::Long;
use Statistics::Descriptive;
use Math::Amoeba qw(MinimiseND);

$|=1;

#declarations
my (%ortho,%options,%distances,%assoc,%fixed,%max,%min,%difOGs,%bmodel,%priormodel,%OGstring)=();
my (%birth,%death,%innovation,%ancestral,%root,@treepath)=();
my (@usedOG,@nodes,%RatesOrder,%RatesOrder_rev,%OGOrder,%anc_cnt,%Q,@ant_param)=();
my ($tree,$OG_analyzed,$version, $error);
my $MIN=0.0000000001;
my $MAX=99999999999999999999999999;
my $U;
$error=0;
my %ERRORS=(
    '0'=>'NO ERRORS',
    '1'=>'Amoeba excedeed the maximum number of iterations. Check convergence (see -start_val option)!',
    '2'=>'Try using a more complex model, or changing the starting values. See the -rmodel, -bmodel and -start_val options',
    );

$version="1.35.00";

$options{'seed'}=int(time()*rand);
$options{'anc'}=0;
$options{'family'}=0;
$options{'outlier'}=0;
$options{'rmodel'}="BDI";
$options{'ep'}="ML";
$options{'help'}=0;
$options{'version'}=0;
$options{'print_ids'}=0;
$options{'priorfile'}="";
$options{'unobs'}=0;
$options{'root_dist'}=1;
$options{'start_val'}=0;
$options{'n_max_int'}=10;

&dieNice if (!defined($ARGV[0]));

if ($ARGV[0]!~/^-/){
    (!&readCtlFile(\%options,$ARGV[0])) && do{die "Unable to read control file\n";};
}else{
    GetOptions(
    'sizefile=s'             => \$options{'sizefile'},        
    'unobs!'                 => \$options{'unobs'},	
    'treefile=s'             => \$options{'treefile'},       
    'root_dist=i'            => \$options{'root_dist'},
    'out:s'                  => \$options{'out'},        
    'seed:i'                 => \$options{'seed'},             
    'bmodel:s'               => \$options{'bmodel'},     
    'family!'                => \$options{'family'},     
    'rmodel:s'               => \$options{'rmodel'},      
    'outlier!'               => \$options{'outlier'},    
    'anc!'                   => \$options{'anc'},        
    'help|h!'                => \$options{'help'},          
    'version!'               => \$options{'version'},    
    'print_ids!'             => \$options{'print_ids'},  
    'priorfile:s'            => \$options{'priorfile'},	
    'start_val:i'            => \$options{'start_val'},	
    'ep:s'                   => \$options{'ep'},	
    'n_max_int:i'            => \$options{'n_max_int'}	
    );
}
die "BadiRate version: $version\n" if ($options{'version'});

#set the seed                                                                    
srand($options{'seed'});

#options checking
$options{'rmodel'}=uc($options{'rmodel'});
$options{'ep'}=uc($options{'ep'});

&dieNice if (!$options{'sizefile'} || !$options{'treefile'} || $options{'help'});
die "model must be BDI, BD, GD, LI or L\n" if ($options{'rmodel'}!~/^(BDI|GD|LI|BD|L)$/);
$options{'anc'}=1 if ($options{'family'} || $options{'outlier'});
die "root_dist option must be between 0 and 4\n" if ($options{'root_dist'} < 0 && $options{'root_dist'}>4);
die "unobs requires root_dist = 1 or root_dist = 3\n" if ($options{'unobs'} && ($options{'root_dist'} == 0 || $options{'root_dist'} == 2 || $options{'root_dist'} == 4));
die "Estimation procedure must be ML, MAP, CML, CMAP, CWP or CSP\n" if ($options{'ep'}!~/^(ML|MAP|CML|CMAP|CWP|CSP)$/);
die "MAP and CMAP estimation procedures require a priorfile defined\n" if ($options{'priorfile'} eq "" && $options{'ep'}=~/MAP/);
die "priorfile defined requires MAP or CMAP estimation procedures\n" if ($options{'priorfile'} ne "" && $options{'ep'}!~/MAP/);
die "-family option must be specified with CML, CMAP, CWP or CSP\n" if ($options{'family'} && $options{'ep'}=~/^(ML|MAP)$/);

#prepare output file
if(defined($options{'out'})) {
    if (uc($options{'out'}) ne "STDOUT"){
	open(ANC_OUT, '>', $options{'out'}) or die "Couldn't open: $!";
    }  else{
	open(ANC_OUT, ">&", \*STDOUT) or die "Couldn't open: $!";
        $options{'out'}="STDOUT";
    }
} else {
    open(ANC_OUT, ">&", \*STDOUT) or die "Couldn't open: $!";
    $options{'out'}="STDOUT";
}

#readortho
(!&readorthotable($options{'sizefile'},\%ortho)) && do{die $!;};

#readnewick
$tree=&readtree($options{'treefile'});
my $root_node=$tree->get_root_node;
my $root_node_id=$root_node->internal_id;

(!&printTreeIDRel(\*ANC_OUT, $tree, \%assoc)) && do{die;};
exit if ($options{'print_ids'});
(!&MakeTreePath(\@treepath)) && do{die;};
(!&GetDistances($tree,\%distances)) && do{die;};
(!&fixedNode) && do{die;};
(!&usedOGs) && do{die;};
die "No valid OG groups or families\n" if (scalar(@usedOG)==0);

foreach my $node ($tree->get_nodes){next unless ($node->ancestor);push(@nodes,$node);}
(!&readmodel($options{'bmodel'})) && do{die;};
if (!$options{'ep'}!~/^(CSP|CWP)$/){
    (!&readprior($options{'priorfile'})) && do{die;};
}

#print input stuff
print ANC_OUT "-"x20,"\n";
print ANC_OUT "INPUT\n";
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);

printf ANC_OUT "\tExecution date: %4d-%02d-%02d %02d:%02d:%02d\n", $year+1900,$mon+1,$mday,$hour,$min,$sec;;
print ANC_OUT "\tVersion: $version\n";
foreach my $option (keys %options){
    if ($option eq "bmodel"){
	print ANC_OUT "\tbmodel= \n";
	foreach my $branch (keys %bmodel){print ANC_OUT "\t\t$branch\t$bmodel{$branch}\n";}
    }elsif ($option ne "bmodel"){
	print ANC_OUT "\t$option = $options{$option}\n";
    }
}
print ANC_OUT "END INPUT\n";

#output
print ANC_OUT "-"x20,"\n";
print ANC_OUT "OUTPUT\n";

print ANC_OUT "\n\t##Family Turnover Rates\n";
(!&getDifferentOGs) && do{die;};
if ($options{'unobs'}){
    foreach my $node ($tree->get_leaf_nodes){
	$ancestral{"UnobservedData"}{$node->internal_id}[0]=0;
	$min{"UnobservedData"}=0;
	$max{"UnobservedData"}=$max{$usedOG[0]}; 
    }
}

(!&AncestralParsimony) && do{die;};
(!&sumAnc) && do{die;};
if ($options{'ep'}=~/(ML|MAP)/){
    (!&Ulimit) && do{die;};
}

#Rates estimation    
if ($options{'ep'} eq 'ML' || $options{'ep'} eq 'MAP'){ #likelihood methods
    (!&FamilyObservedRates(\*ANC_OUT,0)) && do{die;} if ($options{'start_val'} == 0);
    if ($options{'rmodel'} eq "BDI"){
	print ANC_OUT "\t\t#Likelihood: ",&BDI,"\n";
    }elsif ($options{'rmodel'} eq "GD"){
	print ANC_OUT "\t\t#Likelihood: ",&GD,"\n";
    }elsif ($options{'rmodel'} eq "LI"){
	print ANC_OUT "\t\t#Likelihood: ",&LI,"\n";
    }elsif ($options{'rmodel'} eq "BD"){
	print ANC_OUT "\t\t#Likelihood: ",&BD,"\n";
    }elsif ($options{'rmodel'} eq "L"){
	print ANC_OUT "\t\t#Likelihood: ",&L,"\n"
    }

    if ($error!=2){
	warn "\t\t#WARN: ".$ERRORS{$error}."\n" if ($error == 1);
	(!&printRates(\*ANC_OUT)) && do{die;};
    }else{
	warn "\t\t#WARN: ".$ERRORS{$error}."\n";
    }
}elsif ($options{'ep'} eq 'CML' ||$options{'ep'} eq 'CMAP'){ #hybrid methods
    (!&FamilyObservedRates(\*ANC_OUT,0)) && do{die;} if ($options{'start_val'} == 0); 
    if ($options{'rmodel'} eq "BDI"){&BDI;}
    elsif ($options{'rmodel'} eq "GD"){&GD;}
    elsif ($options{'rmodel'} eq "LI"){&LI;}
    elsif ($options{'rmodel'} eq "BD"){&BD;}
    elsif ($options{'rmodel'} eq "L"){&L;}
    if ($error != 2){
	warn "\t\t#WARN: ".$ERRORS{$error}."\n" if ($error == 1);
	my @numAnc=keys %difOGs;
	my $numAncSize=scalar(@numAnc);
	for (my $j=1; $j<=$numAncSize;$j++){
	    $OG_analyzed=$difOGs{$numAnc[$j-1]}[0];
	    my $tmp=&Anc;
	}
	foreach my $OG (@usedOG){ %{$ancestral{$OG}}=%{$ancestral{$difOGs{$OGstring{$OG}}[0]}}; }
	(!&sumAnc) && do{die;};
	(!&FamilyObservedRates(\*ANC_OUT,1)) && do{die;};
	if (scalar(&keysModel)>1){
	    (!&readmodel) && do{die;};
	    (!&FamilyObservedRates(\*ANC_OUT,2)) && do{die;};
	    (!&readmodel($options{'bmodel'})) && do{die;};
	}
    }else{
	warn "\t\t#WARN: ".$ERRORS{$error}."\n";
    }
}elsif ($options{'ep'} eq 'CWP' || $options{'ep'} eq 'CSP'){ #parsimony methods
    (!&FamilyObservedRates(\*ANC_OUT,1)) && do{die;};
}

#ancestral and outlier estimation
if ($options{'anc'} && $error != 2){
    if ($options{'ep'}=~/(ML|MAP)/ && $options{'family'} == 0){
	my @numAnc=keys %difOGs;
	my $numAncSize=scalar(@numAnc);
	for (my $j=1; $j<=$numAncSize;$j++){
	    $OG_analyzed=$difOGs{$numAnc[$j-1]}[0]; 
	    my $tmp=&Anc;
	}
	foreach my $OG (@usedOG){ %{$ancestral{$OG}}=%{$ancestral{$difOGs{$OGstring{$OG}}[0]}}; }
	(!&sumAnc) && do{die;};
    }

    print ANC_OUT "\n\t##Ancestral Family Size\n";
    print ANC_OUT "\t\t#Family\tAncestral Family Size Tree\n";
    foreach my $OG (sort @usedOG){
	my $treecp=$tree->clone;
	foreach my $node ($treecp->get_nodes){
	    if (!$node->is_Leaf()){
	        $node->bootstrap($ancestral{$OG}{$node->internal_id}[0]);
	    }else{
	        $node->id($node->id."_$ancestral{$OG}{$node->internal_id}[0]");
	    }
	}
	print ANC_OUT "\t\t$OG\t",&printTree(\$treecp);
    }

    foreach my $node ($tree->get_nodes){
	if (!$node->is_Leaf()){
	    $node->bootstrap($anc_cnt{$node->internal_id});
	}else{
	    $node->id($node->id."_$anc_cnt{$node->internal_id}");
	}
    }
    print ANC_OUT "\t\tTotal Ancestral Size\t",&printTree(\$tree);
    print ANC_OUT "\n\t##Minimum number of gains and losses per branch\n";
    print ANC_OUT "\t\t#Branch\tGains\tLosses\n";
    foreach my $node (@nodes){
	my $parent=$node->ancestor;
	my $branch=$parent->internal_id."->".$node->internal_id;
	my $branchmodel=$bmodel{$branch};
	my ($gain,$loss)=&calcGainLoss($parent->internal_id,$node->internal_id,0);
	print ANC_OUT "\t\t$branch\t".$gain."\t".$loss."\n";
    }
    
    if ($options{'outlier'}){
	if ($options{'ep'}=~/^(MAP|ML)$/){
	    print ANC_OUT "\n\t##Outlier Families per Branch\n";
	    my (%adjp,@rawp, %pval)=();
	    @rawp=("-5");
	    my $cnt=1;
	    foreach my $OG (@usedOG){
		foreach my $node (@nodes){
		    my $parent=$node->ancestor();
		    my $branch=$parent->internal_id."->".$node->internal_id;
		    my $branchmodel=$bmodel{$branch};
		    my $time=$distances{$parent->internal_id}{$node->internal_id};
		    my $p=&prob($birth{$branchmodel}[0],$death{$branchmodel}[0],$innovation{$branchmodel}[0],$branchmodel."_".$time,$ancestral{$OG}{$parent->internal_id}[0],$ancestral{$OG}{$node->internal_id}[0]);
		    $pval{$OG}{$branch}=$cnt++;
		    push(@rawp,$p);
		}
	    }

	    (!&fdr(\@rawp,\%adjp)) && do{die;};
	    print ANC_OUT "\t\t#Family\tBranch\tP-value\tFDR_P-value\n",;
	    foreach my $OG (@usedOG){
		foreach my $node (@nodes){
		    my $parent=$node->ancestor();
		    my $branch=$parent->internal_id."->".$node->internal_id;
		    print ANC_OUT "\t\t",$OG,"\t",$branch,"\t",&printFloatNice($rawp[$pval{$OG}{$branch}]),"\t",&printFloatNice($adjp{$rawp[$pval{$OG}{$branch}]}),"\n" if ($adjp{$rawp[$pval{$OG}{$branch}]}<0.05);
		}
	    }
	}else{
	    print "To detect outlier families, run ep=ML|MAP\n";
	}
    }
}

print(ANC_OUT "\n\t##Execution time (seconds): ".(time - $^T)."\n");
print(ANC_OUT "END OUTPUT\n");
close(ANC_OUT);
exit(0);

###############################FUNCTIONS########################################
################################################################################
sub readorthotable($$){
#Reads the orthotable. Also an orthogroup id can be selected
    my ($orthofile,$ortho)=@_;
    (!open (ORTHOTABLE, $orthofile)) && do{print "$!\n";return 0;};
    my (@anc,@sp)=();
    my $cnt=0;
    while(<ORTHOTABLE>){
        chomp;
	next if($_ =~ /^#/ || $_!~/\S+/);
        my @orthogroup=split(/\t+|\s+/,$_);
	if ($cnt==0){ #header
	    for (my $j=1; $j<=$#orthogroup; $j++){
		push(@sp,$orthogroup[$j]);
	    }
	}else{
	    return 0 if (scalar(@sp)+1 != scalar(@orthogroup));
	    for (my $j=0; $j<scalar(@sp); $j++){
		push(@{$$ortho{$orthogroup[0]}{'ID'}},{'cnt'=>$orthogroup[$j+1], 'sp' => $sp[$j]});
	    }   
	}
	$cnt++;
    }
    
    close(ORTHOTABLE);
    return 1;
}
########################################################################
sub prior($$){
    my ($x,$hash)=@_;
    if (defined($$hash{'uniform'})){
	return 1/$$hash{'uniform'};
    }elsif (defined($$hash{'shape'})){
	if ($x == 0 && $$hash{'shape'} == 0){
	    return 1;
	}elsif ($x > 0 && $$hash{'shape'} == 0){
	    return 0;
        }else{
	    return &gamma_PDF($x,1,$$hash{'shape'});
	}
    }
}
########################################################################
sub readprior($){
    my ($priorfile)=@_;
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	if ($options{'rmodel'} eq "BDI"){
	    $priormodel{$bmodelstr}{'birth'}{'uniform'}=$priormodel{$bmodelstr}{'death'}{'uniform'}=$priormodel{$bmodelstr}{'innovation'}{'uniform'}=1;
	}elsif ($options{'rmodel'} eq "LI"){
	    $priormodel{$bmodelstr}{'lambda'}{'uniform'}=$priormodel{$bmodelstr}{'innovation'}{'uniform'}=1;
	}elsif ($options{'rmodel'} eq "GD"){
	    $priormodel{$bmodelstr}{'gain'}{'uniform'}=$priormodel{$bmodelstr}{'death'}{'uniform'}=1;
	}elsif ($options{'rmodel'} eq "BD"){
            $priormodel{$bmodelstr}{'birth'}{'uniform'}=$priormodel{$bmodelstr}{'death'}{'uniform'}=1;
	}elsif ($options{'rmodel'} eq "L"){
	    $priormodel{$bmodelstr}{'lambda'}{'uniform'}=1;
	}
    }
    
    return 1 if (!-s $priorfile);

    (!open(FILE,$priorfile)) && do{die;};
    while(<FILE>){
	chomp;
	next if ($_!~/\S+/ || $_=~/^#/);
	my @line=split(/\t+/,$_);
	die "No $line[0] branch class\n" if (scalar(grep(/^$line[0]$/,@keys))==0);
	if ($options{'rmodel'} eq "BDI"){
	    if ($line[1] ne "birth" && $line[1] ne "death" && $line[1] ne "innovation"){
		print STDERR "$line[1] is not birth, death or innovation\n";
		return 0;
	    }
	}elsif ($options{'rmodel'} eq "LI"){
	    if ($line[1] ne "lambda" && $line[1] ne "innovation"){
		print STDERR "$line[1] is not lambda or innovation\n";
		return 0;
	    }
	}elsif ($options{'rmodel'} eq "L"){
            if ($line[1] ne "lambda"){
                print STDERR "$line[1] is not lambda\n";
                return 0;
            }    
	}elsif ($options{'rmodel'} eq "GD"){
	    if ($line[1] ne "gain" && $line[1] ne "death"){
		print STDERR "$line[1] is not gain or death\n";
		return 0;
	    }
	}elsif ($options{'rmodel'} eq "BD"){
            if ($line[1] ne "birth" && $line[1] ne "death"){
                print STDERR "$line[1] is not birth or death\n";
                return 0;
            }
        }

	$priormodel{$line[0]}{$line[1]}=();
	if (uc($line[2]) eq "GAMMA"){
	    if ($line[3]<0){
		die "Wrong gamma parameter\n";
	    }
	    $priormodel{$line[0]}{$line[1]}{'shape'}=$line[3];
	    #$priormodel{$line[0]}{$line[1]}{'scale'}=$line[4];
	}else{
	    die "No Gamma distribution as prior\n";
	}
    }
    close(FILE);
    return 1;
}
########################################################################
sub printTree($){
    my ($ptree)=@_;
    my $tree_string = '';
    open(TMP_OUT, "+<", \$tree_string);
    my $out = new Bio::TreeIO(-fh => \*TMP_OUT, -format => 'newick');
    $out->write_tree($$ptree);
    return $tree_string;
}
########################################################################
sub printRates($$){
    my ($fh,$mode)=@_;
    my @keys=&keysModel;
    $mode = 1 if (!defined($mode));
    if ($mode != 2){
	if ($options{'rmodel'} eq "BDI"){
	    print $fh "\t\t#Branch_Group\tBirth\tDeath\tInnovation\n";
	}elsif ($options{'rmodel'} eq "GD"){
	    print ANC_OUT "\t\t#Branch_Group\tGain\tDeath\tLoss(approx from death)\n";
	}elsif ($options{'rmodel'} eq "LI"){
	    print ANC_OUT "\t\t#Branch_Group\tLambda\tInnovation\n";
	}elsif ($options{'rmodel'} eq "L"){
		print ANC_OUT "\t\t#Branch_Group\tLambda\n";
	}elsif ($options{'rmodel'} eq "BD"){
            print $fh "\t\t#Branch_Group\tBirth\tDeath\n";
	}
    }

    foreach my $bmodelstr (sort {$a <=> $b} @keys){
	my $bmodelstr2=$bmodelstr;
	$bmodelstr2="Global_Rates" if ($mode == 2);
	if ($options{'rmodel'} eq "BDI"){
	    print $fh "\t\t$bmodelstr2\t",&printFloatNice($birth{$bmodelstr}[0]),"\t",&printFloatNice($death{$bmodelstr}[0]),"\t",&printFloatNice($innovation{$bmodelstr}[0]),"\n";
	}elsif ($options{'rmodel'} eq "GD"){
	    my %loss=&death2loss(\%death);
	    print $fh "\t\t$bmodelstr2\t",&printFloatNice($innovation{$bmodelstr}[0]),"\t",&printFloatNice($death{$bmodelstr}[0]),"\t",&printFloatNice($loss{$bmodelstr}[0]),"\n";
	}elsif ($options{'rmodel'} eq "LI"){
	    print $fh "\t\t$bmodelstr2\t",&printFloatNice($birth{$bmodelstr}[0]),"\t",&printFloatNice($innovation{$bmodelstr}[0]),"\n";
	}elsif ($options{'rmodel'} eq "L"){
            print $fh "\t\t$bmodelstr2\t",&printFloatNice($birth{$bmodelstr}[0]),"\n";
	}elsif ($options{'rmodel'} eq "BD"){
            print $fh "\t\t$bmodelstr2\t",&printFloatNice($birth{$bmodelstr}[0]),"\t",&printFloatNice($death{$bmodelstr}[0]),"\n";
	}
    }
    return 1;
}
########################################################################
sub prob($$$$$){
    my ($b,$d,$i,$timeORbranch,$anc,$desc)=@_;
    return $Q{$timeORbranch}[$anc*($U+1)+$desc];
}
########################################################################
sub Anc{

    #joint ML ancestral reconstruction
    my (%L,%C)=();
    foreach my $node ($tree->get_leaf_nodes){ #tips
	my $parent=$node->ancestor;
	my $branch=$parent->internal_id."->".$node->internal_id;
	for (my $i=$min{$OG_analyzed};$i<=$max{$OG_analyzed};$i++){
	    $C{$node->internal_id}{$i}=$ancestral{$OG_analyzed}{$node->internal_id}[0];
	    $L{$node->internal_id}{$i}=&prob(1,1,1,$bmodel{$branch}."_".$distances{$parent->internal_id}{$node->internal_id},$i,$ancestral{$OG_analyzed}{$node->internal_id}[0]);
	}
    }

    for my $node (@treepath){ #internal_nodes
	next if ($node->{'parentid'} eq 'NULL');
	for (my $i=$min{$OG_analyzed};$i<=$max{$OG_analyzed};$i++){
	    my $maxL=-1;
	    my @desc=@{$node->{'desc'}};
	    my $parentid=$node->{'parentid'};
	    my $branch=$parentid."->".$node->{'id'};
	    for (my $j=$min{$OG_analyzed};$j<=$max{$OG_analyzed};$j++){
		my $L_tmp=$L{$desc[0][0]}{$j}*$L{$desc[1][0]}{$j}*&prob(1,1,1,$bmodel{$branch}."_".$distances{$parentid}{$node->{'id'}},$i,$j);
		if ($L_tmp>$maxL){
		    $L{$node->{'id'}}{$i}=$L_tmp;
		    $C{$node->{'id'}}{$i}=$j;
		    $maxL=$L_tmp;
		}
	    }
	}
    }

    for my $node (@treepath){#root
        if ($node->{'parentid'} eq 'NULL'){
	    my @desc=@{$node->{'desc'}};
	    my $maxL=-1;
	    for (my $j=$min{$OG_analyzed};$j<=$max{$OG_analyzed};$j++){
		my $L_tmp=$L{$desc[0][0]}{$j}*$L{$desc[1][0]}{$j}*&rootMLprob($OG_analyzed,$j);
		if ($L_tmp>$maxL){
                    $L{$node->{'id'}}{$j}=$L_tmp;
                    $maxL=$L_tmp;
		    $ancestral{$OG_analyzed}{$node->{'id'}}[0]=$j;
                }
	    }
	}
    }
    
    for my $node (reverse @treepath){
	if ($node->{'parentid'} ne 'NULL'){ #non-root
	    $ancestral{$OG_analyzed}{$node->{'id'}}[0]=$C{$node->{'id'}}{$ancestral{$OG_analyzed}{$node->{'parentid'}}[0]};
	}
    }

    return 1;
}
########################################################################
sub AncestralParsimony{
    my @numAnc=keys %difOGs;
    my $numAncSize=scalar(@numAnc);

    for (my $j=1; $j<=$numAncSize;$j++){
	$OG_analyzed=$difOGs{$numAnc[$j-1]}[0];
	my %dist=();
	foreach my $node ($tree->get_leaf_nodes){$dist{$node->internal_id}{$ancestral{$OG_analyzed}{$node->internal_id}[0]}=1;}
	for my $node (@treepath){
	    if ($options{'ep'} =~/CWP/){
		(!&Wagner(\@{$node->{'desc'}},\%dist, $node->{'id'})) && do{die;};
	    }else{
		(!&Sankoff(\@{$node->{'desc'}},\%dist, $node->{'id'})) && do{die;};
	    }
        }
	for my $node (reverse @treepath){
	    if ($node->{'parentid'} ne 'NULL'){
		$ancestral{$OG_analyzed}{$node->{'id'}}[0]=&MostProb(\%{$dist{$node->{'id'}}},$ancestral{$OG_analyzed}{$node->{'parentid'}}[0]);
	    }else{
		$ancestral{$OG_analyzed}{$node->{'id'}}[0]=&MostProb(\%{$dist{$node->{'id'}}},-1);
	    }
	}
    }
    
    foreach my $OG (@usedOG){ %{$ancestral{$OG}}=%{$ancestral{$difOGs{$OGstring{$OG}}[0]}};}
 
    if ($options{'unobs'}){
	my %dist=();
	$OG_analyzed="UnobservedData";
	foreach my $node ($tree->get_leaf_nodes){$dist{$node->internal_id}{$ancestral{$OG_analyzed}{$node->internal_id}[0]}=1;}
	for my $node (@treepath){
	    if ($options{'ep'} =~/CWP/){
		(!&Wagner(\@{$node->{'desc'}},\%dist, $node->{'id'})) && do{die;};
	    }else{
		(!&Sankoff(\@{$node->{'desc'}},\%dist, $node->{'id'})) && do{die;};
	    }
	}
	for my $node (reverse @treepath){
	    if ($node->{'parentid'} ne 'NULL'){
		$ancestral{$OG_analyzed}{$node->{'id'}}[0]=&MostProb(\%{$dist{$node->{'id'}}},$ancestral{$OG_analyzed}{$node->{'parentid'}}[0]);
	    }else{
		$ancestral{$OG_analyzed}{$node->{'id'}}[0]=&MostProb(\%{$dist{$node->{'id'}}},-1);
	    }
        }
    }

    return 1;
}
########################################################################
sub fdr($$){ #bejamini-hochberg
    my ($rawp,$adjp)=@_;
    my @tmp=@$rawp;
    @$rawp= sort {$a <=> $b} @$rawp;
    my @rawpcp=();
    foreach (@$rawp){push(@rawpcp,$_);}

    $$adjp{$$rawp[-1]}=$$rawp[-1];
    for (my $i=scalar(@$rawp)-2;$i>=1;$i--){
        my @min=($rawpcp[$i+1],(scalar(@$rawp)-1)/$i*$$rawp[$i],1);
        my $minv=&min(\@min);
        $rawpcp[$i]=$minv;
        $$adjp{$$rawp[$i]}=$rawpcp[$i];
    }
    @$rawp=@tmp;
    return 1;
}
########################################################################
sub printFloatNice($){
    my ($float)=@_;
    my $precision=7;
    if ($float eq "NA"){
	return "NA";
    }else{
	return sprintf("%.".$precision."f",$float);
    }
}
########################################################################
sub readmodel($){
    my ($model)=@_;

    %bmodel=();

    if (defined($model)){
	my $cnt=0;
	if (uc($model) eq "FR"){
	    foreach my $node (@nodes){
		my $nodeid=$node->internal_id;
		my $parentid=$node->ancestor->internal_id;
		my $branch=$parentid."->".$nodeid;
		$bmodel{$branch}=$cnt++;
	    }
	}elsif ($model ne ""){
	    my @groups=split(/_/,$model);
            my %branch_ctl=();
            foreach my $group (@groups){
                my @branches=split(/:/,$group);
                foreach my $branch (@branches){
                    $bmodel{$branch}=$cnt;
                    $branch_ctl{$branch}=0 unless ($branch_ctl{$branch});
                    $branch_ctl{$branch}++;
                    die "Wrong branch model specification: $branch repeated\n" if ($branch_ctl{$branch}>1);
                }
                $cnt++;
            }
            foreach my $node (@nodes){
                my $nodeid=$node->internal_id;
                my $parentid=$node->ancestor->internal_id;
                my $branch=$parentid."->".$nodeid;
                next if (grep(/^$branch$/,keys %branch_ctl));
                $bmodel{$branch}=$cnt;
            }
	}else{
	    foreach my $node (@nodes){
		my $nodeid=$node->internal_id;
		my $parentid=$node->ancestor->internal_id;
		my $branch=$parentid."->".$nodeid;
		$bmodel{$branch}=0;
	    }
	}
    }else{
	foreach my $node (@nodes){
	    my $nodeid=$node->internal_id;
	    my $parentid=$node->ancestor->internal_id;
	    my $branch=$parentid."->".$nodeid;
	    $bmodel{$branch}=0;
	}
    }
    
    return 1;
}
########################################################################
sub Sankoff($$$){
    my ($desc,$dist,$nodeid)=@_;
    my $OG=$OG_analyzed;

    my $max_changes=($max{$OG}-$min{$OG})/$distances{$$desc[0][0]}{$nodeid} + ($max{$OG}-$min{$OG})/$distances{$$desc[1][0]}{$nodeid};
    if (scalar(@{$$desc[2]})>0){
        $max_changes+=($max{$OG}-$min{$OG})/$distances{$$desc[2][0]}{$nodeid};
    }
    my $sum=0;

    for (my $i=$min{$OG};$i<=$max{$OG};$i++){
        my $tmp1=0;
        my $tmp2=0;
        my $tmp3=0;
        foreach my $k1 (keys %{$$dist{$$desc[0][0]}}){
            if ($k1!=$i){
                $tmp1+= $max_changes/(abs($k1-$i)/$distances{$$desc[0][0]}{$nodeid})*$$dist{$$desc[0][0]}{$k1};
            }else{
                $tmp1+=$MAX*$$dist{$$desc[0][0]}{$k1};
            }
        }
        
        foreach my $k2 (keys %{$$dist{$$desc[1][0]}}){
            if ($k2!=$i){
                $tmp2+= $max_changes/(abs($k2-$i)/$distances{$$desc[1][0]}{$nodeid})*$$dist{$$desc[1][0]}{$k2};
            }else{
                $tmp2+=$MAX*$$dist{$$desc[1][0]}{$k2};
            }
        }
        $$dist{$nodeid}{$i}=$tmp1*$tmp2;
        if (scalar(@{$$desc[2]})>0){
            my $out_weigth=1/scalar(@{$$desc[2]});
            foreach my $out (sort @{$$desc[2]}){
                foreach my $k3 (keys %{$$dist{$out}}){
                    if ($k3!=$i){
                        $tmp3+= $max_changes/(abs($k3-$i)/$distances{$out}{$nodeid})*$$dist{$out}{$k3}*$out_weigth;
                    }else{
                        $tmp3+=$MAX*$$dist{$out}{$k3}*$out_weigth;
                    }
                }
            }
            $$dist{$nodeid}{$i}*=$tmp3;
        }
        $sum+=$$dist{$nodeid}{$i};
    }
    for (my $i=$min{$OG};$i<=$max{$OG};$i++){
        $$dist{$nodeid}{$i}/=$sum;
    }
    
    return 1;

}
########################################################################
sub Wagner($$$){
    my ($desc,$dist,$nodeid)=@_;
    my @val1=keys %{$$dist{$$desc[0][0]}};
    my @val2=keys %{$$dist{$$desc[1][0]}};
    my @val=(@val1,@val2);
    my $min1=&min(\@val1);
    my $min2=&min(\@val2);
    my $max1=&max(\@val1);
    my $max2=&max(\@val2);
    my $max=$min1;

    if ($min1<$min2){
	$max=$min2;
    }
    my $min=$max1;
    if ($max1>$max2){
        $min=$max2;
    }
    my $start=$max;
    my $end=$min;
    if ($start>$end){
	$start=$min;
	$end=$max;
    }
    for (my $i=$start; $i<=$end;$i++){
	$$dist{$nodeid}{$i}=1;
    }
    return 1;

}
########################################################################
sub keysModel(){

    my %modelval=();

    foreach my $val (values %bmodel){
        $modelval{$val}=0 unless ($modelval{$val});
        $modelval{$val}++;
    }

    return keys %modelval;
}
########################################################################
sub getDifferentOGs(){
    #strategy to increase speed computation: some OGs has the same size distribution in taxa, so reconstructing one the rest are the same 
    foreach my $OG (@usedOG){
        my $string="";
        foreach my $taxa ($tree->get_leaf_nodes()){
            $string.=$ancestral{$OG}{$taxa->internal_id}[0]; #$fixed{$taxa->internal_id}{'v'}{$OG};
        }
        $OGstring{$OG}=$string;
        push(@{$difOGs{$string}},$OG);
    }
    return 1;
}
########################################################################
sub root_dist($$){
    my ($xb,$sm)=@_;
    my $cnt=scalar(@$xb);
    %OGOrder=();
    if ($options{'root_dist'} == 1){
	if ($options{'start_val'}==0){
	    my $sum=0;
	    foreach my $OG (@usedOG){
		$sum+=$ancestral{$OG}{$root_node_id}[0];
	    }
	    push(@$xb,$sum/scalar(@usedOG));
	    push(@$sm,$sum/scalar(@usedOG)/2+$MIN);
	}else{
	    push(@$xb,rand());
            push(@$sm,rand());
	}
        $OGOrder{'ALL'}=$cnt;
    }elsif ($options{'root_dist'} == 2){
        foreach my $string (keys %difOGs){
	    my $OG=$difOGs{$string}[0];
	    $OG_analyzed=$OG;
	    if ($options{'start_val'}==0){
		push(@$xb,$ancestral{$OG_analyzed}{$root_node_id}[0]);
		push(@$sm,$ancestral{$OG_analyzed}{$root_node_id}[0]/2+$MIN);
	    }else{
		push(@$xb,rand());
                push(@$sm,rand());
	    }
            $OGOrder{$OG_analyzed}=$cnt++;
        }
    }elsif ($options{'root_dist'} == 3){
        if ($options{'start_val'} == 0){
            my $sum=0;
            foreach my $OG (@usedOG){
                $sum+=$ancestral{$OG}{$root_node_id}[0];
            }
            my $p=0.5; # p = (variance - mean)/mean
            my $n=$sum/scalar(@usedOG); #n = mean**2/(variance-mean)

            push(@$xb,$n,$p);
            push(@$sm,$n/2+$MIN,$p/2+$MIN);
        }else{
            push(@$xb,rand(),rand());
            push(@$sm,rand(),rand());
        }
        $OGOrder{'ALL'}=$cnt;
    }elsif ($options{'root_dist'} == 4){
        foreach my $string (keys %difOGs){
            my $OG=$difOGs{$string}[0];
            $OG_analyzed=$OG;
            if ($options{'start_val'} == 0){
                my $p=0.5;
                my $n=$ancestral{$OG_analyzed}{$root_node_id}[0];
		push(@$xb,$n,$p);
		push(@$sm,$n/2+$MIN,$p/2+$MIN);
	    }else{
		push(@$xb,rand(),rand());
		push(@$sm,rand(),rand());
	    }
	    $OGOrder{$OG_analyzed}=$cnt++;
	}
    }

    return 1;
}
########################################################################
sub RatesParamBDI($$){
    my ($xb,$sm)=@_;
    %RatesOrder=();
    %RatesOrder_rev=();
    %OGOrder=();
    my $cnt=0;
    my @rates=("birth","death","innovation");
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
    	if ($options{'start_val'}==0){
	    for my $rate (@rates){
		if ($options{'ep'}=~/MAP/ && defined($priormodel{$bmodelstr}{$rate}{'shape'})){			
                        push(@$xb,$priormodel{$bmodelstr}{$rate}{'shape'});
			if ($priormodel{$bmodelstr}{$rate}{'shape'}>0){
			    push(@$sm,$priormodel{$bmodelstr}{$rate}{'shape'}/2 + $MIN);
			}else{
			    push(@$sm,0);
			}
		}else{
		    if ($rate eq "birth"){
			push(@$xb,$birth{$bmodelstr}[0]);
			push(@$sm,$birth{$bmodelstr}[0]/2 + $MIN);
		    }elsif ($rate eq "death"){
			push(@$xb,$death{$bmodelstr}[0]);
                        push(@$sm,$death{$bmodelstr}[0]/2 + $MIN);
		    }elsif ($rate eq "innovation"){
			push(@$xb,$innovation{$bmodelstr}[0]);
                        push(@$sm,$innovation{$bmodelstr}[0]/2 + $MIN);
		    }
		}
	    }
    	}else{
    	    push(@$xb,rand(),rand(),rand());
    	    push(@$sm,rand(),rand(),rand());
    	}
    	$RatesOrder{$bmodelstr}=$cnt++;
    }
    %RatesOrder_rev=reverse %RatesOrder;
    (!&root_dist($xb,$sm)) && do{die;};

    return 1;
}
########################################################################                                                                                                    
sub RatesParamBD($$){
    my ($xb,$sm)=@_;
    %RatesOrder=();
    %RatesOrder_rev=();
    %OGOrder=();
    my $cnt=0;
    my @rates=("birth","death","innovation");
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
        if ($options{'start_val'}==0){
	    for my $rate (@rates){
		if ($rate eq "innovation"){
			push(@$xb,0);
                        push(@$sm,0);	
		}elsif ($options{'ep'}=~/MAP/ && defined($priormodel{$bmodelstr}{$rate}{'shape'})){			
                        push(@$xb,$priormodel{$bmodelstr}{$rate}{'shape'});
			if ($priormodel{$bmodelstr}{$rate}{'shape'}>0){
                            push(@$sm,$priormodel{$bmodelstr}{$rate}{'shape'}/2 + $MIN);
                        }else{
                            push(@$sm,0);
                        }
		}else{
		    if ($rate eq "birth"){
			push(@$xb,$birth{$bmodelstr}[0]);
			push(@$sm,$birth{$bmodelstr}[0]/2 + $MIN);
		    }elsif ($rate eq "death"){
			push(@$xb,$death{$bmodelstr}[0]);
                        push(@$sm,$death{$bmodelstr}[0]/2 + $MIN);
		    }

		}
            }
        }else{
            push(@$xb,rand(),rand(),0);
            push(@$sm,rand(),rand(),0);
        }
        $RatesOrder{$bmodelstr}=$cnt++;
    }
    %RatesOrder_rev=reverse %RatesOrder;
    (!&root_dist($xb,$sm)) && do{die;};

    return 1;
}
########################################################################
sub RatesParamLI($$){
    my ($xb,$sm)=@_;
    %RatesOrder=();
    %RatesOrder_rev=();
    %OGOrder=();
     my $cnt=0;
    my @rates=("lambda","innovation");	 
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	if ($options{'start_val'}==0){
	    for my $rate (@rates){
		if ($options{'ep'}=~/MAP/ && defined($priormodel{$bmodelstr}{$rate}{'shape'})){			
                        push(@$xb,$priormodel{$bmodelstr}{$rate}{'shape'});
			if ($priormodel{$bmodelstr}{$rate}{'shape'}>0){
                            push(@$sm,$priormodel{$bmodelstr}{$rate}{'shape'}/2 + $MIN);
                        }else{
                            push(@$sm,0);
                        }
		}else{
		    if ($rate eq "lambda"){
			push(@$xb,$birth{$bmodelstr}[0]);
			push(@$sm,$birth{$bmodelstr}[0]/2 + $MIN);
		    }elsif ($rate eq "innovation"){
			push(@$xb,$innovation{$bmodelstr}[0]);
                        push(@$sm,$innovation{$bmodelstr}[0]/2 + $MIN);
		    }
		}
	    }
	}else{
	    push(@$xb,rand(),rand());
            push(@$sm,rand(),rand());
	}
	$RatesOrder{$bmodelstr}=$cnt++;
    }
    %RatesOrder_rev=reverse %RatesOrder;
    (!&root_dist($xb,$sm)) && do{die;};

    return 1;
}
########################################################################                                                                                                    
sub RatesParamL($$){
    my ($xb,$sm)=@_;
    %RatesOrder=();
    %RatesOrder_rev=();
    %OGOrder=();
    my $cnt=0;
    my @rates=("lambda","innovation");
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
        if ($options{'start_val'}==0){
	    for my $rate (@rates){
		if ($rate eq "innovation"){
			push(@$xb,0);
                        push(@$sm,0) if ($sm);
		}elsif ($options{'ep'}=~/MAP/ && defined($priormodel{$bmodelstr}{$rate}{'shape'})){			
                        push(@$xb,$priormodel{$bmodelstr}{$rate}{'shape'});
			if ($priormodel{$bmodelstr}{$rate}{'shape'}>0){
                            push(@$sm,$priormodel{$bmodelstr}{$rate}{'shape'}/2 + $MIN);
                        }else{
                            push(@$sm,0);
                        }
		}else{
		    if ($rate eq "lambda"){
			push(@$xb,$birth{$bmodelstr}[0]);
			push(@$sm,$birth{$bmodelstr}[0]/2 + $MIN) if ($sm);
                    }
		}
	    }
        }else{
            push(@$xb,rand(),0);
            push(@$sm,rand(),0);
        }
        $RatesOrder{$bmodelstr}=$cnt++;
    }
    %RatesOrder_rev=reverse %RatesOrder;
    (!&root_dist($xb,$sm)) && do{die;};

    return 1;
}
########################################################################                                                                
sub RatesParamGD($$){
    my ($xb,$sm)=@_;
    %RatesOrder=();
    %OGOrder=();
    %RatesOrder_rev=();
    my $cnt=0;
    my @rates=("gain","death");
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	if ($options{'start_val'}==0){
           for my $rate (@rates){	
		if ($options{'ep'}=~/MAP/ && defined($priormodel{$bmodelstr}{$rate}{'shape'})){			
                        push(@$xb,$priormodel{$bmodelstr}{$rate}{'shape'});
			if ($priormodel{$bmodelstr}{$rate}{'shape'}>0){
                            push(@$sm,$priormodel{$bmodelstr}{$rate}{'shape'}/2 + $MIN);
                        }else{
                            push(@$sm,0);
                        }
		}else{
		    if ($rate eq "gain"){
			push(@$xb,$death{$bmodelstr}[0]);
			push(@$sm,$death{$bmodelstr}[0]/2 + $MIN) if ($sm);
		    }elsif ($rate eq "death"){
			push(@$xb,$innovation{$bmodelstr}[0]);
                        push(@$sm,$innovation{$bmodelstr}[0]/2 + $MIN) if ($sm);
		    }
		}
            }
	}else{
	    push(@$xb,rand(),rand());
            push(@$sm,rand(),rand());
	}
        $RatesOrder{$bmodelstr}=$cnt++;
    }
    %RatesOrder_rev=reverse %RatesOrder;
    (!&root_dist($xb,$sm)) && do{die;};

    return 1;
}
########################################################################
sub usedOGs{
   #ortholog groups are only used if they contain at least one exatnt taxa
     OG: foreach my $OG (keys %ortho){
         foreach my $pair (@{$ortho{$OG}{'ID'}}){
	     if (defined($assoc{$pair->{'sp'}})){
		 if ($options{'unobs'}){
		     if ($ancestral{$OG}{$assoc{$pair->{'sp'}}}[0]>0){
			 push(@usedOG,$OG);
			 next OG;
		     }
		 }else{
		     push(@usedOG,$OG);
		     next OG;
		 }
	     }
	 }
     }
    return 1;
}
########################################################################
sub sumAnc(){
    %anc_cnt=();
    foreach my $OG (@usedOG){
        foreach my $node ($tree->get_nodes){
            $anc_cnt{$node->internal_id}=0 unless ($anc_cnt{$node->internal_id});
            $anc_cnt{$node->internal_id}+=$ancestral{$OG}{$node->internal_id}[0];
        }
    }
    return 1;
}
########################################################################
sub death2loss($){ #conversion from death rate to loss rate (approx)
    my ($death)=@_;
    my %loss=();
    my %cnt=();
    my %sumanc=();

    foreach my $node (@nodes){
	my $parentid=$node->ancestor->internal_id;
        my $nodeid=$node->internal_id;
        my $branch=$parentid."->".$nodeid;
	$sumanc{$bmodel{$branch}} = 0 unless ($sumanc{$bmodel{$branch}});
	$sumanc{$bmodel{$branch}}+=$anc_cnt{$parentid};
	foreach my $OG (@usedOG){
	    $cnt{$bmodel{$branch}} = 0 unless ($cnt{$bmodel{$branch}});
	    $cnt{$bmodel{$branch}}++ if ($ancestral{$OG}{$parentid}[0]>0);
	}
    }

    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	if ($cnt{$bmodelstr} > 0){
	    $loss{$bmodelstr}[0]=$$death{$bmodelstr}[0]*$sumanc{$bmodelstr}/$cnt{$bmodelstr};
	}else{
	    $loss{$bmodelstr}[0]="NA";
	}
    }    
    return %loss;
}
########################################################################
sub calcGainLoss($$$){
    my ($parentid,$nodeid,$mode)=@_;
    my $gain=0;
    my $loss=0;
    my $branch=$parentid."->".$nodeid;
    if ($mode == 0){
	foreach my $OG (@usedOG){
	    if ($ancestral{$OG}{$parentid}[0]<$ancestral{$OG}{$nodeid}[0]){
		$gain+=abs($ancestral{$OG}{$nodeid}[0]-$ancestral{$OG}{$parentid}[0]);
	    }else{
		$loss+=abs($ancestral{$OG}{$parentid}[0]-$ancestral{$OG}{$nodeid}[0]);
	    }
	}
    }
    return ($gain,$loss);
}
########################################################################
sub calcInno($$){
    my ($parentid,$nodeid)=@_;
    my $inno=0;
    if ($options{'family'}){
	if ($anc_cnt{$parentid}==0 && $anc_cnt{$nodeid}>0){
	    $inno=1;
	}
    }else{
	foreach my $OG (@usedOG){
	    if($ancestral{$OG}{$parentid}[0] == 0 && $ancestral{$OG}{$nodeid}[0] > 0){
		$inno++;
	    }
	}
    }
    return $inno;
}
########################################################################
sub FamilyObservedRates($$){

    my ($fh,$mode)=@_;

    my (%distanceA,%distanceB)=();
    my (%birthtmp,%deathtmp,%innovationtmp,%gaintmp, %lambdatmp)=();
    my $family=$options{'family'};
    $options{'family'}=0 if ($mode == 0);

    foreach my $node (@nodes){

        my ($birth,$death,$inno,$gain,$loss,$lambda, $cnt);
        $birth=$death=$inno=$gain=$loss=$lambda=$cnt=0;

        my $parentid=$node->ancestor->internal_id;
        my $nodeid=$node->internal_id;
        my $branch=$parentid."->".$nodeid;

	$inno=&calcInno($parentid,$nodeid);
	($gain,$death)=&calcGainLoss($parentid,$nodeid,0);
	$birth=$gain;

	if ($anc_cnt{$parentid}+$inno>0){
	    my $tmp=$birth-$inno;
	    if ($tmp < 0){
		die "Error $branch has $gain gains and $inno innovations which is impossible\n";
	    }
	    $birthtmp{$bmodel{$branch}}[0]+=($tmp/($anc_cnt{$parentid}+$inno));
	    $deathtmp{$bmodel{$branch}}[0]+=($death/($anc_cnt{$parentid}+$inno));
	    $lambdatmp{$bmodel{$branch}}[0]+=(($death+$tmp)/($anc_cnt{$parentid}+$inno));
	    $distanceA{$bmodel{$branch}}=0 unless ($distanceA{$bmodel{$branch}});
	    $distanceA{$bmodel{$branch}}+=$distances{$parentid}{$nodeid};
	}
	if (!$options{'family'}){
	    $gaintmp{$bmodel{$branch}}[0]+=$gain/scalar(@usedOG);
	    $innovationtmp{$bmodel{$branch}}[0]+=$inno/scalar(@usedOG);
	}else{
	    $gaintmp{$bmodel{$branch}}[0]+=$gain;
	    $innovationtmp{$bmodel{$branch}}[0]+=$inno;
	}
	$distanceB{$bmodel{$branch}}=0 unless ($distanceB{$bmodel{$branch}});
	$distanceB{$bmodel{$branch}}+=$distances{$parentid}{$nodeid};
    }

    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	$innovationtmp{$bmodelstr}[0]=0 if (!defined($innovationtmp{$bmodelstr}[0]));
	$birthtmp{$bmodelstr}[0]=0 if (!defined($birthtmp{$bmodelstr}[0]));
	$deathtmp{$bmodelstr}[0]=0 if (!defined($deathtmp{$bmodelstr}[0]));
	$gaintmp{$bmodelstr}[0]=0 if (!defined($gaintmp{$bmodelstr}[0]));
	$lambdatmp{$bmodelstr}[0]=0 if (!defined($lambdatmp{$bmodelstr}[0]));
	$distanceA{$bmodelstr}=0 unless ($distanceA{$bmodelstr});

	if ($distanceA{$bmodelstr}>0){
	    $lambdatmp{$bmodelstr}[0]/=($distanceA{$bmodelstr}*2);
	    $birthtmp{$bmodelstr}[0]/=$distanceA{$bmodelstr};
	    $deathtmp{$bmodelstr}[0]/=$distanceA{$bmodelstr};
	}else{
	    $lambdatmp{$bmodelstr}[0]="NA";
            $birthtmp{$bmodelstr}[0]="NA";
            $deathtmp{$bmodelstr}[0]="NA";
	}
	$gaintmp{$bmodelstr}[0]/=$distanceB{$bmodelstr};
	$innovationtmp{$bmodelstr}[0]/=$distanceB{$bmodelstr};
    }

    foreach my $bmodelstr (sort {$a <=> $b} @keys){
	if ($birthtmp{$bmodelstr}[0] eq "NA"){
	    $birth{$bmodelstr}[0]=0;
	}else{
	    $birth{$bmodelstr}[0]=$birthtmp{$bmodelstr}[0];
	}
	if ($deathtmp{$bmodelstr}[0] eq "NA"){
	    $death{$bmodelstr}[0]=0;
	}else{
	    $death{$bmodelstr}[0]=$deathtmp{$bmodelstr}[0];
	}
	if ($options{'rmodel'} eq "BDI" || $options{'rmodel'} eq "LI"){
	    $innovation{$bmodelstr}[0]=$innovationtmp{$bmodelstr}[0];
	}else{
	    $innovation{$bmodelstr}[0]=$gaintmp{$bmodelstr}[0];
	}
	if ($options{'rmodel'} eq "LI" || $options{'rmodel'} eq "L"){
	    if ($lambdatmp{$bmodelstr}[0] eq "NA"){
		$birth{$bmodelstr}[0]=0;
	    }else{
		$birth{$bmodelstr}[0]=$lambdatmp{$bmodelstr}[0];
	    }
	}
    }

    if ($mode == 1 ){
	(!&printRates(\*$fh)) && do{die;};
    }elsif ($mode == 2){
	(!&printRates(\*$fh,$mode)) && do{die;};
    }
    $options{'family'}=$family;
    return 1;
}
########################################################################
sub fixedNode{
#fixes leaf nodes 

    my (@idnodes)=();
    foreach my $node ($tree->get_nodes){push(@idnodes,$node->internal_id);}
    
    foreach my $OG (keys %ortho){
	my %tmp=();
        foreach my $node ($tree->get_nodes){
            if ($node->is_Leaf()){
                $fixed{$node->internal_id}{'b'}=1; #all the leaves are fixed
		$ancestral{$OG}{$node->internal_id}[0]=0;
		$tmp{$node->internal_id}=0;
            }else{
                $fixed{$node->internal_id}{'b'}=0;
            }
        }

        foreach my $pair (@{$ortho{$OG}{'ID'}}){ # It gets all species specified in orthotable present as tree leaves
            next unless ($assoc{$pair->{'sp'}}); #only the leaves
	    $ancestral{$OG}{$assoc{$pair->{'sp'}}}[0]=$pair->{'cnt'};
	    $tmp{$assoc{$pair->{'sp'}}}=$pair->{'cnt'};
	}
        
	my @values=values %tmp;
	$min{$OG}=0; 
	$max{$OG}=int(&max(\@values)*2)+$options{'n_max_int'};
    }
    return 1;
}
########################################################################
sub MinMax($$$){
    my ($nodeid,$desc,$likelihood)=@_;
    
    my $OG=$OG_analyzed;
    if ($fixed{$$desc[0][0]}{'b'}==1){
	$$likelihood{$$desc[0][0]}{$ancestral{$OG}{$$desc[0][0]}[0]}=1;
    }
    if ($fixed{$$desc[1][0]}{'b'}==1){
	$$likelihood{$$desc[1][0]}{$ancestral{$OG}{$$desc[1][0]}[0]}=1; 
    }
    return ($min{$OG},$max{$OG});
}
########################################################################                                                
sub Felsenstein($){
    my ($likelihood)=@_;

    foreach my $node (@treepath){
        my $nodeid=$node->{'id'};
	my @desc=@{$node->{'desc'}};
        my ($min,$max)=&MinMax($nodeid,\@desc,$likelihood);
        my $branch=$nodeid."->".$desc[0][0];
        my $branch1=$nodeid."->".$desc[1][0];
	my $key1=$bmodel{$branch}."_".$distances{$nodeid}{$desc[0][0]};
	my $key2=$bmodel{$branch1}."_".$distances{$nodeid}{$desc[1][0]};
	my $hash_ref=likBranch($U,$min,$max,\@{$Q{$key1}},\@{$Q{$key2}},\%{$$likelihood{$desc[0][0]}},\%{$$likelihood{$desc[1][0]}});
	%{$$likelihood{$nodeid}}=%$hash_ref;
    }
    return 1;
}
########################################################################
sub root_dist_p($$$){
    my ($i,$OG,$val)=@_;
    if ($options{'root_dist'} == 0){
	return &poisson($i,$ancestral{$OG}{$root_node_id}[0]);
    }elsif ($options{'root_dist'} == 1){
	return &poisson($i,$$val[$OGOrder{'ALL'}]); 
    }elsif ($options{'root_dist'} == 2){
	return &poisson($i,$$val[$OGOrder{$OG}]);
    }elsif ($options{'root_dist'} == 3){
	return &negbin($i,$$val[$OGOrder{'ALL'}],$$val[$OGOrder{'ALL'}+1]);
    }elsif ($options{'root_dist'} == 4){
	return &negbin($i,$$val[$OGOrder{$OG}],$$val[$OGOrder{$OG}+1]);
    }
}
########################################################################
sub paramChange($$){
    my ($i,$param)=@_;
    if ($options{'rmodel'} eq 'BDI' || $options{'rmodel'} eq 'BD'){
	if (($$param[$i*3] != $ant_param[$i*3]) || ($$param[$i*3+1] != $ant_param[$i*3+1]) || ($$param[$i*3+2] != $ant_param[$i*3+2])){
	    return 1;
	}
    }elsif ($options{'rmodel'} eq "GD" || $options{'rmodel'} eq 'LI' || $options{'rmodel'} eq 'L'){
	if (($$param[$i*2] != $ant_param[$i*2]) || ($$param[$i*2+1] != $ant_param[$i*2+1])){
            return 1;
        }
    }
    return 0;
}
########################################################################
sub Q2do(){
    my @exp=();
    
    if (scalar(@ant_param)>0){
	my @done=();
        for (my $i=0;$i<scalar(values %RatesOrder);$i++){
	    if (&paramChange($i,\@_)){
                my $bmodel_a=$RatesOrder_rev{$i};
                foreach my $node (@treepath){
                    my $nodeid=$node->{'id'};
                    my @desc=@{$node->{'desc'}};
                    my $descid1=$desc[0][0];
                    my $descid2=$desc[1][0];
                    my $branch1=$nodeid."->".$descid1;
                    my $branch2=$nodeid."->".$descid2;
		    my $key1=$bmodel{$branch1}."_".$distances{$descid1}{$nodeid};
		    my $key2=$bmodel{$branch2}."_".$distances{$descid2}{$nodeid};
                    if ($bmodel{$branch1} eq $bmodel_a && !grep(/^$key1$/,@done)){
			push(@exp,{'bmodel'=>$bmodel{$branch1}, 'distance'=>$distances{$descid1}{$nodeid}});
			push(@done,$key1);
                    }
                    if ($bmodel{$branch2} eq $bmodel_a && !grep(/^$key2$/,@done)){
			push(@exp,{'bmodel'=>$bmodel{$branch2}, 'distance'=>$distances{$descid2}{$nodeid}});
			push(@done,$key2);
                    }
                }
            }
        }
    }else{
	my @done=();
        foreach my $node (@treepath){
            my $nodeid=$node->{'id'};
            my @desc=@{$node->{'desc'}};
            my $descid1=$desc[0][0];
            my $descid2=$desc[1][0];
            my $branch1=$nodeid."->".$descid1;
            my $branch2=$nodeid."->".$descid2;
	    my $key1=$bmodel{$branch1}."_".$distances{$descid1}{$nodeid};
	    my $key2=$bmodel{$branch2}."_".$distances{$descid2}{$nodeid};
	    if (!grep(/^$key1$/,@done)){
		push(@exp,{'bmodel'=>$bmodel{$branch1}, 'distance'=>$distances{$descid1}{$nodeid}});
		push(@done,$key1);
            }
            if (!grep(/^$key2$/,@done)){
		push(@exp,{'bmodel'=>$bmodel{$branch2}, 'distance'=>$distances{$descid2}{$nodeid}});
		push(@done,$key2);
	    }
        }
    }
    @ant_param=@_;
    return @exp;
}
########################################################################
sub GDRatesLik{
    my $cnt=1;
    foreach  (@_){
	return $MAX if ($_<0);
        $cnt++;
    }

    if ($options{'root_dist'} == 2){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}]>$max{$OG});
        }
    }elsif ($options{'root_dist'} == 3){
        return $MAX if ($_[$OGOrder{'ALL'}+1]>1);
    }elsif ($options{'root_dist'} == 4){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}+1]>1);
        }
    }

    my @numAnc=keys %difOGs;
    my $numAncSize=scalar(@numAnc);
    my $likelihood=0;
    my $unobslik=0;
    my %likelihood=();
    my @exp=&Q2do(@_);
    foreach my $element (@exp){
	my $key=$element->{'bmodel'}."_".$element->{'distance'};
	my $array_ref=&c_prob($U,0,$_[$RatesOrder{$element->{'bmodel'}}*2],$_[$RatesOrder{$element->{'bmodel'}}*2+1],$element->{'distance'});
	@{$Q{$key}}=@$array_ref;
	
    }

    if ($options{'unobs'}){
	$OG_analyzed="UnobservedData";
	(!&Felsenstein(\%likelihood)) && do{die;};
	for (my $i=$min{$OG_analyzed};$i<=$max{$OG_analyzed};$i++){
	    $unobslik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG_analyzed,\@_);;
	}
	return $MAX if ($unobslik >= 1 || $unobslik eq "nan");
    }

    foreach my $string (keys %difOGs){
	my $OG=$difOGs{$string}[0];
        $OG_analyzed=$OG;
        %likelihood=();
        my $lik=0;
        (!&Felsenstein(\%likelihood)) && do{die;};
        for (my $i=$min{$OG};$i<=$max{$OG};$i++){
	    $lik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG,\@_);
        }
        return $MAX if ($lik == 0 || $lik eq "nan");

        if ($options{'unobs'}){
            $likelihood+=log($lik/(1-$unobslik))*scalar(@{$difOGs{$string}});
        }else{
            $likelihood+=log($lik)*scalar(@{$difOGs{$string}});
        }
    }

    my $tot_prior=0;
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
        my $tmp=&prior($_[$RatesOrder{$bmodelstr}*2],\%{$priormodel{$bmodelstr}{'gain'}});
        my $tmp1=&prior($_[$RatesOrder{$bmodelstr}*2+1],\%{$priormodel{$bmodelstr}{'death'}});
        if (($tmp == 0 || $tmp1 == 0) || ($tmp eq "-inf" || $tmp1 eq "-inf") || ($tmp > 1 || $tmp1 > 1)){
            return $MAX;
        }
        $tot_prior+=log($tmp)+log($tmp1);
    }
    return abs($likelihood+$tot_prior);
}
########################################################################                                                          
sub BDRatesLik(){
    #computes the likelihood of the tree given BDI rates                                                                          
    #constrains                                                                                                                    
    foreach  (@_){
        return $MAX if ($_<0);
    }

    if ($options{'root_dist'} == 2){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}]>$max{$OG});
        }
    }elsif ($options{'root_dist'} == 3){
        return $MAX if ($_[$OGOrder{'ALL'}+1]>1);
    }elsif ($options{'root_dist'} == 4){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}+1]>1);
        }
    }
    my @numAnc=keys %difOGs;
    my $numAncSize=scalar(@numAnc);
    my $likelihood=0;
    my $unobslik=0;
    my %likelihood=();
    my @exp=&Q2do(@_);
    foreach my $element (@exp){
        my $key=$element->{'bmodel'}."_".$element->{'distance'};
        my $array_ref=&c_prob($U,$_[$RatesOrder{$element->{'bmodel'}}*3],$_[$RatesOrder{$element->{'bmodel'}}*3+1],$_[$RatesOrder{$element->{'bmodel'}}*3+2],$element->{'distance'});
	@{$Q{$key}}=@$array_ref;
    }

    if ($options{'unobs'}){
        $OG_analyzed="UnobservedData";
        (!&Felsenstein(\%likelihood)) && do{die;};
        for (my $i=$min{$OG_analyzed};$i<=$max{$OG_analyzed};$i++){
            $unobslik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG_analyzed,\@_);
        }
        return $MAX if ($unobslik >= 1 || $unobslik eq "nan");
    }

    foreach my $string (keys %difOGs){
        my $OG=$difOGs{$string}[0];
        $OG_analyzed=$OG;
        %likelihood=();
        my $lik=0;
        (!&Felsenstein(\%likelihood)) && do{die;};

        for (my $i=$min{$OG};$i<=$max{$OG};$i++){
            $lik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG,\@_);
        }
        return $MAX if ($lik == 0 || $lik eq "nan");
        if ($options{'unobs'}){
            $likelihood+=log($lik/(1-$unobslik))*scalar(@{$difOGs{$string}});
        }else{
            $likelihood+=log($lik)*scalar(@{$difOGs{$string}});
        }
    }

    my $tot_prior=0;
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
        my $tmp=&prior($_[$RatesOrder{$bmodelstr}*3],\%{$priormodel{$bmodelstr}{'birth'}});
        my $tmp1=&prior($_[$RatesOrder{$bmodelstr}*3+1],\%{$priormodel{$bmodelstr}{'death'}});
         if (($tmp == 0 || $tmp1 == 0) || ($tmp eq "-inf" || $tmp1 eq "-inf") || ($tmp > 1 || $tmp1 > 1)){
            return $MAX;
        }
        $tot_prior+=log($tmp)+log($tmp1);
    }
    return abs($likelihood+$tot_prior);
}
########################################################################
sub BDIRatesLik(){
    #computes the likelihood of the tree given BDI rates
    #constrains
    foreach  (@_){
	return $MAX if ($_<0);
    }

    if ($options{'root_dist'} == 2){
	foreach my $OG (keys %OGOrder){
	    return $MAX if ($_[$OGOrder{$OG}]>$max{$OG});
	}
    }elsif ($options{'root_dist'} == 3){
        return $MAX if ($_[$OGOrder{'ALL'}+1]>1);
    }elsif ($options{'root_dist'} == 4){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}+1]>1);
        }
    }

    my @numAnc=keys %difOGs;
    my $numAncSize=scalar(@numAnc);
    my $likelihood=0;
    my $unobslik=0;
    my %likelihood=();
    my @exp=&Q2do(@_);	
    foreach my $element (@exp){
	my $key=$element->{'bmodel'}."_".$element->{'distance'};
	my $array_ref=&c_prob($U,$_[$RatesOrder{$element->{'bmodel'}}*3],$_[$RatesOrder{$element->{'bmodel'}}*3+1],$_[$RatesOrder{$element->{'bmodel'}}*3+2],$element->{'distance'});
	@{$Q{$key}}=@$array_ref;	
    }

    if ($options{'unobs'}){
	$OG_analyzed="UnobservedData";
	(!&Felsenstein(\%likelihood)) && do{die;};
	for (my $i=$min{$OG_analyzed};$i<=$max{$OG_analyzed};$i++){
            $unobslik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG_analyzed,\@_);
        }
        return $MAX if ($unobslik >= 1 || $unobslik eq "nan");
    }

    foreach my $string (keys %difOGs){
	my $OG=$difOGs{$string}[0];
        $OG_analyzed=$OG;
	%likelihood=();
	my $lik=0;
	(!&Felsenstein(\%likelihood)) && do{die;};

	for (my $i=$min{$OG};$i<=$max{$OG};$i++){
	    $lik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG,\@_);
	}
	return $MAX if ($lik == 0 || $lik eq "nan");
	if ($options{'unobs'}){
	    $likelihood+=log($lik/(1-$unobslik))*scalar(@{$difOGs{$string}});
	}else{
	    $likelihood+=log($lik)*scalar(@{$difOGs{$string}});
	}
    }

    my $tot_prior=0;
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	my $tmp=&prior($_[$RatesOrder{$bmodelstr}*3],\%{$priormodel{$bmodelstr}{'birth'}});
	my $tmp1=&prior($_[$RatesOrder{$bmodelstr}*3+1],\%{$priormodel{$bmodelstr}{'death'}});	
	my $tmp2=&prior($_[$RatesOrder{$bmodelstr}*3+2],\%{$priormodel{$bmodelstr}{'innovation'}});
	if (($tmp == 0 || $tmp1 == 0 || $tmp2 == 0) || ($tmp eq "-inf" || $tmp1 eq "-inf" || $tmp2 eq "-inf") || ($tmp > 1 || $tmp1 > 1 || $tmp2 > 1)){
	    return $MAX;
	}
	$tot_prior+=log($tmp)+log($tmp1)+log($tmp2);
    }

    return abs($likelihood+$tot_prior);
}
########################################################################
sub LIRatesLik(){
    #computes the likelihood of the tree given LI rates
    #constrains
    foreach  (@_){
	return $MAX if ($_<0);
    }

    if ($options{'root_dist'} == 2){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}]>$max{$OG});
        }
    }elsif ($options{'root_dist'} == 3){
        return $MAX if ($_[$OGOrder{'ALL'}+1]>1);
    }elsif ($options{'root_dist'} == 4){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}+1]>1);
        }
    }

    my @numAnc=keys %difOGs;
    my $numAncSize=scalar(@numAnc);
    my $likelihood=0;
    my $unobslik=0;
    my %likelihood=();
    my @exp=&Q2do(@_);
    foreach my $element (@exp){
	my $key=$element->{'bmodel'}."_".$element->{'distance'};
	my $array_ref=&c_prob($U,$_[$RatesOrder{$element->{'bmodel'}}*2],$_[$RatesOrder{$element->{'bmodel'}}*2],$_[$RatesOrder{$element->{'bmodel'}}*2+1],$element->{'distance'});
	@{$Q{$key}}=@$array_ref;
	
    }

    if ($options{'unobs'}){
	$OG_analyzed="UnobservedData";
	(!&Felsenstein(\%likelihood)) && do{die;};
	for (my $i=$min{$OG_analyzed};$i<=$max{$OG_analyzed};$i++){
            $unobslik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG_analyzed,\@_);
	}
        return $MAX if ($unobslik >= 1 || $unobslik eq "nan");
    }

    foreach my $string (keys %difOGs){
	my $OG=$difOGs{$string}[0];
        $OG_analyzed=$OG;
        %likelihood=();
        my $lik=0;
        (!&Felsenstein(\%likelihood)) && do{die;};

        for (my $i=$min{$OG};$i<=$max{$OG};$i++){
            $lik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG,\@_); 
        }
        return $MAX if ($lik == 0 || $lik eq "nan");
	if ($options{'unobs'}){
            $likelihood+=log($lik/(1-$unobslik))*scalar(@{$difOGs{$string}});
        }else{
            $likelihood+=log($lik)*scalar(@{$difOGs{$string}});
        }
    }

    my $tot_prior=0;
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
        my $tmp=&prior($_[$RatesOrder{$bmodelstr}*2],\%{$priormodel{$bmodelstr}{'lambda'}});
        my $tmp1=&prior($_[$RatesOrder{$bmodelstr}*2+1],\%{$priormodel{$bmodelstr}{'innovation'}});
        if (($tmp == 0 || $tmp1 == 0) || ($tmp eq "-inf" || $tmp1 eq "-inf") || ($tmp > 1 || $tmp1 > 1)){
            return $MAX;
        }
        $tot_prior+=log($tmp)+log($tmp1);
    }

    return abs($likelihood+$tot_prior);
}
########################################################################                                                                                                    
sub LRatesLik(){
    #computes the likelihood of the tree given LI rates                                                                                                                     
    #constrains                                                                                                                                                             
    foreach  (@_){
        return $MAX if ($_<0);
    }

    if ($options{'root_dist'} == 2){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}]>$max{$OG});
        }
    }elsif ($options{'root_dist'} == 3){
        return $MAX if ($_[$OGOrder{'ALL'}+1]>1);
    }elsif ($options{'root_dist'} == 4){
        foreach my $OG (keys %OGOrder){
            return $MAX if ($_[$OGOrder{$OG}+1]>1);
        }
    }

    my @numAnc=keys %difOGs;
    my $numAncSize=scalar(@numAnc);
    my $likelihood=0;
    my $unobslik=0;
    my %likelihood=();
    my @exp=&Q2do(@_);
    foreach my $element (@exp){
        my $key=$element->{'bmodel'}."_".$element->{'distance'};
        my $array_ref=&c_prob($U,$_[$RatesOrder{$element->{'bmodel'}}*2],$_[$RatesOrder{$element->{'bmodel'}}*2],$_[$RatesOrder{$element->{'bmodel'}}*2+1],$element->{'distance'});
        @{$Q{$key}}=@$array_ref;
    }

    if ($options{'unobs'}){
        $OG_analyzed="UnobservedData";
        (!&Felsenstein(\%likelihood)) && do{die;};
        for (my $i=$min{$OG_analyzed};$i<=$max{$OG_analyzed};$i++){
            $unobslik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG_analyzed,\@_);
        }
        return $MAX if ($unobslik >= 1 || $unobslik eq "nan");
    }

    foreach my $string (keys %difOGs){
        my $OG=$difOGs{$string}[0];
        $OG_analyzed=$OG;
        %likelihood=();
        my $lik=0;
        (!&Felsenstein(\%likelihood)) && do{die;};
        for (my $i=$min{$OG};$i<=$max{$OG};$i++){
            $lik+=$likelihood{$root_node_id}{$i}*&root_dist_p($i,$OG,\@_);
        }
	return $MAX if ($lik == 0 || $lik eq "nan");
	if ($options{'unobs'}){
	    $likelihood+=log($lik/(1-$unobslik))*scalar(@{$difOGs{$string}});
	}else{
	    $likelihood+=log($lik)*scalar(@{$difOGs{$string}});
	}
    }
    
    my $tot_prior=0;
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	my $tmp=&prior($_[$RatesOrder{$bmodelstr}*2],\%{$priormodel{$bmodelstr}{'lambda'}});
	if ($tmp == 0  || $tmp eq "-inf" || $tmp > 1){
	    return $MAX;
	}
	$tot_prior+=log($tmp);
    }
    return abs($likelihood+$tot_prior);
}
########################################################################
sub MostProb($$){
    my ($hash,$anc)=@_;
    my @max=values %$hash;
    my $max=&max(\@max);
    my @putative=();
    foreach my $val (keys %$hash){
        if ($$hash{$val} == $max){
	    push(@putative,$val);
        }
    }
    if (scalar(@putative) == 1){
	return $putative[0];
    }else{
	#the value most similar to the ancestor's one
	@putative= sort {$a <=> $b} @putative;
	if ($anc == -1){
	    return $putative[int(scalar(@putative)/2)];
	}else{
	    my $dif=$MAX;
	    my $returned=$putative[int(scalar(@putative)/2)];
	    foreach my $value (@putative){
		if (abs($value-$anc)<$dif){
		    $returned=$value;
		    $dif=abs($value-$anc);
		}
	    }
	    return $returned;
	}
    }
}
########################################################################
sub rootML{
    my ($vertice)=@_;
    if ($options{'root_dist'} == 0){
        foreach my $OG (@usedOG){
            $root{$OG}=$ancestral{$OG}{$root_node_id}[0];
        }
    }elsif ($options{'root_dist'} == 1){
        $root{'ALL'}=$$vertice->[$OGOrder{'ALL'}];
    }elsif ($options{'root_dist'} == 2){
        foreach my $OG (@usedOG){
            $root{$OG}=$$vertice->[$OGOrder{$OG}];
        }
    }elsif ($options{'root_dist'} == 3){
        $root{'ALLp'}=$$vertice->[$OGOrder{'ALL'}];
        $root{'ALLq'}=$$vertice->[$OGOrder{'ALL'}+1];
    }elsif ($options{'root_dist'} == 4){
        foreach my $OG (@usedOG){
            $root{$OG."p"}=$$vertice->[$OGOrder{$OG}];
            $root{$OG."q"}=$$vertice->[$OGOrder{$OG}+1];
        }
    }
    return 1;
}
########################################################################                                                                                                    
sub rootMLprob($$){
    my ($OG,$k)=(@_);

    if ($options{'root_dist'} == 0 || $options{'root_dist'} == 2){
	return &poisson($k,$root{$OG});
    }elsif ($options{'root_dist'} == 1){
        return &poisson($k,$root{'ALL'});
    }elsif ($options{'root_dist'} == 3){
        return &negbin($k,$root{'ALLp'},$root{'ALLq'});
    }elsif ($options{'root_dist'} == 4){
        return &negbin($k,$root{$OG.'p'},$root{$OG.'q'});
    }
}
########################################################################
sub BD{
    my (@xb,@sm)=();
    (!&RatesParamBD(\@xb,\@sm)) && do{die;};
    my ($vertice,$y,$converged)=MinimiseND(\@xb,\@sm,\&BDRatesLik,1e-8,150000);
    if ($converged == 0){
	$error=1;
    }

    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	$birth{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*3];
	$death{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*3+1];
    }

    (!&rootML(\$vertice)) && do{die;};
    if ($y == $MAX){
	$error=2;
	return "-inf";
    }else{
	return -$y;
    }
}
########################################################################
sub BDI{

    #search the ML rates
    my (@xb,@sm)=();
    (!&RatesParamBDI(\@xb,\@sm)) && do{die;};    
    my ($vertice,$y, $converged)=MinimiseND(\@xb,\@sm,\&BDIRatesLik,1e-8,150000);
    if ($converged == 0){
	$error=1;
    }
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	$birth{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*3];
	$death{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*3+1];
	$innovation{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*3+2];
    }
    (!&rootML(\$vertice)) && do{die;};
    if ($y == $MAX){
	$error=2;
	return "-inf";
    }else{
        return -$y;
    }
}
########################################################################
sub LI{

    #search the ML rates
    my (@xb,@sm)=();
    (!&RatesParamLI(\@xb,\@sm)) && do{die;};
    my ($vertice,$y,$converged)=MinimiseND(\@xb,\@sm,\&LIRatesLik,1e-8,150000);
    if ($converged == 0){
	$error=1;
    }
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	$birth{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*2];
	$innovation{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*2+1];
    }
    (!&rootML(\$vertice)) && do{die;};
    if ($y == $MAX){
	$error=2;
	return "-inf";
    }else{
        return -$y;
    }
}
########################################################################                                                                                                    
sub L{
    #search the ML rates                                                                                                                                                     
    my (@xb,@sm)=();
    (!&RatesParamL(\@xb,\@sm)) && do{die;};
    my ($vertice,$y,$converged)=MinimiseND(\@xb,\@sm,\&LRatesLik,1e-8,150000);
    if ($converged == 0){
	$error=1;
    }
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	$birth{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*2];
    }
    (!&rootML(\$vertice)) && do{die;};
    if ($y == $MAX){
	$error=2;
	return "-inf";
    }else{
        return -$y;
    }
}
########################################################################                                                                
sub GD{
    #search the ML rates                                                                                                                
    my (@xb,@sm)=();
    (!&RatesParamGD(\@xb,\@sm)) && do{die;};

    my ($vertice,$y,$converged)=MinimiseND(\@xb,\@sm,\&GDRatesLik,1e-8,150000);
    if ($converged == 0){
	$error=1;
    }
    my @keys=&keysModel;
    foreach my $bmodelstr (@keys){
	$death{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*2];
	$innovation{$bmodelstr}[0]=$vertice->[$RatesOrder{$bmodelstr}*2+1];
    }
    (!&rootML(\$vertice)) && do{die;};
    if ($y == $MAX){
	$error=2;
	return "-inf";
    }else{
        return -$y;
    }
}
########################################################################
sub GetDistances($$){
    #distance between pair of nodes

    my ($tree, $distances)=@_;
    my @nodes=$tree->get_nodes;
    for (my $i=0;$i<scalar(@nodes)-1;$i++){
        for (my $j=$i;$j<scalar(@nodes);$j++){
            my @distnodes=($nodes[$i],$nodes[$j]);
            $$distances{$nodes[$i]->internal_id}{$nodes[$j]->internal_id}=
		$$distances{$nodes[$j]->internal_id}{$nodes[$i]->internal_id}=
		$tree->distance(\@distnodes);
        }
    }

    return 1;

}
########################################################################
sub printTreeIDRel($$$){ 
    my ($fh,$tree,$assoc)=@_;
    print $fh "##NODES-INTERNAL_ID ASSOCIATION\n";
    my $treecp = $tree->clone();
    foreach my $node ($treecp->get_nodes){
	if (!$node->is_Leaf()){
	    my $parent_tmp=$node->ancestor;
	    $node->bootstrap($node->internal_id);
	}else{
	    $$assoc{$node->id}=$node->internal_id;
	    $node->id($node->id."_".$node->internal_id);
	}
    }

    my $tree_string = '';
    open(TMP_OUT, "+<", \$tree_string);
    my $out = new Bio::TreeIO(-fh => \*TMP_OUT, -format => 'newick');
    $out->write_tree($treecp);
    print $fh $tree_string;

return 1;
}
#######################################################################
sub min($){
    my ($data)=@_;    
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@$data);
    return $stat->min();
}
#######################################################################                               
sub max($){
    my ($data)=@_;
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@$data);
    return $stat->max();
}
#######################################################################                                                           
sub collapse($$){
    my ($tree,$node_col)=@_;
    foreach my $node ($$tree->get_nodes){
        if ($node->internal_id eq $node_col->internal_id){
            foreach my $desc ($node->each_Descendent){
                $node->remove_Descendent($desc);
            }
            last;
        }
    }
    return 1;
}
########################################################################
sub shuffle($) {
    my @a = splice @_;
    for my $i (0 .. $#a) {
        my $j = int rand @a;
        @a[$i, $j] = @a[$j, $i];
    }
    return @a;
}
#######################################################################                                                               
sub getExternalest($$){
    my ($tree)=@_;

    my @nodes = $tree->get_nodes;
    @nodes = &shuffle(@nodes);
  NODE: foreach my $node (@nodes){
      my $id=$node->internal_id;
      my $parent=$node->ancestor;
      next if (!defined($parent));
      if (scalar($node->each_Descendent)==0){
          my $node_sister=&GetTheSister($parent, $node);
          if (scalar($node_sister->each_Descendent)==0){
              return $parent;
          }
      }
  }
    die;
}
########################################################################                                                                
sub GetTheSister($$){
    my ($parent,$node)=@_;
    
    foreach my $desc ($parent->each_Descendent){
        next if ($desc->internal_id eq $node->internal_id);
        return $desc;
    }
    die;
}
################################################################################
sub GetLeaves($){
    my ($node)=@_;
    my @leaves=();
    my @desc=$node->get_all_Descendents;
    foreach my $desc (@desc){
        if ($desc->is_Leaf){
            push(@leaves,$desc->internal_id);
        }
    }
    return @leaves;
}
#######################################################################
sub getDesc($$$){
    
    my ($tree,$node,$descArray)=@_;
    foreach my $desc ($node->each_Descendent){
        push(@$descArray,[$desc->internal_id]);
    }

    my (@all_leaves,@leaves_filt)=();
    @all_leaves=&GetLeaves($root_node);
    foreach my $node_tmp ($tree->get_nodes){
        if ($node_tmp->internal_id eq $node->internal_id){
            @leaves_filt=&GetLeaves($node_tmp);
            my $tmp=join("|",@leaves_filt);
            @leaves_filt=sort {$a <=> $b} grep(!/^($tmp)$/,@all_leaves);
            last;
        }
    }
    push(@$descArray,[@leaves_filt]);
    return 1;
}
#########################################################################
sub readCtlFile($$){
    my ($options,$file)=@_;
    (!open(FILE,$file)) && do{return 0;};
    while(<FILE>){
	next if ($_=~/^#/);
	$_=~s/^(.+)#.+/$1/g;
	next if ($_!~/\S+/);
	chomp;
	my ($option,$value)=($_=~/(\S+)\s+=\s+(\S+)/);
	if (!defined($option) || !defined($value)){
	    close (FILE);
	    warn "Incorrectly defined option in $_\n";
	    return 0;
	}
	if ($option eq 'priorfile' && $value eq "0"){
	    $$options{$option}="";
	}else{
	    $$options{$option}=$value;
	}
    }
    close(FILE);
    return 1;
}
#########################################################################
sub readtree($){
	my ($file)=@_;
	my $treeio = Bio::TreeIO->new('-format' => 'newick', '-file' => $file);
	return $treeio->next_tree;
}
#########################################################################
sub Ulimit(){
    $U=0;
    foreach my $key (keys %max){
	if ($max{$key}>$U){
	    $U=$max{$key};
	}
    }
    $U*=1;
    return 1;
}
#########################################################################
sub MakeTreePath($){
    my ($treepath)=@_;
    my $treecp = $tree->clone();
    while(scalar($treecp->get_nodes)>2){
        my @desc=();
        my $node=&getExternalest($treecp);
        my $nodeid=$node->internal_id;
	my $parentid='NULL';
	if ($node->ancestor){
	    $parentid=$node->ancestor->internal_id;
	}
        (!&getDesc($tree,$node,\@desc)) && do{die;};
	
	push(@treepath,{'id'=>$nodeid, 'desc'=>[@desc], 'parentid'=>$parentid});
	(!&collapse(\$treecp,$node)) && do{die;};
    }
    return 1;
}
#########################################################################
sub dieNice{

    print STDERR "perl BadiRate.pl -treefile NEWICK_FILE -sizefile FAMILY_SIZE_FILE [options] > output.bd
                                                                                                                                      
                            
        -anc         = Reconstruct ancestral family sizes and the minimum number of gains/losses in each lineage                       
        -bmodel      = Run local or free branch models                                                                           
        -sizefile    = Family Size File
        -unobs       = Correct the likelihood for families absent in all extant species                       
        -h|help      = Display this help                                                                                       
        -rmodel      = Family turnover rates to be estimated                                                                           
        -out         = Set the output file                                                                                           
        -outlier     = Report families escaping the estimated stochastic process                                              
        -ep          = Define the estimation procedure                                                                                 
        -n_max_int   = Modify the maximum number of family members in the internal phylogenetic nodes 
        -print_ids   = Display nodes ids in Newick format                                                                        
        -priorfile   = Prior File                                                                                                 
        -root_dist   = Estimation method for the root a priori distribution                                                            
        -seed        = Seed of the pseudo-random number generator                                                                 
        -start_val   = Starting values for the likelihood methods                                                                    
        -treefile    = Phylogenetic tree in Newick format                                                                        
        -version     = Report the BadiRate version                                                                                  
                         
    
                      -----OR------

        perl BadiRate.pl controlfile.bd > output.bd     

";       

    exit(0);
}
############################# C Routines ################################
__DATA__
__C__

   #include <stdlib.h>
   #include <math.h>
 
  double comb (double n, double k)
  {
    double prod=1;
    double cnt=0;
    double i;
    for (i=k;i>=1;i--){
        prod*=(n-cnt)/i;
        if (prod<0){
            abort();
        }

        cnt++;
    }
    return prod;
  }  

  
  double factorial(double x)
  {
    if (x==0)
        return 1;

    double fact,i;
    fact=1;
    for (i=x;i>0;i--){
        fact*=i;
    }
    return fact;
  }
#define M_PI 3.14159265358979323846
#define A 12
double gamma_spouge(double z){
      const int a = A;
      static double c_space[A];
      static double *c = NULL;
      int k;
      double accm;
      
      if ( c == NULL ) {
	  double k1_factrl = 1.0; 
	  c = c_space;
	  c[0] = sqrt(2.0*M_PI);
	  for(k=1; k < a; k++) {
	      c[k] = exp(a-k) * pow(a-k, k-0.5) / k1_factrl;
	      k1_factrl *= -k;
	  }
      }
      accm = c[0];
      for(k=1; k < a; k++) {
	  accm += c[k] / ( z + k );
      }
      accm *= exp(-(z+a)) * pow(z+a, z+0.5);
      return accm/z;
  } 	
  double poisson(double x,double lambda){
      return pow(lambda,x)*exp(-lambda)/factorial(x);
  }
  double binomial (double x, double n, double p){
      return comb(n,x)*pow(p,x)*pow((1-p),(n-x));
  }
  double negbin (double x, double n, double p){
      return comb(x+n-1,x)*pow((1-p),n)*pow(p,x);
  }
  double gamma_PDF(double x, double shape, double scale){
      return pow(x,(shape-1))*(exp(-x/scale)/(gamma_spouge(shape)*pow(scale,shape)));
  }
  SV *c_prob(int U,double birth, double death, double innovation, double time){
      int i;
      int j;
      AV *array=newAV();
      double p,q,r,theta;
      if (birth >0){ /*for BDI processes*/
          theta=innovation/birth;
      }
      if (death >0){ 
	  r=innovation*(1-exp(-death*time))/death;
      }else{
	  r=innovation*time; /*pure gain process*/
      }

      if (birth == death){
	  q=birth*time/(1+birth*time);
          p=q;
      }else{
          p=(death-death*exp(-(death-birth)*time))/(death-birth*exp(-(death-birth)*time));
	  q=(birth-birth*exp(-(death-birth)*time))/(death-birth*exp(-(death-birth)*time)); 
      }
     
      for (j = 0 ;j<=U; j++){
	  double prob;
	  if (birth >0 && innovation > 0){
	      prob=negbin(j,theta,q);
		if (prob != prob || prob == INFINITY){ /*prob is nan or infinity, probably because the birth is so small that it should be a GD process*/
			prob=poisson(j,r);	
		}
	  }else if (innovation > 0 && birth == 0){
	      prob=poisson(j,r);
	  }else{
	      if (j == 0){
		  prob=1;
	      }else{
		  prob=0;
	      }
	  }
	  av_push(array,newSVnv(prob));
      }

      for (i = 1 ;i<=U; i++){
	  for (j = 0 ;j<=U; j++){
	      double prob;	      
	      if (j == 0){
 		  prob=SvNV(*(av_fetch(array,(i-1)*(U+1),0)))*p;
	      }else if (j == 1){
		  prob=SvNV(*(av_fetch(array,(i-1)*(U+1),0)))*(1-p)*(1-q)+SvNV(*(av_fetch(array,(i-1)*(U+1)+1,0)))*p;
	      }else{
		  prob=SvNV(*(av_fetch(array,i*(U+1)+j-1,0)))*q+SvNV(*(av_fetch(array,(i-1)*(U+1)+j-1,0)))*(1-p-q)+SvNV(*(av_fetch(array,(i-1)*(U+1)+j,0)))*p;
	      }
	      av_push(array,newSVnv(prob));
	  }
      }
      return newRV_noinc((SV*)array);
  }


  SV *likBranch (int U, int min, int max, AV *Q1, AV *Q2, SV *lik1_ref, SV *lik2_ref){
      int i=0;
      HV *lik = newHV();
      HV *lik1 = (HV*)SvRV(lik1_ref);
      HV *lik2 = (HV*)SvRV(lik2_ref);
      for (i=min;i<=max;i+=1){
	  double tmp1=0;
	  double tmp2=0;
	  int j=0;
	  int t=0;
	  int num_keys1 = hv_iterinit(lik1);
	  int num_keys2 = hv_iterinit(lik2);
	  for (j=0;j<num_keys1;j+=1){
	      HE *entry1 = hv_iternext(lik1);
	      SV *key = hv_iterkeysv(entry1);
	      SV *val = hv_iterval(lik1, entry1);
	      SV *a=*(av_fetch(Q1,i*(U+1)+SvIV(key),0));
	      tmp1+=SvNV(val)*SvNV(a);
	  }

	  for (t=0;t<num_keys2;t+=1){
	      HE *entry2 = hv_iternext(lik2);
	      SV *key = hv_iterkeysv(entry2);
              SV *val = hv_iterval(lik2, entry2);
	      SV *a=*(av_fetch(Q2,i*(U+1)+SvIV(key),0));
	      tmp2+=SvNV(val)*SvNV(a);
	  }
	  char num[10];
	  sprintf(num,"%d",i);
          hv_store(lik,num,strlen(num),newSVnv(tmp1*tmp2), 0);
      }
      return newRV_noinc((SV*) lik);
   }


      
      
