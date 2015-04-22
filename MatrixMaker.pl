
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# This program will take an NCBI tax id as input and will use Uniprot
# API to obtain proteins for that organism and the align them, find blocks
# and use them to create symmetrical and non symmetrical blosum substitution matrices for that organism

use strict;
use warnings;

# TK Related Imports
require Tk::LabFrame;
require Tk::Text;
require Tk::BrowseEntry;
use Tk;

# Other modules
use Math::Round;
use log2;
use Cwd;

# MartixMaker modules
use giveclusters;
use cluster_number_to_fasta;
use fastatoblocks;


# Main Window TK
my $mw = new MainWindow;
$mw->geometry("400x350+200+200");
#$mw->configure(-background=>'grey');
$mw->title("Matrix Maker");
my $main_menu = $mw->Menu();
$mw->configure(-menu => $main_menu);
$main_menu->command(-label=>"Help", -underline => 0,
                    -command=>sub{`help.txt`});
$main_menu->command(-label=>"About", -underline => 0,
                    -command=>sub{`about.txt`});
my $label = $mw -> Label(-text=>"Enter the following details") -> pack();
my $genus = "eg.Entamoeba";
my $species = "eg.histolytica";
my $NCBI = "eg.5759";
my $Max = "eg.8";
my $Min = "eg.5";
my $CluserNum = "eg.1";
my $frm1 = $mw->Frame(-relief=>'groove',
		      -borderwidth=>3) -> pack(-side=>'top',-fill=>'x');
my $label1 = $frm1 -> Label(-text=>"Enter Genus-") -> pack(-side=>'left');
my $ent1 = $frm1 -> Entry(-textvariable=>\$genus,
			-justify =>'left',
			-font=>'Calibri 10 italic') -> pack(-side=>'right');


my $frm2 = $mw->Frame(-relief=>'groove',
		      -borderwidth=>3) -> pack(-side=>'top',-fill=>'x');
my $label2 = $frm2 -> Label(-text=>"Enter Species-") -> pack(-side=>'left');
my $ent2 = $frm2 -> Entry(-textvariable=>\$species,
			-justify =>'left',
			-font=>'Calibri 10 italic') -> pack(-side=>'right');


my $frm3 = $mw->Frame(-relief=>'groove',
		      -borderwidth=>3) -> pack(-side=>'top',-fill=>'x');
my $label3 = $frm3 -> Label(-text=>"Enter NCBI Tax. ID-") -> pack(-side=>'left');
my $ent3 = $frm3 -> Entry(-textvariable=>\$NCBI,
			-justify =>'left',
			-font=>'Calibri 10 italic') -> pack(-side=>'right');


my $frm4 = $mw->Frame(-relief=>'groove',
		      -borderwidth=>3) -> pack(-side=>'top',-fill=>'x');
my $label4 = $frm4 -> Label(-text=>"Enter Max Per Cluster-") -> pack(-side=>'left');
my $ent4 = $frm4 -> Entry(-textvariable=>\$Max,
			-justify =>'left',
		         -font=>'Calibri 10 italic') -> pack(-side=>'right');

my $frm5 = $mw->Frame(-relief=>'groove',
		      -borderwidth=>3) -> pack(-side=>'top',-fill=>'x');
my $label5 = $frm5 -> Label(-text=>"Enter Min Per Cluster-") -> pack(-side=>'left');
my $ent5 = $frm5 -> Entry(-textvariable=>\$Min,
			-justify =>'left',
			-font=>'Calibri 10 italic') -> pack(-side=>'right');


my $frm6 = $mw->Frame(-relief=>'groove',
		      -borderwidth=>3) -> pack(-side=>'top',-fill=>'x');
my $label6 = $frm6 -> Label(-text=>"Enter Total Number of Clusters-") -> pack(-side=>'left');
my $ent6 = $frm6 -> Entry(-textvariable=>\$CluserNum,
			-justify =>'left',
			-font=>'Calibri 10 italic') -> pack(-side=>'right');

 my $frm7 = $mw->Frame(-relief=>'groove',
		      -borderwidth=>3) -> pack(-side=>'top',-fill=>'x');
my $label7 = $frm7 -> Label(-text=>"Enter Taxons to Select-") -> pack(-side=>'left');
 my $src=$frm7->Text(-height=>2,
               -width=>40,
               -background=>'white',
     )->pack(-side=>'right');

my $Proxy;
my $frame = $mw->LabFrame(
		-label => "Proxy",
		-labelside => 'acrosstop',
		-width => 110,
		-height => 60,
		)->place(-x=>10,-y=>230);

	# Put these values into the frame
	$frame->Radiobutton(
		-variable => \$Proxy,
		-value => '1',
		-text => 'Yes',
		)->place( -x => 10, -y => 0 );
	$frame->Radiobutton(
		-variable => \$Proxy,
		-value => '0',
		-text => 'No',
		)->place( -x => 10, -y => 20 );

my $prioratisation;
my $frame1 = $mw->LabFrame(
		-label => "Prioritisation",
		-labelside => 'acrosstop',
		-width => 110,
		-height => 60,
		)->place(-x=>140,-y=>230);

	# Put these values into the frame
	$frame1->Radiobutton(
		-variable => \$prioratisation,
		-value => '1',
		-text => 'Yes',
		)->place( -x => 10, -y => 0 );
	$frame1->Radiobutton(
		-variable => \$prioratisation,
		-value => '0',
		-text => 'No',
		)->place( -x => 10, -y => 20 );
   
my $unirefno;
my $frame2 = $mw->LabFrame(
		-label => "UnirefNo.",
		-labelside => 'acrosstop',
		-width => 110,
		-height => 60,
		)->place(-x=>270,-y=>230);

	
	$frame2->Radiobutton(
		-variable => \$unirefno,
		-value => 'UniRef90',
		-text => 'UniRef90',
		)->place( -x => 10, -y => 0 );
	$frame2->Radiobutton(
		-variable => \$unirefno,
		-value => 'UniRef50',
		-text => 'UniRef50',
		)->place( -x => 10, -y => 20 );
  

my $button = $mw -> Button(-text => "Submit",
	-command =>sub {
	    my $itext = $src->Contents();
	    
	    my $unirefsend;
	    if($itext !~/([a-z]|[A-Z]|[0-9])/)
	    {
		$itext ="Google";
	    }
	    else
	    {
	    $itext =~s/\n/:/g;
	    if($itext=~/:$/)
	    {
		
		chop($itext);
		chop($itext);
	    }
	    }
	    if($unirefno eq "UniRef100")
	    {
		$unirefsend="1.0";
	    }
	    elsif($unirefno eq "UniRef90")
	    {
		$unirefsend=0.9;
	    }
	    elsif($unirefno eq "UniRef50")
	    {
		$unirefsend=0.5;
	    }
	    my $label1 = $mw -> messageBox(-message=>"The matrices are being computed\nGo,Enjoy a cuppa coffee\nYou will be notified when done");
                commandlineargv_finaldr($genus,$species,$NCBI,$Max,$Min,$prioratisation,$itext,$Proxy,$CluserNum,$unirefsend);
	    

	    sub commandlineargv_finaldr
	    {
	
	    my $mw1 = new MainWindow;
	    $mw1->geometry("250x200+200+200");
	    #$mw->configure(-background=>'grey');
	    $mw1->title("Matrix Results");
	    
	    
	    my $originalcwd = getcwd();
	    `mkdir Results`;
	    chdir("Results");
	    my $permcwd = getcwd();
	    
	    
	    
	    my $genus = shift;
	    my $species = shift;
	    my $ncbiid = shift;
	    my $max =shift;
	    my $min = shift;
	    my $prioratisation=shift;
	    my $restrictcmnds=shift;
	    my $proxy=shift;
	    my $numberofclusers=shift;
	    my $unirefno = shift;
	    
	    my @unirefids = getclusers("$genus","$species",$ncbiid,$proxy,$unirefno);
	    my $restrictions;
	    my @species;
	    if($restrictcmnds !~/Google/)
	    {
	        $restrictions=1;
	     @species= split(":",$restrictcmnds);
	    }
	    else
	    {
	        $restrictions=0;
	    }
	    `mkdir Unirefids`;
	    chdir("Unirefids");
	    open(FA,">","unirefids.txt");
	    
	    foreach(@unirefids)
	    {
	        print FA"$_\n";
	    }
	    close(FA);
	    
	    chdir("$permcwd");
	    
	    
	    
	    
	    my $count=0;
	    my $count1=0;
	    
	    foreach my $unirefid(@unirefids)
	    {
	        my $multiplefasta;
	    if($restrictions==1)
	    {
	    $multiplefasta=getfasta($unirefid,$genus,$species,$max,$min,$permcwd,$restrictions,$prioratisation,$proxy,@species);
	    }
	    elsif($restrictions==0)
	    {
	    $multiplefasta=getfasta($unirefid,$genus,$species,$max,$min,$permcwd,$restrictions,$prioratisation,$proxy);   
	    }
	    chdir("$permcwd");
	    if($multiplefasta !~/Google/i)
	    {
	    my @blocks = blocksretrieve($multiplefasta,$proxy);
	    `mkdir Blocks`;
	    chdir("Blocks");
	    open(FE,">","$unirefid.txt");
	    my $blocktest=0;
	    foreach(@blocks)
	    {
	       if($_ =~/AA$genus/)
	       {
	       print FE "$_\n";
	     
	       $blocktest=1;
	       }
	       
	    }
	    close(FE);
	    chdir("$permcwd");
	    if($blocktest>0)
	    {
	    $count++;
	    }
	    print "+ $count\t";
	    if($count>=$numberofclusers)
	    {
	       last;
	    }
	    
	    }
	    $count1++;
	    print "Cluster Number $count1\n";
	    
	    
	    }
	    
	    
	    chdir("$permcwd");
	    
	    
	    
	    my @scoring_matrix;
	    my @scoringmatrix2;
	    my %aminohash;
	    my $totalnumber=0;
	    my $totalpairs=0;
	    for(my $i=0;$i<23;$i++)
	    {
	        for(my $j=0;$j<23;$j++)
	        {
	            $scoring_matrix[$i][$j]=0;
	            $scoringmatrix2[$i][$j]=0;
	        }
	    }
	    
	    my @aminoarr= qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X);
	    foreach(@aminoarr)
	    {
	        $aminohash{$_}=0;
	    }
	    my @numb = 0..22;
	    my %hash;
	    my $cnt=0;
	    while($cnt<23)
	    {
	    $hash{shift(@aminoarr)}=shift(@numb);
	    $cnt++;
	    }
	    
	    my $multiplefasta;
	    chdir("Blocks");
	     my @files = <*>;
	     foreach my $file (@files) {
	       
	     
	    open(FL,$file);
	    my $string;
	    foreach(<FL>)
	    {
	        if($_=~/^$/)
	        {
	           $string.="@@@"; 
	        }
	        else
	        {
	           $string.=$_; 
	        }
	    }
	    
	    
	    #pop(@returnarray);
	    if(length($string)>0)
	    {
	    my @returnarray = split("@@@",$string);
	    foreach(@returnarray)
	    {
	    
	    oneblock($_,\@scoring_matrix,\@scoringmatrix2,\%hash,\%aminohash,\$totalnumber,\$totalpairs);
	    }
	    
	     }
	    close(FL);
	     }
	    print"";
	    my %matrixhash;
	    my $sum=0;
	    foreach(keys %hash)
	    {
	        $matrixhash{$hash{$_}}=$_;
	        
	    }
	    foreach(keys %aminohash)
	    {
	        
	        $sum+=$aminohash{$_};
	    }
	    $totalpairs= $totalpairs/2;
	    print "";
	    my @scoringmatrix5=newmatrixprint1(\@scoringmatrix2,\%hash,\%matrixhash,\%aminohash,$totalpairs,$sum,$permcwd);
	    print "";
	    
	    sub newmatrixprint1
	    {
	       my $scoringmatrix2 = shift;
	       my $hash = shift;
	       my $matrixhash = shift;
	       my $aminohash = shift;
	       my $totalpairs = shift;
	       my $sum = shift;
	       my $permcwd = shift;
	       my %scoringmatrix3;
	      my  %matrixhash = %{$matrixhash};
	        for(my $i=0;$i<23;$i++)
	    {
	       for(my $j=0;$j<23;$j++)
	        {
	            if(${$scoringmatrix2}[$i][$j] != 0)
	            {
	            if($i==$j)
	            {
	                
	                $scoringmatrix3{$matrixhash{$i}}+=${$scoringmatrix2}[$i][$j];
	            }
	            else
	            {
	                
	                $scoringmatrix3{$matrixhash{$i}}+=(${$scoringmatrix2}[$i][$j]/2);
	            }
	            }
	        
	        }
	        if(exists $scoringmatrix3{$matrixhash{$i}})
	            {
	        $scoringmatrix3{$matrixhash{$i}}/=$totalpairs;
	            }
	        }
	        my @scoringmatrix4;
	        for(my $i=0;$i<23;$i++)
	    {
	       for(my $j=0;$j<23;$j++)
	        {
	            if(${$scoringmatrix2}[$i][$j]!=0&&exists $scoringmatrix3{$matrixhash{$i}}&&exists $scoringmatrix3{$matrixhash{$j}})
	            {
	            if($i==$j)
	            {
	                
	                my $denominator = (($scoringmatrix3{$matrixhash{$i}})*($scoringmatrix3{$matrixhash{$i}}));
	                my $numerator   = ${$scoringmatrix2}[$i][$j]/$totalpairs;
	                $scoringmatrix4[$i][$j]= 2*log2($numerator/$denominator);
	            }
	            else
	            {
	                my $denominator = 2*(($scoringmatrix3{$matrixhash{$i}})*($scoringmatrix3{$matrixhash{$j}}));
	                my $numerator   = ${$scoringmatrix2}[$i][$j]/$totalpairs;
	                $scoringmatrix4[$i][$j]= 2*log2($numerator/$denominator);
	            }
	            }
	        }
	        }
	        my $entropy;
	        for(my $i=0;$i<20;$i++)
	        {
	          for(my $j=0;$j<($i+1);$j++)
	          {
	            my $pehla = ${$scoringmatrix2}[$i][$j]/$totalpairs;
	            my $doosra = $scoringmatrix4[$i][$j];
	            if(defined $pehla&&defined $doosra)
	            {
	            $entropy+=($pehla)*($doosra);
	            }
	          }
	        }
	        print $entropy."\n";
	        my $expect;
	        for(my $i=0;$i<20;$i++)
	        {
	          for(my $j=0;$j<($i+1);$j++)
	          {
	            my $pehla = $scoringmatrix3{$matrixhash{$i}};
	            my $doosra = $scoringmatrix3{$matrixhash{$j}};
	            my $teesra = $scoringmatrix4[$i][$j];
	            if(defined $pehla&&defined $doosra&&defined $teesra)
	            {
	            $expect+=($pehla)*($doosra)*($teesra);
	            }
	          }
	        }
	         print $expect."\n";
	         chdir("$permcwd");
	         `mkdir Matrices`;
	    chdir("Matrices");
	           open(FH,">","SYM-MAT.xls");
	       print FH "\t";
	       my @aminoarr= qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X);
	       my @aminoarr2 = @aminoarr;
	       foreach(@aminoarr)
	    {
	        print FH "$_\t";
	    }
	       print FH "\n";
	       for(my $i=0;$i<23;$i++)
	    {
	        my $tempo = shift(@aminoarr2);
	        print FH $tempo."\t";
	        for(my $j=0;$j<23;$j++)
	        {
	            if(defined $scoringmatrix4[$i][$j])
	            {
	            
	            my $google = nearest(1,$scoringmatrix4[$i][$j]);
	            print FH"$google\t";
	            }
	            else
	            {
	             print FH "0\t";   
	            }
	            
	        }
	        print FH "\n";
	        
	    }
	       close(FH);
	    }
	    
	    
	    
	    sub oneblock
	    {
	       my $block = shift;
	       my $scoringmatrix = shift;
	       my $scoringmatrix2=shift;
	       my $hashtemp = shift;
	       my $aminohash = shift;
	       my $totalnumber = shift;
	       my $totalpairs = shift;
	       my %hashmat = %{$hashtemp};
	       my @blockarr = split("\n",$block);
	       my @coolarr;
	       my %hash = %{$hashtemp};
	       
	      my %aminohash = %{$aminohash};
	       foreach(@blockarr)
	       {
	          if($_!~/^(ID|AC|DE|BL)/)
	          {
	             if($_ =~/\)\s(\D+)/)
	             {
	                my $temp = $1;
	                $temp =~s/\s+//g;
	                push(@coolarr,$temp);
	                
	                
	             }
	          }
	          
	       }
	       my $entam = shift(@coolarr);
	       my $init = 1;
	       while($init<=length($entam))
	       {
	          my $ehist = substr($entam,$init-1,1);
	          my $dastring=$ehist;
	          
	          ${$aminohash}{$ehist}++;
	       foreach(@coolarr)
	                {
	                   
	                   my $other = substr($_,$init-1,1);
	                   ${$scoringmatrix}[$hashmat{$ehist}][$hashmat{$other}]++;
	                   ${$totalnumber}++;
	                   ${$aminohash}{$other}++;
	                   $dastring.=$other;
	                }
	                my %hash1;
	                foreach(split("",$dastring))
	                {
	                    $hash1{$_}++;
	                }
	                foreach my $mainkey(keys %hash1)
	                {
	                    if($hash1{$mainkey}>1)
	                    {
	                        my $limit = $hash1{$mainkey}-1;
	                        for(my $i=$limit;$i>=1;$i--)
	                        {
	                            ${$scoringmatrix2}[$hash{$mainkey}][$hash{$mainkey}]+=$i;
	                            ${$totalpairs}+=$i+$i;
	                        }
	                    }
	                    foreach my $subkey(keys %hash1)
	                    {
	                        if($subkey ne $mainkey)
	                        {
	                            if($hash1{$mainkey}>=$hash1{$subkey})
	                            {
	                               if($hash1{$mainkey}>1&&$hash1{$subkey}>1)
	                               {
	                               ${$scoringmatrix2}[$hash{$mainkey}][$hash{$subkey}]+=$hash1{$mainkey}*$hash1{$subkey}; 
	                               ${$totalpairs}+=$hash1{$mainkey}*$hash1{$subkey};
	                               }
	                               else
	                               {
	                               ${$scoringmatrix2}[$hash{$mainkey}][$hash{$subkey}]+=$hash1{$mainkey};
	                               ${$totalpairs}+=$hash1{$mainkey};
	                               }
	                            }
	                            elsif($hash1{$mainkey}<$hash1{$subkey})
	                            {
	                                if($hash1{$mainkey}>1&&$hash1{$subkey}>1)
	                               {
	                                ${$scoringmatrix2}[$hash{$mainkey}][$hash{$subkey}]+=$hash1{$mainkey}*$hash1{$subkey};
	                                ${$totalpairs}+=$hash1{$mainkey}*$hash1{$subkey};
	                               }
	                               else
	                               {
	                                ${$scoringmatrix2}[$hash{$mainkey}][$hash{$subkey}]+=$hash1{$subkey};
	                               ${$totalpairs}+=$hash1{$subkey};
	                               }
	                            }
	                        }
	                    }
	                    
	                }
	          $init++;
	       }
	       
	       
	    
	    }
	    #______________________________________________________________________________
	    chdir("$permcwd");
	    my @scoring_matrixns;
	    my @scoringmatrix2ns;
	    my %aminohashns;
	    my $totalnumberns=0;
	    my $totalpairsns=0;
	    for(my $i=0;$i<23;$i++)
	    {
	        for(my $j=0;$j<23;$j++)
	        {
	            $scoring_matrixns[$i][$j]=0;
	            $scoringmatrix2ns[$i][$j]=0;
	        }
	    }
	    
	    my @aminoarrns= qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X);
	    foreach(@aminoarrns)
	    {
	        $aminohashns{$_}=0;
	    }
	    my @numbns = 0..22;
	    my %hashns;
	    my $cntns=0;
	    while($cntns<23)
	    {
	    $hashns{shift(@aminoarrns)}=shift(@numbns);
	    $cntns++;
	    }
	    
	    
	    my $multiplefastans;
	    chdir("Blocks");
	     my @filesns = <*>;
	     foreach my $filens (@filesns) {
	       
	     
	    open(FL,$filens);
	    my $stringns;
	    foreach(<FL>)
	    {
	        if($_=~/^$/)
	        {
	           $stringns.="@@@"; 
	        }
	        else
	        {
	           $stringns.=$_; 
	        }
	    }
	    if(length($stringns)>0)
	    {
	    my @returnarrayns = split("@@@",$stringns);
	    #pop(@returnarray);
	    
	    foreach(@returnarrayns)
	    {
	    
	    oneblock1($_,\@scoring_matrixns,\@scoringmatrix2ns,\%hashns,\%aminohashns,\$totalnumberns,\$totalpairsns);
	    }
	    
	    }
	    close(FL);
	     }
	    print"";
	    my %matrixhashns;
	    my $sumns=0;
	    foreach(keys %hashns)
	    {
	        $matrixhashns{$hashns{$_}}=$_;
	        
	    }
	    foreach(keys %aminohashns)
	    {
	        
	        $sumns+=$aminohashns{$_};
	    }
	    $totalpairsns= $totalpairsns/2;
	    print "";
	    my @scoringmatrix5ns=newmatrixprint2(\@scoring_matrixns,\%hashns,\%matrixhashns,\%aminohashns,$totalnumberns,$sumns,$permcwd);
	    print "";
	    
	    sub newmatrixprint2
	    {
	       my $scoringmatrix2 = shift;
	       my $hash = shift;
	       my $matrixhash = shift;
	       my $aminohash = shift;
	       my $totalpairs = shift;
	       my $sum = shift;
	       my $permcwd = shift;
	       my %scoringmatrix3;
	       my %scoringhashh;
	    my %aminohash = %{$aminohash};
	    my %matrixhash = %{$matrixhash};
	        
	        
	        foreach(keys %aminohash)
	        {
	          $scoringmatrix3{$_} = $aminohash{$_}/$sum; 
	        }
	        my @scoringmatrix4;
	        for(my $i=0;$i<23;$i++)
	    {
	       for(my $j=0;$j<23;$j++)
	        {
	            if(${$scoringmatrix2}[$i][$j]!=0&&exists $scoringmatrix3{$matrixhash{$i}}&&exists $scoringmatrix3{$matrixhash{$j}})
	            {
	            if($i==$j)
	            {
	                
	                my $denominator = (($scoringmatrix3{$matrixhash{$i}})*($scoringmatrix3{$matrixhash{$i}}));
	                my $numerator   = ${$scoringmatrix2}[$i][$j]/$totalpairs;
	    	    if($denominator != 0)
	    	    {
	                $scoringmatrix4[$i][$j]= 2*log2($numerator/$denominator);
	    	    }
	            }
	            else
	            {
	                my $denominator = ($scoringmatrix3{$matrixhash{$i}})*($scoringmatrix3{$matrixhash{$j}});
	                my $numerator   = ${$scoringmatrix2}[$i][$j]/$totalpairs;
	    	     if($denominator != 0)
	    	    {
	                $scoringmatrix4[$i][$j]= 2*log2($numerator/$denominator);
	    	    }
	            }
	            }
	        }
	        }
	            chdir("$permcwd");
	      
	    chdir("Matrices");
	        
	           open(FM,">","NSYM-MAT.xls");
	       print FM "\t";
	       my @aminoarr= qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X);
	       my @aminoarr2 = @aminoarr;
	       foreach(@aminoarr)
	    {
	        print FM "$_\t";
	    }
	       print FM "\n";
	       for(my $i=0;$i<23;$i++)
	    {
	        my $tempo = shift(@aminoarr2);
	        print FM $tempo."\t";
	        for(my $j=0;$j<23;$j++)
	        {
	            if(defined $scoringmatrix4[$i][$j])
	            {
	              my $google = nearest(1,$scoringmatrix4[$i][$j]);
	            print FM"$google\t";
	            }
	            else
	            {
	             print FM "0\t";   
	            }
	            
	        }
	        print FM "\n";
	        
	    }
	       close(FM);
	    }
	    
	    
	    
	    sub oneblock1
	    {
	       my $block = shift;
	       my $scoringmatrix = shift;
	       my $scoringmatrix2=shift;
	       my $hashtemp = shift;
	       my $aminohash = shift;
	       my $totalnumber = shift;
	       my $totalpairs = shift;
	       my %hashmat = %{$hashtemp};
	       my @blockarr = split("\n",$block);
	       my @coolarr;
	       
	       foreach(@blockarr)
	       {
	          if($_!~/^(ID|AC|DE|BL)/)
	          {
	             if($_ =~/\)\s(\D+)/)
	             {
	                my $temp = $1;
	                $temp =~s/\s+//g;
	                push(@coolarr,$temp);
	                
	                
	             }
	          }
	          
	       }
	       my $entam = shift(@coolarr);
	       my $init = 1;
	       while($init<=length($entam))
	       {
	          my $ehist = substr($entam,$init-1,1);
	          my $dastring=$ehist;
	          
	          ${$aminohash}{$ehist}++;
	       foreach(@coolarr)
	                {
	                   
	                   my $other = substr($_,$init-1,1);
	                   ${$scoringmatrix}[$hashmat{$ehist}][$hashmat{$other}]++;
	                   ${$totalnumber}++;
	                   ${$aminohash}{$other}++;
	                   $dastring.=$other;
	                }
	                my %hash1;
	                foreach(split("",$dastring))
	                {
	                    $hash1{$_}++;
	                }
	                foreach my $mainkey(keys %hash1)
	                {
	                    if($hash1{$mainkey}>1)
	                    {
	                        my $limit = $hash1{$mainkey}-1;
	                        for(my $i=$limit;$i>=1;$i--)
	                        {
	                            ${$scoringmatrix2}[$hashmat{$mainkey}][$hashmat{$mainkey}]+=$i;
	                            ${$totalpairs}+=$i+$i;
	                        }
	                    }
	                    foreach my $subkey(keys %hash1)
	                    {
	                        if($subkey ne $mainkey)
	                        {
	                            if($hash1{$mainkey}>=$hash1{$subkey})
	                            {
	                               if($hash1{$mainkey}>1&&$hash1{$subkey}>1)
	                               {
	                               ${$scoringmatrix2}[$hashmat{$mainkey}][$hashmat{$subkey}]+=$hash1{$mainkey}*$hash1{$subkey}; 
	                               ${$totalpairs}+=$hash1{$mainkey}*$hash1{$subkey};
	                               }
	                               else
	                               {
	                               ${$scoringmatrix2}[$hashmat{$mainkey}][$hashmat{$subkey}]+=$hash1{$mainkey};
	                               ${$totalpairs}+=$hash1{$mainkey};
	                               }
	                            }
	                            elsif($hash1{$mainkey}<$hash1{$subkey})
	                            {
	                                if($hash1{$mainkey}>1&&$hash1{$subkey}>1)
	                               {
	                                ${$scoringmatrix2}[$hashmat{$mainkey}][$hashmat{$subkey}]+=$hash1{$mainkey}*$hash1{$subkey};
	                                ${$totalpairs}+=$hash1{$mainkey}*$hash1{$subkey};
	                               }
	                               else
	                               {
	                                ${$scoringmatrix2}[$hashmat{$mainkey}][$hashmat{$subkey}]+=$hash1{$subkey};
	                               ${$totalpairs}+=$hash1{$subkey};
	                               }
	                            }
	                        }
	                    }
	                    
	                }
	          $init++;
	       }
	       
	       
	       
	    }#___________________________________________________________________________
	    chdir("$permcwd");
	    chdir("Matrices");
	    my @files1 = <*>;
	    my $counter=1;
	    my @arr = qw(NSYM-MAT SYM-MAT);
	     foreach my $file (@files1) {
	        if($file =~/xls/)
	        {
	    open(AI,"$file");
	    #------------The output file
	     
	    my $temp = shift(@arr);
	    open(AO,">","$temp.txt");
	     
	    
	    
	    #-------VARIABLE DECLARATIONS
	    #----------------Scalars
	    #----------------Arrays
	    my(@filecontents);
	    #------------------@filecontents_______ for storing the input file into and array
	    #------------------This is a better tradeoff than doing thousand iffs
	    
	    
	    #--------BODY
	    print AO"#  Matrix made by MatrixMaker for $genus $species($ncbiid)\n";
	    print AO"#  Type:$temp Cluster Percentage: >= 50\n";
	    print AO"#  Parameters Given Max:$max Min:$min Total Number of Clusters:$numberofclusers\n";
	    print AO"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n";
	    #taking file into an array
	    @filecontents = <AI>;
	    #removing the first line of amino acids
	    shift(@filecontents);
	    #iterating over each line of the file
	    foreach my $fileline(@filecontents)
	    {
	        #taking control of newline in my hands
	        chomp($fileline);
	        #iterating over elements of each line
	        foreach my $indelement(split(/\t/,$fileline))
	        {
	            if($indelement =~/\d/)
	            {
	                if(($indelement<0&&$indelement>-10)||$indelement>=10)
	                {
	                    print AO " $indelement";
	                }
	                elsif($indelement<0&&$indelement<=-10)
	                {
	                    print AO "$indelement";
	                }
	                else
	                {
	                    print AO "  $indelement";
	                }
	            }
	            else
	            {
	                print AO "$indelement";
	            }
	        }
	        print AO " -5\n";
	    }
	    print AO "* -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1\n";
	    close(AO);
	    $counter++;
	    }
	    
	     }
	    
	    
	    
	    
	    
	    
	    chdir("$originalcwd");
	    my $google = strftime('%d/%b/%Y-%H:%M',localtime );
	    my $resultzip =`7za a -tzip results.zip "Results"`;
	    
	    if($resultzip =~/Everything is Ok/)
	    {
	        `rd /s /q "Results"`;
	        my $label = $mw1 -> Label(-text=>"Huff...\nYour Results are Here\nclick Below to retireve\n") -> pack();
	        my $button = $mw1 -> Button(-text => "Submit",
	    		-command =>sub {
	                        `results.zip`;
	                        exit;
	    		    })
	    	-> pack( -side=>'bottom');
}
	    }
                exit;})
-> pack( -side=>'bottom');

MainLoop;





sub exitProgam {
    my $stuff = shift;
 	$mw->messageBox(-message=>"$stuff");
	exit;
}
sub Help
{
    `help.txt`;
}

sub About
{
    `about.txt`;
}