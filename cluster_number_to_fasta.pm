require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

use strict;
use warnings;
use LWP::UserAgent;
my $ua = LWP::UserAgent->new();

sub getfasta
{
   my $unirefid = shift;
   my $genus = shift;
   my $speices = shift;
   my $max = shift;
   my $min = shift;
   my @restrictions = @_;
   my $finalcount=0;   
   my $resp = $ua->get("http://www.uniprot.org/uniprot/?query=cluster%3a%28$unirefid%29&format=txt");
   my @temparr=%{$resp};
   my @temparr2 = split(/\n/,$temparr[3]);
   my $ossil=1;
   my $ossil2=0;
   my ($ocline,$osline,$sequence);
   my $multiplefasta;
   my %hash;
   my %hash2;
   my $filelength = scalar(@temparr2);
   my $checker=0;
   my $prints=0;
   my %hash3;
   
   foreach(@temparr2)
   {
       $checker++;
       print $filelength."-----".$checker."\n";
           if($checker==$filelength)
           {
               if(exists $hash{"$genus $speices"})
               {
                  $hash3{"$genus $speices"}=$hash{"$genus $speices"};
                  delete $hash{"$genus $speices"};
                  delete $hash2{"$genus $speices"};
               if(scalar(keys %hash2)>$min)
               {
                   my %storagehash;
                   my $level=0;
                   
                   my $comparer=0;
                   my $letterscore=0;
                   while($comparer==0)
                   {
                   foreach(keys %hash2)
                   {
                      my @temparr = split(/;/,$hash2{$_});
                      my ($one,$two) = split("\\s+",$_);
                      if(defined $temparr[$level])
                      {
                      if($temparr[$level]!~/$one/)
                      {
                      $letterscore=1;
                     if(exists $storagehash{$temparr[$level]})
                      {
                       
                      }
                       else
                      {
                       $hash3{$_}++;
                       $storagehash{$temparr[$level]}++;
                       
                      }
                      }
                     }
                  if(scalar(keys %hash3)>=$max)
                  {
                     $letterscore=0;
                     last;
                  }
                   }
               
                     if($letterscore==1)
                     {
                     $level++;
                     $letterscore=0;
                     }
                    
                    else
                    {
                     if(scalar(keys %hash3)>=$min)
                     {
                     $prints=1;
                     }
                     $comparer=1;
                    }
                      
                   }
               }
           }
                 }
       chomp($_);
      if($_ =~/^\/\//)
       {
           if($ossil2==1)
           {
            $sequence =~s/\s+//g;
            $hash{$osline} = "$sequence";
            $hash2{$osline} = $ocline;
            $ossil2=0;
           }
            $ocline="";
            $osline="";
            $sequence="";
       }
       
       if($ossil==1)
       {
         if($_ =~/^OS\s{3}(.*)/)
         {
            my($one,$two) = split("\\s+",$1);
            if($two eq "")
            {
              $two = "spp";
            }
            $osline ="$one $two";
            $osline =~s/\.//g;
            $ossil=0;
         }
       }
       if($_ =~/^OC\s{3}(.*)/)
       {
         $ocline.=$1;
       }
       if($_ =~/SQ   SEQUENCE/)
       {
         $ossil=1;  
         chomp($_);
         $hash{$osline}++;
         $ossil2=1;
         next; 
       }
            
           
           
       
       if($ossil2==1)
       {
         $sequence.=$_;
       }
     
       
       
       
       
   }
   
   if($prints==1)
   {
      
      $multiplefasta.=">AA$genus $speices\n".$hash3{"$genus $speices"}."\n\n";
      delete $hash3{"$genus $speices"};
      my $gmail=0;
      foreach(keys %hash3)
      {
         $gmail++;
         $multiplefasta.=">$_\n".$hash{$_}."\n\n";
         
      }
   
     return $multiplefasta; 
   }
   else
   {
      my $returnvalue = "Google";
      return $returnvalue;
   }

}

1;