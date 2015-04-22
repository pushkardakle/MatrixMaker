require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
use WWW::Mechanize;
use LWP::UserAgent;
my $mech = WWW::Mechanize->new();
my $ua = LWP::UserAgent->new();


sub blocksretrieve
{
my $multiplefasta = shift;
my $id;
$mech->get("http://blocks.fhcrc.org/blocks/make_blocks.html");
$mech->set_fields(desc =>'gmail',sequences =>"$multiplefasta");
$mech->submit();

if($mech->content =~/\?(.*)\"/)
{
   $id=$1;
}

my $args ="../tmp/bm/$id/$id.mblks";
my $resp = $ua->get("http://blocks.fhcrc.org/blocks-bin/catfile.sh?../tmp/bm/$id/$id.mblks");
my @temparr=%{$resp};
$temparr[3]=~s/(<PRE>\n|<\/PRE>\n)//g;
my @returnarray = split("//\n",$temparr[3]);
return @returnarray;
}
1;