require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
use strict;
use warnings;
use LWP::UserAgent;
my $ua = LWP::UserAgent->new();

sub getclusers
{
my($genus,$species,$ncbiid)= @_;
my $resp = $ua->get("http://www.uniprot.org/uniref/?query=identity%3a0.5+AND+taxonomy%3a%22$genus+$species+[$ncbiid]%22+AND+count%3a[10+TO+*]&format=list");
my @idlist=%{$resp};
return $idlist[3];
}
1;


