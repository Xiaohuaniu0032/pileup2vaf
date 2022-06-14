use strict;
use warnings;

my ($infile,$cutoff,$outfile) = @ARGV;

open O, ">$outfile" or die;

open IN, "$infile" or die;
my $header = <IN>;
print O "$header";
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $var_freq = $arr[-1];
	if ($var_freq >= $cutoff){
		print O "$_\n";
	}
}
close IN;
close O;
