
=head1 NAME

contig_loci2_matching.pl

=head1 SYNOPSIS

This creates an easy-to-read file that determines which GenBank mRNA and/or unigene sequences were used to calculate each tentative consensus sequence.

=head1 AUTHOR

Jean Hooi Lee

=cut

use strict;
use warnings;

use File::Slurp;
use SGN::Context;

#---connect to database---

my $c = SGN::Context->new();
my $phenome_schema = $c->dbic_schema('CXGN::Phenome::Schema');

#---input and output file names to be used---

my $input  = 'loci2_info_contigreadlist.txt';

my $output = 'contig_loci2_matches2.txt';

#---read in loci_info_contigreadlist.txt file---

my @lines = read_file($input);


#---header info---
write_file($output, 
"#
# This is a list of tentative consensus sequences
# and their component sequences for tomato loci.
# ");


#---determine which loci were used to create each contig--- 

my $start_contig = "";

for(my $i = 0; $i <= $#lines; $i++){

### do not print out any lines that are commented in the text file
    if (!($lines[$i] =~ /#/)){
	
### print out contig that is currently being evaluated
	my $current_contig;
	if($lines[$i] =~ /(\w+)/){
	    $current_contig = $1;
	}
	if($current_contig ne $start_contig){
	    $start_contig = $current_contig;
	    append_file($output, "\n".$start_contig."\n");
	}
	
#---determine locus_ids---
	my $locus_id;
	
### check if it is a locus - if so, locus_id = digits immediately after "SGN-L"
### otherwise, must be a unigene - locus_id found at end of name
	if ($lines[$i] =~ /SGN-L(\d+)/ || $lines[$i] =~ /SGN-U\d+-(\d+)/){
	    $locus_id = $1;
	    append_file($output, "locus_id: ".$locus_id."\t\t");
	}
	
### print out the full name
	if ($lines[$i] =~ /(SGN-L[\w-]+)/ || $lines[$i] =~ /(SGN-U[\w-]+)/){
	    append_file($output, "name: ".$1."    \t");
	}
	
### print out the locus_name
	my $locus_name = $phenome_schema->resultset('Locus')
	                                ->find(locus_id => $locus_id)
	                                ->get_column('locus_name');
	
	append_file($output, "locus_name: ".$locus_name."\n");
    }
}

print "Output writing complete.\n";
