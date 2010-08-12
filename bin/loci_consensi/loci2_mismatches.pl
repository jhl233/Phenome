
=head1 NAME

loci2_mismatches.pl

=head1 SYNOPSIS

This determines which consensus sequences have errors (i.e. where one sequence matches more than one locus) and prints them to a new file. It can also create another file in a tabular format, which can be easily imported into Excel.

=head1 AUTHOR

Jean Hooi Lee <jhl233@cornell.edu>

=cut

use strict;
use warnings;

use File::Slurp;

#---input and output files to be used---

my $input = 'contig_loci2_matches.txt';

my $output = 'loci2_mismatches.txt';

my $output2 = 'loci2_mismatches_tab.txt';

#---read input file---

my @lines = read_file($input);

#---information---

write_file($output,
"#
#
# This file shows consensi that have more than one tomato sequence associated with them
# (which is a problem, because each sequence should match only one locus). This could
# be due to a mistake in the database or because one of the sequences was actually an
# allele. Please see the bottom of the file to see the number of consensi affected.
#
#");

# write_file($output2,
# "#
# #
# # This file shows consensi that have more than one tomato sequence associated with them
# # (which is a problem, because each sequence should match only one locus). This could
# # be due to a mistake in the database or because one of the sequences was actually an
# # allele. Please see the bottom of the file to see the number of consensi affected.
# #
# #\n");

# append_file($output2, "consensi_id\t locus_id\t       name\t\t     locus_name\n");

#---Determine if all loci in each contig are the same - if not, print to the output file---

my $consensus;
my @loci_list;
my $num_of_indices = 0;
my $errors = 0;

for(my $i = 0; $i <= $#lines; $i++){

### ignore commented lines
    if (!($lines[$i] =~ /#/)){

	### read in consensus id
	if($lines[$i] =~ /(\w+\d_\w+_\w\d+)/ || $lines[$i] =~ /(\w+\d_\w\d+)/){
	    $consensus = $1;
	} 
	
	### read in locus_id
	elsif($lines[$i] =~ /locus_id: (\d+)/){
		push(@loci_list, $1);
		$num_of_indices++;
	}
	
	### compare the locus_ids, determine if they should be printed to the file, then reset variables
	elsif($lines[$i] =~ /(\s+)/){
	    my $start_index = $i - $num_of_indices;
	    my $end_index = $i - 1;
	    
	    for(my $j = 0; $j <= $#loci_list; $j++){
		if($loci_list[0] != $loci_list[$j]){
		    $errors++;
		    
		    #########################################################################################
		    #---USED TO MAKE FIRST OUTPUT - COMMENT OUT IF SECOND OUTPUT DESIRED---
		    
		    append_file($output, "\n".$consensus."\n");
  	            append_file($output, @lines[$start_index..$end_index]);
		    
		    #########################################################################################
		    #---USED TO MAKE SECOND OUTPUT IN TABULAR FORMAT - COMMENT OUT IF FIRST OUTPUT DESIRED---
		    
		   #  ###parse lines to only get data
# 		    my ($locus_id, $name, $locus_name);
		    
# 		    for(my $k = $start_index; $k <= $end_index; $k++){
			
#        			### add locus_id
# 			if ($lines[$k] =~ /locus_id: (\d+)/){
# 			    $locus_id = $1;
# 			}
			
# 			### store name
# 			if($lines[$k] =~ /(SGN-L[\w-]+)/ || $lines[$k] =~ /(SGN-U[\w-]+)/){
# 			    $name = $1;
# 			}
			
# 			###store locus_name
# 			if($lines[$k] =~ /locus_name: ([\d\D^\n]+)/){
# 			    $locus_name = $1;
# 			}
# 			append_file($output2, $consensus."\t   ".$locus_id."\t\t".$name."   \t".$locus_name);
# 		    }
		    ########################################################################################
		    
		    last;
		}
	    }
	    @loci_list = ();
	    $num_of_indices = 0;
	}
    }
}
#---show number of consensi with more than one loci---

append_file($output , "#\n"."#\n"."# The number of consensi printed above: ".$errors."\n#\n#");
#append_file($output2, "#\n"."#\n"."# The number of consensi printed above: ".$errors."\n#\n#");

