
=head1 NAME

loci_compilation.pl

=head1 SYNOPSIS

This perl script extracts the GenBank mRNA sequences and unigene sequences from the SGN database and stores it in a FASTA file. It also creates an error file that stores the dbxref_ids which have no corresponding features from the feature table.

=head1 AUTHOR

Jean Hooi Lee <jhl233@cornell.edu>

=cut

use strict;
use warnings;

use Bio::PrimarySeq;
use Bio::SeqIO;

use CXGN::Transcript::Unigene;
use SGN::Context;

#---create error file to store all errors along the way---

open ERROR, ">loci_error_file.txt";

#---connect to schemas - i.e. public, phenome, and sgn schemas---

my $c = SGN::Context->new();

my $public_schema  = $c->dbic_schema('Bio::Chado::Schema');
my $phenome_schema = $c->dbic_schema('CXGN::Phenome::Schema');
my $sgn_schema     = $c->dbic_schema('SGN::Schema');


#---find all GenBank loci---

### find the dbxrefs which are linked to DB:GenBank_GI

my @dbxref = $public_schema->resultset('General::Db')
                           ->find({name => 'DB:GenBank_GI'})
                           ->search_related('dbxrefs', {});

my @dbxref_ids;
foreach(@dbxref){
    push (@dbxref_ids, $_->dbxref_id);
}

### find common_name_id for Tomato

my $common_name = $sgn_schema->resultset('CommonName')
                             ->find({common_name => 'Tomato'});

### find tomato locus_dbxrefs

my @locus_dbxref = $phenome_schema->resultset('Locus')
                                  ->search({common_name_id => $common_name->common_name_id})
                                  ->search_related('locus_dbxrefs', {dbxref_id => {'in' => [@dbxref_ids]}});

### remove obsolete locus_dbxrefs

for(my $i = 0; $i <= $#locus_dbxref; $i++){
   if($locus_dbxref[$i]->obsolete){
       splice(@locus_dbxref, $i, 1);
       $i--;  #prevents loop from skipping indexes after previous index deleted
   }
}

### remove obsolete loci

for(my $i = 0; $i <= $#locus_dbxref; $i++){
    my $obsolete_locus = $phenome_schema->resultset('Locus')
	                                ->find(locus_id => $locus_dbxref[$i]->get_column('locus_id'))
					->get_column('obsolete');

    if($obsolete_locus == 1){
	print $i.") locus_id: ".$locus_dbxref[$i]->get_column('locus_id')."\n";
	splice(@locus_dbxref, $i, 1);
	$i--;
    }
}

### make sure features exist, otherwise remove dbxref_id

print ERROR "dbxref_ids with no features: \n";
print ERROR "index) dbxref_id \n";

for(my $i = 0; $i <= $#locus_dbxref; $i++){
    my $feature = $public_schema->resultset('Sequence::Feature')
	                        ->find({dbxref_id => $locus_dbxref[$i]->dbxref_id});

    if(! defined $feature){
	print ERROR $i.") ".$locus_dbxref[$i]->dbxref_id."\n";
	splice(@locus_dbxref, $i, 1);
	$i--;
    }
}

print ERROR "\n\n";

### get list of loci and features

my @loci;
my @features;
for(my $i = 0; $i <= $#locus_dbxref; $i++){
    $loci[$i] = $phenome_schema->resultset('Locus')
	                       ->find({locus_id => $locus_dbxref[$i]->get_column('locus_id')});

    $features[$i] = $public_schema->resultset('Sequence::Feature')
	                          ->find({dbxref_id => $locus_dbxref[$i]->dbxref_id});
}

#---find all unigenes---

### find tomato locus_unigenes

my @locus_unigene = $phenome_schema->resultset('Locus')
                                   ->search({common_name_id => $common_name->common_name_id})
                                   ->search_related('locus_unigenes', {});

### remove obsolete locus_unigenes

for(my $i = 0; $i <= $#locus_unigene; $i++){
    if($locus_unigene[$i]->obsolete){
	splice(@locus_unigene, $i, 1);
	$i--;
    }
}

### remove obsolete loci related to unigenes

for(my $i = 0; $i <= $#locus_unigene; $i++){
    my $obsolete_locus = $phenome_schema->resultset('Locus')
	                                ->find(locus_id => $locus_unigene[$i]->get_column('locus_id'))
					->get_column('obsolete');

    if($obsolete_locus == 1){
	splice(@locus_unigene, $i, 1);
	$i--;
    }
}

### find the current build for tomatoes

my $current_build = $sgn_schema->resultset('Group')
                               ->find({comment => {'like' => '%Tomato%'}})
                               ->find_related('unigene_builds', {status => 'C'})
                               ->get_column('unigene_build_id');

### keep locus_unigenes with the latest tomato build
### (i.e. delete the other locus_unigenes with lower builds)

for(my $i = 0; $i <= $#locus_unigene; $i++){
    my $unigene = $sgn_schema->resultset('Unigene')
	                     ->find({unigene_id => $locus_unigene[$i]->unigene_id});

    if($unigene->unigene_build_id != $current_build){
	splice(@locus_unigene, $i, 1);
	$i--;
    }
}

### get loci information

my @unigene_loci;
for(my $i = 0; $i <= $#locus_unigene; $i++){
    $unigene_loci[$i] = $phenome_schema->resultset('Locus')
	                               ->find({locus_id => $locus_unigene[$i]->get_column('locus_id')});
}

### find unigene sequences

my @unigene_seq;
for(my $i = 0; $i <= $#locus_unigene; $i++){
    my $id = $locus_unigene[$i]->unigene_id;
    my $seq = CXGN::Transcript::Unigene->new($sgn_schema, $id);
    $unigene_seq[$i] = $seq->get_sequence();
}


#---Change to fasta format and write loci and unigenes to file---

my $output = Bio::SeqIO->new(-file => ">loci_sequences.fasta",
			     -format => 'Fasta');

### write loci to file

for(my $i = 0; $i <= $#locus_dbxref; $i++){

    my $locus = new Bio::PrimarySeq->
	new(-seq          => $features[$i]->residues,
	    -display_id   => "SGN-L".$loci[$i]->locus_id."-".
	                     $features[$i]->uniquename."|",     	   
	    -description  => $features[$i]->uniquename."|".
	                     $loci[$i]->locus_symbol."|".
	                     $loci[$i]->locus_name);

    $output->write_seq($locus);
}

### write unigenes to file

for(my $i = 0; $i <= $#locus_unigene; $i++){

    my $unigene = new Bio::PrimarySeq->
	new(-seq          => $unigene_seq[$i],
	    -display_id   => "SGN-U".$locus_unigene[$i]->unigene_id."-".
	                     $unigene_loci[$i]->locus_id."|",
	    -description  => $unigene_loci[$i]->locus_symbol."|".
	                     $unigene_loci[$i]->locus_name);

    $output->write_seq($unigene);
}

$output->close();

print "Output writing to loci_sequences.fasta complete.\n";

close ERROR;
