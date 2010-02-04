
# script to load IL lines genotype data. 
# Lukas, 2006-12-09
#

use strict;
use Getopt::Std;

use CXGN::DB::InsertDBH;
use CXGN::Phenome::Genotype;
use CXGN::Phenome::PolymorphicFragment;
use CXGN::Marker::Tools;

=head1 NAME

load_mapping_data.pl - loads mapping data for a population in the phenome database from a mapmaker file.

=head1 USAGE

load_mapping_data.pl [-s -t ] -p 9  mapmaker_inputfile

=head1 DESCRIPTION

Command line options:

=over 5

=item -s

create the schema and exit.

=item -t

trial mode. Do not make any changes to the database. Good for debugging the program!

=item -p

the population_id for which we want to load these data.

=back

=head1 AUTHOR(S)

Lukas Mueller (lam87@cornell.edu)

=cut

our ($opt_s, $opt_t, $opt_p);

getopts('stp');



my $mode = shift;
my $file = shift;

my $MAP_ID = 5; # the map id for the map that was used to obtain genome data.
my $SP_PERSON_ID=222; # the owner of these data.
my $ZYGOCITY = "homozygous"; # either homozygous or heterozygous.
my $BACKGROUND_ACCESSION_ID = 100; 
my $TYPE = "inbred"; # either "inbred" or "mutant". used to adapt display.
if (!$file && !$opt_s) { 
    usage();
}

# mapmaker codes
# 
# 1 = parent1 homozygous
# 2 = heterozygous
# 3 = parent2 homozygous
# 4 = parent1, zygosity unknown (1 or 2)
# 5 = paretn2, zygosity unknown (3 or 2)


my $dbhost = "";
my $dbname = "";
my $dbschema = "";

if ($mode eq "cxgn") { 
    $dbhost = "hyoscine";
    $dbname = "cxgn";
    $dbschema = "public";
}
elsif ($mode eq "sandbox" || $mode eq "trial") { 
    $dbhost = "scopolamine";
    $dbname = "sandbox";
    $dbschema = "public";
}
else { 
    print "\nERROR! NEED A DATABASE NAME! (either sandbox or cxgn)\n";
    usage();
}

print "PLEASE VERIFY: Loading will occur to dbhost: $dbhost. DB name: $dbname. \n";
print "CONTINUE? ";

my $a = (<STDIN>);
if ($a !~ /[yY]/) { exit(); }

my $dbh = CXGN::DB::InsertDBH::connect( {
    dbhost =>$dbhost,
    dbname =>$dbname,
    dbschema => "",
});

if ($opt_s) { 
    eval { 
	my $g = CXGN::Phenome::Genotype->new($dbh);
	$g->create_schema();
	my $p = CXGN::Phenome::PolymorphicFragment->new($dbh);
	$p->create_schema();
    };
    if ($@ || $opt_t) { 
	print "An error occurred or trial mode (-t): $@.\n"; 
	exit(-1);
    }
    else { 
	$dbh->commit();
    }
    exit();
}

open (F, $file) || die "Can't open $file... :( \n";

my $IL_name ="";
my $old_IL_name = "";

my $marker_query = "SELECT marker_id FROM marker_alias WHERE alias=?";
my $marker_query_h = $dbh->prepare($marker_query);

eval { 
    while (<F>) { 
	chomp;
	my ($marker_name, $individual_scores) = split /\s+/;

	# put the individual_scores in an array - which is zero-based,
	# thus we need to add 1 to get the corresponding name of the 
        # individual!
	#
	my @individual_scores = split //, $individual_scores;

	# remove leading asterisk that is sometimes present.
	#
	my $marker_name =~ s/\*(.*)/$1/;

	

	print STDERR "READ: $marker_name, $individual_scores\n";

	$marker = CXGN::Marker::Tools::clean_marker_name($marker_name);

	$marker_query_h->execute($marker_name);
	
	my ($marker_id) = $marker_query_h->fetchrow_array();

	if (!$marker_id) { 
	    print STDERR "Marker $marker was not found. \n";
	    next();
	}

	

	my $genotype;

	foreach my $n ($IL_name, $heterozygous_name) { 
	    if ($n =~ /H/) { $ZYGOCITY="heterozygous"; }
	    else { $ZYGOCITY="homozygous"; }
	    # find individual with name IL_name
	    my $q = "SELECT individual_id FROM phenome.individual WHERE name=?";
	    my $h = $dbh->prepare($q);
	    $h->execute($n);
	    	    
	    my ($individual_id)=$h->fetchrow_array();
	    if (!$individual_id) { 
		die "$n not found in the database! It's like so over!\n";
	    }

	    $genotype = CXGN::Phenome::Genotype->new($dbh);
	    $genotype->set_experiment_name("IL-lines mapped using ExPEN1992 map");
	    $genotype->set_reference_map_id($MAP_ID);
	    $genotype->set_background_accession_id($BACKGROUND_ACCESSION_ID);
	    $genotype->set_preferred(1);
	    $genotype->set_sp_person_id($SP_PERSON_ID);
	    $genotype->set_individual_id($individual_id);
	    $genotype_id = $genotype->store();
	
	
	    # get the marker_ids
	    print STDERR "Getting the marker ids...\n";
	    my $mq = "SELECT marker_id FROM sgn.marker_alias 
                  JOIN sgn.marker_experiment using (marker_id)
                  JOIN sgn.marker_location using (location_id) 
                  JOIN sgn.map_version using (map_version_id)
                  WHERE alias=? and map_id=?";
	    my $mh = $dbh->prepare($mq);
	    print STDERR "Marker1...($marker1, $MAP_ID)\n";
	    $mh -> execute($marker1, $MAP_ID);
	    my ($marker1_id) = $mh->fetchrow_array();
	    print STDERR "Marker2...($marker2, $MAP_ID)\n";
	    $mh -> execute($marker2, $MAP_ID);
	    my ($marker2_id) = $mh->fetchrow_array();
	    
	    if (!$marker1_id || !$marker2_id) { die "marker $marker1 or $marker2 not found in db"; }
	    
	    my $linkage_group = $n;
	    $linkage_group =~ s/[A-Za-z]+(\d+?)\-.*/$1/i;
	    
	    # insert the polymorphic_fragment data
	    #
	    print STDERR "inserting polymorphic fragment $marker1_id, $marker2_id\n";
	    my $polyfrag = CXGN::Phenome::PolymorphicFragment->new($dbh);
	    $polyfrag->set_zygocity($ZYGOCITY);
	    $polyfrag->set_flanking_marker1_id($marker1_id);
	    $polyfrag->set_flanking_marker2_id($marker2_id);
	    $polyfrag->set_linkage_group($linkage_group);
	    $polyfrag->set_type($TYPE);
	    $polyfrag->set_genotype_id($genotype_id);
	    $polyfrag->store();
	    
	    #$old_IL_name= $IL_name;
	
	}
    }
};

if ($@ || $opt_t ) { 
    $dbh->rollback();
    print "Everything has been rolled back due to an error or because of trial mode (-t). ($@)\n";
}
else { 
   $dbh->commit();
   print STDERR "Committed transaction!  :-) \n";
}

print STDERR "We're done! Yeah!\n";


sub usage { 
    print <<USAGE;

Usage: load_IL_lines.pl {-t} {-s} [ sandbox | cxgn ] file
  options:    -t trial. rolls back any transactions initiated.
	      -s create necessary database schemas and then exit.
  parameters: database [ sandbox or cxgn ]
              file [ tab delimited: ( chromosome | IL-name | marker1 | marker2 )]
	  
USAGE

    exit();
}
