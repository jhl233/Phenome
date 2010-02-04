


=head1 NAME

CXGN::Phenome::Population - a class that deals with populations, such as mapping populations, mutant populations, inbred lines etc in the SGN database.

=head1 DESCRIPTION

This class inherits from CXGN::DB::ModifiableI and can therefore be used to implement user-modifiable pages easily. It also inherits (wow, multiple inheritance) from CXGN::Phenome::Main, which handles the database connection.

=head1 AUTHOR(S)

Lukas Mueller (lam87@cornell.edu), Isaak Tecle (iyt2@cornell.edu), Naama Menda<nm249@cornell.edu>

=head1 FUNCTIONS

This class implements the following functions:

=cut

package CXGN::Phenome::Population;

use strict;
use CXGN::DB::Connection;
use CXGN::People::Person;
use CXGN::Chado::Cvterm;
use CXGN::Phenome::PopulationDbxref;
use List::Compare;
use Cache::File;


use base qw / CXGN::DB::ModifiableI  /;


=head2 function new

  Synopsis:	my $p = CXGN::Phenome::Population->new($dbh, $population_id)
  Arguments:	a database handle and a population id
  Returns:	a population object
  Side effects:	if $population_id is omitted, an empty project is created. 
                if an illegal $population_id is supplied, undef is returned.
  Description:	

=cut

sub new {
    my $class = shift;
    my $dbh = shift;
    my $population_id = shift;
    
    if (!$dbh->isa("CXGN::DB::Connection")) { 
	die "First argument to CXGN::Phenome::Population constructor needs to be a database handle.";
    }
    my $self=$class->SUPER::new($dbh);


    if ($population_id) { 
	$self->set_population_id($population_id);
	$population_id = $self->fetch();
	if (!$population_id) { return undef; }
    }
    return $self;
    
}

=head2 function new_with_name

  Synopsis:	my $population = CXGN::Phenome::Population->new_with_name($dbh, $population_name)
  Arguments:	a database handle and a population name
  Returns:	a population object
  Side effects:	if a non - existing  $population_name is supplied, undef is returned.
  Description:	

=cut

sub new_with_name {
    my $self = shift;
    my $dbh=shift;
    my $name=shift;
    my $population=undef;
    my $query = "SELECT population_id FROM phenome.population WHERE name ilike ?";
    my $sth=$dbh->prepare($query);
    $sth->execute($name);
    my $pop_id = $sth->fetchrow_array;
    if ($pop_id) {   $population= $self->new($dbh,$pop_id); }
    return $population;
}


sub fetch { 
    my $self= shift;
    my $query = "SELECT population.population_id, population.name, population.description, population.background_accession_id, population.sp_person_id, population.create_date, population.modified_date, population.obsolete FROM phenome.population 
LEFT JOIN sgn.accession ON (population.background_accession_id= sgn.accession.accession_id)
LEFT JOIN phenome.individual USING (population_id)
LEFT JOIN sgn.common_name ON (individual.common_name_id = common_name.common_name_id)
WHERE population_id=? and population.obsolete='f'";
  
    my $sth = $self->get_dbh()->prepare($query);
    $sth->execute($self->get_population_id());

    my ($population_id, $name, $description, $background_accession_id, $sp_person_id, $create_date, $modified_date, $obsolete) = $sth->fetchrow_array();

    $self->set_population_id($population_id);
    $self->set_name($name);
    $self->set_description($description);
    $self->set_background_accession_id($background_accession_id);
    $self->set_sp_person_id($sp_person_id);
    $self->set_create_date($create_date);
    $self->set_modification_date($modified_date);
    $self->set_obsolete($obsolete);
    
    return $population_id;

}

=head2 function store

  Synopsis: $self->store()
  Arguments:	none
  Returns:	database id
  Side effects:	Update an existing population
  Description:	

=cut

sub store {
    my $self = shift;
    if ($self->get_population_id()) { 
	my $query = "UPDATE phenome.population SET
                       name = ?,
                       description = ?,
                       background_accession_id=?,
                       sp_person_id=?,
                       modified_date = now()
                     WHERE
                       population_id = ?
                     ";
	my $sth = $self->get_dbh()->prepare($query);
	$sth->execute($self->get_name(), $self->get_description(), $self->get_background_accession_id(), $self->get_sp_person_id(), $self->get_population_id());
	return $self->get_population_id();
    }
    else { 
	my $query = "INSERT INTO phenome.population
                      (name, description, background_accession_id, sp_person_id, create_date)
                     VALUES 
                      (?, ?, ?, ?, now())";
	my $sth = $self->get_dbh()->prepare($query);
	$sth->execute($self->get_name(), $self->get_description(), $self->get_background_accession_id(), 
                      $self->get_sp_person_id());
	return $self->get_dbh()->last_insert_id("population", "phenome");
    }
}


=head2 get_all_populations

 Usage:        my $names_ref = CXGN::Phenome::Population::get_all_populations($dbh);
 Desc:         This is a static function. Selects the distinct population names and their IDs from phenome.population.
               Useful for populating a unique drop-down menu with only the population names that exist in the table.
 Ret:          Returns a 2D arrayrefs with all the
               population names  and their corresponding database ids
 Args:         a database handle
 Side Effects: none
 Example:

=cut

sub get_all_populations {
    my $dbh= shift;
    my $query = "SELECT distinct(name), population_id FROM phenome.population 
                    WHERE obsolete = 'f'";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my @aoa = ();
    while (my($name, $id) = $sth->fetchrow_array()) { 
	push @aoa, [$id, $name];
    }
    return \@aoa;
}




=head2 function get_individuals

  Synopsis:  $self->get_individuals()
  Arguments: none
  Returns:	a list of individual objects
  Side effects:	none
  Description:	Finds all individuals of this population

=cut

sub get_individuals {
    my $self = shift;
    
    my $query = "SELECT individual_id FROM phenome.individual WHERE population_id=? ORDER BY individual_id";
    my $sth = $self->get_dbh()->prepare($query);
    $sth->execute($self->get_population_id());
    my $individual;
    my @individuals =();
    while (my ($individual_id) = $sth->fetchrow_array()) { 
	$individual = CXGN::Phenome::Individual->new($self->get_dbh(), $individual_id);
	push @individuals, $individual;
    }
    return @individuals;
}

=head2 function get_all_individual_ids

  Synopsis:	$self->get_all_individual_ids()
  Arguments:	none 
  Returns:	arrayref
  Side effects:	none
  Description:  Get a list of individual ids from this population	

=cut

sub get_all_individual_ids {
    my $self = shift;
    
    my $query = "SELECT individual_id FROM phenome.individual WHERE population_id=? ORDER BY individual_id";
    my $sth = $self->get_dbh()->prepare($query);
    $sth->execute($self->get_population_id());
    my @individual_ids =();
    while (my ($individual_id) = $sth->fetchrow_array()) { 
	push @individual_ids, $individual_id;
    }
    return \@individual_ids;
}


=head2 accessors available in this class (get/set_accessor_name(
    
    population_id
    name
    description
    background_accession_id

=cut

sub get_population_id { 
    my $self=shift;
    return $self->{population_id};
}

sub set_population_id { 
    my $self=shift;
    $self->{population_id}=shift;
}

sub get_name { 
    my $self=shift;
    return $self->{name};
}

sub set_name { 
    my $self=shift;
    $self->{name}=shift;
}

sub get_description { 
    my $self=shift;
    return $self->{description};
}

sub set_description { 
    my $self=shift;
    $self->{description}=shift;
}

sub get_background_accession_id { 
    my $self=shift;
    return $self->{background_accession_id};
}


sub set_background_accession_id { 
    my $self=shift;
    $self->{background_accession_id}=shift;
}


sub get_sp_person_id {
  my $self=shift;
  return $self->{sp_person_id};

}

sub set_sp_person_id {
  my $self=shift;
  $self->{sp_person_id}=shift;
}
 

=head2 get_owners

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

   
sub get_owners {
    my $self=shift;
    my $query = "SELECT sp_person_id FROM phenome.population
                 WHERE population_id = ? AND obsolete = 'f'";
    my $sth=$self->get_dbh()->prepare($query);
    $sth->execute($self->get_population_id());
    my $person;
    my @owners = ();
    while (my ($sp_person_id) = $sth->fetchrow_array()) { 
        $person = CXGN::People::Person->new($self->get_dbh(), $sp_person_id);
	push @owners, $sp_person_id;
    }
    return @owners;
}


=head2 get_pop_data_summary

 Usage: my ($min, $max, $avg, $std, $count) =$pop->get_pop_data_summary($cvterm_id)
 Desc:  returns the minimum, maximum, average, and standard deviation values of phenotype data for a cvterm and number of individuals phenotyped for  a cvterm in a population
 Ret:   $min (| 0.0) , $max ( | 0.0) , $avg, $std, $count
 Args:   cvterm id
 Side Effects: accesses database
 Example:

=cut

sub get_pop_data_summary {
    my $self=shift;
    my $cvterm_id=shift;
    my $population_id=$self->get_population_id();   
    my $query = "SELECT  MIN(cast(value as numeric)), MAX(cast(value as numeric)), ROUND(AVG(cast(value as numeric)), 2), ROUND(STDDEV(cast(value as numeric)), 2), count(distinct individual_id)
                      FROM public.phenotype 
                      LEFT JOIN phenome.individual USING (individual_id)              
                      LEFT JOIN public.cvterm ON (phenotype.observable_id = cvterm.cvterm_id)  
                      WHERE individual.population_id =? AND cvterm_id =? AND cast(value as numeric) is not null
                      GROUP BY population_id";
    my $sth = $self->get_dbh()->prepare($query);
    
    $sth->execute($population_id, $cvterm_id);
   
    my ($min, $max, $ave, $std, $count) =$sth->fetchrow_array() ;
    
    if ($min == 0) {$min = '0.0';}
    if ($max == 0) {$min = '0.0';}
    #if ($ave == 0) {$ave = '0.0';}
    #push @min, $min;
    #push @max, $max;
    #push @ave, $ave;
	

    return $min,  $max, $ave, $std, $count;
}


=head2 get_cvterms

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub get_cvterms {
    my $self=shift;
    my $population_id= $self->get_population_id();
    my $query = "SELECT distinct(observable_id), cvterm.name  FROM public.phenotype 
                 JOIN phenome.individual USING (individual_id) JOIN phenome.population using (population_id) JOIN public.cvterm ON (observable_id = cvterm_id)
                 WHERE population_id = ?
                 ORDER BY cvterm.name";
    my $sth=$self->get_dbh()->prepare($query);
    $sth->execute($population_id);
    my @cvterms;
    

    while (my ($cvterm_id, $cvterm_name) = $sth->fetchrow_array() ) {
	my $cvterm= CXGN::Chado::Cvterm->new($self->get_dbh(), $cvterm_id);

	push @cvterms, $cvterm;

	
    }
    
    return @cvterms;
}


=head2 get_pop_raw_data

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut


sub get_pop_raw_data {
    my $self=shift;
    my $population_id=shift;
   
    my $query = "SELECT individual.population_id, population.name, individual.individual_id, individual.name, observable_id, cvterm.name, cvterm.definition, phenotype.value  
                      FROM public.phenotype 
                      LEFT JOIN phenome.individual USING (individual_id)  
                      LEFT JOIN phenome.population USING (population_id)
                      LEFT JOIN public.cvterm ON (phenotype.observable_id = cvterm.cvterm_id)  
                      WHERE individual.population_id =? 
                      ORDER BY individual.name, cvterm.name";

    my $sth = $self->get_dbh()->prepare($query);
   
    $sth->execute($population_id);
    
    my (@pop_id, @pop_name, @ind_id, @ind_name, @obs_id, @cvterm, @definition, @value);
   
    while (my ($pop_id, $pop_name, $ind_id, $ind_name, $obs_id, $cvterm, $definition, $value) =$sth->fetchrow_array()) {

	if ($value == 0) {$value = '0.0';}
	elsif (!defined($value)) {$value = 'NA';}
	push @pop_id, $pop_id;
	push @pop_name, $pop_name;
	push @ind_id, $ind_id;
	push @ind_name, $ind_name;
	push @obs_id, $obs_id;
	push @cvterm, $cvterm;
        push @definition, $definition;
	push @value, $value;

    }
    return \@pop_id, \@pop_name, \@ind_id, \@ind_name, \@obs_id, \@cvterm, \@definition, \@value;
}

=head2 add_population_dbxref

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub add_population_dbxref {
    my $self=shift;
    my $dbxref=shift; #dbxref object
    my $pop_dbxref_id=shift;
    my $sp_person_id=shift;
    
    my $pop_dbxref=CXGN::Phenome::PopulationDbxref->new($self->get_dbh(), $pop_dbxref_id );
    $pop_dbxref->set_population_id($self->get_population_id() );
    $pop_dbxref->set_dbxref_id($dbxref->get_dbxref_id() );
    $pop_dbxref->set_sp_person_id($sp_person_id);
    $pop_dbxref->store();
    
    return $self->{pop_dbxref} ;
}


=head2 get_population_dbxref

 Usage: $population->get_population_dbxref($dbxref)
 Desc:  access population_dbxref object for a given population and its dbxref object
 Ret:    an PopulationDbxref object
 Args:   dbxref object
 Side Effects:
 Example:

=cut

sub get_population_dbxref {
    my $self=shift;
    my $dbxref=shift; # my dbxref object..
    my $query=$self->get_dbh()->prepare("SELECT population_dbxref_id from phenome.population_dbxref
                                        WHERE population_id=? AND dbxref_id=? ");
    $query->execute($self->get_population_id(), $dbxref->get_dbxref_id() );
    my ($pop_dbxref_id) = $query->fetchrow_array();
    my $pop_dbxref= CXGN::Phenome::PopulationDbxref->new($self->get_dbh(), $pop_dbxref_id);
    return $pop_dbxref;
}


sub get_all_population_dbxrefs {   #get an array of dbxref objects (from public dbxref):
    my $self=shift;
    my @pop_dbxrefs;
    my $dbxref_query=$self->get_dbh()->prepare("SELECT dbxref_id from phenome.population_dbxref JOIN public.dbxref USING (dbxref_id) WHERE population_id=? ORDER BY public.dbxref.accession"); 
    my @dbxref_id;
    $dbxref_query->execute($self->get_population_id() );
    while (my ($d) = $dbxref_query->fetchrow_array() ) {
	push @dbxref_id, $d;
    }
    #an array of dbxref objects
    foreach my $id (@dbxref_id) {
	my $dbxref_obj= CXGN::Chado::Dbxref->new($self->get_dbh(), $id);
	#$self->add_dbxref($dbxref_obj);
	push @pop_dbxrefs, $dbxref_obj;
    }
    return @pop_dbxrefs;
}

sub get_population_publications {
    my $self=shift;
    my $query = $self->get_dbh()->prepare("SELECT pub_id FROM pub_dbxref 
                                           JOIN dbxref USING (dbxref_id)
                                           JOIN phenome.population_dbxref USING (dbxref_id)
                                           WHERE population_id = ?");
    $query->execute($self->get_population_id());
    my $publication;
    my @publications;
    while (my ($pub_id) = $query->fetchrow_array()) { 
	$publication = CXGN::Chado::Publication->new($self->get_dbh(), $pub_id);
	push @publications, $publication;
    }
    return @publications;

}

=head2 get_all_indls_cvterm

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub get_all_indls_cvterm {
    my $self = shift;
    my $pop_id = $self->get_population_id();
    my $observable_id = shift;

    my $query = "SELECT phenotype.individual_id, individual.name, phenotype.value 
                 FROM public.phenotype 
                 LEFT JOIN phenome.individual USING (individual_id)
                 LEFT JOIN public.cvterm ON (observable_id = cvterm_id)
                 WHERE individual.population_id =? AND observable_id = ? 
                 ORDER BY CAST(value AS NUMERIC)";
    
    my $sth = $self->get_dbh()->prepare($query);
    $sth->execute($pop_id, $observable_id);

    my (@indl_id, @indl_name, @value,);
    while (my ($indl_id,  $indl_name, $value ) = $sth->fetchrow_array()) {

	if ($value eq 'null') {$value = 'NA';}
	if ($value == 0) { $value = '0.0';}
	#push @obs_id, $obs_id;
	
	push @indl_name, $indl_name;
	push @value, $value;
	#push @cvterm, $cvterm;
        push @indl_id, $indl_id;
    }
    
  return \@indl_id, \@indl_name, \@value;  
    

}
=head2 plot_cvterm

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub plot_cvterm {
    my $self = shift;
    my $pop_id = $self->get_population_id();
    my $observable_id = shift;

    my $query = "SELECT phenotype.individual_id, individual.name, ROUND(CAST(value AS NUMERIC), 2)
                 FROM public.phenotype 
                 LEFT JOIN phenome.individual USING (individual_id)
                 LEFT JOIN public.cvterm ON (observable_id = cvterm_id)
                 WHERE individual.population_id =? AND observable_id = ?";
    
    my $sth = $self->get_dbh()->prepare($query);
    $sth->execute($pop_id, $observable_id);

    my (@indl_id, @indl_name, @value,);
    while (my ($indl_id,  $indl_name, $value ) = $sth->fetchrow_array()) {

	if ($value eq 'null') {$value = undef;}
	if ($value == 0) { $value = '0.0';}
	#push @obs_id, $obs_id;
	
	push @indl_name, $indl_name;
	push @value, $value;
	#push @cvterm, $cvterm;
        push @indl_id, $indl_id;
    }
    
  return \@indl_id, \@indl_name, \@value;  


}
   

=head2 indls_range_cvterm

 Usage: ($indl_id, $indl_name, $value)= $population->indls_range_cvterm($cvterm_id, $lower, $upper);
 Desc: returns individuals with phenotypic values within a given range 
 Ret: array refs for individual id, individual name, and the corresponding phenotypic values
 Args: cvterm_id, lower and upper limits of the range
 Side Effects:
 Example:

=cut


sub indls_range_cvterm {
my $self=shift;
my $cvterm_id = shift;
my $lower = shift;
my $upper = shift;
my $pop_id = $self->get_population_id();
my $query;
if ($lower == 0) {
    #my $lower_lowest = $lower;

    $query = "SELECT individual_id, individual.name, value 
                     FROM public.phenotype 
                     LEFT JOIN phenome.individual USING (individual_id) 
                     WHERE individual.population_id = $pop_id AND observable_id = $cvterm_id 
                           AND (CAST (value AS NUMERIC) <=$upper) 
                     ORDER BY CAST(value AS NUMERIC)";
} else {
    $query = "SELECT individual_id, individual.name, value 
                     FROM public.phenotype 
                     LEFT JOIN phenome.individual USING (individual_id) 
                     WHERE individual.population_id = $pop_id AND observable_id = $cvterm_id 
                           AND (CAST(value AS NUMERIC) > $lower AND CAST(value AS NUMERIC) <= $upper) 
                     ORDER BY CAST (value AS NUMERIC)";
}

my $sth = $self->get_dbh()->prepare($query);
$sth->execute();

my ($indl_id, $indl_name, $value);
my (@indl_id, @indl_name, @value);

while (my ($indl_id, $indl_name, $value ) = $sth->fetchrow_array()) {
    if ($value eq 'null') {$value = 'NA';}
    if ($value == 0 ) {$value = '0.0';}

    push @indl_id, $indl_id;
    push @indl_name, $indl_name;
    push @value, $value;
	
    }
    
return \@indl_id,  \@indl_name, \@value;  


}

=head2 get_cvterm_acronyms

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut
sub get_cvterm_acronyms {
    my $self=shift;
    my $population_id=$self->get_population_id();
    
    my $query = "SELECT DISTINCT(observable_id), cvterm.name  
                      FROM public.phenotype
                      LEFT JOIN phenome.individual USING (individual_id)  
                      LEFT JOIN phenome.population USING (population_id)
                      LEFT JOIN public.cvterm ON (phenotype.observable_id = cvterm.cvterm_id)  
                      WHERE individual.population_id =?
                      ORDER BY cvterm.name";

    my $sth = $self->get_dbh()->prepare($query);
   
    $sth->execute($population_id);
    
    my (@cvterms, @cvterm_acronyms);

   
    while (my ($observable_id, $cvterm) =$sth->fetchrow_array()) {
	my @words = split(/\s/, $cvterm);
	my $acronym;
	
	if (scalar(@words)== 1) {
	    my $word =shift(@words);
		my $l = substr($word,0,2,q{}); 
		$acronym .= $l;
	     $acronym = uc($acronym);
	 }  else {
	    foreach my $word (@words) {
		if ($word=~/^\D/){
		    my $l = substr($word,0,1,q{});
		    
		    $acronym .= $l;
		} else {$acronym .= $word;}
		$acronym = uc($acronym);
		$acronym =~/(\w+)/;
		$acronym = $1;
	    }
	    
	}
	push @cvterm_acronyms, $acronym;
	push @cvterms, $cvterm;
    }

    return \@cvterms, \@cvterm_acronyms;
}

=head2 cvterm_acronym

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut
sub cvterm_acronym {
    my $self=shift;
    my $cvterm = shift;
   
   
   
    my @words = split(/\s/, $cvterm);
    my $acronym;
	
    if (scalar(@words)== 1) {
	my $word =shift(@words);
	my $l = substr($word,0,2,q{}); 
	$acronym .= $l;
	$acronym = uc($acronym);
    }  else {
	foreach my $word (@words) {
	    if ($word=~/^\D/){
		my $l = substr($word,0,1,q{}); 
		$acronym .= $l;
	    } else {$acronym .= $word;}
	    $acronym = uc($acronym);
	    $acronym =~/(\w+)/;
	    $acronym = $1;
	}
	   
    }
    return $acronym;

}


=head2 get_genotype_data

 Usage:my ($ind_name, $ind_id, $marker, $marker_id, $map_version, $genotype, $lg_name, $position)=$pop_obj->get_genotype_data()
 Desc:retrieves genetic markers, their genotype values, linkage group and map location for individuals in a population study. 
 Ret: array refs for individual name, individual id, marker alias, marker_id, map_version,  genotype value, linkage name, map position
 Args: none
 Side Effects:
 Example:

=cut
 
sub get_genotype_data {
    my $self = shift;
    my $pop_id = $self->get_population_id();
    my $query = "SELECT individual.name, genotype.individual_id, marker_alias.alias, marker_alias.marker_id, map_version.map_version_id, map_version.map_id, genotype_region.zygocity_code, lg_name, position 
                         FROM phenome.genotype 
                         JOIN phenome.individual ON (genotype.individual_id = individual.individual_id)
                         JOIN phenome.genotype_region ON (genotype.genotype_id = genotype_region.genotype_id) 
                         JOIN sgn.marker_alias ON (marker_alias.marker_id=marker_id_nn) 
                         JOIN sgn.marker_experiment ON (marker_alias.marker_id=marker_experiment.marker_id) 
                         JOIN sgn.marker_location USING(location_id) 
                         JOIN sgn.linkage_group ON (linkage_group.lg_id=marker_location.lg_id) 
                         JOIN sgn.map_version ON (marker_location.map_version_id=map_version.map_version_id) 
                         WHERE genotype_experiment_id = (SELECT DISTINCT(genotype_experiment_id) 
                                    FROM phenome.genotype_experiment 
                                    JOIN phenome.genotype USING (genotype_experiment_id) 
                                    JOIN phenome.individual ON (genotype.individual_id = individual.individual_id) 
                                    WHERE individual.population_id = $pop_id)  
                         AND map_version.map_version_id = (SELECT DISTINCT(map_version_id) 
                                    FROM sgn.map_version 
                                    JOIN phenome.genotype_experiment ON (map_version.map_id = 
                                                        genotype_experiment.reference_map_id) 
                                    JOIN phenome.genotype USING (genotype_experiment_id) 
                                    JOIN phenome.individual USING (individual_id) 
                                    WHERE map_version.current_version = 't' and individual.population_id = $pop_id)
                         AND individual.population_id = $pop_id
                         ORDER BY individual.name, marker_alias.alias";

    my $sth = $self->get_dbh()->prepare($query);
    $sth->execute();
    my (@ind_name, @ind_id, @marker, @marker_id, @map_version, @map_id,  @genotype, @lg_name, @position) =();
   
 
    while (my ($ind_name, $ind_id, $marker, $marker_id, $map_version, $map_id, $genotype, $lg_name, $position) = $sth->fetchrow_array()) {
	
	push  @ind_name, $ind_name;
	push  @ind_id, $ind_id;
	push  @marker, $marker;
	push  @marker_id, $marker_id;	
	push  @map_version, $map_version;
	push  @map_id, $map_id;
	if ($genotype eq "a") {$genotype = 1;}
	elsif ($genotype eq "h") {$genotype = 2;}
	elsif ($genotype eq "b") {$genotype = 3;}
	elsif ($genotype eq "d") {$genotype = 4;}
	elsif ($genotype eq "c") {$genotype = 5;}
	push  @genotype, $genotype;
	push  @lg_name, $lg_name;
	push  @position, $position;
    }
    


    return \@ind_name, \@ind_id, \@marker, \@marker_id, \@map_version, \@map_id, \@genotype, \@lg_name, \@position;

}

=head2 get_all_markers

 Usage:my ($marker_id, $marker_alias)=$pop_obj->get_all_markers()
 Desc: useful for retrieving markers assayed on all individual accessions in a population. Not all markers are genotyped on every individual. 
 Ret: array references for marker ids and aliases.
 Args:none
 Side Effects:
 Example:

=cut



sub get_all_markers {
    my $self = shift;
    my $pop_id = $self->get_population_id();
    my $query = "SELECT DISTINCT(marker_alias.marker_id), marker_alias.alias
                        FROM phenome.genotype 
                        JOIN phenome.individual ON (genotype.individual_id = individual.individual_id)
                        JOIN phenome.genotype_region ON (genotype.genotype_id = genotype_region.genotype_id) 
                        JOIN sgn.marker_alias ON (marker_alias.marker_id=marker_id_nn) 
                        JOIN sgn.marker_experiment ON (marker_alias.marker_id=marker_experiment.marker_id) 
                        JOIN sgn.marker_location USING(location_id) 
                        JOIN sgn.linkage_group ON (linkage_group.lg_id=marker_location.lg_id) 
                        JOIN sgn.map_version ON (marker_location.map_version_id=map_version.map_version_id) 
                        WHERE genotype_experiment_id = (SELECT DISTINCT(genotype_experiment_id) 
                                   FROM phenome.genotype_experiment 
                                   JOIN phenome.genotype USING (genotype_experiment_id) 
                                   JOIN phenome.individual ON (genotype.individual_id = individual.individual_id) 
                                   WHERE individual.population_id = $pop_id)  
                        AND map_version.map_version_id = (SELECT DISTINCT(map_version_id) FROM sgn.map_version 
                                   JOIN phenome.genotype_experiment ON (map_version.map_id = genotype_experiment.reference_map_id)                                      JOIN phenome.genotype USING (genotype_experiment_id) 
                                   JOIN phenome.individual USING (individual_id) 
                                   WHERE map_version.current_version = 't' and individual.population_id = $pop_id)
                        ORDER BY marker_alias.alias";

  
  

                    
    my $sth = $self->get_dbh()->prepare($query);
    
    $sth->execute();
    my (@marker_id, @marker_alias);

    while (my ($marker_id, $marker_alias) =$sth->fetchrow_array()) {

	
        push @marker_id, $marker_id;
	push @marker_alias, $marker_alias;
    }
   
    return \@marker_id, \@marker_alias;
}

=head2 get_genotyped_indls

 Usage:my ($indls_id, $indl_name)=$pop_obj->get_genotyped_indls()
 Desc: useful for retrieving markers all individual accessions genotyped for one or more markers in a population. Not all markers are genotyped on every individual. 
 Ret: array references for individual ids and names.
 Args:none
 Side Effects:
 Example:

=cut



sub get_genotyped_indls {
    my $self = shift;
    my $pop_id = $self->get_population_id();
    my $query = "SELECT DISTINCT(genotype.individual_id), individual.name
                        FROM phenome.genotype 
                        JOIN phenome.individual ON (genotype.individual_id = individual.individual_id)
                        JOIN phenome.genotype_region ON (genotype.genotype_id = genotype_region.genotype_id) 
                        JOIN sgn.marker_alias ON (marker_alias.marker_id=marker_id_nn) 
                        JOIN sgn.marker_experiment ON (marker_alias.marker_id=marker_experiment.marker_id) 
                        JOIN sgn.marker_location USING(location_id) 
                        JOIN sgn.linkage_group ON (linkage_group.lg_id=marker_location.lg_id) 
                        JOIN sgn.map_version ON (marker_location.map_version_id=map_version.map_version_id) 
                        WHERE genotype_experiment_id = (SELECT DISTINCT(genotype_experiment_id) 
                                   FROM phenome.genotype_experiment 
                                   JOIN phenome.genotype USING (genotype_experiment_id) 
                                   JOIN phenome.individual ON (genotype.individual_id = individual.individual_id) 
                                   WHERE individual.population_id = $pop_id)  
                        AND map_version.map_version_id = (SELECT DISTINCT(map_version_id) FROM sgn.map_version 
                                   JOIN phenome.genotype_experiment ON (map_version.map_id = genotype_experiment.reference_map_id)                      
                                   JOIN phenome.genotype USING (genotype_experiment_id) 
                                   JOIN phenome.individual USING (individual_id) 
                                   WHERE map_version.current_version = 't' and individual.population_id = $pop_id)
                        ORDER BY individual.name";

  
  

                    
    my $sth = $self->get_dbh()->prepare($query);
    
    $sth->execute();
    my (@indl_id, @indl_name);

    while (my ($indl_id, $indl_name) =$sth->fetchrow_array()) {

	
        push @indl_id, $indl_id;
	push @indl_name, $indl_name;
    }
   
    return \@indl_id, \@indl_name;
}



sub get_ind_marker_genotype {
    my $self = shift;
    my $pop_id = $self->get_population_id();
    my $ind_id = shift;
    my $marker_id = shift;
    my ( $marker, $map_version, $genotype);
    my $query = "SELECT genotype.individual_id, marker_alias.marker_id, map_version.map_version_id,
                        genotype_region.zygocity_code
                        FROM phenome.genotype 
                        JOIN phenome.individual ON (genotype.individual_id = individual.individual_id)
                        JOIN phenome.genotype_region ON (genotype.genotype_id = genotype_region.genotype_id) 
                        JOIN sgn.marker_alias ON (marker_alias.marker_id=genotype_region.marker_id_nn) 
                        JOIN sgn.marker_experiment ON (marker_alias.marker_id=marker_experiment.marker_id) 
                        JOIN sgn.marker_location USING(location_id) 
                        JOIN sgn.linkage_group ON (linkage_group.lg_id=marker_location.lg_id) 
                        JOIN sgn.map_version ON (marker_location.map_version_id=map_version.map_version_id) 
                        WHERE genotype_experiment_id = (SELECT DISTINCT(genotype_experiment_id) 
                                   FROM phenome.genotype_experiment 
                                   JOIN phenome.genotype USING (genotype_experiment_id) 
                                   JOIN phenome.individual ON (genotype.individual_id = individual.individual_id) 
                                   WHERE individual.population_id = $pop_id)  
                        AND map_version.map_version_id = (SELECT DISTINCT(map_version_id) 
                                   FROM sgn.map_version 
                                   JOIN phenome.genotype_experiment ON (map_version.map_id = genotype_experiment.reference_map_id)                                            JOIN phenome.genotype USING (genotype_experiment_id) 
                                   JOIN phenome.individual USING (individual_id) 
                                   WHERE map_version.current_version = 't' and individual.population_id = $pop_id) 
                       AND genotype.individual_id = $ind_id
                       AND marker_alias.marker_id = $marker_id";


    my  $sth = $self->get_dbh()->prepare($query);
    $sth->execute();



    ($ind_id, $marker, $map_version, $genotype) = $sth->fetchrow_array();
	
#	push @ind_id, $ind_id;
#	push @marker_id, $marker_id;	
#	push @map_version, $map_version;
#	if (!defined($genotype)) {
#	    $genotype = 'NA';
#	}
#	    push @genotype, $genotype;
#	}else {push @genotype, $genotype;}
#	push @lg_name, $lg_name;
#	push @position, $position;
#    }

    return $ind_id, $marker, $map_version, $genotype;

}
# sub get_marker_genotype {
#     my $self = shift;
#     my $pop_id = $self->get_population_id();
#     my $marker_id = shift;
#     my $query = "SELECT genotype.individual_id, marker_alias.marker_id, map_version.map_version_id,
#                         genotype_region.zygocity_code, lg_name, position  
#                         FROM phenome.genotype 
#                         JOIN phenome.individual ON (genotype.individual_id = individual.individual_id)
#                         JOIN phenome.genotype_region ON (genotype.genotype_id = genotype_region.genotype_id) 
#                         JOIN sgn.marker_alias ON (marker_alias.marker_id=marker_id_nn) 
#                         JOIN sgn.marker_experiment ON (marker_alias.marker_id=marker_experiment.marker_id) 
#                         JOIN sgn.marker_location USING(location_id) 
#                         JOIN sgn.linkage_group ON (linkage_group.lg_id=marker_location.lg_id) 
#                         JOIN sgn.map_version ON (marker_location.map_version_id=map_version.map_version_id) 
#                         WHERE genotype_experiment_id = (SELECT DISTINCT(genotype_experiment_id) 
#                                    FROM phenome.genotype_experiment 
#                                    JOIN phenome.genotype USING (genotype_experiment_id) 
#                                    JOIN phenome.individual ON (genotype.individual_id = individual.individual_id) 
#                                    WHERE individual.population_id = $pop_id)  
#                         AND map_version.map_version_id = (SELECT DISTINCT(map_version_id) FROM sgn.map_version 
#                                    JOIN phenome.genotype_experiment ON (map_version.map_id = genotype_experiment.reference_map_id)                                      JOIN phenome.genotype USING (genotype_experiment_id) 
#                                    JOIN phenome.individual USING (individual_id) 
#                                    WHERE map_version.current_version = 't' and individual.population_id = $pop_id) 
#                        AND marker_alias.marker_id = ? 
#                        ";


#     my  $sth = $self->get_dbh()->prepare($query);
#     $sth->execute($marker_id);

#     my (@ind_id, @marker_id, @map_version, @genotype, @lg_name, @position);

#     while (my ($ind_id, $marker,$map_version, $genotype, $lg_name, $position) = $sth->fetchrow_array()) {
	
# 	push @ind_id, $ind_id;
# 	push @marker_id, $marker_id;	
# 	push @map_version, $map_version;
# 	if ($genotype == 'null') {
# 	    $genotype = 'NA';
# 	}
# 	push @genotype, $genotype;
# 	push @lg_name, $lg_name;
# 	push @position, $position;
#     }

#     return \@ind_id, \@marker_id, \@map_version, \@genotype, \@lg_name, \@position;

# }

sub get_marker_chr_position {
    my $self = shift;
    my $pop_id = $self->get_population_id();
    my $marker_id = shift;
    
    my $query = "SELECT marker_alias.alias, map_version.map_version_id, lg_name, position
                        FROM phenome.genotype 
                        JOIN phenome.individual ON (genotype.individual_id = individual.individual_id)
                        JOIN phenome.genotype_region ON (genotype.genotype_id = genotype_region.genotype_id) 
                        JOIN sgn.marker_alias ON (marker_alias.marker_id=marker_id_nn) 
                        JOIN sgn.marker_experiment ON (marker_alias.marker_id=marker_experiment.marker_id) 
                        JOIN sgn.marker_location USING(location_id) 
                        JOIN sgn.linkage_group ON (linkage_group.lg_id=marker_location.lg_id) 
                        JOIN sgn.map_version ON (marker_location.map_version_id=map_version.map_version_id) 
                        WHERE genotype_experiment_id = (SELECT DISTINCT(genotype_experiment_id) 
                                   FROM phenome.genotype_experiment 
                                   JOIN phenome.genotype USING (genotype_experiment_id) 
                                   JOIN phenome.individual ON (genotype.individual_id = individual.individual_id) 
                                   WHERE individual.population_id = $pop_id)  
                        AND map_version.map_version_id = (SELECT DISTINCT(map_version_id) FROM sgn.map_version 
                                   JOIN phenome.genotype_experiment ON 
                                          (map_version.map_id = genotype_experiment.reference_map_id)   
                                   JOIN phenome.genotype USING (genotype_experiment_id) 
                                   JOIN phenome.individual USING (individual_id) 
                                   WHERE map_version.current_version = 't' and individual.population_id = $pop_id) 
                       AND marker_alias.marker_id = ?
                       ORDER BY marker_alias.alias";
                        
                      


    my  $sth = $self->get_dbh()->prepare($query);
    $sth->execute($marker_id);

    my (@marker_name, @map_version, @lg_name, @position);

    while (my ($marker_name, $map_version, $lg_name, $position) = $sth->fetchrow_array()) {
	
	
	push @marker_name, $marker_name;	
	push @map_version, $map_version;	
	push @lg_name, $lg_name;
	push @position, $position;
    }

    return  \@marker_name, \@map_version, \@lg_name, \@position;

}


sub get_marker_position {
    my $self = shift;
    my $map_version_id = shift;
    my $marker_name = shift;
    
    my $query = "SELECT  position FROM sgn.marker_alias LEFT JOIN sgn.marker_experiment USING (marker_id) LEFT JOIN sgn.marker_location USING (location_id) WHERE map_version_id = ?  and marker_alias.alias LIKE ?";
                        
    my  $sth = $self->get_dbh()->prepare($query);
    $sth->execute($map_version_id, $marker_name);

    my ($position) = $sth->fetchrow_array(); 
	

    return  $position;

}


sub phenotype_dataset {
    my $self = shift;
    my $dbh = $self->get_dbh();
    my $pop_id = $self->get_population_id();

    my $phe_dataset;

    my @individuals = $self->get_individuals();
    my $individual_obj = $individuals[1];
    #my $ind_id = $individual_obj->get_individual_id();

    my @cvterms = $individual_obj->get_unique_cvterms( $individual_obj->get_individual_id());
    my ($pop_name, $ind_id, $ind_name, $obs_id, $cvterm, $definition, $value);
    ($pop_id, $pop_name, $ind_id, $ind_name, $obs_id, $cvterm, $definition, $value) = $self->get_pop_raw_data($pop_id);
    
    my ($cvterms2, $cvterm_acronyms) = $self->get_cvterm_acronyms();
    my @cvterms2 = @{$cvterms2};
    my @cvterm_acronyms = @{$cvterm_acronyms};
    
    for (my $i = 0; $i<@cvterms; $i++){
	if  ($cvterms[$i] eq $cvterms2[$i]) {
	 #print "cvterms match \n"; 
	 #do nothing
	} else { print  "There is a mismatch between $cvterms[$i] and $cvterms2[$i]\n";
		   exit();
		}
    }
    $phe_dataset .= "ID" . ",";
    $phe_dataset .= join (",", @cvterm_acronyms);   
    my $old_ind_id = "";
    my @pheno_indls=();
 
    for (my $i=0; $i<@$pop_id; $i++) { 
	
	if ($old_ind_id != $ind_id->[$i]) {
	    if ($old_ind_id) {
		$phe_dataset = substr $phe_dataset, 0, -1;
	    }	
	    $phe_dataset .= "\n$ind_name->[$i]" . ",";		
	}
	    
	foreach my $t (@cvterms) { 
	    my $term = $cvterm->[$i];
	    if ($t =~ /^$term$/i) { 
		$phe_dataset .="$value->[$i]" . ",";
	    }	    
	}
	$old_ind_id=$ind_id->[$i];
	   
    }
    
    $phe_dataset = substr $phe_dataset, 0, -1;

    return \$phe_dataset;

}
=head2 genotype_dataset

 Usage:my ($genodataset)=$pop_obj->genotype_dataset()
 Desc:  returns a genotype dataset for the relevant population (csv file). 
        It is  formatted as per RQTL requirement. 
 Ret: reference
 Args:none
 Side Effects:
 Example:

=cut

sub genotype_dataset {
    my $self = shift;
    my $dbh = $self->get_dbh();
    
    my $pop_id = $self->get_population_id();
    my ($pop_name, $pheno_ind_id, $pheno_ind_name, $obs_id, $cvterm, $definition, $value);
    
    ($pop_id, $pop_name, $pheno_ind_id, $pheno_ind_name, $obs_id, $cvterm, $definition, $value) = $self->get_pop_raw_data($pop_id);
    
    my ($ind_name, $ind_id, $marker, $g_marker_id, $map_version, $map_id, $genotype, $lg_name, $position) = $self->get_genotype_data();
    my ($marker_ids, $marker_alias) = $self->get_all_markers();
    my ($genotyped_indls_id, $genotyped_indls_name) = $self->get_genotyped_indls();
    
    my @genotyped_indls_id = @{$genotyped_indls_id};
    my @genotyped_indls_name = @{$genotyped_indls_name};
    my @pop_marker_ids = @{$marker_ids};
    my @pop_marker_alias = @{$marker_alias};
 #  my @g_marker_id = @{$g_marker_id};
    

    my $gen_dataset = "ID" . ",";

    for (my $i=0; $i<@pop_marker_ids; $i++) {
	$gen_dataset .= "$marker_alias->[$i]" . ",";
    }
    
    $gen_dataset = substr $gen_dataset, 0, -1;
    $gen_dataset .= "\n";
        
    my $m_marker_name;
    foreach my $m_id (@pop_marker_ids) {
	($m_marker_name, $map_version, $lg_name, $position) = $self->get_marker_chr_position($m_id);
 	$gen_dataset .= ",";
 	$gen_dataset .= $lg_name->[1];
    }
    
    $gen_dataset .= "\n";
   
    foreach my $m_id (@pop_marker_ids) {	
 	($m_marker_name, $map_version, $lg_name, $position) = $self->get_marker_chr_position($m_id);
 	$gen_dataset .= ",";		 
 	$gen_dataset .= sprintf('%0.2f',$position->[1]);
	       
    }
    my $old_ind_id = " ";
    my @ind_marker_list=();
    my %marker_genotype={};
    
    my $same_ind_id = $genotyped_indls_id->[0];
   
    my %ind_marker_genotype=();
    my $ind_name_i;
    my @ind_name_list;
    
    for (my $i=0; $i<@$ind_id; $i++) { 	 
	if ($old_ind_id != $ind_id->[$i]) {
	    $ind_name_i= "$ind_name->[$i]";
	    my $ind_name_list = $ind_name->[$i];
	    push @ind_name_list, $ind_name_list;
	}
	my ($ind_markers, $ind_genotype);
	my $j;	
	if ($i > 2) {
	    $j = $i - 1;
	}	
	if ($i == 0)  {
	    %ind_marker_genotype = ($ind_name_i => {$marker->[$i] => $genotype->[$i]});

	}elsif ($ind_name->[$j] eq $ind_name_i) {
	    $ind_marker_genotype{$ind_name_i}->{$marker->[$i]} = $genotype->[$i];
	} elsif ($ind_name->[$j] ne $ind_name_i) { 
	    $ind_marker_genotype{$ind_name_i}={$marker->[$i] => $genotype->[$i]}; 
	}

	$old_ind_id = $ind_id->[$i];
    
    }

    foreach my $genotyped_ind (@genotyped_indls_name) {    
	my @ind_markers = ( sort (keys %{$ind_marker_genotype{$genotyped_ind}}));  
	my $compare = List::Compare->new(\@ind_markers, \@pop_marker_alias);
	my @not_genotyped_m = $compare->get_complement();

	if ( @not_genotyped_m ) {
	    foreach my $not_g_m (@not_genotyped_m) {		    
		$ind_marker_genotype{$genotyped_ind}->{$not_g_m} = 'NA';		   
	    }
	}
    	
    }

    foreach my $ind_i ( sort (keys %ind_marker_genotype)) {
	$gen_dataset .= "\n$ind_i" . ",";
	foreach my $marker_i (@pop_marker_alias) {
	    my $marker_i_genotype = $ind_marker_genotype{$ind_i}{$marker_i};
	    $gen_dataset .= "$marker_i_genotype" . ",";
	}
	$gen_dataset = substr $gen_dataset, 0, -1;	
    } 
  
    return \$gen_dataset;

}


=head2 has_qtl_data

 Usage: my @pop_objs = $pop_obj->has_qtl_data();
 Desc: returns a list of population (objects)  with genetic and phenotypic data (qtl data). 
       The assumption is if a trait genetic and phenotype data, it is from qtl study.
 Ret: an array of population objects
 Args: none
 Side Effects: accesses the database
 Example:

=cut

sub has_qtl_data {
    my $self = shift;    
    my $dbh = CXGN::DB::Connection->new();
    my $query = "SELECT DISTINCT (population_id) FROM public.phenotype LEFT JOIN phenome.individual USING (individual_id)";

    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    my (@pop_objs, @pop_ids)=();
    
    while (my ($pop_id) = $sth->fetchrow_array()) {
	push @pop_ids, $pop_id;
    }

    foreach my $pop_id2 (@pop_ids) {

	my $query2 = "SELECT DISTINCT (population_id) FROM phenome.genotype LEFT JOIN phenome.individual USING (individual_id) WHERE individual.population_id = ?";
	my $sth2 = $dbh->prepare($query2);
	$sth2->execute($pop_id2);
	
	my ($qtl_pop_id) = $sth2->fetchrow_array();
	
	if ($qtl_pop_id) {
	    my $pop_obj = CXGN::Phenome::Population->new($dbh, $qtl_pop_id);
	

	    push @pop_objs, $pop_obj;
	}
    }

        return  @pop_objs; 
 
    

}


=head2 mapversion_id

 Usage: my $map_version_id = $pop->mapversion_id();
 Desc: a quick way to get the map version id of a genetic map for markers genotyped in a population
 Ret: map version id
 Args: none
 Side Effects: accesses database
 Example:

=cut


sub mapversion_id {
     my $self = shift;
     my $dbh = $self->get_dbh();
     my $pop_id = $self->get_population_id();
    
     my $query = "SELECT DISTINCT (map_version_id) 
                          FROM sgn.map_version 
                          JOIN phenome.genotype_experiment ON (map_version.map_id = genotype_experiment.reference_map_id) 
                          JOIN phenome.genotype USING (genotype_experiment_id) 
                          JOIN phenome.individual USING (individual_id) 
                          WHERE map_version.current_version = 't' AND individual.population_id = ?"; 

      my $sth = $dbh->prepare($query);
      $sth->execute($pop_id);

     my ( $map_version_id ) = $sth->fetchrow_array();
    
     
     return $map_version_id;

}


=head2 linkage_groups

  Usage: my @lg = $pop->linkage_groups();
 Desc: useful for getting a list of the linkage groups in a genetic map 
 Ret: an array of the linkage group names
 Args: none
 Side Effects:
 Example:

=cut

sub linkage_groups {
    my $self=shift;
    my $mapversion_id = $self->mapversion_id();
    
    my $query = "SELECT lg_name FROM sgn.linkage_group WHERE map_version_id = ?";
    
    my $sth = $self->get_dbh()->prepare($query);
    $sth->execute($mapversion_id);
    
    my @lg_names;
   
    while (my $lg_name = $sth->fetchrow_array()) {
	push @lg_names, $lg_name;
    }

    return @lg_names;
	
}




############# 
return  1;
############
