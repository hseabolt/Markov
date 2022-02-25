#!/usr/bin/perl

# markov_driver.pl
# Driver and testing code for Markov.pm and MarkovNode.pm modules

use strict;
use warnings;

#use lib '/media/hunter/Data/scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';
use Markov;

# Begin driver code
print "\n\n************************\n\n";

# Generate a set of states for the model
my @weather_states = ("sunny", "snowy", "rainy");
my $probs = [
	[ 0.7, 0.1, 0.2 ],
	[ 0.05, 0.7, 0.25 ],
	[ 0.25, 0.2, 0.55 ],
];

# Instantiate a new Markov object
my $model = Markov->new();

# Add the states to the model and set their names
for ( my $u=0; $u < scalar @weather_states; $u++ )	{
	print "$weather_states[$u] -- Probs[u]: ", join(", ", @{$probs->[$u]}), "\n";
	$model->add_state( $weather_states[$u], $probs->[$u]  );
}

# Print the adjacency matrix
my @Adj = $model->get_adj;
print "\n\nAdj: \n";
print join(", ", @{$model->get_states}), "\n";
foreach my $u ( @Adj )	{
	foreach my $v ( @{$u} )	{
		print join(", ", @{$v}), "\n";
	}
}

print "\n\n";

# Add a 4th state
$model->add_state( "squanch", [ 0.2, 0.5, 0.1, 0.2 ]  );

# Print the adjacency matrix
@Adj = $model->get_adj;
print "\n\nUpdated Adj: \n";
print join(", ", @{$model->get_states}), "\n";
foreach my $u ( @Adj )	{
	foreach my $v ( @{$u} )	{
		print join(", ", @{$v}), "\n";
	}
}

print "\n\n";


# Get a random starting state
my $first_state = $model->random_state();
print "Randomly chosen starting first state: ", $first_state->get_name(), "\n";

# Generate the next state of the model
my $second_state = $model->next_state( $first_state );
print "Second state: ", $first_state->get_name(), "\n";

# Generate the next 10 states of the model
my @future_states = @{ $model->generate_states( $second_state, 10 ) };
print "Future state: ", $_->get_name(), "\n" foreach ( @future_states );

print "\n\n************************\n\n";












exit;
