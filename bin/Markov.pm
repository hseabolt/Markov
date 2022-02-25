#!/usr/bin/perl

# package Markov.pm
# A Perl package for Markov modeling (Markov chain processes)

# Author: MH Seabolt
# Last Updated: 11-14-2019	

# Current implementation is primarily for discrete-time Markov chain processes (represented as directed graph adjacency matrices)
# Future plans to include HMM algorithms

package Markov;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(max min sum lg log10 pow round);			 #Import from other packages

use strict;
use warnings;
use experimental 'smartmatch';
use List::Util qw( sum first shuffle );
use Scalar::Util;
use Storable qw(dclone);
use Data::Dumper;
use Carp;
our $AUTOLOAD;

use lib '/media/hunter/Data/scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';
use MarkovNode;


# ESSENTIAL CPAN module for Sets
use Set::Scalar;

################################################################################## 

# Class data and methods --> this section doesnt really differ much from a basic Graph object implentation, but with different terminology tailored to Markov modellings and processes.
# Attributes 
{
	my %_attribute_properties = (
		_states		=> [ ],								# An anonymous array containing the list of states in the model, where each state will be a MarkovNode object
		_adj		=> [ ],								# An anonymous 2D AoA structure that represents the transition probability matrix
	);
	
	# Global variable counter
	my $_count = 0;
	
	# Return a list of all attributes
	sub _all_attributes	{
		keys %_attribute_properties;
	}
	
	# Return the default value for a given attribute
    	sub _attribute_default 	{
     	my( $self, $attribute ) = @_;
        	$_attribute_properties{$attribute};
    	}
    
	# Manage the count of existing objects
	sub get_count	{
		$_count;
	}
	sub _incr_count	{
		++$_count;
	}
	sub _decr_count	{
		--$_count;
	}	
}


############################################################
#                       CONSTRUCTORS                       #
############################################################

# The contructor method
# Construct a new graph (my $node = Markov->new() );
# Returns a scalar reference to a
sub new				{
	my ( $class, %arg ) = @_;
	
	# Create the new object
	my $self = bless {}, $class;

	foreach my $attribute ( $self->_all_attributes() ) {
        	# E.g. attribute = "_name",  argument = "name"
        	my ($argument) = ( $attribute =~ /^_(.*)/ );
        	# If explicitly given
        	if (exists $arg{$argument}) 	{
            	$self->{$attribute} = $arg{$argument};
        	}
        	else	{
            	$self->{$attribute} = $self->_attribute_default($attribute);
        	}
   	}
    	$class->_incr_count();
	return $self;
}

# The clone method
# All attributes are copied from the calling object, unless specifically overriden
# Called from an existing object ( Syntax: $cloned_obj = $obj->clone(); )
sub clone	{
	my ( $caller, %arg ) = @_;
	# Extract the class name from the calling object
	my $class =ref($caller);
		
	# Create a new object
	my $self = bless {}, $class;
		
	foreach my $attribute ( $self->_all_attributes() )	{
		my ($argument) = ( $attribute =~ /^_(.*)/ );
			
		# If explicitly given
		if ( exists $arg{$argument} )	{
			$self->{$attribute} = $arg{$argument};
		}
			
		# Otherwise, copy attribute of new object from the calling object
		else	{
			$self->{$attribute} = $caller->{$attribute};
		}
	}
	$self->_incr_count();
	return $self;
}

#######################################
# Autoload getters and setters
sub AUTOLOAD {
    	my ( $self, $newvalue ) = @_;
    	my ( $operation, $attribute ) = ( $AUTOLOAD =~ /(get|set)(_\w+)$/ );
    
    	# Is this a legal method name?
    	unless( $operation && $attribute ) {
        	croak "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
    	}
    	unless( exists $self->{$attribute} ) {
        	croak "No such attribute $attribute exists in the class ", ref($self);
    	}

    	# Turn off strict references to enable "magic" AUTOLOAD speedup
    	no strict 'refs';

    	# AUTOLOAD accessors
    	if( $operation eq 'get' ) 	{
        	# Install this accessor definition in the symbol table
        	*{$AUTOLOAD} = sub {
            	my ($self) = @_;
          	$self->{$attribute};
     	};
    }
    # AUTOLOAD mutators
    elsif( $operation eq 'set' ) 	{
		# Set the attribute value
        	$self->{$attribute} = $newvalue;

        	# Install this mutator definition in the symbol table
        	*{$AUTOLOAD} = sub {
        		my ($self, $newvalue) = @_;
            	$self->{$attribute} = $newvalue;
        	};
    }

    # Turn strict references back on
    use strict 'refs';

    # Return the attribute value
    return $self->{$attribute};
}


############################################################
#                BASIC UTILITY SUBROUTINES                 #
############################################################


# Returns an empty MarkovNode object
sub _getMarkovNode	{
	return MarkovNode->new();
}

# This hash is associated with the _stateToIndex() subroutine, but it must be outside the subroutine scope to be maintainable and accessible.
my $p = 0;			# The initial index in the prebuilt hash of states below, we will increment it if needed
my %Index = ();
my %ReverseIndex = ();		# Associated with the _indexToState() subroutine and automatically updated by _stateToIndex()
my %Names = ();			
my %NamesIndex = ();
my %NamesState = ();

# Converts key current character into index
sub _stateToIndex	{
	my ( $self, $state )	= @_;
	return if ( not $state );		# Sanity check
	my $index;
	
	if ( exists $Index{$state}  ) 	{
		$index = $Index{$state};
	}
	# If the index doesnt exist in the Index hash, then add it at the end and increment the value
	# This should be fine for multiple letter k-mer style indices and some symbols
	else		{
		$Index{$state} = $p;				# Be careful here, as odd symbols may cause errors
		$ReverseIndex{$p} = $state;
		my $name = $state->{_name};
		$NamesIndex{$p} = $name;
		$Names{$name} = $p;
		$NamesState{$name} = $state;
		$p++;
		$index = $Index{$state};
	}

	return $index;
}

# Converts key current index into character
sub _indexToState	{
	my ( $self, $index )	= @_;
	my $state;
	if ( exists( $ReverseIndex{$index} ) ) 	{
		$state = $ReverseIndex{$index};
	}
	# If the index doesnt exist in the ReverseIndex hash, then there is nothing we can do...
	
	return $state;
}

# Converts from an index to a state NAME
sub _indexToName	{
	my ( $self, $index ) = @_;
	my $name;
	if ( exists $NamesIndex{$index} )	{
		$name = $NamesIndex{$index};
	}
	return $name;
}

# Converts from the name of a state to the index
sub _nameToIndex	{
	my ( $self, $name ) = @_;
	my $index;
	if ( exists $Names{$name} )	{
		$index = $Names{$name};
	}
	return $index;
}

sub _nameToState	{
	my ( $self, $name ) = @_;
	my $state;
	if ( exists $Names{$name} )	{
		$state = $ReverseIndex{ $Names{$name} };
	}
	return $state;
}


# Chooses a random state from the list of ALL possible states in the graph, without regard to any weighting
# Works by getting all the shuffling the array to get a random one.
sub random_state		{
	my ( $self ) = @_;
	my @all_states = shuffle @{ $self->get_states() };	
	return $all_states[0]->get_name;		# Return the first element of the randomly shuffled array
}

# Advances the model to the next state from a given state by randomly sampling from a weighted distribution
sub next_state		{
	my ( $self, $name ) = @_;	
	my $current_state = $self->_nameToState( $name );	
	my $index = $self->_stateToIndex( $current_state );
	my @trPr = @{$self->{_adj}->[$index]};		# The ROW in _adj 
	my $next_state = $current_state->weighted_random_sample( \@trPr );
	return $next_state->get_name;
} 

# Generates the next k states of the system
# Returns an array REFERENCE containing the k future states of the model
sub generate_states		{
	my ( $self, $current_state, $k ) = @_;
	my @future_states;
	
	my $i = 0;
	while ( $i < $k )		{
		my $next_state = $self->next_state( $current_state );
		push @future_states, $next_state;
		$current_state = $next_state;
		$i++;
	}
	return \@future_states;	
}

# Checks if a node is a leaf (terminal) node or not
sub isLeaf	{
	my ( $self ) = @_;
	( $self->{_isFinalState} >= 1 )? return 1 : return 0;
}

# Checks if the sum of transition probabilities exiting a state/node sum to 1.0 (not strictly required, but in theory, should be enforced).
sub _validTransitions		{
	my ( $self ) = @_;
	( sum($self->get_children) == 1.0 )? return 1 : return 0;
}

# Calculates the degree of a specified node
# Returns a scalar integer value of the degree of a node.
sub degree		{
	my ( $self, $u ) = @_;
	my $degree = 0;	
	my $r = $self->_nameToIndex( $u );	
	my @children = @{ $self->{_adj}->[$r] };
	
	for ( my $v = 0; $v < scalar @children; $v++ )	{
		next if ( $children[$v] == 0 || $children[$v] eq "inf" || $children[$v] eq "nan" );	# Increment degree if the edge exists
		$degree++;
		$degree++ if ( $self->_stateToName($self->{_states}->[$v]) ~~ $u );			# Increment degree AGAIN if we have a self loop.
	}
	return $degree;
}

# Calculates the outdegree of a specified node (directed edges originating from this node, to another node) --> the row in the Adj
# Returns a scalar integer value of the degree of a node.
# Only use for a directed graph structure
sub outdegree		{
	my ( $self, $u ) = @_;
	my $degree = 0;	
	my $r = $self->_nameToIndex( $u );	
	my @children = @{ $self->{_adj}->[$r] };
	
	for ( my $v = 0; $v < scalar @children; $v++ )	{
		$degree++ if ( $children[$v] != 0 || $children[$v] ne "inf" || $children[$v] ne "nan" );	# Increment degree if the edge exists

	}
	return $degree;
}

# Calculates the outdegree of a specified node (directed edges originating from this node, to another node) --> the column in the Adj
# Returns a scalar integer value of the degree of a node.
# Only use for a directed graph structure
sub indegree		{
	my ( $self, $u ) = @_;
	my $degree = 0;
	my @children;
	my $c = $self->_nameToIndex( $u );	
	foreach my $row ( @{ $self->{_adj}	} )	{
		push @children, $self->{_adj}->[$row]->[$c];
	}
	
	for ( my $v = 0; $v < scalar @children; $v++ )	{
		$degree++ if ( $children[$v] != 0 || $children[$v] ne "inf" || $children[$v] ne "nan" );	# Increment degree if the edge exists

	}
	return $degree;
}

# Brandes algorithm for betweeness-centrality
# Computes the shortest-path betweenness centrality for all nodes
# Betweenness centrality of a node V is the sum of the fraction of all-pairs shortest paths that pass through V
# NOTE: For weighted graphs the edge weights must be greater than zero.   Zero edge weights can produce an infinite number of equal length paths between pairs of nodes. 
# Returns a hash with KEYS = nodes, VALUES= betweenness centrality as the value
sub betweenness_centrality		{
	my ( $self ) = @_;
	my @names;
	
	# Initalize the hash of values with KEYS= vertices and VALUES= betweeness centrality value
	my %C = ();
	my %P = ();
	my %G = ();
	my %D = ();
	my %E = ();
	foreach my $state ( @{$self->{_states}} )	{
		my $name = $self->_stateToName($state);
		push @names, $name;
		$C{$name} = 0;
	}
	
	foreach my $s ( @names )	{
		# Initialize various structures for each node in the graph
		my @S;
		my @q;
		foreach ( @names )	{
			$P{$_} = [];
			$G{$_} = 0;
			$D{$_} = -1;
			$E{$_} = 0;
		}
		$G{$s} = 1;
		$D{$s} = 0;
		push @q, $s;
		
		# Loop
		while ( scalar @q > 0 )	{
			my $v = shift(@q);
			my $vname = $self->_nameToIndex($v);
			push @S, $v;
			foreach my $w ( @names )		{
				my $wname = $self->_nameToIndex($w);
				next if ( $self->{_adj}->[$v]->[$w] == 0 || $self->{_adj}->[$v]->[$w] eq "inf" || $self->{_adj}->[$v]->[$w] eq "nan" );
				if ( $D{$w} < 0 )	{
					push @q, $w;
					$D{$w} = $D{$v} + 1;
				}
				if ( $D{$w} == $D{$v} + 1 )	{
					$G{$w} = $G{$w} + $G{$v};
					push @{$P{$w}}, $v;
				}
			}
		}
		
		while ( scalar @S > 0 )	{
			my $w = pop @S;
			foreach my $v ( @{ $P{$w} } )	{
				$E{$v} = $E{$v} + ($G{$v}/$G{$w}) * (1+$E{$w});
			}
			if ( $w !~ $s )	{
				$C{$w} = $C{$w} + $E{$w};
			}
		}
	}
	return \%C;
}

# Utility subroutine to make a deep copy of a Markov object (similar to the clone() subroutine)
sub copy_model	{
	my ( $self ) = @_;
	my $clone = dclone $self;
	return $clone;
}

# Utility subroutine to return an EMPTY adjacency matrix with the same structure (states) as the input.
sub empty	{
	my ( $self ) = @_;
	my @Adj = @{$self->get_adj};
	for ( my $u=0; $u < scalar @{ $self->{_states} }; $u++ )	{
		my $state = $self->_indexToState( $u );	 
		
		# This loop just resizes the array and initializes zeros where needed
		for ( my $v=0; $v < scalar @{ $self->{_states} }; $v++ )	{
			my $child = $self->_indexToState( $v );	
			$child->{_children}->[$v] = $child if ( $Adj[$u]->[$v] && $Adj[$u]->[$v] > 0 );
			next if ( $Adj[$u]->[$v] );
			$Adj[$u]->[$v] = 0;			
		}	
	}
	$self->set_adj( \@Adj );
}



# When an object is no longer being used, garbage collect it and adjust count of existing objects
sub DESTROY	{
	my ( $self ) = @_;
	$self->_decr_count();
}

############################################################
#          INSERTION AND DELETION SUBROUTINES              #
############################################################

# Adds a vertex (a "state") to the Markov chain graph
# Updates the _states array and the %Index /%ReverseIndex hashes
# Transition probabilties should be a hash ref of probabilities with keys as the state objects and the values as the probability weights
# Operates DIRECTLY on the Markov object.
sub add_state		{
	my ( $self, $new_state, $transitions ) = @_;

	# Initialize a new node object
	$new_state = ( $new_state )? $new_state : $p;		# ALL states must have a name! Use $p as a placeholder if no name is supplied.
	my $node = MarkovNode->new( "name" => "$new_state" );
	
	# Update the requisite hashes and lists with the new node
	push @{ $self->{_states} }, $node;
	$self->_stateToIndex( $node );
		
	# Update _adj and set the children of the new state
	my @Adj = @{$self->get_adj};
	for ( my $u=0; $u < scalar @{ $self->{_states} }; $u++ )	{
		my $state = $self->_indexToState( $u );	 
		
		# This loop just resizes the array and initializes zeros where needed
		for ( my $v=0; $v < scalar @{ $self->{_states} }; $v++ )	{
			my $child = $self->_indexToState( $v );	
			$child->{_children}->[$v] = $child if ( $Adj[$u]->[$v] && $Adj[$u]->[$v] > 0 );
			next if ( $Adj[$u]->[$v] );
			$Adj[$u]->[$v] = 0;			
		}	
	}
	
	if ( $transitions )	{
		for ( my $w=0; $w < scalar @{$transitions}; $w++ )	{
			my $state = $self->_indexToState( $w );
			$node->{_children}->[$w] = $state if ( $transitions->[$w] > 0 );
			$Adj[-1]->[$w] = $transitions->[$w];
		}
	}
	$self->set_adj( \@Adj );
}

# Deletes a state from both the list of states and adj
sub delete_state	{
	my ( $self, $name ) = @_;
	my $kill_state = $self->_nameToState( $name );
	
	# Update _states
	my $index = $self->_stateToIndex( $kill_state );
	
	# Update _adj and the children of the affected nodes
	for ( my $u=0; $u < scalar @{ $self->{_adj} }; $u++ )	{
		# Delete the column of the dead node
		splice( @{ $self->{_adj}->[$u] }, $index, 1 );
		my $state = $self->{_states}->[$u];
		splice( @{ $state->{_children} }, $index, 1 );
	}
	# Finally, delete the entire row for the node
	splice( @{ $self->{_adj} }, $index, 1 );
	splice( @{ $self->{_states} }, $index, 1 );
	$kill_state->DESTROY;
	
	# Update the organizational hashes
	# Under construction
}

# Adds an edge with weight $weight (required) between source state $u and target state $v
# If $weight is not given as an arg, then nothing will happen
# If an edge already exists between $u --> $v, then this will overwrite it
sub add_edge		{
	my ( $self, $uname, $vname, $weight ) = @_;
	my $u = $self->_nameToState( $uname );
	my $v = $self->_nameToState( $vname );
	
	# Sanity checks
	if ( not $u || not $v || not $weight )	{
		warn " --- Markov::add_edge() ERROR:  missing subroutine argument!\n";
		return;
	}
	elsif ( not exists $Names{$u} || not exists $Names{$v} )	{
		warn " --- Markov::add_edge() ERROR:  cannot add edges between nodes that don't exist!\n";
		return;
	}
	
	# Update the edge in _adj
	$self->{_adj}->[ $self->_stateToIndex($u) ]->[ $self->_stateToIndex($v) ] = $weight;
	
	# Update the children of the $u (the only node with a new child
	my @children = @{ $self->{_adj}->[$u] };
	$u->set_children( \@children );
	
}

# Deletes an edge (here, we are actually just resetting the Pr to 0) between source $u and target $v.
sub delete_edge	{
	my ( $self, $uname, $vname ) = @_;
	my $u = $self->_nameToState( $uname );
	my $v = $self->_nameToState( $vname );
	
	# Sanity checks
	if ( not exists $Index{$u} || not exists $Index{$v} )	{
		warn " --- Markov::delete_edge() ERROR:  cannot delete edges between nodes that don't exist!\n";
		return;
	}
	
	# Update the edge Pr to 0
	$self->{_adj}->[ _stateToIndex($u) ]->[ _stateToIndex($v) ] = 0;
	
	# Update the children of the $u (the only node with a new child
	my @children = @{ $self->{_adj}->[$u] };
	$u->set_children( \@children );
}

# Adds an edge with weight $weight (required) between source state $u and target state $v
# If $weight is not given as an arg, then nothing will happen
# If an edge already exists between $u --> $v, then this will overwrite it
# This code is the same as add_edge(), but functions as an alias.
sub update_edge		{
	my ( $self, $uname, $vname, $weight ) = @_;
	my $u = $self->_nameToState( $uname );
	my $v = $self->_nameToState( $vname );
	
	# Sanity checks
	if ( not $u || not $v || not $weight )	{
		warn " --- Markov::update_edge() ERROR:  missing subroutine argument!\n";
		return;
	}
	elsif ( not exists $Index{$u} || not exists $Index{$v} )	{
		warn " --- Markov::update_edge() ERROR:  cannot alter edges between nodes that don't exist!\n";
		return;
	}
	
	# Update the edge in _adj
	$self->{_adj}->[ $self->_stateToIndex($u) ]->[ $self->_stateToIndex($v) ] = $weight;
	
	# Update the children of the $u (the only node with a new child
	my @children = @{ $self->{_adj}->[$u] };
	$u->set_children( \@children );
}







##########################################################################
# DISCRETE-TIME MARKOV CHAIN PROCESSES
##########################################################################

# helper function that executes a DFS search and returns a boolean true or false if node $v is reachable from node $u
sub _isAccessible	{
	my ( $self, $u, $v ) = @_;
	my $preorder = $self->DFS( $u );
	my %Elements = map { $_ => 1 } @{$preorder};
	if ( exists $Elements{$v} )	{
		return 1;
	}
	else		{
		return 0;
	}
}

# Helper function that returns 1 if a Markov model cannot be reduced (if all states are reachable from any other state)
# Returns 0 if not.
sub _isIrreducible	{
	my ( $self ) = @_;
	my %AllStates = ();
	foreach my $state ( @{$self->{_states}} )	{
		my $name = $self->_stateToName($state);
		$AllStates{$name} = 1;
	}
	
	for ( my $u=0; $u < scalar @{ $self->{_states} }; $u++ )	{
		my $uname = $self->_stateToName( $self->{_states}->[$u] );
		my @reachable_nodes = @{ $self->DFS($uname) };
		
		# Check if @reachable nodes contains the right number of nodes, 
		# Then confirm that each node is in the hash
		if ( scalar @reachable_nodes == scalar @{ $self->{_states} } )	{
			for ( my $v=0; $v < scalar @reachable_nodes; $v++ )	{
				return 0 if ( not exists($AllStates{$reachable_nodes[$v]}) );
			}
		}
		else		{
			return 0;
		}
	}
	return 1;
}


# Depth-first search of the markov chain
# Returns an array reference containing either the pre- or post-ordered names of the nodes (tree-edges only, not cross/forward/back edges)
sub DFS		{
	my ( $self, $src, $type ) = @_;
	my @queue = ($src);
	my %Seen = ();
	my @preorder;
	my @postorder;
	$type = ( $type && $type =~ /post/ )? "post" : "pre";
	
	while( scalar( @queue ) > 0 )	{
		my $vertex = pop( @queue ); 
		push @preorder, $vertex if ( not exists($Seen{$vertex}) );			

		# Get the neighbors of $vertex --> remember: neighbors is a Set!
		my $vstate = $self->_nameToState( $vertex );
		my $neighbors = $vstate->get_neighbors;
			 
		# Add the neighbors to the queue and update some information about the current node
		foreach ( $neighbors->members() )	{
			if ( exists $Seen{$_} )	{
				# Can identify cycles, cross/back/forward edges here!
				next;
			}
			push @queue, $_;
		}
		
		# Add the current vertex to %Seen
		$Seen{$vertex} = 1;
	}
	
	# Return the ordering requested
	( $type eq "post" )? return \@postorder : return \@preorder; 			 # NOTE: 10-10-2019 MHS: post-order return currently not implemented, so only use pre-order
}

# Shamelessly jacked and refactored from somewhere on Google
# Subroutine for a BREADTH first search of a graph structure
sub BFS		{
	my ( $self, $src ) = @_;
	my @queue = ($src);
	my @bfs;
	my %Seen = ();
	
	while( scalar( @queue ) > 0 )	{
		my $vertex = shift( @queue ); 
		push @bfs, $vertex if ( not exists($Seen{$vertex}) );
			
		# Get the neighbors of $vertex --> remember: neighbors is a Set!
		my $vstate = $self->_nameToState( $vertex );
		my $neighbors = $vstate->get_neighbors;
		foreach ( $neighbors->members() )	{
			next if ( exists $Seen{$_} );
			push @queue, $_;
		}
			
		# Update the current vertex information
		$Seen{$vertex} = 1;
	}
	return \@bfs;
}

# Floyd-Warshall algorithm for computing all pairs shortest paths
sub APSP_floyd_warshall		{
	my ( $self ) = @_;
	
	# Get the corrected diagonals _adj for this algorithm
	my @M = @{$self->get_adj};
	my $n = scalar @M;
	
	# Initialize an empty HoH to store the length of the shortest paths
	# AND an empty matrix to hold traceback data
	my %FloydShortestPaths = ();
	my %Traceback = ();
	
	# Save off an enumeration of the nodes so that we can map them back later
	my @order = @{ $self->get_states };
	my %Indices = ();
	for ( my $i = 0; $i < scalar @order; $i++ )	{
		my $node = $self->_stateToName($order[$i]);
		$Indices{$i} = $node;
	}
	
	# Set all cells in the traceback matrix to 0
	for ( my $i=0; $i < $n; $i++ )	{
		for ( my $j=0; $j < $n; $j++ )	{
			$Traceback{$Indices{$i}}{$Indices{$j}} = 0;
		}
	}
	
	foreach my $i ( 0..$n-1 )	{
		foreach my $j ( 0..$n-1 )	{
			$Traceback{$Indices{$i}}{$Indices{$j}} = $i;
			
			if ( $i != $j && $M[$i][$j] == 0 )	{
				$Traceback{$Indices{$i}}{$Indices{$j}} = "inf";
				$M[$i][$j] = "inf";
			}
		}
	}
	
	# Set up the Floyd-Warshall nested loops:
	# Use the indices we set up in the enumeration to act as the keys for the ShortestPaths HoH
	foreach my $k ( 0..$n-1 )	{
		my $kX = $Indices{$k};
		foreach my $i ( 0..$n-1 )	{
			my $iX = $Indices{$i};
			foreach my $j ( 0..$n-1 )	{
				my $jX = $Indices{$j};
				if ( $M[$i][$j] > $M[$i][$k] + $M[$k][$j] )	{
					$FloydShortestPaths{$iX}{$jX} = $M[$i][$k] + $M[$k][$j];
					$Traceback{$iX}{$jX} = $Traceback{$kX}{$jX};
					$M[$i][$j] = $M[$i][$k] + $M[$k][$j];
				}
				else		{
					$FloydShortestPaths{$iX}{$jX} = $M[$i][$j];
					$Traceback{$iX}{$jX} = $Traceback{$iX}{$jX};
				}				
			}
		}
	}
	return ( \%FloydShortestPaths, \%Traceback );
}

############################## FLOYD WARSHALL UTILITY SUBROUTINES #####################################################
# Traces back the sequence of nodes from source U --> target V using the Traceback matrix built in the Floyd-Warshall subroutine
# Input args are ( @Traceback AoAref rom FW algorithm, vertex U, vertex V, and a path arrayref, which wil typically be empty when originally calling this routine )
sub get_path	{
	my ( $traceback, $i, $j, $path ) = @_;	
	my @order = sort keys %{$traceback};

	if ( $i eq $j )	{
		push @{$path}, $i;
	}
	elsif ( $traceback->{$i}->{$j} eq "inf" ) 	{
		push @{$path}, "-";
	}
	else		{
		my $k = $order[ $traceback->{$i}->{$j} ];
		get_path( $traceback, $i, $k, $path );
		push @{$path}, $j;
	}
	return $path;
}

# Returns the period of a given state in the Markov chain
sub get_period	{
	my ( $self, $state ) = @_;
	my $bellmanford = $self->bellman_ford($state);
	my @paths = values %{$bellmanford};
	return find_gcd( \@paths );
}

# Find the GCD of a list of numbers
sub find_gcd 	{
	my $list = @_;
	my $gcd = gcd_euclid( $list->[0], $list->[1] );		# The first two elements
	if ( scalar @{$list} > 2 )	{
		for ( my $i=2; $i < scalar @{$list}; $i++ )	{
			next if ( $list->[$i] eq "inf" || $list->[$i] eq "nan" );
			$gcd = gcd_euclid( $gcd, $list->[$i] );
		}
	}
	return $gcd;
}


# Helper function for computing periodicity
sub gcd_euclid 	{
  	my ($a, $b) = @_;
  	($a,$b) = ($b,$a) if $a > $b;
  	while ($a) {
    		($a, $b) = ($b % $a, $a);
  	}
  	return $b;
}

# Returns boolean true if the Markov chain is aperiodic, and false if it is periodic
sub _isAperiodic	{
	my ( $self ) = @_;
	
	my @periods;
	for ( my $u=0; $u < scalar @{ $self->{_states} }; $u++ ) 	{
		my $name = $self->_stateToName( $self->{_states}->[$u] );
		my $period = $self->get_period($name);
		return 0 if ( $period != 1 );
	}
	return 1;
}

# Checks if a given state is transient or not, returns boolean true if yes, false if no.
# If a state is transient, then by definition, once the system leaves this state, it will not be able to come back.
sub _isTransient	{
	my ( $self, $state ) = @_;
	my $index = $self->_nameToIndex($state);
	my $sum = 0;
	for ( my $c=0; $c < scalar @{$self->{_adj}}; $c++ )	{
		$sum += $self->{_adj}->[$c]->[$index];
	}
	( $sum < 1 )? return 1 : return 0;
	
}

# Checks if a given state in the Markov chain is absorbing, meaning 
# that once we arrive at this state, it is not possible to transition out of it
sub _isAbsorbing	{
	my ( $self, $state ) = @_;
	my $index = $self->_nameToIndex($state);
	my @row = @{ $self->{_adj}->[$index] };
	
	for ( my $r=0; $r < scalar @row; $r++ )	{
		if ( $row[$r] > 0  )	{
			if ( $r == $index )		{
				next;
			}
			else		{
				return 0; 
			}
		}
		else		{
			return 1;
		}
	}
}

# An implementation of the Bellman-Ford algorithm for finding single source shortest paths
sub bellman_ford		{
	my ( $self, $src ) = @_;	
	my $graph = $self->AoA2HoH;
	my @vertices = keys %{$graph};
	my @path; 
	
	# Initialize the hash to store the distances (init as "inf") and a Traceback matrix with null placeholders
	my %Distance = ();
	my %Traceback = ();
	foreach my $vertex ( @vertices )	{
		$Distance{$vertex} = "inf";
		foreach my $v ( @vertices ) 	{
			$Traceback{$vertex}{$v} = 0;
		}
	}
	$Distance{$src} = 0;		# Distance to the original source is zero
	
	# Relax all the edges repeatedly
	foreach my $i ( scalar @vertices - 1 )	{
		foreach my $u ( keys @vertices )	{
			foreach my $v ( @vertices )	{
				if ( $Distance{$u} + $graph->{$u}->{$v} < $Distance{$v} )	{
					$Distance{$v} = $Distance{$u} + $graph->{$u}->{$v};
					$Traceback{$u}{$v} = first { $vertices[$_] ~~ $u } 0..$#vertices;
				}
			}
		}
	}
	
	# Check for negative weight cycles
	foreach my $u ( keys @vertices )	{
		foreach my $v ( @vertices )	{
			if ( $Distance{$u} + $graph->{$u}->{$v} < $Distance{$v} )	{
				print STDERR "ERROR:  Bellman-Ford subroutine has identified a negative-weight cycle in the graph!\n";
				return;
			}
		}
	}
	
	return ( \%Distance, \%Traceback );		 # 10-11-2019 MHS: Traceback currently not working well, but it will run without error.  Final output is incorrect.
}

# Utility subroutine
# Convert the _adj matrix built into our model into an HoH, since some subroutines here operate more easily on HoH (and were originally written as such)
# Returns a reference to a HoH
sub AoA2HoH	{
	my ( $self ) = @_;
	my %HoH = ();
	for ( my $i=0; $i < scalar @{$self->{_states}}; $i++ )	{
		my $name1 = $self->_stateToName( $self->{_states}->[$i] );
		for ( my $j=0; $j < scalar @{$self->{_states}}; $j++ )	{
			my $name2 = $self->_stateToName( $self->{_states}->[$j] );
			$HoH{$name1}{$name2} = $self->{_adj}->[$i]->[$j];
		}
	}
	return \%HoH;
}

# Print a matrix in a nice format --> does not assume a square matrix
sub print_matrix	{  
	my ( $matrix, $sep ) = @_;
	$sep = ( $sep )? $sep : "\t";
	for ( my $u=0; $u < scalar @{$matrix}; $u++ )	{
		print join("$sep", @{$matrix->[$u]}), "\n";
	}
}




























# Last line in a package must always be 1 !
1;
