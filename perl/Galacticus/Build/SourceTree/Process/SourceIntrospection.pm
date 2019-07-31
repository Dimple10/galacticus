# Contains a Perl module which implements source code introspection.

package Galacticus::Build::SourceTree::Process::SourceIntrospection;
use strict;
use warnings;
use File::Slurp qw(slurp);
use Data::Dumper;
use XML::Simple;
use Encode;
    
# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'sourceIntrospection'} = \&Process_Source_Introspection;

sub ReadFile {
    # Read a file and add instrumentation to assist in source introspection.
    my $fileName = shift();
    my $code     = slurp($fileName);
    $code        = &Instrument($code);
    return $code;
}

sub Instrument {
    # Add instrumentation to assist in source introspection.
    my $codeIn     = shift();
    my $codeOut             ;
    my $lineNumber = 0      ;
    open(my $code,"<",\$codeIn);
    while ( my $line = <$code> ) {
	++$lineNumber;
	if ( $line =~ m/\{introspection:location(:compact)??\}/ ) {
	    $line =~ s/\{(introspection:location(:compact)??)\}/{$1:$lineNumber}/g;
	}
	$codeOut .= $line;
    }
    return $codeOut;

}

sub Process_Source_Introspection {
    # Get the tree.
    my $tree = shift();
    # Walk the tree.
    my $node  = $tree;
    my $depth = 0;    
    while ( $node ) {
	if ( $node->{'type'} eq "code" ) {
	    die("Galacticus::Build::SourceTree::Process::SourceIntrospection::Process_Source_Introspection(): content does not exist")
		unless ( exists($node->{'content'}) && defined($node->{'content'}) );
	    my $newCode;
	    my $content = encode(q{utf8},$node->{'content'});
	    open(my $code,q{<:utf8},\$content);
	    while ( my $line = <$code> ) {
		if ( $line =~ m/\{introspection:location(:compact)??:(\d+)\}/ ) {
		    my $compact    = defined($1);
		    my $lineNumber = $2;
		    my $location   = &Location($node,$lineNumber,compact => $compact);
		    $line =~ s/\{introspection:location(:compact)??:\d+\}/$location/g;
		}
		$newCode .= $line;
	    }
	    close($code);				    
	    $node->{'content'} = $newCode;	    
	}
	# Move on to the next node.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

sub Location {
    # Get the node and line number.
    my $node       = shift();
    my $lineNumber = shift();
    my %options    = @_;
    $options{'compact'} = 0
	unless ( exists($options{'compact'}) );
    # Construct the location.
    my $location   = $options{'compact'} ? "''" : "char(10)//' Occurred at:'";
    my $branch     = $node;
    while ( $branch ) {
	if (
	    $branch->{'type'} eq "file"
	    ||
	    $branch->{'type'} eq "module"
	    ||
	    $branch->{'type'} eq "function"
	    ||
	    $branch->{'type'} eq "subroutine"
	    ) {
	    if ( $options{'compact'} ) {
		$location .= "//'".$branch->{'type'}."(".$branch->{'name'}.")'";
	    } else {
		$location .= "//char(10)//'   ".(" " x (10-length($branch->{'type'}))).$branch->{'type'}.":".$branch->{'name'}."'";
	    }
	}
	$location .= "//char(10)//'    directive:".$branch->{'type'}."'"
	    if ( exists($branch->{'directive'}) && ! $options{'compact'} );
	$branch = $branch->{'parent'};
    }
    if ( $options{'compact'} ) {
	$location .= "//':".$lineNumber."'";
    } else {
	$location .= "//'   [line ".$lineNumber."]'";
    }
    return $location;
}

1;
