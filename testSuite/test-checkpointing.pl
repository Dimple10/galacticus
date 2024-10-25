#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;

# Run a set of Galacticus models to test checkpointing functionality.
# Andrew Benson (05-July-2023)

# Ensure output directory exists.
system("mkdir -p outputs");

# # Run model without checkpointing.
system("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/checkpointingNoCheckpoints.xml"  );
die("FAILED: failed to run model without checkpointing")
    unless ( $? == 0 );

# # Run model with checkpointing - interrupted.
system("export OMP_NUM_THREADS=1; cd ..; rm -f checkpoint.chk; ./Galacticus.exe testSuite/parameters/checkpointingCheckpoints.xml");
die("FAILED: failed to produce a checkpoint file")
    unless ( -e "../checkpoint.chk" );

# # Run model with checkpointing - resuming.
system("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/checkpointingResume.xml");
die("FAILED: failed to run model resuming from checkpoint")
    unless ( $? == 0 );

# Read model data.
my $data;
foreach my $modelName ( "checkpointingNoCheckpoints", "checkpointingCheckpoints" ) {
    my $model = new PDL::IO::HDF5("outputs/".$modelName.".hdf5");
    my $nodes = $model->group("Outputs/Output1/nodeData");
    $data->{$modelName}->{$_} = $nodes->dataset($_)->get()
	foreach ( "nodeIndex", "basicMass" );
}

# Compare results.
my $status            = "SUCCESS";
my $toleranceRelative = pdl 0.01;
for(my $pass=0;$pass<2;++$pass) {
    for(my $i=0;$i<nelem($data->{'checkpointingNoCheckpoints'}->{'nodeIndex'});++$i) {
	my $match = which($data->{'checkpointingCheckpoints'}->{'nodeIndex'} == $data->{'checkpointingNoCheckpoints'}->{'nodeIndex'}->(($i)));
	next
	    unless ( nelem($match) == 1 );
	my $massNoCheckpoints = $data->{'checkpointingNoCheckpoints'}->{'basicMass'}          ->(($i));
	my $massCheckpoints   = $data->{'checkpointingCheckpoints'  }->{'basicMass'}->($match)->(( 0));
	my $agrees            = abs($massNoCheckpoints-$massCheckpoints) < $toleranceRelative*$massNoCheckpoints;
	$status = "FAIL"
	    unless ( $agrees );
	print "\t(".$i.") ".$massNoCheckpoints." == ".$massCheckpoints." ? : ".($agrees ? "success" : "FAILURE")."\n"
	    if ( $status eq "FAIL" && $pass == 1 );
    }
}
print $status.": resume from checkpoint file\n";

exit;
