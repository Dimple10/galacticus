#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use File::Copy;
use Data::Dumper;
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Galacticus::Constraints::Parameters;

# Construct halo mass function data from a variety of cosmological N-body simulations.
# Andrew Benson (14-October-2020)

# Get command line options.
my %options =
    (
     submitSleepDuration =>  5,
     waitSleepDuration   => 30,
     pbsJobMaximum       => 64,
     slurmJobMaximum     => 64,
     ompThreads          => "max"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Validate required parameters are present.
die('simulationDataPath is required but is not present')
    unless ( exists($options{'simulationDataPath'}) );

# Define simulations to process.
my @simulations =
    (
     {
         label               => "VSMDPL",
         subpath             => "CosmoSim",
         description         => "Halo mass function for non-backsplash halos from the VSMDPL simulation.",
         simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
         simulationURL       => "https://www.cosmosim.org/cms/simulations/vsmdpl/",
         hubbleConstant      => 0.6777,
         massParticle        => 6.2e6,
         subvolumes          => 10,
         expansionFactors    => [ 1.00000, 0.67110, 0.50250, 0.33000, 0.24750 ]
     },
     {
         label               => "SMDPL",
         subpath             => "CosmoSim",
         description         => "Halo mass function for non-backsplash halos from the SMDPL simulation.",
         simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
         simulationURL       => "https://www.cosmosim.org/cms/simulations/smdpl/",
         hubbleConstant      => 0.6777,
         massParticle        => 9.63e7,
         subvolumes          => 10,
         expansionFactors    => [ 1.00000, 0.66430, 0.50000, 0.33100, 0.24800 ]
     },
     {
         label               => "MDPL2",
         subpath             => "CosmoSim",
         description         => "Halo mass function for non-backsplash halos from the MDPL2 simulation.",
         simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
         simulationURL       => "https://www.cosmosim.org/cms/simulations/mdpl2/",
         hubbleConstant      => 0.6777,
         massParticle        => 1.51e9,
         subvolumes          => 10,
         expansionFactors    => [ 1.00000, 0.67120, 0.50320, 0.33030, 0.24230 ]
     },
     {
         label               => "BigMDPL",
         subpath             => "CosmoSim",
         description         => "Halo mass function for non-backsplash halos from the BigMDPL simulation.",
         simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
         simulationURL       => "https://www.cosmosim.org/cms/simulations/bigmdpl/",
         hubbleConstant      => 0.6777,
         massParticle        => 2.359e10,
         subvolumes          => 10,
         expansionFactors    => [ 1.00000, 0.67040, 0.50000, 0.31800, 0.25700 ]
     },
     {
         label               => "HugeMDPL",
         subpath             => "CosmoSim",
         description         => "Halo mass function for non-backsplash halos from the HugeMDPL simulation.",
         simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
         simulationURL       => "https://www.cosmosim.org/cms/simulations/hugemdpl/",
         hubbleConstant      => 0.6777,
         massParticle        => 7.9e10,
         subvolumes          => 10,
         expansionFactors    => [ 1.00000, 0.67120, 0.50320, 0.33030, 0.24770 ]
     },
     {
     	 # Milky Way zoom-in simulations from the Symphony suite. Host halo IDs were found (by Ethan Nadler) by matching the MAH to
     	 # that of the target halo in the originasl cosmological box simulation.
     	 label                   => "MilkyWay",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo023", "Halo088", "Halo119", "Halo188", "Halo247", "Halo268", "Halo270", "Halo288", "Halo327", "Halo349", "Halo364", "Halo374", "Halo414", "Halo415", "Halo416", "Halo440", "Halo460", "Halo469", "Halo490", "Halo530", "Halo558", "Halo567", "Halo570", "Halo606", "Halo628", "Halo641", "Halo675", "Halo718", "Halo738", "Halo749", "Halo797", "Halo800", "Halo825", "Halo829", "Halo852", "Halo878", "Halo881", "Halo925", "Halo926", "Halo937", "Halo939", "Halo967", "Halo9749", "Halo9829", "Halo990", "Halo004", "Halo113" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM3",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 3keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 2.81981e5,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM4",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 4keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 2.81981e5,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM5",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 5keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 2.81981e5,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM6",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 6keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 2.81981e5,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM6.5",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 6.5keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 2.81981e5,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM10",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 10keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 2.81981e5,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1GeV_envelope",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1GeV (envelope) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1GeV_halfmode",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo023", "Halo113" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1GeV (halfmode) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1e-2GeV_envelope",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1e-2GeV (envelope) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1e-2GeV_halfmode",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1e-2GeV (halfmode) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1e-4GeV_envelope",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1e-4GeV (envelope) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1e-4GeV_halfmode",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1e-4GeV (halfmode) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_25.9e-22eV",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 25.9-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_69.4e-22eV",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 69.4e-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_113e-22eV",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 113e-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_151e-22eV",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 151e-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_185e-22eV",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 185-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_490e-22eV",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004", "Halo113", "Halo023" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 490e-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 2.81981e5,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM3_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 3keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 3.5247625e4,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM4_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 4keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 3.5247625e4,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM5_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 5keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 3.5247625e4,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM6_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 6keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 3.5247625e4,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM6.5_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 6.5keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 3.5247625e4,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_WDM10_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 10keV WDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 massParticle            => 3.5247625e4,
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 processExtract          => \&zoomInsProcessExtract,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1GeV_envelope_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1GeV (envelope) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1GeV_halfmode_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1GeV (halfmode) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1e-2GeV_envelope_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1e-2GeV (envelope) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1e-2GeV_halfmode_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1e-2GeV (halfmode) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1e-4GeV_envelope_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1e-4GeV (envelope) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_IDM1e-4GeV_halfmode_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 1e-4GeV (halfmode) IDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_25.9e-22eV_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 25.9-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_69.4e-22eV_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 69.4e-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_113e-22eV_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 113e-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_151e-22eV_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 151e-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_185e-22eV_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 185-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 label                   => "MilkyWay_fdm_490e-22eV_hires",
     	 subpath                 => "ZoomIns",
     	 realizations            => [ "Halo004" ],
     	 description             => "Halo mass function for non-backsplash halos from Milky Way, 490e-22eV FDM zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.5247625e4,
     	 massHostLogMin          => 12.0, 
     	 massHostLogMax          => 12.3, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       => [
     	     \&zoomInsPreProcessExtractLocate,
     	     \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     },
     {
     	 # LMC zoom-in simulations from the Symphony suite. Host halo IDs were found (by Ethan Nadler) by matching the MAH to that of
     	 # the target halo in the originasl cosmological box simulation.
     	 label                   => "LMC",
     	 subpath                 => "ZoomIns",
     	 # Halo479 and Halo888 have been excluded due to contamination by low-res particles 
     	 realizations            => [ "Halo032", "Halo059", "Halo0662", "Halo083", "Halo088", "Halo097", "Halo104", "Halo110", "Halo202", "Halo208", "Halo218", "Halo296", "Halo301", "Halo303", "Halo340", "Halo374", "Halo380", "Halo391", "Halo405", "Halo440", "Halo463", "Halo4662", "Halo511", "Halo524", "Halo539", "Halo567", "Halo575", "Halo602", "Halo697", "Halo711", "Halo721", "Halo767", "Halo802", "Halo824", "Halo850", "Halo853", "Halo914", "Halo932", "Halo933" ],
     	 description             => "Halo mass function for non-backsplash halos from LMC zoom-in simulations.",
     	 simulationReference     => "Nadler et al.",
     	 simulationURL           => "https://www",
     	 hubbleConstant          => 0.7,
     	 massParticle            => 3.52476e4,
     	 massHostLogMin          => 10.95, 
     	 massHostLogMax          => 11.50, 
     	 subvolumes              => 1,
     	 expansionFactors        => [   1.0000,   0.66503,   0.50239,   0.32987,   0.20064 ],
     	 snapshots               => [ 235     , 203      , 181      , 148      , 109       ],
     	 processIdentify         => \&zoomInsProcessIdentify,
     	 preprocessExtract       =>
     	     [
     	      \&zoomInsPreProcessExtractLocate,
     	      \&zoomInsPreProcessExtractUncontaminated,
     	     ],
     	 processExtract          => \&zoomInsProcessExtract,
     	 postprocessExtract      =>
     	     [
     	      \&zoomInsPostprocessSelectInSphere    ,
     	      \&zoomInsPostprocessExtractSelectedIDs,
     	      \&zoomInsPostprocessSelectInICs       ,
     	      \&zoomInsPostprocessAnalyze           ,
     	      \&zoomInsPostprocessSetVolume
     	     ],
     	 postprocessMassFunction =>
     	     [
     	      \&zoomInsPostProcessMassFunction
     	     ]
     }
    );

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Determine number of OpenMP threads to use.
my $ompThreads = $options{'ompThreads'} eq "max" ? $queueConfig->{'ppn'} : $options{'ompThreads'};

# Get an XML object.
my $xml = new XML::Simple();

# Iterate over simulations to identify always-isolated halos.
{
    my @jobsIdentify;
    foreach my $simulation ( @simulations ) {
	# Convert particle mass to Solar masses.
	$simulation->{'massParticle'} /= $simulation->{'hubbleConstant'};
	# Determine minimum and maximum masses for the mass function.
	$simulation->{'massMinimum'} = 10.0**(int(log($simulation->{'massParticle'})/log(10.0))+2);
	$simulation->{'massMaximum'} = 1.0e16;
	# Construct the simulation path.
	$simulation->{'path'} = $options{'simulationDataPath'};
	$simulation->{'path'} .= "/"
	    unless ( $simulation->{'path'} =~ m/\/$/ );
	$simulation->{'path'} .= $simulation->{'subpath'}."/".$simulation->{'label'}."/";
	# Add a single, null realization is the simulation has none.
	$simulation->{'realizations'} = [ "" ]
	    unless ( exists($simulation->{'realizations'}) );
	# Identify always isolated halos at z=0 and export these to IRATE format files.
	## Iterate over realizations.
	foreach my $realization ( @{$simulation->{'realizations'}} ) {
	    my $pathName = $simulation->{'path'}.($realization eq "" ? "" : $realization."/");
	    ## Iterate over subvolumes.
	    for(my $i=0;$i<$simulation->{'subvolumes'};++$i) {
		for(my $j=0;$j<$simulation->{'subvolumes'};++$j) {
		    for(my $k=0;$k<$simulation->{'subvolumes'};++$k) {
		        # Parse the base parameters.
			my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/haloMassFunctionIdentifyAlwaysIsolated.xml");
		        # Modify file names.
			$parameters->{'nbodyImporter'}                        ->{'fileName'}->{'value'} = $pathName."tree_"                   .$i."_".$j."_".$k.".dat" ;
			$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[3]->{'fileName'}->{'value'} = $pathName."alwaysIsolated_subVolume".$i."_".$j."_".$k.".hdf5";
		        # If a custom process function is defined, call it.
			&{$simulation->{'processIdentify'}}($simulation,$realization,$pathName,$parameters)
			    if ( exists($simulation->{'processIdentify'}) );
		        # Write parmeter file.
			my $parameterFileName = $pathName."identifyAlwaysIsolated_".$i."_".$j."_".$k.".xml";
			open(my $outputFile,">",$parameterFileName);
			print $outputFile $xml->XMLout($parameters, RootName => "parameters");
			close($outputFile);
		        # Skip if the file exists.
		        next
			    if ( -e $pathName."alwaysIsolated_subVolume".$i."_".$j."_".$k.".hdf5" );
		        # Generate a job.
			my $job;
			$job->{'command'   } =
			    "Galacticus.exe ".$parameterFileName;
			$job->{'launchFile'} = $pathName."identifyAlwaysIsolated_".$i."_".$j."_".$k.".sh" ;
			$job->{'logFile'   } = $pathName."identifyAlwaysIsolated_".$i."_".$j."_".$k.".log";
			$job->{'label'     } =           "identifyAlwaysIsolated_".$i."_".$j."_".$k       ;
			$job->{'ppn'       } = $ompThreads;
			$job->{'ompThreads'} = $ompThreads;
			$job->{'nodes'     } = 1 ;
			$job->{'mem'       } = "8G";
			$job->{'walltime'  } = "8:00:00";
			$job->{'mpi'       } = "no";
			push(@jobsIdentify,$job);
		    }
		}
	    }
	}
    }
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobsIdentify)
	if ( scalar(@jobsIdentify) > 0 );
}

# Perform any preprocessing.
{
    my $workDone  =  1;
    my $iteration = -1;
    while ( $workDone ) {
	my @preprocessingJobs;
	++$iteration;
	$workDone = 0;
	foreach my $simulation ( @simulations ) {
	    # Iterate over realizations.
	    foreach my $realization ( @{$simulation->{'realizations'}} ) {
		my $pathName = $simulation->{'path'}.($realization eq "" ? "" : $realization."/");
		# If a custom preprocess function is defined, call it.
		if ( exists($simulation->{'preprocessExtract'}) && scalar(@{$simulation->{'preprocessExtract'}}) > $iteration ) {
		    $workDone = 1;
		    &{$simulation->{'preprocessExtract'}->[$iteration]}($simulation,$realization,$pathName,\@preprocessingJobs);
		}
	    }
	}
	&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@preprocessingJobs)
	    if ( scalar(@preprocessingJobs) > 0 );
    }
}

# Iterate over simulations to extract snapshots
{
    my @jobsExtract;
    foreach my $simulation ( @simulations ) {
	# Identify always isolated halos at z=0 and export these to IRATE format files.
	## Iterate over realizations.
	foreach my $realization ( @{$simulation->{'realizations'}} ) {
	    my $pathName = $simulation->{'path'}.($realization eq "" ? "" : $realization."/");
	    ## Iterate over subvolumes.
	    for(my $i=0;$i<$simulation->{'subvolumes'};++$i) {
		for(my $j=0;$j<$simulation->{'subvolumes'};++$j) {
		    for(my $k=0;$k<$simulation->{'subvolumes'};++$k) {
			# Iterate over expansion factors.
			foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
			    my $redshift            =  1.0        /$expansionFactor-1.0;
			    my $expansionFactorLow  = (1.0-5.0e-4)*$expansionFactor    ;
			    my $expansionFactorHigh = (1.0+5.0e-4)*$expansionFactor    ;
			    my $redshiftLabel       = sprintf("z%5.3f",$redshift)      ;
			    # Always-Isolated halos.
			    {
				# Parse the base parameters.
				my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/haloMassFunctionExtractSnapshot.xml");
				# Modify file names.
				$parameters->{'nbodyImporter'}                        ->{'fileName'     }->{'value'} = $pathName."alwaysIsolated_subVolume"                   .$i."_".$j."_".$k.".hdf5";
				$parameters->{'nbodyImporter'}                        ->{'properties'   }->{'value'} = "particleID alwaysIsolated expansionFactor massVirial";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'rangeLow'     }->{'value'} = "1 ".$expansionFactorLow ;
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'rangeHigh'    }->{'value'} = "1 ".$expansionFactorHigh;
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'propertyNames'}->{'value'} = "alwaysIsolated expansionFactor";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'fileName'     }->{'value'} = $pathName."alwaysIsolated_".$redshiftLabel."_subVolume".$i."_".$j."_".$k.".hdf5";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'redshift'     }->{'value'} =                             $redshift                                           ;
				# If a custom process function is defined, call it.
				&{$simulation->{'processExtract'}}($simulation,$realization,$pathName,$parameters,$expansionFactor)
				    if ( exists($simulation->{'processExtract'}) );
				# Write parameter file.
				my $parameterFileName = $pathName."identifyAlwaysIsolated_".$redshiftLabel."_".$i."_".$j."_".$k.".xml";
				open(my $outputFile,">",$parameterFileName);
				print $outputFile $xml->XMLout($parameters, RootName => "parameters");
				close($outputFile);
				# Skip if the file exists.
				next
				    if ( -e $pathName."alwaysIsolated_".$redshiftLabel."_subVolume".$i."_".$j."_".$k.".hdf5" );
				# Generate a job.
				my $job;
				$job->{'command'   } =
				    "Galacticus.exe ".$parameterFileName;
				$job->{'launchFile'} = $pathName."identifyAlwaysIsolated_".$redshiftLabel."_".$i."_".$j."_".$k.".sh" ;
				$job->{'logFile'   } = $pathName."identifyAlwaysIsolated_".$redshiftLabel."_".$i."_".$j."_".$k.".log";
				$job->{'label'     } =           "identifyAlwaysIsolated_".$redshiftLabel."_".$i."_".$j."_".$k       ;
				$job->{'ppn'       } = 1;
				$job->{'ompThreads'} = 1;
				$job->{'nodes'     } = 1;
				$job->{'mem'       } = "8G";
				$job->{'walltime'  } = "8:00:00";
				$job->{'mpi'       } = "yes";
				push(@jobsExtract,$job);
			    }
			    # Non-flyby halos.
			    {
				# Parse the base parameters.
				my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/haloMassFunctionExtractSnapshot.xml");
				# Modify file names.
				$parameters->{'nbodyImporter'}                        ->{'fileName'     }->{'value'} = $pathName."alwaysIsolated_subVolume"             .$i."_".$j."_".$k.".hdf5";
				$parameters->{'nbodyImporter'}                        ->{'properties'   }->{'value'} = "particleID isFlyby expansionFactor massVirial";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'propertyNames'}->{'value'} = "isFlyby expansionFactor";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'rangeLow'     }->{'value'} = "0 ".$expansionFactorLow ;
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'rangeHigh'    }->{'value'} = "0 ".$expansionFactorHigh;
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'propertyNames'}->{'value'} = "isFlyby expansionFactor";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'fileName'     }->{'value'} = $pathName."nonFlyby_".$redshiftLabel."_subVolume".$i."_".$j."_".$k.".hdf5";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'redshift'     }->{'value'} =                       $redshift                                           ;
				# If a custom process function is defined, call it.
				&{$simulation->{'processExtract'}}($simulation,$realization,$pathName,$parameters,$expansionFactor)
				    if ( exists($simulation->{'processExtract'}) );
				# Write parameter file.
				my $parameterFileName = $pathName."identifyNonFlyby_".$redshiftLabel."_".$i."_".$j."_".$k.".xml";
				open(my $outputFile,">",$parameterFileName);
				print $outputFile $xml->XMLout($parameters, RootName => "parameters");
				close($outputFile);
				# Skip if the file exists.
				next
				    if ( -e $pathName."nonFlyby_".$redshiftLabel."_subVolume".$i."_".$j."_".$k.".hdf5" );
				# Generate a job.
				my $job;
				$job->{'command'   } =
				    "Galacticus.exe ".$parameterFileName;
				$job->{'launchFile'} = $pathName."identifyNonFlyby_".$redshiftLabel."_".$i."_".$j."_".$k.".sh" ;
				$job->{'logFile'   } = $pathName."identifyNonFlyby_".$redshiftLabel."_".$i."_".$j."_".$k.".log";
				$job->{'label'     } =           "identifyNonFlyby_".$redshiftLabel."_".$i."_".$j."_".$k       ;
				$job->{'ppn'       } = 1;
				$job->{'ompThreads'} = 1;
				$job->{'nodes'     } = 1;
				$job->{'mem'       } = "8G";
				$job->{'walltime'  } = "8:00:00";
				$job->{'mpi'       } = "yes";
				push(@jobsExtract,$job);
			    }
			    # All halos.
			    {
				# Parse the base parameters.
				my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/haloMassFunctionExtractSnapshot.xml");
				# Modify file names.
				$parameters->{'nbodyImporter'}                        ->{'fileName'     }->{'value'} = $pathName."alwaysIsolated_subVolume"        .$i."_".$j."_".$k.".hdf5";
				$parameters->{'nbodyImporter'}                        ->{'properties'   }->{'value'} = "particleID expansionFactor massVirial";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'propertyNames'}->{'value'} = "expansionFactor";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'rangeLow'     }->{'value'} = $expansionFactorLow ;
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'rangeHigh'    }->{'value'} = $expansionFactorHigh;
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'propertyNames'}->{'value'} = "expansionFactor";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'fileName'     }->{'value'} = $pathName."all_".$redshiftLabel."_subVolume".$i."_".$j."_".$k.".hdf5";
				$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'redshift'     }->{'value'} =                  $redshift                                           ;
				# If a custom process function is defined, call it.
				&{$simulation->{'processExtract'}}($simulation,$realization,$pathName,$parameters,$expansionFactor)
				    if ( exists($simulation->{'processExtract'}) );
				# Write parameter file.
				my $parameterFileName = $pathName."identifyAll_".$redshiftLabel."_".$i."_".$j."_".$k.".xml";
				open(my $outputFile,">",$parameterFileName);
				print $outputFile $xml->XMLout($parameters, RootName => "parameters");
				close($outputFile);
				# Skip if the file exists.
				next
				    if ( -e $pathName."all_".$redshiftLabel."_subVolume".$i."_".$j."_".$k.".hdf5" );
				# Generate a job.
				my $job;
				$job->{'command'   } =
				    "Galacticus.exe ".$parameterFileName;
				$job->{'launchFile'} = $pathName."identifyAll_".$redshiftLabel."_".$i."_".$j."_".$k.".sh" ;
				$job->{'logFile'   } = $pathName."identifyAll_".$redshiftLabel."_".$i."_".$j."_".$k.".log";
				$job->{'label'     } =           "identifyAll_".$redshiftLabel."_".$i."_".$j."_".$k       ;
				$job->{'ppn'       } = 1;
				$job->{'ompThreads'} = 1;
				$job->{'nodes'     } = 1;
				$job->{'mem'       } = "8G";
				$job->{'walltime'  } = "8:00:00";
				$job->{'mpi'       } = "yes";
				push(@jobsExtract,$job);
			    }
			}
		    }
		}
	    }
	}
    }
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobsExtract)
	if ( scalar(@jobsExtract) > 0 );
}

# Perform any postprocessing.
{
    my $workDone  =  1;
    my $iteration = -1;
    while ( $workDone ) {
	my @postprocessingJobs;
	++$iteration;
	$workDone = 0;
	foreach my $simulation ( @simulations ) {
	    # Iterate over realizations.
	    foreach my $realization ( @{$simulation->{'realizations'}} ) {
		my $pathName = $simulation->{'path'}.($realization eq "" ? "" : $realization."/");
		# If a custom postprocess function is defined, call it.
		if ( exists($simulation->{'postprocessExtract'}) && scalar(@{$simulation->{'postprocessExtract'}}) > $iteration ) {
		    $workDone = 1;
		    &{$simulation->{'postprocessExtract'}->[$iteration]}($simulation,$realization,$pathName,\@postprocessingJobs);
		}
	    }
	}
	&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@postprocessingJobs)
	    if ( scalar(@postprocessingJobs) > 0 );
    }
}

# Iterate over simulations to construct the mass functions.
{
    my @massFunctionJobs;
    foreach my $simulation ( @simulations ) {
	## Iterate over expansion factors.
	foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	    my $redshift            =  1.0/$expansionFactor-1.0;
	    my $redshiftLabel       = sprintf("z%5.3f",$redshift);
	    # Iterate over realizations.
	    foreach my $realization ( @{$simulation->{'realizations'}} ) {
		my $pathName = $simulation->{'path'}.($realization eq "" ? "" : $realization."/");
		# Iterate over halo types.
		foreach my $haloType ( "alwaysIsolated", "nonFlyby", "all" ) {
		    # Iterate over subvolumes.
		    my @nbodyImporters;
		    for(my $i=0;$i<$simulation->{'subvolumes'};++$i) {
			for(my $j=0;$j<$simulation->{'subvolumes'};++$j) {
			    for(my $k=0;$k<$simulation->{'subvolumes'};++$k) {
				# Add an importer for this subvolume.
				push(
				    @nbodyImporters,
				    {
					value      => "IRATE"                                                                                ,
					fileName   => {value => $pathName.$haloType."_".$redshiftLabel."_subVolume".$i."_".$j."_".$k.".hdf5"},
					snapshot   => {value => "1"},
					properties => {value => "massVirial"}
				    }
				    );
			    }
			}
		    }
		    ## Compute the mass function.
		    unless ( -e $pathName."haloMassFunction_".$haloType."_".$redshiftLabel.":MPI0000.hdf5" ) {
			## Parse the base parameters.
			my $massFunctionParameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/haloMassFunctionCompute.xml");
			## Modify parameters.
			@{$massFunctionParameters->{'nbodyImporter' }->{'nbodyImporter'}}                                          = @nbodyImporters;
			$massFunctionParameters  ->{'outputFileName'}                                                  ->{'value'} = $pathName."haloMassFunction_".$haloType."_".$redshiftLabel.".hdf5";
			$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[0]->{'values'             }->{'value'} = $simulation->{'massParticle'       };
			$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'description'        }->{'value'} = $simulation->{'description'        };
			$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'simulationReference'}->{'value'} = $simulation->{'simulationReference'};
			$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'simulationURL'      }->{'value'} = $simulation->{'simulationURL'      };
			$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'massMinimum'        }->{'value'} = $simulation->{'massMinimum'        };
			$massFunctionParameters  ->{'nbodyOperator' }->{'nbodyOperator'} ->[1]->{'massMaximum'        }->{'value'} = $simulation->{'massMaximum'        };
			## Write the parameter file.
			my $parameterFileName = $pathName."haloMassFunction_".$haloType."_".$redshiftLabel.".xml";
			open(my $outputFile,">",$parameterFileName);
			print $outputFile $xml->XMLout($massFunctionParameters, RootName => "parameters");
			close($outputFile);
			## Construct the job.
			my $job;
			$job->{'command'   } =
			    "Galacticus.exe ".$parameterFileName;
			$job->{'launchFile'} = $pathName."haloMassFunction_".$haloType."_".$redshiftLabel.".sh" ;
			$job->{'logFile'   } = $pathName."haloMassFunction_".$haloType."_".$redshiftLabel.".log";
			$job->{'label'     } =           "haloMassFunction_".$haloType."_".$redshiftLabel       ;
			$job->{'ppn'       } = $ompThreads;
			$job->{'ompThreads'} = $ompThreads;
			$job->{'nodes'     } = 1;
			$job->{'mem'       } = "8G";
			$job->{'walltime'  } = "8:00:00";
			$job->{'mpi'       } = "no";
			push(@massFunctionJobs,$job);
		    }
		}
	    }
	}
    }
    &{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@massFunctionJobs)
	if ( scalar(@massFunctionJobs) > 0 );
}

# Perform any postprocessing.
{
    my $workDone  =  1;
    my $iteration = -1;
    while ( $workDone ) {
	my @postprocessingJobs;
	++$iteration;
	$workDone = 0;
	foreach my $simulation ( @simulations ) {
	    # Iterate over realizations.
	    foreach my $realization ( @{$simulation->{'realizations'}} ) {
		my $pathName = $simulation->{'path'}.($realization eq "" ? "" : $realization."/");
		# Iterate over halo types.
		foreach my $haloType ( "alwaysIsolated", "nonFlyby", "all" ) {
		    # If a custom postprocess function is defined, call it.
		    if ( exists($simulation->{'postprocessMassFunction'}) && scalar(@{$simulation->{'postprocessMassFunction'}}) > $iteration ) {
			$workDone = 1;
			&{$simulation->{'postprocessMassFunction'}->[$iteration]}($simulation,$realization,$haloType,$pathName,\@postprocessingJobs);
		    }
		}
	    }
	}
	&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@postprocessingJobs)
	    if ( scalar(@postprocessingJobs) > 0 );
    }
}

# Move the resulting mass functions.
foreach my $simulation ( @simulations ) {
    # Iterate over realizations.
    foreach my $realization ( @{$simulation->{'realizations'}} ) {
	my $pathName = $simulation->{'path'}.($realization eq "" ? "" : $realization."/");
	# Iterate over expansion factors.
	foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	    my $redshift            =  1.0/$expansionFactor-1.0;
	    my $redshiftLabel       = sprintf("z%5.3f",$redshift);
	    copy($pathName."haloMassFunction_nonFlyby_".$redshiftLabel.":MPI0000.hdf5",$ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$simulation->{'label'}.($realization eq "" ? "" : "_".$realization)."_".$redshiftLabel.".hdf5");
	}
    }
}

exit 0;

sub zoomInsProcessIdentify {
    # Set the appropriate cosmology.
    my $simulation  = shift();
    my $realization = shift();
    my $pathName    = shift();
    my $parameters  = shift();
    my $expansionFactor = shift();
    # Find the host halo IDs for this simulation.
    my $xml         = new XML::Simple();
    my $hostHaloIDs = $xml->XMLin("constraints/pipelines/darkMatter/symphonyZoomInHostHaloIDs.xml");
    @{$simulation->{'hostHaloIDs'}} = map {$hostHaloIDs->{$simulation->{'label'}}->{$_}} @{$simulation->{'realizations'}};
    # Set the relevant cosmology.
    $parameters->{'cosmologyParameters'}->{'HubbleConstant' }->{'value'} = 70.000;
    $parameters->{'cosmologyParameters'}->{'OmegaMatter'    }->{'value'} =  0.286;
    $parameters->{'cosmologyParameters'}->{'OmegaDarkEnergy'}->{'value'} =  0.714;
    $parameters->{'cosmologyParameters'}->{'OmegaBaryon'    }->{'value'} =  0.047;
    # Add read of additional columns.
    my @propertiesImport = split(" ",$parameters->{'nbodyImporter'}->{'readColumns'}->{'value'});
    foreach my $property ( "X", "Y", "Z", "Rvir", "rs" ) {
	$parameters->{'nbodyImporter'}->{'readColumns'}->{'value'} .= " ".$property
	    unless ( grep {$_ eq $property} @propertiesImport );
    }
    # Remove any delete of these properties.
    my @propertiesDeleted = split(" ",$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'propertyNames'}->{'value'});
    my @propertiesToDelete;
    foreach my $property ( @propertiesDeleted ) {
	push(@propertiesToDelete,$property)
	    unless
	    (
	     $property eq "position"
	     ||
	     $property eq "radiusScale"
	     ||
	     $property eq "radiusVirial"
	    );
    }
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[2]->{'propertyNames'}->{'value'} = join(" ",@propertiesToDelete);
    # Add physical to comoving conversion.
    splice(
	@{$parameters->{'nbodyOperator'}->{'nbodyOperator'}},
	2,0,
	{
	    value => "physicalToComoving"
	}
	);
}

sub zoomInsPreProcessExtractLocate {
    # Identify the primary progenitor.
    my $simulation  = shift();
    my $realization = shift();
    my $pathName    = shift();
    my $jobs        = shift();
    # Find the host halo ID for this realization.
    my $hostHaloID;
    if ( exists($simulation->{'hostHaloIDs'}) ) {
	($hostHaloID) = map {$simulation->{'realizations'}->[$_] eq $realization ? $simulation->{'hostHaloIDs'}->[$_] : ()} 0..$#{$simulation->{'realizations'}};
    } else {
	$hostHaloID = -1;
    }
    # Iterate over expansion factors.
    my $job;
    foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	# Look for an existing record of the primary halo.
	my $redshift               =  1.0/$expansionFactor-1.0;
	my $redshiftLabel          = sprintf("z%5.3f",$redshift);
	my $primaryHaloFileName    = $pathName."primaryHalo_".$redshiftLabel.".xml";
	my $xml                    = new XML::Simple();
	unless ( -e $primaryHaloFileName ) {
	    # Generate a job to extract the primary halo.
	    unless ( defined($job) ) {
		$job->{'command'   } = "";
		$job->{'launchFile'} = $pathName."zoomInExtract_".$simulation->{'label'}.".sh" ;
		$job->{'logFile'   } = $pathName."zoomInExtract_".$simulation->{'label'}.".log";
		$job->{'label'     } =           "zoomInExtract_".$simulation->{'label'}       ;
		$job->{'ppn'       } = 1;
		$job->{'ompThreads'} = 1;
		$job->{'nodes'     } = 1;
		$job->{'mem'       } = "32G";
		$job->{'walltime'  } = "8:00:00";
		$job->{'mpi'       } = "no";
	    }
	    $job->{'command'} .= $ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/haloMassFunctionZoomInExtract.pl ".$pathName." ".$primaryHaloFileName." ".$expansionFactor." ".$simulation->{'hubbleConstant'}." ".$simulation->{'massParticle'}." ".$simulation->{'massHostLogMin'}." ".$simulation->{'massHostLogMax'}." ".$hostHaloID."\n";
	}
    }   
    push(@{$jobs},$job)
	if ( defined($job) );
}

sub zoomInsPreProcessExtractUncontaminated {
    # Find the uncontaminated region.
    my $simulation  = shift();
    my $realization = shift();
    my $pathName    = shift();
    my $jobs        = shift();
    # Iterate over expansion factors.
    foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	# Look for an existing record of the primary halo.
	my $redshift             =  1.0/$expansionFactor-1.0;
	my $redshiftLabel        = sprintf("z%5.3f",$redshift);
	my $primaryHaloFileName  = $pathName."primaryHalo_".$redshiftLabel.".xml";
	# Find properties of the central halo.
	my $xml                  = new XML::Simple();
	my $primaryHaloData      = $xml->XMLin($primaryHaloFileName);
	my $expansionFactorLabel = sprintf("sphericalOrigin:a%5.3f",$expansionFactor);
	$simulation->{$realization}->{$expansionFactorLabel} =
	    $primaryHaloData->{'x'}." ".
	    $primaryHaloData->{'y'}." ".
	    $primaryHaloData->{'z'}    ;
	$simulation->{'massCentral'} = pdl $primaryHaloData->{'mc'};
	# Generate job to find the uncontaminated region.
	my $uncontaminatedFileName = $pathName."uncontaminated_".$redshiftLabel.".hdf5";
	unless ( -e $uncontaminatedFileName ) {
	    my $parametersUncontaminated = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/zoomInSelectUncontaminated.xml");
	    # Modify file names.
	    (my $snapshot) = map {$simulation->{'expansionFactors'}->[$_] == $expansionFactor ? $simulation->{'snapshots'}->[$_] : ()} 0..$#{$simulation->{'expansionFactors'}};
	    $parametersUncontaminated->{'nbodyImporter'}                        ->{'fileName'      }->{'value'} = $pathName."snapshots/snapshot_".$snapshot;
	    $parametersUncontaminated->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'point'         }->{'value'} = $simulation->{$realization}->{$expansionFactorLabel};
	    $parametersUncontaminated->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'fileName'      }->{'value'} = $uncontaminatedFileName;
	    $parametersUncontaminated                                           ->{'outputFileName'}->{'value'} = $pathName."uncontaminatedExtract_".$redshiftLabel.".hdf5";
	    # Write parameter file.
	    my $parameterFileName = $pathName."uncontaminated_".$redshiftLabel.".xml";
	    open(my $outputFile,">",$parameterFileName);
	    print $outputFile $xml->XMLout($parametersUncontaminated, RootName => "parameters");
	    close($outputFile);
	    # Generate a job to extract the uncontaminated region around the primary halo.
	    my $job;
	    $job->{'command'   } = $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$pathName."uncontaminated_".$redshiftLabel.".xml";
	    $job->{'launchFile'} = $pathName."uncontaminatedExtract_".$simulation->{'label'}."_".$redshiftLabel.".sh" ;
	    $job->{'logFile'   } = $pathName."uncontaminatedExtract_".$simulation->{'label'}."_".$redshiftLabel.".log";
	    $job->{'label'     } =           "uncontaminatedExtract_".$simulation->{'label'}."_".$redshiftLabel       ;
	    $job->{'ppn'       } = $ompThreads;
	    $job->{'ompThreads'} = $ompThreads;
	    $job->{'nodes'     } = 1;
	    $job->{'mem'       } = "16G";
	    $job->{'walltime'  } = "8:00:00";
	    $job->{'mpi'       } = "no";
	    push(@{$jobs},$job)
	}
    }
}

sub zoomInsProcessExtract {
    # Determine the central point, and extent of the high-resolution region in a zoom in simulation.
    my $simulation      = shift();
    my $realization     = shift();
    my $pathName        = shift();
    my $parameters      = shift();
    my $expansionFactor = shift();
    print "Processing extraction for ".$simulation->{'label'}." : ".$realization."\n";
    # Set the relevant cosmology.
    $parameters->{'cosmologyParameters'}->{'HubbleConstant' }->{'value'} = 70.000;
    $parameters->{'cosmologyParameters'}->{'OmegaMatter'    }->{'value'} =  0.286;
    $parameters->{'cosmologyParameters'}->{'OmegaDarkEnergy'}->{'value'} =  0.714;
    $parameters->{'cosmologyParameters'}->{'OmegaBaryon'    }->{'value'} =  0.047;
    # Set path names.
    my $redshift               =  1.0/$expansionFactor-1.0;
    my $redshiftLabel          = sprintf("z%5.3f",$redshift);
    my $expansionFactorLabel   = sprintf("sphericalOrigin:a%5.3f",$expansionFactor);
    my $primaryHaloFileName    = $pathName."primaryHalo_".$redshiftLabel.".xml";
    my $uncontaminatedFileName = $pathName."uncontaminated_".$redshiftLabel.".hdf5";
    # Find properties of the central halo.
    my $xml                    = new XML::Simple();
    my $primaryHaloData        = $xml->XMLin($primaryHaloFileName);
    # Find uncontaminated region.
    my $uncontaminatedFile      = new PDL::IO::HDF5($uncontaminatedFileName);
    my $uncontaminatedParticles = $uncontaminatedFile     ->group  ('Snapshot00001/HaloCatalog');
    (my $radiusUncontaminated)  = $uncontaminatedParticles->attrGet('radiusUncontaminated'     );
    $primaryHaloData->{'ur'} = $radiusUncontaminated->sclr();
    open(my $primaryHaloDataFile,">",$primaryHaloFileName);
    print $primaryHaloDataFile $xml->XMLout($primaryHaloData, RootName => "primaryHalo");
    close($primaryHaloDataFile);
    # Add read of (x,y,z) coordinate columns, and subsequent delete.
    $parameters->{'nbodyImporter'}                        ->{'properties'   }->{'value'} .= " position"         ;
    $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'propertyNames'}->{'value'} .= " distanceFromPoint";
    # Add calculation of distance from primary halo.
    $simulation->{$realization}->{$redshiftLabel}->{'sphericalRadiusMinimum'} = 0.0;
    $simulation->{$realization}->{$redshiftLabel}->{'sphericalRadiusMaximum'} = $radiusUncontaminated->sclr();
    splice(
	@{$parameters->{'nbodyOperator'}->{'nbodyOperator'}},
	1,0,
	{
	    value         => "distanceFromPoint",
	    point         => {value => $simulation->{$realization}->{$expansionFactorLabel  }}
	},
	{
	    value         =>           "filterProperties"                                                        ,
	    propertyNames => {value => "distanceFromPoint"                                                      },
	    rangeLow      => {value => $simulation->{$realization}->{$redshiftLabel}->{'sphericalRadiusMinimum'}},
	    rangeHigh     => {value => $simulation->{$realization}->{$redshiftLabel}->{'sphericalRadiusMaximum'}},
	},
	{
	    value         => "setBoxSize",
	    boxSize       => {value => $primaryHaloData->{'l'}}
	}
	);
}

sub zoomInsPostprocessSelectInSphere {
    # Select all particles within the sphere of interest.
    my $simulation  = shift();
    my $realization = shift();
    my $pathName    = shift();
    my $jobs        = shift();
    # Generate a parameter file to extract particles within the spherical volume of interest.
    my $xml        = new XML::Simple();
    my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/zoomInSelectInSphere.xml");
    # Iterate over expansion factors.
    foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	my $redshift             =  1.0/$expansionFactor-1.0;
	my $redshiftLabel        = sprintf("z%5.3f"                ,$redshift       );
	my $expansionFactorLabel = sprintf("sphericalOrigin:a%5.3f",$expansionFactor);
	(my $snapshot)           = map {$simulation->{'expansionFactors'}->[$_] == $expansionFactor ? $simulation->{'snapshots'}->[$_] : ()} 0..$#{$simulation->{'expansionFactors'}};
	$parameters->{'nbodyImporter'}                        ->{'fileName' }->{'value'} = $pathName."snapshots/snapshot_".$snapshot;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'fileName' }->{'value'} = $pathName."selectedParticles_".$redshiftLabel.".hdf5";
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'redshift' }->{'value'} =                                $redshift             ;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'point'    }->{'value'} = $simulation->{$realization}->{$expansionFactorLabel}                            ;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'rangeLow' }->{'value'} = $simulation->{$realization}->{$redshiftLabel       }->{'sphericalRadiusMinimum'};
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'rangeHigh'}->{'value'} = $simulation->{$realization}->{$redshiftLabel       }->{'sphericalRadiusMaximum'};
	my $parameterFileName = $pathName."selectInSphere_".$redshiftLabel.".xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($parameters, RootName => "parameters");
	close($outputFile);
	my $job;
	$job->{'command'   } =
	    "Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $pathName."selectInSphere_".$redshiftLabel.".sh" ;
	$job->{'logFile'   } = $pathName."selectInSphere_".$redshiftLabel.".log";
	$job->{'label'     } =           "selectInSphere_".$redshiftLabel       ;
	$job->{'ppn'       } = $ompThreads;
	$job->{'ompThreads'} = $ompThreads;
	$job->{'nodes'     } =  1;
	$job->{'mem'       } = "16G";
	$job->{'walltime'  } = "8:00:00";
	$job->{'mpi'       } = "no";
	push(@{$jobs},$job)
	    unless ( -e $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[4]->{'fileName' }->{'value'} );
    }
}

sub zoomInsPostprocessExtractSelectedIDs {
    # Extract the IDs of selected particles.
    my $simulation  = shift();
    my $realization = shift();
    my $pathName    = shift();
    my $jobs        = shift();
    # Iterate over expansion factors.
    foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	my $redshift            =  1.0/$expansionFactor-1.0;
	my $redshiftLabel       = sprintf("z%5.3f",$redshift);
	# Read particle IDs.
	next
	    if ( -e $pathName."selectedParticlesIDs_".$redshiftLabel.".hdf5" );
	my $selectedParticleFile = new PDL::IO::HDF5(    $pathName."selectedParticles_"   .$redshiftLabel.".hdf5");
	my $selectedIDFile       = new PDL::IO::HDF5(">".$pathName."selectedParticlesIDs_".$redshiftLabel.".hdf5");
	my $selectedIDs          = $selectedParticleFile->group('Snapshot00001')->group('HaloCatalog')->dataset('HaloID')->get();
	$selectedIDFile->dataset('id')->set($selectedIDs);
    }
}

sub zoomInsPostprocessSelectInICs {
    # Identify selected particles in the initial conditions.
    my $simulation  = shift();
    my $realization = shift();
    my $pathName    = shift();
    my $jobs        = shift();
    # Find the redshift of the ICs.
    my $musicFileName = $pathName."music.conf";
    die("missing file: ".$musicFileName)
	unless ( -e $musicFileName );
    open(my $musicFile,$musicFileName);
    while ( my $line = <$musicFile> ) {
	if ( $line =~ m/^zstart\s*=\s*([\d\.]+)/ ) {
	    $simulation->{'redshiftICs'} = $1;
	}
    }
    close($musicFile);
    # Generate a parameter file to extract selected particles in the ICs.
    my $xml        = new XML::Simple();
    my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/zoomInSelectInICs.xml");
    # Iterate over expansion factors.
    foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	my $redshift            =  1.0/$expansionFactor-1.0;
	my $redshiftLabel       = sprintf("z%5.3f",$redshift);
	$parameters->{'nbodyImporter'}                        ->{'fileName'           }->{'value'} =                $pathName                   ."ic/ic_gadget_dist"                                ;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[0]->{'idSelectionFileName'}->{'value'} =                $pathName                   ."selectedParticlesIDs_".$redshiftLabel.    ".hdf5" ;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'fileName'           }->{'value'} =                $pathName                   ."selectedParticles_"   .$redshiftLabel."_ICs.hdf5" ;
	$parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'redshift'           }->{'value'} = sprintf("%.3f",$simulation->{'redshiftICs'}                                                   );
	$parameters                                           ->{'outputFileName'     }->{'value'} =                $pathName                   ."selectInICs_"         .$redshiftLabel.    ".hdf5" ;
	my $parameterFileName = $pathName."selectInICs_".$redshiftLabel.".xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($parameters, RootName => "parameters");
	close($outputFile);
	my $job;
	$job->{'command'   } =
	    "Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $pathName."selectInICs_".$redshiftLabel.".sh" ;
	$job->{'logFile'   } = $pathName."selectInICs_".$redshiftLabel.".log";
	$job->{'label'     } =           "selectInICs_".$redshiftLabel.""    ;
	$job->{'ppn'       } = $ompThreads;
	$job->{'ompThreads'} = $ompThreads;
	$job->{'nodes'     } =  1;
	$job->{'mem'       } = "16G";
	$job->{'walltime'  } = "8:00:00";
	$job->{'mpi'       } = "no";
	push(@{$jobs},$job)
	    unless ( -e $parameters->{'nbodyOperator'}->{'nbodyOperator'}->[1]->{'fileName' }->{'value'} );
    }
}

sub zoomInsPostprocessAnalyze {
    # Measure mass and overdensity of the selected region.
    my $simulation  = shift();
    my $realization = shift();
    my $pathName    = shift();
    my $jobs        = shift();
    # Generate a parameter file to measure mass and overdensity.
    my $xml        = new XML::Simple();
    my $parameters = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/constraints/pipelines/darkMatter/zoomInAnalyze.xml");
    # Iterate over expansion factors.
    foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	my $redshift            =  1.0/$expansionFactor-1.0;
	my $redshiftLabel       = sprintf("z%5.3f",$redshift);
	$parameters->{'outputFileName'}                                                         ->{'value'} =                $pathName                   ."environment_"      .$redshiftLabel.".hdf5"     ;
	$parameters->{'nbodyImporter' }->{'nbodyImporter'}->[0]                   ->{'fileName'}->{'value'} =                $pathName                   ."selectedParticles_".$redshiftLabel."_ICs.hdf5" ;
	$parameters->{'nbodyImporter' }->{'nbodyImporter'}->[1]                   ->{'fileName'}->{'value'} =                $pathName                   ."ic/ic_gadget_dist"                             ;
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[0]                   ->{'values'  }->{'value'} = sprintf("%.3f",$simulation->{'redshiftICs'}                                                );
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[5]->{'nbodyOperator'}->{'fileName'}->{'value'} =                $pathName                   ."allParticles_"     .$redshiftLabel."_ICs.hdf5" ;
	$parameters->{'nbodyOperator' }->{'nbodyOperator'}->[5]->{'nbodyOperator'}->{'redshift'}->{'value'} = sprintf("%.3f",$simulation->{'redshiftICs'}                                                );
	my $parameterFileName = $pathName."analyze_".$redshiftLabel.".xml";
	open(my $outputFile,">",$parameterFileName);
	print $outputFile $xml->XMLout($parameters, RootName => "parameters");
	close($outputFile);
	my $job;
	$job->{'command'   } =
	    "Galacticus.exe ".$parameterFileName;
	$job->{'launchFile'} = $pathName."analyze_".$redshiftLabel.".sh" ;
	$job->{'logFile'   } = $pathName."analyze_".$redshiftLabel.".log";
	$job->{'label'     } =           "analyze_".$redshiftLabel       ;
	$job->{'ppn'       } = $ompThreads;
	$job->{'ompThreads'} = $ompThreads;
	$job->{'nodes'     } =  1;
	$job->{'mem'       } = "32G";
	$job->{'walltime'  } = "8:00:00";
	$job->{'mpi'       } = "no";
	push(@{$jobs},$job)
	    unless ( -e $pathName."environment_".$redshiftLabel.":MPI0000.hdf5" );
    }
}

sub zoomInsPostprocessSetVolume {
    # Set the effective Lagrangian volume for the selected region.
    my $simulation              = shift();
    my $realization             = shift();
    my $pathName                = shift();
    my $jobs                    = shift();
    # Set cosmological parameters.
    my $HubbleConstant          = pdl 70.000;
    my $OmegaMatter             = pdl  0.286;
    # Iterate over expansion factors.
    foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	my $redshift            =  1.0/$expansionFactor-1.0;
	my $redshiftLabel       = sprintf("z%5.3f",$redshift);
	# Extract the mass of the region.
	my $analysisFile            = new PDL::IO::HDF5($pathName."selectedParticles_".$redshiftLabel.".hdf5");
	die("error: 'massTotal' attribute is missing from file '".$pathName."selectedParticles_".$redshiftLabel.".hdf5'")
	    unless ( grep {$_ eq "massTotal"} $analysisFile->group('Snapshot00001')->group('HaloCatalog')->attrs() );
	(my $mass)                  = $analysisFile->group('Snapshot00001')->group('HaloCatalog')->attrGet('massTotal');
	# Compute the mean density of the universe.
	my $gravitationalConstant   = pdl 4.3011827419096073e-9;
	my $densityMean             = 3.0*$OmegaMatter*$HubbleConstant**2/8.0/PI/$gravitationalConstant;
	# Compute the volume of a cube containing the mass of our selected region.
	my $boxSize                                                     = (    $mass/$densityMean       )**(1.0/3.0);
	$simulation->{$realization}->{$redshiftLabel}->{'radiusRegion'} = (3.0*$mass/$densityMean/4.0/PI)**(1.0/3.0);
	$simulation->{$realization}->{$redshiftLabel}->{'massRegion'  } =      $mass                                ;
	# Iterate over halo types.
	foreach my $haloType ( "alwaysIsolated", "nonFlyby", "all" ) {
	    # Set the box size.
	    my $halosFile = new PDL::IO::HDF5(">".$pathName.$haloType."_".$redshiftLabel."_subVolume0_0_0.hdf5");
	    $halosFile->group('Snapshot00001'       )->group('HaloCatalog')->attrSet(boxSize => $boxSize);
	    $halosFile->group('SimulationProperties')                      ->attrSet(boxSize => $boxSize);
	}
    }
}

sub zoomInsPostProcessMassFunction {
    # Store the region mass and overdensity in the mass function file.
    my $simulation           = shift();
    my $realization          = shift();
    my $haloType             = shift();
    my $pathName             = shift();
    my $jobs                 = shift();
    # Iterate over expansion factors.
    foreach my $expansionFactor ( @{$simulation->{'expansionFactors'}} ) {
	my $redshift            =  1.0/$expansionFactor-1.0;
	my $redshiftLabel       = sprintf("z%5.3f",$redshift);
	# Extract overdensity from the analysis file.
	my $analysisFile         = new PDL::IO::HDF5($pathName."environment_".$redshiftLabel.":MPI0000.hdf5");
	die("error: 'simulation0002' group is missing from file '".$pathName."environment_".$redshiftLabel.":MPI0000.hdf5'")
	    unless ( grep {$_ eq "simulation0002"} $analysisFile->groups() );
	foreach my $attributeName ( 'massTotal', 'convexHullOverdensity' ) {
	    die("error: '".$attributeName."' attribute is missing from file '".$pathName."environment_".$redshiftLabel.":MPI0000.hdf5'")
		unless ( grep {$_ eq $attributeName} $analysisFile->group('simulation0002')->attrs() );
	}
	(my $mass)               = $analysisFile->group('simulation0002')->attrGet('massTotal'            );
	(my $overdensity)        = $analysisFile->group('simulation0002')->attrGet('convexHullOverdensity');
	# Store these to the halo mass function.
	my $haloMassFunctionFile = new PDL::IO::HDF5(">".$pathName."haloMassFunction_".$haloType."_".$redshiftLabel.":MPI0000.hdf5");
	my $simulationGroup      = $haloMassFunctionFile->group('simulation0001')                                     ;
	$simulationGroup->attrSet(massRegion             => $simulation->{$realization}->{$redshiftLabel}->{'massRegion'  });
	$simulationGroup->attrSet(radiusRegion           => $simulation->{$realization}->{$redshiftLabel}->{'radiusRegion'});
	$simulationGroup->attrSet(overdensityEnvironment => $overdensity                                                   );
	$simulationGroup->attrSet(massEnvironment        => $mass                                                          );
    }
}
