# Contains a Perl module which handles building the base types of the class hierarchy during build.

package Galacticus::Build::Components::Components;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::Components::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     baseTypes => 
     {
	 types =>
	     [
	      \&Build_Node_Component_Class
	     ]
     }
    );

sub Build_Node_Component_Class {
    # Build the "nodeComponent" base class from which all other component classes and implementations inherit.
    my $build = shift();
    # Define type-bound functions.
    my @typeBoundFunctions = 
	(
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "type"                                                                                                 ,
	     function    => "Node_Component_Generic_Type"                                                                          ,
	     description => "Return the type of this object."                                                                      ,
	     returnType  => "\\textcolor{red}{\\textless type(varying\\_string)\\textgreater}"                                     ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "host"                                                                                                 ,
	     function    => "Node_Component_Host_Node"                                                                             ,
	     description => "Return a pointer to the host {\\normalfont \\ttfamily treeNode} object."                              ,
	     returnType  => "\\textcolor{red}{\\textless *type(treeNode)\\textgreater}"                                            ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "destroy"                                                                                              ,
	     function    => "Node_Component_Generic_Destroy"                                                                       ,
	     description => "Destroy the object."                                                                                  ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "addMetaProperty"                                                                                      ,
	     function    => "Node_Component_Generic_Add_Meta_Property"                                                             ,
	     description => "Add a meta-property to this class."                                                                   ,
	     returnType  => "\\intzero"                                                                                            ,
	     arguments   => "\\textcolor{red}{\\textless type(varying\_string)\\textgreater label, \\textcolor{red}{\\textless character(len=*)\\textgreater name, \\logicalzero\ [isEvolvable], \\logicalzero\ [isCreator]"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "addRank1MetaProperty"                                                                                      ,
	     function    => "Node_Component_Generic_Add_Rank1_Meta_Property"                                                             ,
	     description => "Add a rank-1 meta-property to this class."                                                                   ,
	     returnType  => "\\intzero"                                                                                            ,
	     arguments   => "\\textcolor{red}{\\textless type(varying\_string)\\textgreater label, \\textcolor{red}{\\textless character(len=*)\\textgreater name, \\logicalzero\ [isEvolvable], \\logicalzero\ [isCreator]"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "addIntegerMetaProperty"                                                                               ,
	     function    => "Node_Component_Generic_Add_Integer_Meta_Property"                                                     ,
	     description => "Add an integer meta-property to this class."                                                          ,
	     returnType  => "\\intzero"                                                                                            ,
	     arguments   => "\\textcolor{red}{\\textless type(varying\_string)\\textgreater label, \\textcolor{red}{\\textless character(len=*)\\textgreater name, \\logicalzero\ [isCreator]"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeCount"                                                                                       ,
	     function    => "Node_Component_Serialize_Count_Zero"                                                                  ,
	     description => "Return a count of the number of evolvable quantities to be evolved."                                  ,
	     returnType  => "\\intzero"                                                                                            ,
	     arguments   => "\\intzero\\ propertyType\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializationOffsets"                                                                                 ,
	     function    => "Node_Component_Serialization_Offsets"                                                                 ,
	     description => "Set offsets into serialization arrays."                                                               ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeValues"                                                                                      ,
	     function    => "Node_Component_Serialize_Null"                                                                        ,
	     description => "Serialize the evolvable quantities to an array."                                                      ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\doubleone\\ array\\argout, \\intzero\\ propertyType\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "deserializeRaw"                                                                                       ,
	     function    => "Node_Component_Read_Raw_Null"                                                                         ,
	     description => "Read properties from raw file."                                                                       ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ fileHandle\\argin"
	 },	 
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "deserializeValues"                                                                                    ,
	     function    => "Node_Component_Deserialize_Null"                                                                      ,
	     description => "Deserialize the evolvable quantities from an array."                                                  ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\doubleone\\ array\\argin, \\intzero\\ propertyType\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "odeStepRatesInitialize"                                                                               ,
	     function    => "Node_Component_ODE_Step_Initialize_Null"                                                              ,
	     description => "Initialize rates for evolvable properties."                                                           ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "odeStepScalesInitialize"                                                                              ,
	     function    => "Node_Component_ODE_Step_Initialize_Null"                                                              ,
	     description => "Initialize scales for evolvable properties."                                                          ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""          
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeASCII"                                                                                       ,
	     function    => "Node_Component_Dump_Null"                                                                             ,
	     description => "Generate an ASCII dump of all properties."                                                            ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeXML"                                                                                         ,
	     function    => "Node_Component_Dump_XML_Null"                                                                         ,
	     description => "Generate an XML dump of all properties."                                                              ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => ""
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "serializeRaw"                                                                                         ,
	     function    => "Node_Component_Dump_Raw_Null"                                                                         ,
	     description => "Generate a binary dump of all properties."                                                            ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ fileHandle\\argin"                   
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "outputCount"                                                                                          ,
	     function    => "Node_Component_Output_Count_Null"                                                                     ,
	     description => "Compute a count of outputtable properties."                                                           ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ integerPropertyCount\\arginout, \\intzero\\ doublePropertyCount\\arginout, \\doublezero\\ time\\argin, \\intzero\\ instance\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "outputNames"                                                                                          ,
	     function    => "Node_Component_Output_Names_Null"                                                                     ,
	     description => "Generate names of outputtable properties."                                                            ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ integerProperty\\arginout, \\textcolor{red}{\\textless type(outputPropertyInteger)(:)\\textgreater} integerProperties\\arginout, \\intzero\\ doubleProperty\\arginout, \\textcolor{red}{\\textless type(otuputPropertyDouble)(:)\\textgreater} doubleProperties\\arginout, \\doublezero\\ time\\argin, \\intzero\\ instance\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "output"                                                                                               ,
	     function    => "Node_Component_Output_Null"                                                                           ,
	     description => "Generate values of outputtable properties."                                                           ,
	     returnType  => "\\void"                                                                                               ,
	     arguments   => "\\intzero\\ integerProperty\\arginout, \\intzero\\ integerBufferCount\\arginout, \\textcolor{red}{\\textless type(outputPropertyInteger)(:)\\textgreater} integerProperties\\arginout, \\intzero doubleProperty\\arginout, \\intzero\\ doubleBufferCount\\arginout, \\textcolor{red}{\\textless type(otuputPropertyDouble)(:)\\textgreater} doubleProperties\\arginout, \\doublezero\\ time\\argin, \\intzero\\ instance\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "enclosedMass"                                                                                         ,
	     function    => "Node_Component_Enclosed_Mass_Null"                                                                    ,
	     description => "Compute the mass enclosed within a radius."                                                           ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "acceleration"                                                                                         ,
	     function    => "Node_Component_Acceleration_Null"                                                                     ,
	     description => "Compute the gravitational acceleration at a point."                                                   ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doubleone"                                                                                          ,
	     arguments   => "\\doubleone\\ position\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "chandrasekharIntegral"                                                                                ,
	     function    => "Node_Component_Chandrasekhar_Integral_Null"                                                           ,
	     description => "Compute the Chandrasekhar integral for a given position and velocity."                                ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doubleone"                                                                                          ,
	     arguments   => "\\doubleone\\ position\\argin, \\doubleone\\ velocity\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "tidalTensor"                                                                                          ,
	     function    => "Node_Component_Tidal_Tensor_Null"                                                                     ,
	     description => "Compute the gravitational tidal tensor at a point."                                                   ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\textcolor{red}{\\textless type(tensorRank2Dimension3Symmetric)\\textgreater}"                       ,
	     arguments   => "\\doubleone\\ position\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "density"                                                                                              ,
	     function    => "Node_Component_Density_Null"                                                                          ,
	     description => "Compute the density."                                                                                 ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\textcolor{red}{\\textless double(3)\\textgreater} positionSpherical\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "surfaceDensity"                                                                                       ,
	     function    => "Node_Component_Surface_Density_Null"                                                                  ,
	     description => "Compute the surface density."                                                                         ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\textcolor{red}{\\textless double(3)\\textgreater} positionCylindrical\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin, \\enumWeightBy\\ [weightBy]\\argin, \\intzero\\ [weightIndex]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "potential"                                                                                            ,
	     function    => "Node_Component_Potential_Null"                                                                        ,
	     description => "Compute the gravitational potential."                                                                 ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "rotationCurve"                                                                                        ,
	     function    => "Node_Component_Rotation_Curve_Null"                                                                   ,
	     description => "Compute the rotation curve."                                                                          ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin"
	 },
	 {
	     type        => "procedure"                                                                                            ,
	     name        => "rotationCurveGradient"                                                                                ,
	     function    => "Node_Component_Rotation_Curve_Gradient_Null"                                                          ,
	     description => "Compute the rotation curve gradient."                                                                 ,
	     mappable    => "summation"                                                                                            ,
	     returnType  => "\\doublezero"                                                                                         ,
	     arguments   => "\\doublezero\\ radius\\argin, \\enumComponentType\\ [componentType]\\argin, \\enumMassType\\ [massType]\\argin"
	 }
	);
    # Specify the data content.
    my @dataContent =
	(
	 {
	     intrinsic  =>   "type",
	     type       =>   "treeNode"            ,
	     attributes => [ "pointer" , "public" ],
	     variables  => [ "hostNode => null()" ]
	 },
	 {
	     intrinsic  =>   "double precision",
	     attributes => [ "allocatable" , "dimension(:)" ],
	     variables  => [ "metaProperties" ]
	 },
	 {
	     intrinsic  =>   "type",
	     type       =>   "rank1MetaProperty",
	     attributes => [ "allocatable" , "dimension(:)" ],
	     variables  => [ "rank1MetaProperties" ]
	 },
	 {
	     intrinsic  =>   "integer(kind_int8)",
	     attributes => [ "allocatable" , "dimension(:)" ],
	     variables  => [ "integerMetaProperties" ]
	 }
	);
    # Create the nodeComponent class.
    $build->{'types'}->{'nodeComponent'} = {
	name           => "nodeComponent"                           ,
	comment        => "A class for components in \\glspl{node}.",
	isPublic       => 1                                         ,
	boundFunctions => \@typeBoundFunctions,
	dataContent    => \@dataContent
    };
}

1;
