!--> 1 "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()"
!--> 1 "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()"
submodule (Merger_Tree_Branching) mergerTreeBranchingProbabilityParkinsonColeHelly_
!--> 1 "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()"
!--> 1 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!!{
!< Implements a merger tree branching probability class using the algorithm of \cite{parkinson_generating_2008}.
!!}

!--> 24 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 101 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 109 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

  ! Module-scope pointer to self used for root-finding.
!--> 111 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  class(mergerTreeBranchingProbabilityParkinsonColeHelly), pointer :: self_
!--> 112 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  !$omp threadprivate(self_)

  ! Branching probability integrand integration tolerance.
!--> 115 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  double precision, parameter :: toleranceintegralrelative=1.0d-3
!--> 116 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

  ! Limit on α for use in effective γ parameters.
!--> 119 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  double precision, parameter :: alphaminimum=5.0d-3
!--> 118 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

contains
!--> 1 "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()"
!--> 121 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 122 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyConstructorParameters
!--> 1 "Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses()"
!--> 1 "Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses()"
use :: Input_Parameters
use :: ISO_Varying_String
use :: MPI_Utilities             , only : mpiSelf
use :: Error
use :: Cosmological_Density_Field, only : cosmologicalMassVariance, cosmologicalMassVarianceClass, criticalOverdensity, criticalOverdensityClass
!--> 123 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Constructor for the ``parkinsonColeHelly'' merger tree branching probability class which reads parameters from a provided
!<     parameter list.
    !!}
!--> 127 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  class(cosmologicalMassVarianceClass), pointer :: cosmologicalmassvariance_
  class(criticalOverdensityClass), pointer :: criticaloverdensity_
  double precision :: gamma1, gamma2, g0, accuracyfirstorder, precisionhypergeometric
  logical :: hypergeometrictabulate, cdmassumptions
  type(varying_string), dimension(:), allocatable :: allowedParameterNames_
  logical, save :: warnObjectBuilder0__
 !$omp threadprivate(warnObjectBuilder0__)
  type(inputParameters), pointer :: parametersCurrent
  type(inputParameter), pointer :: parameterNode
  class(*), pointer :: genericObject
  logical, save :: warnObjectBuilder1__
 !$omp threadprivate(warnObjectBuilder1__)
  integer, save :: referenceCount_
 !$omp threadprivate(referenceCount_)
!--> 136 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    ! Check and read parameters.
!--> 139 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 139 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <inputParameter>
!<       <name>G0</name>
!<       <defaultValue>0.57d0</defaultValue>
!<       <description>The parameter $G_0$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.</description>
!<       <source>parameters</source>
!<     </inputParameter>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::InputParameter::Process_InputParameters()"
  ! Auto-generated input parameter
  call parameters%value('G0',G0,defaultValue=0.57d0)
  ! End auto-generated input parameter

!--> 144 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 144 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <inputParameter>
!<       <name>gamma1</name>
!<       <defaultValue>0.38d0</defaultValue>
!<       <description>The parameter $\gamma_1$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.</description>
!<       <source>parameters</source>
!<     </inputParameter>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::InputParameter::Process_InputParameters()"
  ! Auto-generated input parameter
  call parameters%value('gamma1',gamma1,defaultValue=0.38d0)
  ! End auto-generated input parameter

!--> 150 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 150 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <inputParameter>
!<       <name>gamma2</name>
!<       <defaultValue>-0.01d0</defaultValue>
!<       <description>The parameter $\gamma_2$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.</description>
!<       <source>parameters</source>
!<     </inputParameter>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::InputParameter::Process_InputParameters()"
  ! Auto-generated input parameter
  call parameters%value('gamma2',gamma2,defaultValue=-0.01d0)
  ! End auto-generated input parameter

!--> 156 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 156 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <inputParameter>
!<       <name>accuracyFirstOrder</name>
!<       <defaultValue>0.1d0</defaultValue>
!<       <description>Limits the step in $\delta_\mathrm{crit}$ when constructing merger trees using the \cite{parkinson_generating_2008}
!<          algorithm, so that it never exceeds {\normalfont \ttfamily accuracyFirstOrder}$\sqrt{2[\sigma^2(M_2/2)-\sigma^2(M_2)]}$.</description>
!<       <source>parameters</source>
!<     </inputParameter>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::InputParameter::Process_InputParameters()"
  ! Auto-generated input parameter
  call parameters%value('accuracyFirstOrder',accuracyFirstOrder,defaultValue=0.1d0)
  ! End auto-generated input parameter

!--> 163 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 163 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <inputParameter>
!<       <name>precisionHypergeometric</name>
!<       <defaultValue>1.0d-6</defaultValue>
!<       <description>The fractional precision required in evaluates of hypergeometric functions in the modified Press-Schechter tree branching calculations.</description>
!<       <source>parameters</source>
!<     </inputParameter>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::InputParameter::Process_InputParameters()"
  ! Auto-generated input parameter
  call parameters%value('precisionHypergeometric',precisionHypergeometric,defaultValue=1.0d-6)
  ! End auto-generated input parameter

!--> 169 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 169 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <inputParameter>
!<       <name>hypergeometricTabulate</name>
!<       <defaultValue>.true.</defaultValue>
!<       <description>Specifies whether hypergeometric factors should be precomputed and tabulated in modified Press-Schechter tree branching functions.</description>
!<       <source>parameters</source>
!<     </inputParameter>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::InputParameter::Process_InputParameters()"
  ! Auto-generated input parameter
  call parameters%value('hypergeometricTabulate',hypergeometricTabulate,defaultValue=.true.)
  ! End auto-generated input parameter

!--> 175 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 175 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <inputParameter>
!<       <name>cdmAssumptions</name>
!<       <defaultValue>.false.</defaultValue>
!<       <description>If true, assume that $\alpha(=-\mathrm{d}\log \sigma/\mathrm{d}\log M)&gt;0$ and $\mathrm{d}\alpha/\mathrm{d}M&gt;0$ (as is true in the case of \gls{cdm}) when constructing merger trees using the \cite{parkinson_generating_2008}.</description>
!<       <source>parameters</source>
!<     </inputParameter>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::InputParameter::Process_InputParameters()"
  ! Auto-generated input parameter
  call parameters%value('cdmAssumptions',cdmAssumptions,defaultValue=.false.)
  ! End auto-generated input parameter

!--> 181 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 181 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::ObjectBuilder::Process_ObjectBuilder()"
   ! Determine where to build+store or point to the required object....
   parametersCurrent => parameters
   do while (.not.parametersCurrent%isPresent('cosmologicalMassVariance').and.associated(parametersCurrent%parent))
      parametersCurrent => parametersCurrent%parent
   end do
   if (parametersCurrent%isPresent('cosmologicalMassVariance')) then
      ! Object should belong to the parameter node. Get the node and test whether the object has already been created in it.
      parameterNode => parametersCurrent%node('cosmologicalMassVariance')
      if (parameterNode%objectCreated()) then
         ! Object already exists - simply get a pointer to it. Increment the reference counter as this is a new reference to an existing object.
         genericObject => parameterNode%objectGet()
         select type (genericObject)
         class is (cosmologicalMassVarianceClass)
            cosmologicalMassVariance_ => genericObject
            call cosmologicalMassVariance_%referenceCountIncrement()
         class default
            call Error_Report('parameter-stored object is not of [cosmologicalMassVariance] class'//char(10)//' Occurred at:'//char(10)//'    directive:objectBuilder'//char(10)//'     function:parkinsonColeHellyConstructorParameters'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 181]')
         end select
      else
         ! Object does not yet exist - build it and store in the parameter node. Increment reference counter here as this is a newly constructed object.
         cosmologicalMassVariance_ => cosmologicalMassVariance(parametersCurrent)
            call cosmologicalMassVariance_%referenceCountIncrement()
         call parameterNode%objectSet(cosmologicalMassVariance_)
         call cosmologicalMassVariance_%autoHook()
      end if
   else
      ! Object is not explicitly defined. Cause a default object of the class to be added to the parameters. Increment the reference count here as this is a new object.
      cosmologicalMassVariance_ => cosmologicalMassVariance(parametersCurrent)
      call cosmologicalMassVariance_%referenceCountIncrement()
      call cosmologicalMassVariance_%autoHook()
      if (mpiSelf%isMaster() .and. .not.warnObjectBuilder0__) then
         call Warn('Using default class for parameter ''['//char(parametersCurrent%path())//'cosmologicalMassVariance]''')
         warnObjectBuilder0__=.true.
      end if
   end if
!--> 182 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 182 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::ObjectBuilder::Process_ObjectBuilder()"
   ! Determine where to build+store or point to the required object....
   parametersCurrent => parameters
   do while (.not.parametersCurrent%isPresent('criticalOverdensity').and.associated(parametersCurrent%parent))
      parametersCurrent => parametersCurrent%parent
   end do
   if (parametersCurrent%isPresent('criticalOverdensity')) then
      ! Object should belong to the parameter node. Get the node and test whether the object has already been created in it.
      parameterNode => parametersCurrent%node('criticalOverdensity')
      if (parameterNode%objectCreated()) then
         ! Object already exists - simply get a pointer to it. Increment the reference counter as this is a new reference to an existing object.
         genericObject => parameterNode%objectGet()
         select type (genericObject)
         class is (criticalOverdensityClass)
            criticalOverdensity_ => genericObject
            call criticalOverdensity_%referenceCountIncrement()
         class default
            call Error_Report('parameter-stored object is not of [criticalOverdensity] class'//char(10)//' Occurred at:'//char(10)//'    directive:objectBuilder'//char(10)//'     function:parkinsonColeHellyConstructorParameters'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 182]')
         end select
      else
         ! Object does not yet exist - build it and store in the parameter node. Increment reference counter here as this is a newly constructed object.
         criticalOverdensity_ => criticalOverdensity(parametersCurrent)
            call criticalOverdensity_%referenceCountIncrement()
         call parameterNode%objectSet(criticalOverdensity_)
         call criticalOverdensity_%autoHook()
      end if
   else
      ! Object is not explicitly defined. Cause a default object of the class to be added to the parameters. Increment the reference count here as this is a new object.
      criticalOverdensity_ => criticalOverdensity(parametersCurrent)
      call criticalOverdensity_%referenceCountIncrement()
      call criticalOverdensity_%autoHook()
      if (mpiSelf%isMaster() .and. .not.warnObjectBuilder1__) then
         call Warn('Using default class for parameter ''['//char(parametersCurrent%path())//'criticalOverdensity]''')
         warnObjectBuilder1__=.true.
      end if
   end if
!--> 183 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    self=mergerTreeBranchingProbabilityParkinsonColeHelly(G0,gamma1,gamma2,accuracyFirstOrder,precisionHypergeometric,hypergeometricTabulate,cdmAssumptions,cosmologicalMassVariance_,criticalOverdensity_)
!--> 187 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 187 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <inputParametersValidate source="parameters"/>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::InputParametersValidate::Process_InputParametersValidate()"
   call self%allowedParameters(allowedParameterNames_,'parameters',.false.)
   if (.not.mergerTreeBranchingProbabilityDsblVldtn) call parameters%checkParameters(allowedParameterNames_)
   if (allocated(allowedParameterNames_)) deallocate(allowedParameterNames_)
!--> 187 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 187 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::ObjectBuilder::Process_ObjectBuilder()"
if (associated(cosmologicalMassVariance_)) then
   ! Decrement the reference count, and decide if this object can be destroyed.
   referenceCount_=cosmologicalMassVariance_%referenceCountDecrement()
   if (referenceCount_ == 0) then
      ! Deallocate the pointer.
      deallocate(cosmologicalMassVariance_)
   else if (referenceCount_ < 0) then
      ! Negative counter - should not happen.
      call Error_Report('negative reference counter in object "cosmologicalMassVariance_"'//char(10)//' Occurred at:'//char(10)//'    directive:objectDestructor'//char(10)//'     function:parkinsonColeHellyConstructorParameters'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 187]')
   else
      ! Nullify the pointer.
      nullify(cosmologicalMassVariance_)
   end if
end if
!--> 188 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 188 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <objectDestructor name="criticalOverdensity_"     />
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::ObjectBuilder::Process_ObjectBuilder()"
if (associated(criticalOverdensity_)) then
   ! Decrement the reference count, and decide if this object can be destroyed.
   referenceCount_=criticalOverdensity_%referenceCountDecrement()
   if (referenceCount_ == 0) then
      ! Deallocate the pointer.
      deallocate(criticalOverdensity_)
   else if (referenceCount_ < 0) then
      ! Negative counter - should not happen.
      call Error_Report('negative reference counter in object "criticalOverdensity_"'//char(10)//' Occurred at:'//char(10)//'    directive:objectDestructor'//char(10)//'     function:parkinsonColeHellyConstructorParameters'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 188]')
   else
      ! Nullify the pointer.
      nullify(criticalOverdensity_)
   end if
end if
!--> 189 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    return
end procedure parkinsonColeHellyConstructorParameters
!--> 193 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 194 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyConstructorInternal
!--> 195 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Internal constructor for the ``parkinsonColeHelly'' merger tree branching probability class.
    !!}
!--> 198 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    use :: Error                , only : Error_Report
    use :: Numerical_Integration, only : GSL_Integ_Gauss15
!--> 200 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
!--> 209 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 209 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <constructorAssign variables="G0, gamma1, gamma2, accuracyFirstOrder, precisionHypergeometric, hypergeometricTabulate, cdmAssumptions, *cosmologicalMassVariance_, *criticalOverdensity_"/>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::Constructor::Process_Constructors()"
  ! Auto-generated constructor assignment
   self%G0=G0
   self%gamma1=gamma1
   self%gamma2=gamma2
   self%accuracyFirstOrder=accuracyFirstOrder
   self%precisionHypergeometric=precisionHypergeometric
   self%hypergeometricTabulate=hypergeometricTabulate
   self%cdmAssumptions=cdmAssumptions
   self%cosmologicalMassVariance_ => cosmologicalMassVariance_
   if (associated(self%cosmologicalMassVariance_))  call self%cosmologicalMassVariance_%referenceCountIncrement()
   self%criticalOverdensity_ => criticalOverdensity_
   if (associated(self%criticalOverdensity_))  call self%criticalOverdensity_%referenceCountIncrement()
  ! End auto-generated constructor assignment.

!--> 209 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    ! Validate inputs.
    if (gamma1 == 1.0d0) call Error_Report('γ₁=1 leads to divergent integrals'//char(10)//' Occurred at:'//char(10)//'     function:parkinsonColeHellyConstructorInternal'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 213]')
    ! Initialize.
    self%subresolutionHypergeometricInitialized=.false.
    self%upperBoundHypergeometricInitialized   =.false.
    self%massResolutionTabulated               =-1.0d0
    self%haloMassPrevious                      =-1.0d0
    self%deltaCriticalPrevious                 =-1.0d0
    self%massResolutionPrevious                =-1.0d0
    self%probabilityPrevious                   =-1.0d0
    self%integrator_                           =integrator(                                                                     &
            &                                                                parkinsonColeHellyProbabilityIntegrandLogarithmic, &
            &                                              toleranceRelative=toleranceIntegralRelative                        , &
            &                                              integrationRule  =GSL_Integ_Gauss15                                  &
            &                                             )
    return
end procedure parkinsonColeHellyConstructorInternal
!--> 229 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 230 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyDestructor
!--> 1 "Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses()"
!--> 1 "Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses()"
use :: Error
!--> 231 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  integer, save :: referenceCount_
 !$omp threadprivate(referenceCount_)
!--> 233 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    if (self%subresolutionHypergeometricInitialized) call self%subresolutionHypergeometric%destroy()
    if (self%upperBoundHypergeometricInitialized   ) call self%upperBoundHypergeometric   %destroy()
!--> 237 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 237 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <objectDestructor name="self%criticalOverdensity_"     />
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::ObjectBuilder::Process_ObjectBuilder()"
if (associated(self%criticalOverdensity_)) then
   ! Decrement the reference count, and decide if this object can be destroyed.
   referenceCount_=self%criticalOverdensity_%referenceCountDecrement()
   if (referenceCount_ == 0) then
      ! Deallocate the pointer.
      deallocate(self%criticalOverdensity_)
   else if (referenceCount_ < 0) then
      ! Negative counter - should not happen.
      call Error_Report('negative reference counter in object "self%criticalOverdensity_"'//char(10)//' Occurred at:'//char(10)//'    directive:objectDestructor'//char(10)//'   subroutine:parkinsonColeHellyDestructor'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 237]')
   else
      ! Nullify the pointer.
      nullify(self%criticalOverdensity_)
   end if
end if
!--> 237 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
!--> 237 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !![
!<     <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
!--> 1 "Galacticus::Build::SourceTree::Process::ObjectBuilder::Process_ObjectBuilder()"
if (associated(self%cosmologicalMassVariance_)) then
   ! Decrement the reference count, and decide if this object can be destroyed.
   referenceCount_=self%cosmologicalMassVariance_%referenceCountDecrement()
   if (referenceCount_ == 0) then
      ! Deallocate the pointer.
      deallocate(self%cosmologicalMassVariance_)
   else if (referenceCount_ < 0) then
      ! Negative counter - should not happen.
      call Error_Report('negative reference counter in object "self%cosmologicalMassVariance_"'//char(10)//' Occurred at:'//char(10)//'    directive:objectDestructor'//char(10)//'   subroutine:parkinsonColeHellyDestructor'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 237]')
   else
      ! Nullify the pointer.
      nullify(self%cosmologicalMassVariance_)
   end if
end if
!--> 238 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    return
end procedure parkinsonColeHellyDestructor
!--> 242 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 243 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyMassBranch
!--> 244 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     A merger tree branch split mass function.
    !!}
!--> 247 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  double precision :: b, mu, beta, halfmassalpha, halfmasssigma, eta, halfmassv, massfractionresolution, halfpowereta, x, massfraction, resolutionsigma, massfractionresolutionpowereta
!--> 261 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !$GLC attributes unused :: node

    ! Simply branch to the relevant function.
    if (self%cdmAssumptions) then
       parkinsonColeHellyMassBranch=massBranchCDMAssumptions()
    else
       parkinsonColeHellyMassBranch=massBranchGeneric       ()
    end if
    return

!--> 271 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  contains
!--> 272 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 273 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    double precision function massBranchCDMAssumptions()
!--> 274 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      !!{
!<       A merger tree branch split mass function which assumes a \gls{cdm}-like power spectrum. With these assumptions, it can
!<       employ the mass sampling algorithm of \cite{parkinson_generating_2008}. One difference with respect to the algorithm of
!<       \cite{parkinson_generating_2008} is that here the normalization of their function $S(q)$ (eqn. A2) is irrelevant, since a
!<       branch split has already been decided to have occurred---all that remains necessary is to determine its mass. Variable and
!<       function names follow \cite{parkinson_generating_2008}.
      !!}
!--> 281 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      implicit none
      logical :: reject
!--> 283 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

      ! Get parent and half-mass σ and α.
      self%timeParent        =self%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=     deltaCritical,mass=haloMass,node=node)
      self%sigmaParentSquared=self%cosmologicalMassVariance_%rootVariance  (time               =self%timeParent   ,mass=haloMass          )**2
      call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*haloMass,self%timeParent,halfMassSigma,halfMassAlpha)
      ! Compute parameters β, μ, and B.
      massFractionResolution=+massResolution                               &
           &                 /haloMass
      halfMassV             =+self%V(0.5d0,haloMass)
      beta                  =+log(                                         &
           &                      +self%V(massFractionResolution,haloMass) &
           &                      /halfMassV                               &
           &                     )                                         &
           &                 /log(                                         &
           &                      +massFractionResolution                  &
           &                      /0.5d0                                   &
           &                     )
      B                     =+halfMassV                                    &
           &                 *2.0d0    **beta
      if (self%gamma1 >= 0.0d0) then
         mu                 =-halfMassAlpha
      else
         resolutionSigma    =+self%cosmologicalMassVariance_%rootVariance(massResolution,self%timeParent)
         mu                 =-log(                        &
              &                   +resolutionSigma        &
              &                   /halfMassSigma          &
              &                  )                        &
              &              /log(                        &
              &                   +massFractionResolution &
              &                   /0.5d0                  &
              &                  )
      end if
      eta                           =+beta                        &
           &                         -1.0d0                       &
           &                         -mu                          &
           &                         *self%gamma1
      massFractionResolutionPowerEta=+massFractionResolution**eta
      halfPowerEta                  =+0.5d0                 **eta
      ! Sample from S(q), using rejection sampling on R(q) to decide whether to keep/reject the
      ! proposed q.
      reject=.true.
      do while (reject)
         ! Draw a random q from S(q).
         x           =randomNumberGenerator_%uniformSample()
         massFraction=(                                                 &
              &        +              massFractionResolutionPowerEta    &
              &        +(halfPowerEta-massFractionResolutionPowerEta)*x &
              &       )**(1.0d0/eta)
         x           =randomNumberGenerator_%uniformSample()
         reject=x > R(massFraction)
      end do
      massBranchCDMAssumptions=massFraction*haloMass
      return
    end function massBranchCDMAssumptions
!--> 337 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 338 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    double precision function R(massFraction)
!--> 339 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      !!{
!<       The function $R(q)$ from \cite[][eqn. A3]{parkinson_generating_2008}.
      !!}
!--> 342 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      implicit none
      double precision, intent(in   ) :: massFraction
      double precision                :: massFractionSigma, massFractionAlpha
!--> 345 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

      call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massFraction*haloMass,self%timeParent,massFractionSigma,massFractionAlpha)
      R      =+(                             &
           &    +massFractionAlpha           &
           &    /halfMassAlpha               &
           &   )                             &
           &  *self%V(massFraction,haloMass) &
           &  /B                             &
           &  /massFraction**beta            &
           &  *(                             &
           &    +(                           &
           &      +2.0d0                     &
           &      *massFraction              &
           &     )**mu                       &
           &    *massFractionSigma           &
           &    /halfMassSigma               &
           &   )**self%gamma1
      return
    end function R
!--> 364 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 365 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    double precision function massBranchGeneric()
!--> 366 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      !!{
!<       Determine the mass of one of the halos to which the given halo branches, given the branching probability, {\normalfont
!<       \ttfamily probability}. Typically, {\normalfont \ttfamily probabilityFraction} is found by multiplying {\normalfont \ttfamily probability}
!<       by a random variable drawn in the interval 0--1 if a halo branches. This routine then finds the progenitor mass
!<       corresponding to this value.
      !!}
!--> 372 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      use :: Root_Finder, only : GSL_Root_fSolver_Brent, rootFinder
!--> 373 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      implicit none
      double precision            , parameter :: toleranceAbsolute=0.0d0  , toleranceRelative=1.0d-9
      type            (rootFinder), save      :: finder
      logical                     , save      :: finderConstructed=.false.
!--> 377 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      !$omp threadprivate(finder,finderConstructed)
!--> 378 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
      double precision                        :: logMassMinimum           , logMassMaximum
!--> 379 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

      ! Initialize global variables.
      call self%computeCommonFactors(deltaCritical,time,haloMass,node)
      self_                            => self
      self  %probabilityMinimumMass    =            massResolution
      self  %probabilityMinimumMassLog =  log(      massResolution)
      self  %probabilityMaximumMassLog =  log(0.5d0*haloMass      )
      self  %probabilitySeek           =  probabilityFraction
      ! Check the sign of the root function at half the halo mass.
      if (parkinsonColeHellyMassBranchRoot(self%probabilityMaximumMassLog) >= 0.0d0) then
         ! The root function is zero, or very close to it (which can happen due to rounding errors
         ! occasionally). Therefore we have an almost perfect binary split.
         massBranchGeneric=0.5d0*haloMass
      else
         ! Initialize our root finder.
         if (.not.finderConstructed) then
            finder           =rootFinder(                                                    &
                 &                       rootFunction     =parkinsonColeHellyMassBranchRoot, &
                 &                       toleranceAbsolute=toleranceAbsolute               , &
                 &                       toleranceRelative=toleranceRelative,                &
                 &                       solverType       =GSL_Root_fSolver_Brent            &
                 &                      )
            finderConstructed=.true.
         end if
         ! Split is not binary - seek the actual mass of the smaller progenitor.
         logMassMinimum                  =log(      massResolution)
         logMassMaximum                  =log(0.5d0*haloMass      )
         self_%probabilityGradientMinimum=parkinsonColeHellyMassBranchRootDerivative(logMassMinimum)
         self_%probabilityGradientMaximum=parkinsonColeHellyMassBranchRootDerivative(logMassMaximum)
         self_%probabilityMaximum        =parkinsonColeHellyMassBranchRoot          (logMassMaximum)
         massBranchGeneric=exp(finder%find(rootRange=[logMassMinimum,logMassMaximum]))
      end if
      return
    end function massBranchGeneric
!--> 413 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

end procedure parkinsonColeHellyMassBranch
!--> 415 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 416 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyV
!--> 417 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     The function $V(q)$ from \cite[][eqn. A4]{parkinson_generating_2008}.
    !!}
!--> 420 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  double precision :: childsigmasquared
!--> 424 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    childSigmaSquared  =+self%cosmologicalMassVariance_%rootVariance(massFraction*haloMass,self%timeParent)**2
    parkinsonColeHellyV=+       childSigmaSquared  &
         &              /(                         &
         &                +     childSigmaSquared  &
         &                -self%sigmaParentSquared &
         &               )**1.5d0
    return
end procedure parkinsonColeHellyV
!--> 433 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 434 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  double precision function parkinsonColeHellyMassBranchRoot(logMassMaximum)
!--> 435 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Used to find the mass of a merger tree branching event.
    !!}
!--> 438 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    implicit none
    double precision, intent(in   ) :: logMassMaximum
    double precision                :: integral      , massMaximum
!--> 441 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    if      (logMassMaximum < self_%probabilityMinimumMassLog) then
       parkinsonColeHellyMassBranchRoot=self_%probabilitySeek   +self_%probabilityGradientMinimum*(logMassMaximum-self_%probabilityMinimumMassLog)
    else if (logMassMaximum > self_%probabilityMaximumMassLog) then
       parkinsonColeHellyMassBranchRoot=self_%probabilityMaximum+self_%probabilityGradientMaximum*(logMassMaximum-self_%probabilityMaximumMassLog)
    else
       massMaximum=+exp(logMassMaximum)
       integral   =+self_%branchingProbabilityPreFactor                                                           &
            &      *self_%integrator_                  %integrate(self_%probabilityMinimumMassLog,logMassMaximum)
       parkinsonColeHellyMassBranchRoot=self_%probabilitySeek-integral
    end if
    return
  end function parkinsonColeHellyMassBranchRoot
!--> 454 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 455 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  double precision function parkinsonColeHellyMassBranchRootDerivative(logMassMaximum)
!--> 456 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Used to find the mass of a merger tree branching event.
    !!}
!--> 459 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    implicit none
    double precision, intent(in   ) :: logMassMaximum
    double precision                :: integral
!--> 462 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    integral=+self_%branchingProbabilityPreFactor                                                    &
         &   *parkinsonColeHellyProbabilityIntegrandLogarithmic(                                     &
         &                                                      max(                                 &
         &                                                          logMassMaximum                 , &
         &                                                          self_%probabilityMinimumMassLog  &
         &                                                         )                                 &
         &                                                     )
    parkinsonColeHellyMassBranchRootDerivative=-integral
    return
  end function parkinsonColeHellyMassBranchRootDerivative
!--> 473 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 474 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyStepMaximum
!--> 475 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Return the maximum allowed step in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time {\normalfont \ttfamily
!<     deltaCritical} should be allowed to take.
    !!}
!--> 479 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  double precision, parameter :: largestep=1.0d10
  double precision :: parenthalfmasssigma, parentsigma
!--> 485 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !$GLC attributes unused :: deltaCritical, time

    ! Get σ and δ_critical for the parent halo.
    if (haloMass > 2.0d0*massResolution) then
       parentSigma                  =+self%cosmologicalMassVariance_%rootVariance(      haloMass,self%timeParent)
       parentHalfMassSigma          =+self%cosmologicalMassVariance_%rootVariance(0.5d0*haloMass,self%timeParent)
       parkinsonColeHellyStepMaximum=+self%accuracyFirstOrder        &
            &                        *sqrt(                          &
            &                              +2.0d0                    &
            &                              *(                        &
            &                                +parentHalfMassSigma**2 &
            &                                -parentSigma        **2 &
            &                               )                        &
            &                             )
    else
       parkinsonColeHellyStepMaximum=largeStep
    end if
    return
end procedure parkinsonColeHellyStepMaximum
!--> 504 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 505 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyRate
!--> 506 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Return the rate per unit mass and per unit change in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time
!<     {\normalfont \ttfamily deltaCritical} will undergo a branching to progenitors with mass {\normalfont \ttfamily massBranch}.
    !!}
!--> 510 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  double precision :: massbranch_
!--> 516 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    
    ! Always use the rate from the lower half of the mass range.
    if (massBranch > 0.5d0*mass) then
       massBranch_=+mass-massBranch
    else
       massBranch_=     +massBranch
    end if
    call self%computeCommonFactors(deltaCritical,time,mass,node)
    self_                  =>  self
    parkinsonColeHellyRate =  +self%branchingProbabilityPreFactor                                  &
         &                    *parkinsonColeHellyProbabilityIntegrandLogarithmic(log(massBranch_)) &
         &                    /massBranch_
    return
end procedure parkinsonColeHellyRate
!--> 530 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 531 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyProbability
!--> 532 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Return the probability per unit change in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time
!<     {\normalfont \ttfamily deltaCritical} will undergo a branching to progenitors with mass greater than {\normalfont \ttfamily massResolution}.
    !!}
!--> 536 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  double precision :: massmaximum, massminimum
!--> 542 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !$GLC attributes unused :: node

    ! Recompute branching probability if necessary.
    if     (                                               &
         &   haloMass       /= self%haloMassPrevious       &
         &  .or.                                           &
         &   deltaCritical  /= self%deltaCriticalPrevious  &
         &  .or.                                           &
         &   massResolution /= self%massResolutionPrevious &
         & ) then
       self_                        => self
       self %haloMassPrevious       =  haloMass
       self %deltaCriticalPrevious  =  deltaCritical
       self %massResolutionPrevious =  massResolution
       ! Get σ and δ_critical for the parent halo.
       if (haloMass > 2.0d0*massResolution) then
          call self%computeCommonFactors(deltaCritical,time,haloMass,node)
          massMinimum             =+           massResolution
          massMaximum             =+0.5d0*self%massHaloParent
          self%probabilityPrevious=+self%branchingProbabilityPreFactor                             &
               &                   *self%integrator_                  %integrate(                  &
               &                                                                 log(massMinimum), &
               &                                                                 log(massMaximum)  &
               &                                                                )
       else
          self%probabilityPrevious=0.0d0
       end if
    end if
    parkinsonColeHellyProbability=self%probabilityPrevious
    return
end procedure parkinsonColeHellyProbability
!--> 573 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 574 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  double precision function parkinsonColeHellyProbabilityIntegrandLogarithmic(logChildHaloMass)
!--> 575 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Integrand for the branching probability.
    !!}
!--> 578 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    implicit none
    double precision, intent(in   ) :: logChildHaloMass
    double precision                :: childAlpha      , childSigma, &
         &                             childHaloMass
!--> 582 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    childHaloMass=exp(logChildHaloMass)
    call self_%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(childHaloMass,self_%timeParent,childSigma,childAlpha)
    parkinsonColeHellyProbabilityIntegrandLogarithmic=+parkinsonColeHellyProgenitorMassFunction(childHaloMass,childSigma,childAlpha)&
         &                                            *                                         childHaloMass
    return
  end function parkinsonColeHellyProbabilityIntegrandLogarithmic
!--> 589 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 590 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  double precision function parkinsonColeHellyProgenitorMassFunction(childHaloMass,childSigma,childAlpha)
!--> 591 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Progenitor mass function from Press-Schechter. The constant factor of the parent halo mass is not included here---instead
!<     it is included in a multiplicative prefactor by which integrals over this function are multiplied.
    !!}
!--> 595 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    implicit none
    double precision, intent(in   ) :: childAlpha, childHaloMass, childSigma
!--> 597 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    parkinsonColeHellyProgenitorMassFunction=+      parkinsonColeHellyMergingRate(childSigma   ,childAlpha)    &
         &                                   *self_%modifier                     (childSigma              )    &
         &                                   /                                    childHaloMass            **2
    return
  end function parkinsonColeHellyProgenitorMassFunction
!--> 603 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 604 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  double precision function parkinsonColeHellyMergingRate(childSigma,childAlpha)
!--> 605 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Merging rate from Press-Schechter. The constant factor of $\sqrt{2/\pi}$ not included here---instead it is included in a
!<     multiplicative prefactor by which integrals over this function are multiplied.
    !!}
!--> 609 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    implicit none
    double precision, intent(in   ) :: childAlpha       , childSigma
    double precision                :: childSigmaSquared
!--> 612 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    childSigmaSquared=childSigma**2
    if (childSigmaSquared > self_%sigmaParentSquared .and. childAlpha < 0.0d0) then
       parkinsonColeHellyMergingRate=(childSigmaSquared/((childSigmaSquared-self_%sigmaParentSquared)**1.5d0))*abs(childAlpha)
    else
       parkinsonColeHellyMergingRate=0.0d0
    end if
    return
  end function parkinsonColeHellyMergingRate
!--> 621 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 622 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyModifier
!--> 623 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Empirical modification of the progenitor mass function from \cite{parkinson_generating_2008}. The constant factors of
!<     $G_0 (\delta_\mathrm{p}/\sigma_\mathrm{p})^{\gamma_2}$ and $1/\sigma_\mathrm{p}^{\gamma_1}$ are not included
!<     here---instead they are included in a multiplicative prefactor by which integrals over this function are multiplied.
    !!}
!--> 628 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
!--> 631 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    parkinsonColeHellyModifier=childSigma**self%gamma1
    return
end procedure parkinsonColeHellyModifier
!--> 635 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 636 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyHypergeometricA
!--> 637 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Compute the $a$ parameter of the hypergeometric function.
    !!}
!--> 640 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
!--> 644 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    
    a=[1.5d0,0.5d0-0.5d0*gamma]
    return
end procedure parkinsonColeHellyHypergeometricA
!--> 648 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  
!--> 649 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyProbabilityBound
!--> 650 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Return a bound on the probability per unit change in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily
!<     haloMass} at time {\normalfont \ttfamily deltaCritical} will undergo a branching to progenitors with mass greater than
!<     {\normalfont \ttfamily massResolution}.
    !!}
!--> 655 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    use            :: Display                 , only : displayMessage    , verbosityLevelWarn, displayMagenta, displayReset
    use            :: Error                   , only : Error_Report
    use            :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use, intrinsic :: ISO_C_Binding           , only : c_int
    use            :: Interface_GSL           , only : GSL_Success
    use            :: Numerical_Constants_Math, only : Pi
!--> 661 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  double precision, parameter :: sqrttwooverpi=sqrt(2.0d0/pi)
  double precision :: probabilityintegrandlower, probabilityintegrandupper, halfparentsigma, halfparentalpha, gammaeffective
  double precision :: hypergeometricfactorlower, hypergeometricfactorupper, resolutionsigmaoverparentsigma
  integer(c_int) :: statuslower, statusupper
  logical :: usingcdmassumptions
  integer :: ibound
!--> 676 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !$GLC attributes unused :: node

    ! Get σ and δ_critical for the parent halo.
    if (haloMass > 2.0d0*massResolution) then
       call self%computeCommonFactors(deltaCritical,time,haloMass,node)
       call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massResolution,self%timeParent,self%resolutionSigma,self%resolutionAlpha)
       if (massResolution /= self%massResolutionTabulated .or. self%cosmologicalMassVariance_%growthIsMassDependent()) then
          ! Resolution changed - recompute σ and α at resolution limit. Also reset the hypergeometric factor tables since
          ! these depend on resolution.
          self%upperBoundHypergeometricInitialized=.false.
       end if
       resolutionSigmaOverParentSigma=self%resolutionSigma/self%sigmaParent
       ! Estimate probability.
       if (resolutionSigmaOverParentSigma > 1.0d0) then
          ! Compute relevant σ and α.
          call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*self%massHaloParent,self%timeParent,halfParentSigma,halfParentAlpha)
          ! Iterative over available bounds.
          parkinsonColeHellyProbabilityBound=0.0d0
          do iBound=1,2
             ! Determine if CDM assumptions can be used. Do this only is these have been explicitly allowed, if this is our first
             ! pass through the bounds evaluation, and if both αs are sufficiently large. (This last condition is required
             ! since we raise quantities to the power of 1/α which can cause problems for very small α.)
             usingCDMAssumptions= self%cdmAssumptions                       &
                  &              .and.                                      &
                  &               iBound                    == 1            &
                  &              .and.                                      &
                  &               abs(self%resolutionAlpha) >  alphaMinimum &
                  &              .and.                                      &
                  &               abs(     halfParentAlpha) >  alphaMinimum
             ! Compute the effective value of γ.
             gammaEffective=self%gamma1
             if (usingCDMAssumptions) then
                select case (bound%ID)
                case (mergerTreeBranchingBoundLower%ID)
                   gammaEffective=gammaEffective-1.0d0/self%resolutionAlpha
                case (mergerTreeBranchingBoundUpper%ID)
                   gammaEffective=gammaEffective-1.0d0/     halfParentAlpha
                end select
             end if
             ! Compute probability factors. The logic here becomes complicated, as we use various optimizations and tabulations to
             ! speed up calculation.
             !
             ! Tabulations will only be used if self%tabulateHypergeometric is true.
             !
             ! Set status to success by default.
             statusLower=GSL_Success
             statusUpper=GSL_Success
             ! First, check if CDM assumptions are not being used and we're allowed to tabulate hypergeometric factors,
             if (.not.usingCDMAssumptions.and.self%hypergeometricTabulate) then
                ! CDM assumptions are not being used. In this case we can use the same table of hypergeometric factors as the
                ! subresolution merger fraction.
                call parkinsonColeHellySubresolutionHypergeometricTabulate(self,resolutionSigmaOverParentSigma)
                call parkinsonColeHellySubresolutionHypergeometricTabulate(self,halfParentSigma   /self%sigmaParent)
                probabilityIntegrandLower=+self%factorG0Gamma2*self%subresolutionHypergeometric%interpolate(+resolutionSigmaOverParentSigma  -1.0d0)/self%sigmaParent
                probabilityIntegrandUpper=+self%factorG0Gamma2*self%subresolutionHypergeometric%interpolate(+halfParentSigma/self%sigmaParent-1.0d0)/self%sigmaParent
             else
                ! Next, check if CDM assumptions are being used, we're allowed to tabulate hypergeometric factors, and the bound
                ! requested is the upper bound.
                if     ( usingCDMAssumptions                    &
                     &  .and.                                   &
                     &   self%hypergeometricTabulate            &
                     &  .and.                                   &
                     &   bound == mergerTreeBranchingBoundUpper &
                     & ) then
                   ! Use a tabulation of the hypergeometric functions for the upper bound, made using CDM assumptions. Since the
                   ! tables already include the difference between the upper and lower integrand, we simply set the lower
                   ! integrand to zero here.
                   call parkinsonColeHellyUpperBoundHypergeometricTabulate(self,self%massHaloParent,massResolution)
                   probabilityIntegrandUpper=self%factorG0Gamma2*self%upperBoundHypergeometric%interpolate(self%massHaloParent)/self%resolutionSigma
                   probabilityIntegrandLower=0.0d0
                else
                   ! Use a direct calculation of the hypergeometric factors in this case.
                   hyperGeometricFactorLower=Hypergeometric_2F1(                                                           &
                        &                                                         self%hypergeometricA(gammaEffective)   , &
                        &                                                         [      1.5d0-0.5d0*gammaEffective]     , &
                        &                                                         1.0d0/resolutionSigmaOverParentSigma**2, &
                        &                                       toleranceRelative=self%precisionHypergeometric           , &
                        &                                       status           =statusLower                              &
                        &                                      )
                   if (statusLower /= GSL_Success) then
                      if (usingCDMAssumptions) then
                         if (.not.self%hypergeometricFailureWarned) then
                            self%hypergeometricFailureWarned=.true.
                            call displayMessage(                                                                                                                      &
                                 &              displayMagenta()//'WARNING:'//displayReset()//' hypergeometric function evaluation failed when computing'//char(10)// &
                                 &              'merger tree branching probability bounds - will revert to more'                                         //char(10)// &
                                 &              'robust (but less stringent) bound in this and future cases'                                                       ,  &
                                 &              verbosityLevelWarn                                                                                                    &
                                 &             )
                         end if
                         cycle
                      else
                         parkinsonColeHellyProbabilityBound=0.0d0
                         call Error_Report('hypergeometric function evaluation failed'//char(10)//' Occurred at:'//char(10)//'     function:parkinsonColeHellyProbabilityBound'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 769]')
                      end if
                   end if
                   probabilityIntegrandLower=+sqrtTwoOverPi                                            &
                        &                    *(self%factorG0Gamma2/self%sigmaParent)                   &
                        &                    *(resolutionSigmaOverParentSigma**(gammaEffective-1.0d0)) &
                        &                    /(1.0d0-gammaEffective)                                   &
                        &                    *hyperGeometricFactorLower
                   ! Check if we can use a table to compute the upper factor.
                   hyperGeometricFactorUpper=Hypergeometric_2F1(                                                          &
                        &                                                         self%hypergeometricA(gammaEffective)  , &
                        &                                                         [      1.5d0-0.5d0*gammaEffective]    , &
                        &                                                         self%sigmaParent**2/halfParentSigma**2, &
                        &                                       toleranceRelative=self%precisionHypergeometric          , &
                        &                                       status           =statusUpper                             &
                        &                                      )
                   if (statusUpper /= GSL_Success) then
                      if (usingCDMAssumptions) then
                         if (.not.self%hypergeometricFailureWarned) then
                            self%hypergeometricFailureWarned=.true.
                            call displayMessage(                                                                                                                      &
                                 &              displayMagenta()//'WARNING:'//displayReset()//' hypergeometric function evaluation failed when computing'//char(10)// &
                                 &              'merger tree branching probability bounds - will revert to more'                                         //char(10)// &
                                 &              'robust (but less stringent) bound in this and future cases'                                                       ,  &
                                 &              verbosityLevelWarn                                                                                                    &
                                 &             )
                         end if
                         cycle
                      else
                         parkinsonColeHellyProbabilityBound=0.0d0
                         call Error_Report('hypergeometric function evaluation failed'//char(10)//' Occurred at:'//char(10)//'     function:parkinsonColeHellyProbabilityBound'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 799]')
                      end if
                   end if
                   probabilityIntegrandUpper=+sqrtTwoOverPi                                                &
                        &                    *(self%factorG0Gamma2/self%sigmaParent)                       &
                        &                    *((halfParentSigma/self%sigmaParent)**(gammaEffective-1.0d0)) &
                        &                    /(1.0d0-gammaEffective)                                       &
                        &                    *hyperGeometricFactorUpper
                end if
             end if
             ! Compute the bound.
             select case (bound%ID)
             case (mergerTreeBranchingBoundLower%ID)
                if (usingCDMAssumptions) then
                   parkinsonColeHellyProbabilityBound=+(                               &
                        &                               +probabilityIntegrandUpper     &
                        &                               -probabilityIntegrandLower     &
                        &                              )                               &
                        &                             *self%massHaloParent             &
                        &                             /massResolution                  &
                        &                             *(                               &
                        &                               +self%resolutionSigma          &
                        &                               /self%sigmaParent              &
                        &                              )**(1.0d0/self%resolutionAlpha)
                else
                   parkinsonColeHellyProbabilityBound=+(                               &
                        &                               +probabilityIntegrandUpper     &
                        &                               -probabilityIntegrandLower     &
                        &                              )                               &
                        &                             *       self%massHaloParent      &
                        &                             /(0.5d0*self%massHaloParent)
                end if
             case (mergerTreeBranchingBoundUpper%ID)
                if (usingCDMAssumptions) then
                   parkinsonColeHellyProbabilityBound=+(                               &
                        &                               +probabilityIntegrandUpper     &
                        &                               -probabilityIntegrandLower     &
                        &                              )                               &
                        &                             *self%massHaloParent             &
                        &                             /massResolution                  &
                        &                             *(                               &
                        &                               +self%resolutionSigma          &
                        &                               /self%sigmaParent              &
                        &                              )**(1.0d0/halfParentAlpha)
                else
                   parkinsonColeHellyProbabilityBound=+(                               &
                        &                               +probabilityIntegrandUpper     &
                        &                               -probabilityIntegrandLower     &
                        &                              )                               &
                        &                             *self%massHaloParent             &
                        &                             /massResolution
                end if
             case default
                parkinsonColeHellyProbabilityBound=-1.0d0
                call Error_Report('unknown bound type'//char(10)//' Occurred at:'//char(10)//'     function:parkinsonColeHellyProbabilityBound'//char(10)//'         file:merger_trees.branching_probability.Parkinson_Cole_Helly.F90'//'   [line 853]')
             end select
             if (statusUpper == GSL_Success .and. statusLower == GSL_Success) exit
          end do
       else
          parkinsonColeHellyProbabilityBound=-1.0d0
       end if
    else
       parkinsonColeHellyProbabilityBound=0.0d0
    end if
    return
end procedure parkinsonColeHellyProbabilityBound
!--> 865 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 866 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyFractionSubresolution
!--> 867 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Return the fraction of mass accreted in subresolution halos, i.e. those below {\normalfont \ttfamily massResolution}, per unit change in
!<     $\delta_\mathrm{crit}$ for a halo of mass {\normalfont \ttfamily haloMass} at time {\normalfont \ttfamily deltaCritical}. The integral is computed analytically in
!<     terms of the $_2F_1$ hypergeometric function.
    !!}
!--> 872 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Constants_Math, only : Pi
!--> 874 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  double precision, parameter :: sqrttwooverpi=sqrt(2.0d0/pi)
  double precision :: hypergeometricfactor, resolutionsigmaoverparentsigma, resolutionsigma
!--> 882 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !$GLC attributes unused :: node

    ! Get σ and δ_critical for the parent halo.
    call self%computeCommonFactors(deltaCritical,time,haloMass,node)
    resolutionSigma               =self%cosmologicalMassVariance_%rootVariance(massResolution,self%timeParent)
    resolutionSigmaOverParentSigma=resolutionSigma/self%sigmaParent
    if (resolutionSigmaOverParentSigma > 1.0d0) then
       if (self%hypergeometricTabulate) then
          ! Use tabulation of hypergeometric factors.
          call parkinsonColeHellySubresolutionHypergeometricTabulate(self,resolutionSigmaOverParentSigma)
          parkinsonColeHellyFractionSubresolution=+self%factorG0Gamma2                                                          &
               &                                  *self%subresolutionHypergeometric%interpolate(                                &
               &                                                                                +resolutionSigmaOverParentSigma &
               &                                                                                -1.0d0                          &
               &                                                                               )                                &
               &                                  /self%sigmaParent
       else
          ! Compute hypergeometric factors directly.
          hyperGeometricFactor=Hypergeometric_2F1(                                                           &
               &                                                    self%hypergeometricA(self%gamma1)      , &
               &                                                    [      1.5d0-0.5d0*self%gamma1]        , &
               &                                                    1.0d0/resolutionSigmaOverParentSigma**2, &
               &                                  toleranceRelative=self%precisionHypergeometric             &
               &                                 )
          parkinsonColeHellyFractionSubresolution=+sqrtTwoOverPi                                             &
               &                                  *self%factorG0Gamma2                                       &
               &                                  /self%sigmaParent                                          &
               &                                  *resolutionSigmaOverParentSigma**(+self%gamma1-1.0d0)      &
               &                                  /                                (-self%gamma1+1.0d0)      &
               &                                  *hyperGeometricFactor
       end if
    else
       parkinsonColeHellyFractionSubresolution=-1.0d0
    end if
    return
end procedure parkinsonColeHellyFractionSubresolution
!--> 918 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 919 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
module procedure parkinsonColeHellyComputeCommonFactors
!--> 920 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Precomputes some useful factors that are used in the modified Press-Schechter branching integrals.
    !!}
!--> 923 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    use :: Numerical_Constants_Math, only : Pi
!--> 924 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
implicit none
  double precision, parameter :: sqrttwooverpi=sqrt(2.0d0/pi)
!--> 930 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    self%deltaParent                  =                                                                       deltaParent
    self%massHaloParent               =                                                                                             massHaloParent
    self%timeParent                   =                                                                       time
    self%sigmaParent                  =self%cosmologicalMassVariance_%rootVariance  (time               =self%timeParent ,mass=self%massHaloParent          )
    self%sigmaParentSquared           =self%sigmaParent**2
    self%factorG0Gamma2               =self%G0*((max(self%criticalOverdensity_%value(time=self%timeParent,mass=self%massHaloParent,node=node),0.0d0)/self%sigmaParent)**self%gamma2)
    self%branchingProbabilityPreFactor=sqrtTwoOverPi*self%massHaloParent*self%factorG0Gamma2/self%sigmaParent**self%gamma1
    return
end procedure parkinsonColeHellyComputeCommonFactors
!--> 940 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 941 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  subroutine parkinsonColeHellySubresolutionHypergeometricTabulate(self,x,xMinimumIn,xMaximumIn)
!--> 942 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Tabulate the hypergeometric term appearing in the subresolution merger fraction expression.
    !!}
!--> 945 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Constants_Math, only : Pi
    use :: Table_Labels            , only : extrapolationTypeAbort
!--> 948 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout)           :: self
    double precision                                                  , intent(in   )           :: x
    double precision                                                  , intent(in   ), optional :: xMinimumIn                    , xMaximumIn
    integer                                                           , parameter               :: xCountPerDecade=10
    double precision                                                  , parameter               :: sqrtTwoOverPi  =sqrt(2.0d0/Pi)
    double precision                                                                            :: xMinimum                      , xMaximum
    integer                                                                                     :: xCount                        , i
    logical                                                                                     :: tabulate
!--> 957 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    tabulate=.false.
    if (.not.self%subresolutionHypergeometricInitialized) then
       tabulate=.true.
       if (present(xMinimumIn)) then
          xMinimum=xMinimumIn
       else
          xMinimum=min( 1.0d-9 ,     (x-1.0d0))
       end if
       if (present(xMaximumIn)) then
          xMaximum=xMaximumIn
       else
          xMaximum=max(12.5d+0,2.0d0*(x-1.0d0))
       end if
    else
       if     (                                                    &
            &   (x-1.0d0) < self%subresolutionHypergeometric%x(+1) &
            &  .or.                                                &
            &   (x-1.0d0) > self%subresolutionHypergeometric%x(-1) &
            & ) then
          tabulate=.true.
          xMinimum=min(self%subresolutionHypergeometric%x(+1),      (x-1.0d0))
          xMaximum=max(self%subresolutionHypergeometric%x(-1),2.0d0*(x-1.0d0))
       end if
    end if
    if (tabulate) then
       xCount=int(log10(xMaximum/xMinimum)*dble(xCountPerDecade))+1
       if (.not.self%subresolutionHypergeometricInitialized) call self%subresolutionHypergeometric%destroy()
       call self%subresolutionHypergeometric%create(xMinimum,xMaximum,xCount,1,extrapolationType=spread(extrapolationTypeAbort,1,2))
       do i=1,xCount
          call self%subresolutionHypergeometric%populate(                                                                                         &
               &                                    +sqrtTwoOverPi                                                                                &
               &                                    *(self%subresolutionHypergeometric%x(i)+1.0d0)**(+self%gamma1-1.0d0)                          &
               &                                    /                                               (-self%gamma1+1.0d0)                          &
               &                                    *Hypergeometric_2F1(                                                                          &
               &                                                                          self%hypergeometricA(self%gamma1)                     , &
               &                                                                          [      1.5d0-0.5d0*self%gamma1]                       , &
               &                                                                          1.0d0/(self%subresolutionHypergeometric%x(i)+1.0d0)**2, &
               &                                                        toleranceRelative=self%precisionHypergeometric                            &
               &                                                       )                                                                        , &
               &                                    i                                                                                             &
               &                                   )
       end do
       self%subresolutionHypergeometricInitialized=.true.
    end if
    return
  end subroutine parkinsonColeHellySubresolutionHypergeometricTabulate
!--> 1004 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

!--> 1005 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
  subroutine parkinsonColeHellyUpperBoundHypergeometricTabulate(self,mass,massResolution,massMinimumIn,massMaximumIn)
!--> 1006 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    !!{
!<     Tabulate the hypergeometric term appearing in the upper bound branching probability rate expression.
    !!}
!--> 1009 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    use :: Hypergeometric_Functions, only : Hypergeometric_2F1
    use :: Numerical_Constants_Math, only : Pi
    use :: Table_Labels            , only : extrapolationTypeAbort
!--> 1012 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"
    implicit none
    class           (mergerTreeBranchingProbabilityParkinsonColeHelly), intent(inout)           :: self
    double precision                                                  , intent(in   )           :: mass                              , massResolution
    double precision                                                  , intent(in   ), optional :: massMinimumIn                     , massMaximumIn
    integer                                                           , parameter               :: massCountPerDecade =30
    double precision                                                  , parameter               :: sqrtTwoOverPi      =sqrt(2.0d0/Pi)
    double precision                                                                            :: massMinimum                       , massMaximum
    integer                                                                                     :: massCount                         , i
    logical                                                                                     :: tabulate
    double precision                                                                            :: massSigma                         , gammaEffective     , &
         &                                                                                         halfMassSigma                     , halfMassAlpha      , &
         &                                                                                         resolutionMassSigma               , resolutionMassAlpha
!--> 1024 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"

    tabulate=.false.
    if (.not.self%upperBoundHypergeometricInitialized) then
       tabulate=.true.
       if (present(massMinimumIn)) then
          massMinimum=massMinimumIn
       else
          massMinimum=           2.0d0*massResolution
       end if
       if (present(massMaximumIn)) then
          massMaximum=massMaximumIn
       else
          massMaximum=max(1.0d16,2.0d0*mass          )
       end if
    else
       if     (                                            &
            &   mass < self%upperBoundHypergeometric%x(+1) &
            &  .or.                                        &
            &   mass > self%upperBoundHypergeometric%x(-1) &
            & ) then
          tabulate=.true.
          massMinimum=                                        2.0d0*massResolution
          massMaximum=max(self%upperBoundHypergeometric%x(-1),2.0d0*mass          )
       end if
    end if
    if (tabulate) then
       self%massResolutionTabulated=massResolution
       massCount=int(log10(massMaximum/massMinimum)*dble(massCountPerDecade))+1
       if (.not.self%upperBoundHypergeometricInitialized) call self%upperBoundHypergeometric%destroy()
       call self%upperBoundHypergeometric%create(massMinimum,massMaximum,massCount,1,extrapolationType=spread(extrapolationTypeAbort,1,2))
       ! Evaluate σ and α at the mass resolution.
       call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(massResolution,self%timeParent,resolutionMassSigma,resolutionMassAlpha)
       do i=1,massCount
          ! Evaluate σ and α.
          call           self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(0.5d0*self%upperBoundHypergeometric%x(i),self%timeParent,halfMassSigma,halfMassAlpha)
          massSigma     =self%cosmologicalMassVariance_%rootVariance                      (      self%upperBoundHypergeometric%x(i),self%timeParent                            )
          gammaEffective=self%gamma1-1.0d0/halfMassAlpha
          call self%upperBoundHypergeometric%populate(                                                                             &
               &                                      +sqrtTwoOverPi                                                               &
               &                                      *resolutionMassSigma                                                         &
               &                                      /massSigma                                                                   &
               &                                      *(                                                                           &
               &                                        +(halfMassSigma/massSigma)**(+gammaEffective-1.0d0)                        &
               &                                        /                           (-gammaEffective+1.0d0)                        &
               &                                        *Hypergeometric_2F1(                                                       &
               &                                                                             self%hypergeometricA(gammaEffective), &
               &                                                                             [      1.5d0-0.5d0*gammaEffective]  , &
               &                                                                             (massSigma/halfMassSigma)**2        , &
               &                                                           toleranceRelative=self%precisionHypergeometric          &
               &                                                          )                                                        &
               &                                        -(resolutionMassSigma/massSigma)**(+gammaEffective-1.0d0)                  &
               &                                        /                                 (-gammaEffective+1.0d0)                  &
               &                                        *Hypergeometric_2F1(                                                       &
               &                                                                             self%hypergeometricA(gammaEffective), &
               &                                                                             [      1.5d0-0.5d0*gammaEffective]  , &
               &                                                                             (massSigma/resolutionMassSigma)**2  , &
               &                                                           toleranceRelative=self%precisionHypergeometric          &
               &                                                          )                                                        &
               &                                       )                                                                         , &
               &                                      i                                                                            &
               &                                     )
       end do
       self%upperBoundHypergeometricInitialized=.true.
    end if
    return
  end subroutine parkinsonColeHellyUpperBoundHypergeometricTabulate
!--> 1090 "/home/dsarnaaik/Galacticus/forked_galacticus/source/merger_trees.branching_probability.Parkinson_Cole_Helly.F90"


end submodule mergerTreeBranchingProbabilityParkinsonColeHelly_
