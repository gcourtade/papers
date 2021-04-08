# YASARA MACRO
# This macro is a modification of md_run.mcr by Elmar Krieger, combined with parts of a macro from Reinhard Wimmer.
# It runs a short (100 ps) restrained MD simulation on user-provided NMR-derived distance restraints and coordinate files for the protein BlLPMO10A. It was used to generate the PDB ID 6TWE.
# It also incorporates copper restraints from the following citation:
# Bissaro et al. (2018) "How a lytic polysaccharide monooxygenase binds crystalline chitin" Biochemistry 57:1893-1906
# 
# 
# LICENSE: GPL
#
# Gaston Courtade - Jan 2020
# ====================================================================================


filename='BlAA10-Cu'
restrainfile='(filename).tbl'
startfile='5lw4'

# The structure to simulate must be present with a .pdb or .sce extension.
# If a .sce (=YASARA scene) file is present, the cell must have been added.
# You can either set the target structure by clicking on Options > Macro > Set target,
# or by uncommenting the line below and specifying it directly.

# pH at which the simulation should be run, by default physiological pH
ph=5.5

# The ion concentration as a mass fraction, here we use 0.9% NaCl (physiological solution)
ions='Na,Cl,0.9'

# Simulation temperature
# If you run at a temperature that differs from 298K, you also need
# to adapt the pressure control below, look in the PressureCtrl documentation.
temperature='298K'

# Pressure control mode
# Default: Rescale the cell such that residues named HOH reach a density of 0.997 g/l. 
# For solvents other than water, you have to create your own solvent box
# as described in the FillCellObj documentation and save it as .._solvent.sce.
density=0.997
pressurectrl='SolventProbe,Name=HOH,Density=(density)'

# Alternative: Uncomment below to calculate the pressure from the virial and
# rescale the cell to reach a pressure of 1 bar. Use this method if you do not
# know the correct density. 
#pressurectrl='Manometer,Pressure=1' 

# Alternative: Do not control pressure
#pressurectrl='Off'

# The simulation speed, either 'normal' (conservative) or 'fast' (maximize performance)
# Do not use 'fast' if you simulate incorrect molecules (that would not be stable in reality) 
# 'if !count speed' simply checks if variable 'speed' as been defined previously (e.g. by the macro md_runfast) 
if !count speed
  speed='normal'

# The format used to save the trajectories: 'sim' or 'xtc'. If you choose the latter, a single
# *.sim restart file will be saved too, since XTC does not contain velocities, only positions
format='xtc'

# Duration of the simulation, alternatively use e.g. duration=5000 to simulate for 5000 picoseconds
duration=100

# Extension of the cell on each side of the protein
# '10' means that the cell will be 20 A larger than the protein.
# Cell settings only apply if you do not provide your own cell in a *.sce file.
extension=5

# Flag to use a cubic simulation cell. This makes sure that also elongated
# molecules can rotate freely during very long simulations. If only a short
# simulation is planned, it can be speeded up by setting the flag to 0,
# creating a rectangular cell that fits the solute more tightly.
cubic=1

# Forcefield to use (these are all YASARA commands, so no '=' used)
ForceField YASARA


# Cutoff
Cutoff 7.86

# Cell boundary
Boundary periodic

# Use longrange coulomb forces (particle-mesh Ewald)
Longrange Coulomb

# Number of simulation steps per screen and pairlist update
SimSteps 1,1

  # The normal restraining function (see RestrainPot command)
defaultpot='SoftSquare, SqConstant=1.00, SqOffset=0.0, SqExponent=2, rSwitch=1.00, SoExponent=1, Asymptote=2.00'

  # The normal restraining parameters (see RestrainPar command)
defaultpar='Average=Sum, Ceil=9999, Monomers=1, JoinDis=-1'

  # The maximum allowed distance restraint violation in A, 9999 to ignore
violdismax=0.25
  

  # The maximum allowed dihedral angle restraint violation in degrees, 9999 to ignore
violdihmax=30

  # Flag if cis-peptide bonds before prolines should be corrected.
  # If you do not want to potentially miss a cis-proline in your structure,
  # set this flag to 'No'. If the lowest energy structures in the ensemble
  # then all have a certain cis-proline, fix it in the cis-conformation.
  # In any case, create a new ensemble with the flag set to 'Yes' to avoid
  # a random collection of cis-prolines.
CorrectCis On,Proline='No'
CorrectIso On
CorrectConv On



defaultscale='Distance=25, Dihedral=2'

# Normally no change required below this point
# ============================================

RequireVersion 14.5.1

# Keep the solute from diffusing around and crossing periodic boundaries
CorrectDrift On

# Treat all simulation warnings as errors that stop the macro
WarnIsError Off

# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

Console off
logfile='(MacroTarget).log'
LogAs (logfile),append=Yes, print 'MD refinement with Cu FF'
Clear
LoadPDB (MacroTarget)
#Correction of naming deviations, residues in excess, etc. Highly individually different from protein to protein
AddTer
CleanAll
#Build the copper atom
SwapRes 1,Hid,Isomer=L
SwapRes 90,Hie,Isomer=L
BuildAtom Element=Copper,Copies=1,1 19 1322
JoinObj 2,1
SavePDB 1,(MacroTarget)-Cu.pdb

#Prepare the copper
CleanAll
#SaveSce (MacroTarget)_Cu.sce
Cell Auto,Extension=5
FillCellWater
SwapRes 1 and obj 1 and name his,Hie,Isomer=L
SwapRes 90 and obj 1 and name his,Hid,Isomer=L
#BuildAtom Element=Copper,Copies=1,1 13 1328
#JoinObj 2,1
EnergyUnit kcal/mol
AddSpring ND1 and obj 1 and obj 1 and res 1,Cu,Len=1.986,SFC=121.21
AddSpring NE2 and obj 1 and res 90,Cu,Len=1.986,SFC=121.21
AddSpring N and obj 1 and res 1,Cu,Len=2.078,SFC=64.50
AddAngle CE1 and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.9,BFC=82.8
AddAngle CE1 and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.9,BFC=82.8
AddAngle ND1 and obj 1 and res 1,Cu,NE2 and obj 1 and res 90,Min=174.9,BFC=52.5
AddAngle CG and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.3,BFC=84.4
AddAngle CG and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.3,BFC=84.4
AddAngle CD2 and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.2,BFC=83.7
AddAngle CD2 and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.2,BFC=83.7
AddAngle CA and obj 1 and res 1,N and obj 1 and res 1,Cu,Min=120.8,BFC=83.6
AddAngle H2 and obj 1 and res 1,N and obj 1 and res 1,Cu,Min=105.5,BFC=48.7
AddAngle N and obj 1 and res 1,Cu,ND1 and obj 1 and res 1,Min=91.0,BFC=74.1
AddAngle N and obj 1 and res 1,Cu,NE2 and obj 1 and res 90,Min=91.0,BFC=74.1
EnergyUnit kJ/mol
ChargeAtom Cu,e=0.1104
ChargeAtom N and obj 1 and res 1,e=-0.4740
ChargeAtom H3 and obj 1 and res 1,e=0.2246
ChargeAtom H2 and obj 1 and res 1 ,e=0.2246
ChargeAtom CA and obj 1 and res 1,e=0.1924
ChargeAtom HA and obj 1 and res 1 ,e=0.0045
ChargeAtom CB and obj 1 and res 1,e=-0.0999
ChargeAtom HB2 and obj 1 and res 1,e=0.0602
ChargeAtom HB3 and obj 1 and res 1,e=0.0602
ChargeAtom CG and obj 1 and res 1,e=0.0691
ChargeAtom ND1 and obj 1 and res 1,e=-0.0344
ChargeAtom CE1 and obj 1 and res 1,e=-0.0973
ChargeAtom HE1 and obj 1 and res 1,e=0.1949
ChargeAtom NE2 and obj 1 and res 1,e=-0.1678
ChargeAtom HE2 and obj 1 and res 1,e=0.3449
ChargeAtom CD2 and obj 1 and res 1,e=-0.2031
ChargeAtom HD2 and obj 1 and res 1,e=0.2025
ChargeAtom C and obj 1 and res 1,e=0.6123
ChargeAtom CB and obj 1 and res 90,e=-0.0462
ChargeAtom HB2 and obj 1 and res 90,e=0.01075
ChargeAtom HB3 and obj 1 and res 90,e=0.01075
ChargeAtom CG and obj 1 and res 90,e=0.2893
ChargeAtom ND1 and obj 1 and res 90,e=-0.2749
ChargeAtom HD1 and obj 1 and res 90,e=0.3609
ChargeAtom CE1 and obj 1 and res 90,e=-0.0936
ChargeAtom HE1 and obj 1 and res 90,e=0.1916
ChargeAtom NE2 and obj 1 and res 90,e=0.0317
ChargeAtom CD2 and obj 1 and res 90,e=-0.3208
ChargeAtom HD2 and obj 1 and res 90,e=0.1952
ChargeAtom CA and obj 1 and res 90,e=0.0188
Experiment Neutralization
  WaterDensity (density)
  pH (ph)
  Ions (ions)
  pKaFile (MacroTarget).pka
  Speed Fast
  Experiment On
  Wait ExpEnd
# Save scene with water

SaveSce (MacroTarget)_water.sce

# Choose timestep and activate constraints
if speed=='normal'
  # Normal simulation speed
  # Remove any constraints
  FreeBond all,all
  FreeAngle all,all,all
  # Choose a multiple timestep of 2*1.25 = 2.5 fs
  tslist=2,1.25
  # Save simulation snapshots every 25 fs
  saveinterval=25
ts=tslist1*tslist2
# Snapshots are saved every 'savesteps'
savesteps=saveinterval/ts
# Set final simulation parameters
TimeStep (tslist)
Temp (temperature)
# Make sure all atoms are free to move
FreeAll
# Uncomment to completely fix some atoms
# FixAtom Backbone Obj 1
# Alread a snapshot/trajectory present?
i=00000
if format=='sim'
  trajectfilename='(MacroTarget)(i).sim'
else  
  trajectfilename='(MacroTarget).xtc'
  restartfilename='(MacroTarget).sim'  
  # Backwards compatibility: Starting with YASARA version 12.8.1, XTC trajectories no longer contain a number in the filename
  old = FileSize (MacroTarget)(i).xtc
  if old
    RenameFile (MacroTarget)(i).xtc,(trajectfilename)
running = FileSize (trajectfilename)
if not running
  SwapRes 1 and obj 1 and name his,Hie,Isomer=L
  SwapRes 90 and obj 1 and name his,Hid,Isomer=L
  #BuildAtom Element=Copper,Copies=1,1 13 1328
  #JoinObj 2,1
  EnergyUnit kcal/mol
  AddSpring ND1 and obj 1 and obj 1 and res 1,Cu,Len=1.986,SFC=121.21
  AddSpring NE2 and obj 1 and res 90,Cu,Len=1.986,SFC=121.21
  AddSpring N and obj 1 and res 1,Cu,Len=2.078,SFC=64.50
  AddAngle CE1 and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.9,BFC=82.8
  AddAngle CE1 and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.9,BFC=82.8
  AddAngle ND1 and obj 1 and res 1,Cu,NE2 and obj 1 and res 90,Min=174.9,BFC=52.5
  AddAngle CG and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.3,BFC=84.4
  AddAngle CG and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.3,BFC=84.4
  AddAngle CD2 and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.2,BFC=83.7
  AddAngle CD2 and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.2,BFC=83.7
  AddAngle CA and obj 1 and res 1,N and obj 1 and res 1,Cu,Min=120.8,BFC=83.6
  AddAngle H2 and obj 1 and res 1,N and obj 1 and res 1,Cu,Min=105.5,BFC=48.7
  AddAngle N and obj 1 and res 1,Cu,ND1 and obj 1 and res 1,Min=91.0,BFC=74.1
  AddAngle N and obj 1 and res 1,Cu,NE2 and obj 1 and res 90,Min=91.0,BFC=74.1
  EnergyUnit kJ/mol
  ChargeAtom Cu,e=0.1104
  ChargeAtom N and obj 1 and res 1,e=-0.4740
  ChargeAtom H3 and obj 1 and res 1,e=0.2246
  ChargeAtom H2 and obj 1 and res 1 ,e=0.2246
  ChargeAtom CA and obj 1 and res 1,e=0.1924
  ChargeAtom HA and obj 1 and res 1 ,e=0.0045
  ChargeAtom CB and obj 1 and res 1,e=-0.0999
  ChargeAtom HB2 and obj 1 and res 1,e=0.0602
  ChargeAtom HB3 and obj 1 and res 1,e=0.0602
  ChargeAtom CG and obj 1 and res 1,e=0.0691
  ChargeAtom ND1 and obj 1 and res 1,e=-0.0344
  ChargeAtom CE1 and obj 1 and res 1,e=-0.0973
  ChargeAtom HE1 and obj 1 and res 1,e=0.1949
  ChargeAtom NE2 and obj 1 and res 1,e=-0.1678
  ChargeAtom HE2 and obj 1 and res 1,e=0.3449
  ChargeAtom CD2 and obj 1 and res 1,e=-0.2031
  ChargeAtom HD2 and obj 1 and res 1,e=0.2025
  ChargeAtom C and obj 1 and res 1,e=0.6123
  ChargeAtom CB and obj 1 and res 90,e=-0.0462
  ChargeAtom HB2 and obj 1 and res 90,e=0.01075
  ChargeAtom HB3 and obj 1 and res 90,e=0.01075
  ChargeAtom CG and obj 1 and res 90,e=0.2893
  ChargeAtom ND1 and obj 1 and res 90,e=-0.2749
  ChargeAtom HD1 and obj 1 and res 90,e=0.3609
  ChargeAtom CE1 and obj 1 and res 90,e=-0.0936
  ChargeAtom HE1 and obj 1 and res 90,e=0.1916
  ChargeAtom NE2 and obj 1 and res 90,e=0.0317
  ChargeAtom CD2 and obj 1 and res 90,e=-0.3208
  ChargeAtom HD2 and obj 1 and res 90,e=0.1952
  ChargeAtom CA and obj 1 and res 90,e=0.0188
  # Perform energy minimization
  Experiment Minimization
  Experiment On
  Wait ExpEnd

  # And now start the real simulation
  SwapRes 1 and obj 1 and name his,Hie,Isomer=L
  SwapRes 90 and obj 1 and name his,Hid,Isomer=L
  #BuildAtom Element=Copper,Copies=1,1 13 1328
  #JoinObj 2,1
  EnergyUnit kcal/mol
  AddSpring ND1 and obj 1 and and res 1,Cu,Len=1.986,SFC=121.21
  AddSpring NE2 and obj 1 and res 90,Cu,Len=1.986,SFC=121.21
  AddSpring N and obj 1 and res 1,Cu,Len=2.078,SFC=64.50
  AddAngle CE1 and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.9,BFC=82.8
  AddAngle CE1 and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.9,BFC=82.8
  AddAngle ND1 and obj 1 and res 1,Cu,NE2 and obj 1 and res 90,Min=174.9,BFC=52.5
  AddAngle CG and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.3,BFC=84.4
  AddAngle CG and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.3,BFC=84.4
  AddAngle CD2 and obj 1 and res 90,NE2 and obj 1 and res 90,Cu,Min=126.2,BFC=83.7
  AddAngle CD2 and obj 1 and res 1,ND1 and obj 1 and res 1,Cu,Min=126.2,BFC=83.7
  AddAngle CA and obj 1 and res 1,N and obj 1 and res 1,Cu,Min=120.8,BFC=83.6
  AddAngle H2 and obj 1 and res 1,N and obj 1 and res 1,Cu,Min=105.5,BFC=48.7
  AddAngle N and obj 1 and res 1,Cu,ND1 and obj 1 and res 1,Min=91.0,BFC=74.1
  AddAngle N and obj 1 and res 1,Cu,NE2 and obj 1 and res 90,Min=91.0,BFC=74.1
  EnergyUnit kJ/mol
  ChargeAtom Cu,e=0.1104
  ChargeAtom N and obj 1 and res 1,e=-0.4740
  ChargeAtom H3 and obj 1 and res 1,e=0.2246
  ChargeAtom H2 and obj 1 and res 1 ,e=0.2246
  ChargeAtom CA and obj 1 and res 1,e=0.1924
  ChargeAtom HA and obj 1 and res 1 ,e=0.0045
  ChargeAtom CB and obj 1 and res 1,e=-0.0999
  ChargeAtom HB2 and obj 1 and res 1,e=0.0602
  ChargeAtom HB3 and obj 1 and res 1,e=0.0602
  ChargeAtom CG and obj 1 and res 1,e=0.0691
  ChargeAtom ND1 and obj 1 and res 1,e=-0.0344
  ChargeAtom CE1 and obj 1 and res 1,e=-0.0973
  ChargeAtom HE1 and obj 1 and res 1,e=0.1949
  ChargeAtom NE2 and obj 1 and res 1,e=-0.1678
  ChargeAtom HE2 and obj 1 and res 1,e=0.3449
  ChargeAtom CD2 and obj 1 and res 1,e=-0.2031
  ChargeAtom HD2 and obj 1 and res 1,e=0.2025
  ChargeAtom C and obj 1 and res 1,e=0.6123
  ChargeAtom CB and obj 1 and res 90,e=-0.0462
  ChargeAtom HB2 and obj 1 and res 90,e=0.01075
  ChargeAtom HB3 and obj 1 and res 90,e=0.01075
  ChargeAtom CG and obj 1 and res 90,e=0.2893
  ChargeAtom ND1 and obj 1 and res 90,e=-0.2749
  ChargeAtom HD1 and obj 1 and res 90,e=0.3609
  ChargeAtom CE1 and obj 1 and res 90,e=-0.0936
  ChargeAtom HE1 and obj 1 and res 90,e=0.1916
  ChargeAtom NE2 and obj 1 and res 90,e=0.0317
  ChargeAtom CD2 and obj 1 and res 90,e=-0.3208
  ChargeAtom HD2 and obj 1 and res 90,e=0.1952
  ChargeAtom CA and obj 1 and res 90,e=0.0188
  Sim On
HideMessage
  
# Set temperature and pressure control
TempCtrl Rescale
PressureCtrl (pressurectrl)

# Uncomment to add distance constraints
# Load the restraints
RestrainPar (defaultpar)
RestrainPot (defaultpot)
LoadTbl (filename)-upl.tbl,1
LoadTbl (filename)-aco.tbl,1
SaveTbl 1,(restrainfile)
ScaleRest (defaultscale)

# And finally, make sure that future snapshots are saved
#Save(format) (trajectfilename),(savesteps)
#if format=='xtc'
  # We save an XTC trajectory plus a single SIM restart file with velocities
  #SaveSim (restartfilename),(savesteps),Number=no

if duration=='forever'
  Console On
  Wait forever
else
  Console Off
  # Wait for given number of picoseconds
  do
    Wait 10
    t = Time
  while t<(duration)*1000+1
  Sim Off
  DelAtom Element H
  AddBond N and obj 1 and res 1,Cu,1.00
  Clean
  SavePDB 1,(MacroTarget)_(duration)ps.pdb

  # Exit YASARA if this macro was provided as command line argument in console mode
  if runWithMacro and ConsoleMode
    Exit

  