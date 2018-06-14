// Work started 08/14/2012
// New code cleaned on June 10, 2013. 
// For past cpmd codes refer to versions before June 5, 2013
// NB joined fun 05/01/2017.
// This is main.
// This is (Car Parrinello) molecular dynamics for electrostatics on membranes

// VJ Lessons Learned, other notes (pre-May 2017):
/*
Lessons:

1. Do not use exit(1) to exit the program. This leads to memory leaks.

Other Notes:
Aug 22 2012 - NOTE THREADSIZE different in different regions of code.
for precal it is 8, for fmd forces, energies it is 4.
Testing parameters: echo 3 4 100 35 0.0001 300 0.8 0 1 1
*/

//  NB Changes (May 2017 - Present):
/* Changes :

2017.05.12:  NB has edited main.cpp to incorporate a single (1) type of fractionally ionizable group with parameter (alpha). This includes ability to echo (input) this new parameter (alpha) below and to call function "new assign_random_q_values" (found in Interface.cpp) using this parameter. [Preliminary studies reporting multiple groups have assumed the groups to have the same alpha values, i.e. one group is basic the other acidic and you are at pH = Mean[pKa1, pKa2].]

2017.06.05:  NB has added working area constraint functions (SHAKE_for_area & RATTLE_for_area).  Area is conserved to within a fraction of a percent with default masses = 1.0, etc.  

                2017.09.20:  SHAKE_for_area & RATTLE_for_area found to be in error, math mistake.  Only coincide

2017.06.13:  NB has added a pre-MD line below to dump the initial state in the movie file as well; it calls this timestep "-1" in the movie file due to the way the function "interface_movie" is defined.

2017.06.26:  NB is working to make the optimization dynamics a real-time study.  Because the temperature uniformly started at (.5T), an excess factor of 2 in the denominator of the standard deviation of the initial velocity distribution "p_sigma" was removed (within "functions.cpp").  The temperature now begins sufficiently close to the intended temperature (within 10% @ t = 1 step).  Energy conservation is improved after this correction. [Adding a factor of 2 in the denominator of all DOF's supplied to the thermostat also had this effect, but it seems correcting "p_sigma" - not thermostat DOF's - is the correct solution.  Some sources imply DOF's should be 3N - 1 due to constraint, however energy conservation seems to work marginally better without subtracting the constraint DOF, so subtraction may be wrong.]

2017.06.28:  NB has finished adapting the code to a real-time study (sans barostat).  The parameters that optimally conserve energy for discretizations up to (3,8) are {C = 5,  T = 1.0, Q = 0.01, M = 1.0, dt = .00001, DOF = 3N}.

2017.07.04:  NB has added a for statement allowing "assign_random_q_values" to create 5 divisions.

2017.07.08:  NB made a few trivial changes such as 1) making "total_time" now specify the number of steps to run itself rather than a function of (dt) because it would have an error if the ratio became to large and 2) eliminating the unused variable M introduced to specify vertex mass (we must keep as 1.0).  Additionally, a character flag ("annealFlag") has been added (that must be echoed after all numerical parameters) that signifies if annealing should occur ('y' or 'n').

2017.07.09:  NB is refining functionality to continue a simulation from a pre-existing off-file ("Optimized.off", must have this name). This involved adding a second character flag (that also must be echoed, after the numerical parameters and annealing flag) signifying if said off-file is to be loaded. If you use Nicholas' automated-PBS script generator, it will automatically denote the values of these flags in the cluster job name (and warn you if you are submitting a RT job that anneals).

2017.08.08:  NB has randomized the charge distribution assignment scheme (rather than using permutations as it introduces significant stripes and polar deficits in charge density & potential, though it produces similar shapes for the case tested.)  Also, I am changing the function for generation of axially distributed patches of charge (spherical caps).

2017.09.07:  NB has verified all mathematical components of the volume and area constraints.  All are correct, a factor of (1/2) is missing in the face area gradient within the "compute_neg_VGrads" function; however, adding it here causes major issues, meaning it is probably within those functions already (must check).  It is only used in bending energy calculations and (if area constrained) the "SHAKE/RATTLE_for_area" functions.  Adding (1/2) factors in the appropriate place within the "SHAKE/RATTLE_for_area" seems to be fine.

*/

#include "utility.h"
#include "interface.h"
#include "vertex.h"
#include "edge.h"
#include "face.h"
#include "control.h"
#include "functions.h"
#include "thermostat.h"
#include <boost/filesystem/operations.hpp>  // 2017.09.13 NB added:  For filesystem manipulation (currently just deletes outfiles, even if inhabited), then recreates an empty one.  // Working in CLion, must check on HPC.  My Boost version:
#include <boost/program_options.hpp>

using namespace boost::program_options;

// Declaring function to initiate and propagate the MD (NB added char = constraint flag):
void md_interface(INTERFACE&, vector<THERMOSTAT>&, CONTROL&, char, char);

int main(int argc, const char* argv[])
{
  // (2017.09.14 NB added.)  Delete the existing outfiles directory and any files it contains from a previous run, then recreate it empty.
  if(boost::filesystem::remove_all("outfiles") == 0) perror("Error deleting outfiles directory");
  else cout << "Pre-existing outfiles directory successfully, replaced with an empty directory." << endl;
  boost::filesystem::create_directory("outfiles");

  // Declare variables to make the system:
  double ein = epsilon_water; 		// permittivity of inside medium
  double eout = epsilon_water; 		// permittivity of outside medium
  double T, Q, dt;				        // Temperature of the system, mass of thermostate (0 if inactive), timestep.
  unsigned int chain_length_real;	// length of nose-hoover thermostat chain, 1 minimum, 1 means no thermostat
  
  INTERFACE boundary;			        // interface (modelled as soft LJ wall)
  unsigned int disc1, disc2;
  double lambda_a, lambda_v;		  // lambdas measuring strength of area & volume constraints  (unused)

  double q_strength, alpha, conc_out, z_out; // net charge (if all charged), fractional q-occupancy, salt conc (MOLAR), salt valency
  char annealFlag, offFlag, geomConstraint, constraintForm;
  
  int total_time, anneal_freq;

  // Specify variables via command line (-X x):
  options_description desc("Usage:\nrandom_mesh <options>");
  desc.add_options()
          ("help,h", "print usage message")
          ("netCharge,q", value<double>(&q_strength)->default_value(900),
           "Net charge on the particle if fully occupied (elementary charges).")    // enter in nanometers
          ("saltConc,c", value<double>(&conc_out)->default_value(0.02),
           "The bulk salt concentration (Molar).")
          ("tensSigma,t", value<double>(&boundary.sigma_a)->default_value(3),
          "The surface tension constant (dynes/cm).")
          ("bending,b", value<double>(&boundary.bkappa)->default_value(30),
          "The bending modulus of the particle.")
          ("stretch,s", value<double>(&boundary.sconstant)->default_value(30),
          "The stretching modulus of the particle.")
          ("totTime,D", value<int>(&total_time)->default_value(250000),
          "The duration of the simulation (total timesteps).")
          ("initTemp,T", value<double>(&T)->default_value(0.01),
           "The target initial temperature (prior to annealing).")
          ("annealFlag,a", value<char>(&annealFlag)->default_value('y'),
           "If annealing occurs or not ('y' for yes).");
  variables_map vm;
  store(parse_command_line(argc, argv, desc), vm);
  notify(vm);

  //cin >> q_strength >> conc_out >> boundary.sigma_a >> boundary.bkappa >> boundary.sconstant >> total_time >> T >> annealFlag;

  disc1 = 3;                                          //  Discretization parameters (h, k).
  disc2 = 8;
  //q_strength = 800;                                   //  Bare membrane net charge (if fully occupied).
  alpha = 1.0;                                        //  Fractional charge occupancy.
  //conc_out = 0.06;                                    //  Salt concentration in Molar.
  //boundary.sigma_a = 3;                               //  Desired tension in (dyne/cm) (multiplier enforces this below).
  //boundary.bkappa = 30;                               //  Reduced bending modulus.
  //boundary.sconstant = 30;                            //  Reduced stretching modulus.
  //total_time = 250000;                              // Net number of timesteps.
  chain_length_real = 5;                            //  Number of global entities (C = thermostat count + 1).
  //T = 0.0005;                                         // Target temperature of thermostat (if C != 1).
  if(chain_length_real == 1) Q = 0;                 //  Reduced mass of the thermostat particle(s).
  else Q = 10*T;                                     // Thermostat mass must be set to zero for safety if (C = 1).
  dt = 0.001;                                       //  Timestep coefficient.
  //annealFlag = 'y';
  offFlag = 'n';

  boundary.sigma_a = (24.3053 * boundary.sigma_a);    //  Coefficient is 1 (dyne/cm) in (kB T_room/(10nm^2)).

  // Control of the dynamics:
  CONTROL mdremote;
  mdremote.timestep = dt;				    // NOTE adjust time step for faster simulation ... testing purposes
  mdremote.steps = total_time;      //int(total_time/mdremote.timestep);  NB changed, error (negative) if ratio too great.
  mdremote.anneal = annealFlag;	    // If set to 'n', no annealing, simulation proceeds at fixed temperature.
  mdremote.annealfreq = 50000;//anneal_freq;	      // Interval number of steps between which annealing commences.
  mdremote.annealDuration = 1;      // Number of steps over which to reduce (T = fT*T) & (Q = fQ*Q), same as before if 1.
  mdremote.TAnnealFac = 1.25;        // Fold reduction to temperature at initial annealing call.  Was 10 before (and constant throughout).
  mdremote.QAnnealFac = 1.00 + (mdremote.TAnnealFac / 10.0);  // Maintain previous scaling, same as before (fQ = 2) if fT = 10.

  // Constrain the dynamics:
  geomConstraint = 'V';
  constraintForm = 'L';

  // ### Simulation setup: ###
    // Unused quantities for antiquated constraints & energy functionals:
  lambda_v = 0;	                    // hardwiring volume constraint. Unused.
  lambda_a = 0;	                    // hardwiring area constraint. Unused.
  boundary.lambda_l = 0;            // forgot what was this?!  NB:  I think its an edge line tension force. Unused.
  boundary.lj_length = 0.1;         // Reduced LJ diameter associated with vertex; safety reasons, rare to come near as they repel.

    // Electrostatics setup (membrane, implicit salt ions):
  z_out = 1;                        // Valency of implicit (salt) ions.
  boundary.ein = ein;               // Permittivity inside the membrane.
  boundary.eout = eout;             // Permittivity outside the membrane.  (Both equal by default.)
  boundary.set_up();
      // Assigning the Debye length outside (in reduced units):
  if (conc_out == 0)
    boundary.inv_kappa_out = 100000000;  // Effectively infinite, no screening.
  else
    boundary.inv_kappa_out = (0.257 / sqrt(z_out * boundary.lB_out * unit_radius_sphere * conc_out)) / unit_radius_sphere;
  
    // Generate the membrane:
  boundary.discretize(disc1, disc2);			// discretize the interface
  if (disc1 != 0 || disc2 != 0)
    boundary.assign_random_q_values(1,q_strength,alpha);		// number of components input
    
  boundary.dressup(lambda_a, lambda_v);			// dress the interface with normals,... 
  
  boundary.ref_area = boundary.total_area; // NB changed from (4 * pi), discrete membrane differs.
  boundary.ref_volume = boundary.total_volume;
  boundary.ref_Area_Vertices = boundary.total_Area_Vertices;    // (2017.09.19 NB added.)  Initial area by vertices.

  // stretching parameters
  // boundary.l0 = boundary.avg_edge_length;
  
  // LJ length
  boundary.lj_length = boundary.lj_length * boundary.avg_edge_length;
  boundary.elj = 1.0;

  // Thermostat for MD of the interface:
  vector<THERMOSTAT> real_bath;		// thermostat for real system (interface + ions) ; no ions yet.
  if (chain_length_real == 1)
    real_bath.push_back(THERMOSTAT(0,T,3*boundary.V.size(),0.0,0,0));		
  else
  {
    real_bath.push_back(THERMOSTAT(Q, T, 3*boundary.V.size(), 0, 0, 0));
    while (real_bath.size() != chain_length_real - 1)
      real_bath.push_back(THERMOSTAT(Q, T, 1, 0, 0, 0));				
    real_bath.push_back(THERMOSTAT(0,T,3*boundary.V.size(),0.0,0,0));	// dummy bath (has 0 mass)
  }
  
  // Assign masses to vertices
  for (unsigned int i = 0; i < boundary.V.size(); i++)
    boundary.V[i].m = 1.0;      // NB:  Energy conservation fails invariably if not = 1.0 .
  
  // NB added following lines to load off-file (does not warn if the file is not found, yet):
  if(offFlag == 'y') boundary.load_configuration("380000.off");

  // Define an output stream to provide the same things as the (*.out) file on the HPC cluster (minus run outputs):
    // Useful in automated analysis when console output is not saved by default on local PC.
  ofstream list_out("outfiles/simulation.out", ios::app);

  // ### Output to the console (or HPC *.out file) & backup "simulation.out" file) the simulation parameters: ###
  cout << "\n===================\nMembrane\n===================\n";
  cout << "Unit radius of the sphere (in nm): " << unit_radius_sphere << endl; //NB added in observing size change.
  cout << "Number of vertices: " << boundary.V.size() << endl;
  cout << "Number of edges: " << boundary.E.size() << endl;
  cout << "Number of faces: " << boundary.F.size() << endl;
  cout << "Lambda_v: " << lambda_v << endl;
  cout << "Lambda_a: " << lambda_a << endl;
  cout << "Lambda_l: " << boundary.lambda_l << endl;
  cout << "Total intial volume: " << boundary.total_volume << "  Sphere volume: " << (4./3.) * pi << "  Ref volume: " << boundary.ref_volume << endl;
  cout << "Total intial face area: " << boundary.total_area << "  Sphere area: " << 4 * pi << "  Ref area: " << boundary.ref_area << endl;
  cout << "Total intial vertex area: " << boundary.total_Area_Vertices << "  Sphere area: " << 4 * pi << "  Ref area: " << boundary.ref_Area_Vertices << endl;
  cout << "Bending rigidity: " << boundary.bkappa << endl;
  cout << "Stretching constant: " << boundary.sconstant << endl;
  cout << "Surface Tension constant: " << (boundary.sigma_a / 24.3053) << endl;
  cout << "Unstretched edge length: " << boundary.avg_edge_length << endl; // NB uncommented & replaced the output value.
  cout << "LJ strength: " << boundary.elj << endl;
  cout << "LJ distance cutoff: " << boundary.lj_length << endl;
  
  cout << "\n===================\nElectrostatics\n===================\n";
  cout << "Total charge on membrane: " << q_strength << endl;
  cout << "Fractional charge occupancy (alpha): " << alpha << endl;
  //cout << "The number of patches (divisions) is: " << num_divisions << endl;
  cout << "Salt concentration: " << conc_out << endl;
  cout << "Debye length: " << boundary.inv_kappa_out << endl;
  cout << "Bjerrum length: " << boundary.lB_out << endl;
 
  cout << "\n===================\nThermodynamics\n===================\n";
  cout << "Initial Temperature (T): " << T << endl;
  cout << "Thermostat mass (Q): " << real_bath[0].Q << endl;
  cout << "Number of thermostats (C): " << real_bath.size() << endl;
  cout << "Timestep (dt): " << mdremote.timestep << "  &  Total steps: " << mdremote.steps << "." << endl;
  cout << "Annealing: " << mdremote.anneal << endl;
  cout << "Anneal frequency: " << mdremote.annealfreq << endl;
  cout << "Anneal duration: " << mdremote.annealDuration << endl;
  cout << "Temperature annealing decrement: " << mdremote.TAnnealFac << endl;
  cout << "Thermostat mass decrement: " << mdremote.QAnnealFac << endl;
  cout << "Off-file loaded: " << offFlag << endl;
    // Same thing in a file even if you're not on an HPC that can store the console output automatically on request:
  list_out << "===================\nMembrane\n===================\n";
  list_out << "Unit radius of the sphere (in nm): " << unit_radius_sphere << endl; //NB added in observing size change.
  list_out << "Number of vertices: " << boundary.V.size() << endl;
  list_out << "Number of edges: " << boundary.E.size() << endl;
  list_out << "Number of faces: " << boundary.F.size() << endl;
  list_out << "Lambda_v: " << lambda_v << endl;
  list_out << "Lambda_a: " << lambda_a << endl;
  list_out << "Lambda_l: " << boundary.lambda_l << endl;
  list_out << "Total intial volume: " << boundary.total_volume << "  Sphere volume: " << (4./3.) * pi << "  Ref volume: " << boundary.ref_volume << endl;
  list_out << "Total intial face area: " << boundary.total_area << "  Sphere area: " << 4 * pi << "  Ref area: " << boundary.ref_area << endl;
  list_out << "Total intial vertex area: " << boundary.total_Area_Vertices << "  Sphere area: " << 4 * pi << "  Ref area: " << boundary.ref_Area_Vertices << endl;
  list_out << "Bending rigidity: " << boundary.bkappa << endl;
  list_out << "Stretching constant: " << boundary.sconstant << endl;
  list_out << "Surface Tension constant: " << (boundary.sigma_a / 24.3053) << endl;
  list_out << "Unstretched edge length: " << boundary.avg_edge_length << endl; // NB uncommented & replaced the output value.
  list_out << "LJ strength: " << boundary.elj << endl;
  list_out << "LJ distance cutoff: " << boundary.lj_length << endl;

  list_out << "\n===================\nElectrostatics\n===================\n";
  list_out << "Total charge on membrane: " << q_strength << endl;
  list_out << "Fractional charge occupancy (alpha): " << alpha << endl;
  //list_out << "The number of patches (divisions) is: " << num_divisions << endl;
  list_out << "Salt concentration: " << conc_out << endl;
  list_out << "Debye length: " << boundary.inv_kappa_out << endl;
  list_out << "Bjerrum length: " << boundary.lB_out << endl;

  list_out << "\n===================\nThermodynamics\n===================\n";
  list_out << "Initial Temperature (T): " << T << endl;
  list_out << "Thermostat mass (Q): " << real_bath[0].Q << endl;
  list_out << "Number of thermostats (C): " << real_bath.size() << endl;
  list_out << "Timestep (dt): " << mdremote.timestep << "  &  Total steps: " << mdremote.steps << "." << endl;
  list_out << "Annealing: " << mdremote.anneal << endl;
  list_out << "Anneal frequency: " << mdremote.annealfreq << endl;
  list_out << "Anneal duration: " << mdremote.annealDuration << endl;
  list_out << "Temperature annealing decrement: " << mdremote.TAnnealFac << endl;
  list_out << "Thermostat mass decrement: " << mdremote.QAnnealFac << endl;
  list_out << "Off-file loaded: " << offFlag << endl;

  // (NB added) Print the initial state in the movie file, prior to MD:
  interface_movie(0, boundary.V, boundary);  // NB added.  Note, calls this "time step -1" in movie file.

  // Initiate MD of the boundary/membrane:
  boundary.dressup(lambda_a, lambda_v);
  boundary.output_configuration();    // NB added.  More information on the initial state.
  boundary.compute_local_energies();  // NB added to output energy profiles pre-MD & rename the files to avoid conflict.
  rename("outfiles/local_electrostatic_E.off","outfiles/local_electrostatic_E_initial.off");
  rename("outfiles/local_elastic_E.off","outfiles/local_elastic_E_initial.off");
  md_interface(boundary,real_bath, mdremote, geomConstraint, constraintForm);
  boundary.compute_local_energies();
  
  return 0;
} 
// End of main

