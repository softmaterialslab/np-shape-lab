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
#include "particle.h"
#include "control.h"
#include "functions.h"
#include "thermostat.h"
#include <boost/filesystem/operations.hpp>
#include <boost/program_options.hpp>

using namespace boost::program_options;

// Declaring function to initiate and propagate the MD:
void md_interface(INTERFACE &, vector<PARTICLE> &, vector<THERMOSTAT> &, vector<THERMOSTAT> &, CONTROL &, char, char, char, const double scalefactor, double box_radius);

//MPI boundary parameters
unsigned int lowerBoundMesh;
unsigned int upperBoundMesh;
unsigned int sizFVecMesh;
//unsigned int extraElementsMesh;
int lowerBoundIons;
int upperBoundIons;
int sizFVecIons;
unsigned int extraElementsIons;
mpi::environment env;
mpi::communicator world;

int main(int argc, const char *argv[]) {
    // (2017.09.14 NB added.)  Delete existing outfiles directory if it exists, then recreate it empty.
    if (world.rank() == 0) {
        if (boost::filesystem::remove_all("outfiles") == 0) perror("Error deleting outfiles directory");
        else cout << "Pre-existing outfiles directory successfully, replaced with an empty directory." << endl;
        boost::filesystem::create_directory("outfiles");
    }
    // Declare variables to make the system:
    double ein = epsilon_water;        // permittivity of inside medium
    double eout = epsilon_water;        // permittivity of outside medium
    double T, Q_mesh, Q_ion;                        // Temperature of the system, mass of thermostate (0 if inactive), timestep.
    unsigned int chain_length;    // length of nose-hoover thermostat chain, 1 minimum, 1 means no thermostat

    INTERFACE boundary;                    // interface (modelled as soft LJ wall)
    unsigned int disc1, disc2;
    double lambda_a, lambda_v;          // lambdas measuring strength of area & volume constraints  (unused)

    double unit_radius_sphere, youngsModulus, q_strength, alpha, conc_out, z_out; // radius (in nm), net charge (if all charged), fractional q-occupancy, salt conc (MOLAR), salt valency
    int numPatches;
    double fracChargedPatch;
    char randomFlag, offFlag, geomConstraint, bucklingFlag, constraintForm, counterionFlag, functionFlag;  // functionFlag indicates different pattern initializations, -y for yingyang pattern, -c for cube formation.
    string externalPattern;

    // Counterion Initializations:
    double packing_fraction, box_radius;
    int counterion_valency;
    double counterion_diameter = 0.6; // Counterion diameter in nanometers.
    vector<PARTICLE> counterions;

    // Control of the dynamics:
    CONTROL mdremote;

    // Specify variables via command line (-X x):
    options_description desc("Usage:\nrandom_mesh <options>");
    desc.add_options()
            ("help,h", "print usage message")
                 // Physical Parameters:
            ("unitRadius,R", value<double>(&unit_radius_sphere)->default_value(10),
             "Radius of the initial sphere & simulation unit of length (in nanometers).")
            ("netCharge,q", value<double>(&q_strength)->default_value(600),
             "Net NP charge, if fully occupied (elementary charges).")
            ("saltConc,c", value<double>(&conc_out)->default_value(0.005),
             "Salt concentration (Molar).")
            ("tensSigma,t", value<double>(&boundary.sigma_a)->default_value(1),
             "Surface tension constant (dynes/cm).")
            ("volTensSigma,v", value<double>(&boundary.sigma_v)->default_value(1),
             "Volume tension constant (atm/nm^3).")
            ("bending,b", value<double>(&boundary.bkappa)->default_value(40),
             "The bending modulus of the particle (kB*T).")
            ("Stretching,s", value<double>(&boundary.sconstant)->default_value(40),
             "Reduced stretching modulus of the particle (kB*T/R0^2).")
            ("GeomConstraint,G", value<char>(&geomConstraint)->default_value('N'),
             "Specification of rigid geometric constraints, 'V' for volume.")
            ("functionFlag,H", value<char>(&functionFlag)->default_value('n'),
             "specification of function, 'y' for yinyang function.")
            ("bucklingFlag,B", value<char>(&bucklingFlag)->default_value('n'),
             "Specification of stretching form, if spontaneous buckling should occur.")
                 // Physical Parameters for patterned (Janus, Striped, Polyhedral) particles:
            ("numPatches,N", value<int>(&numPatches)->default_value(1),
             "The number of distinct charge patches (of tunable size if N = 2).")
            ("fracChargePatch,p", value<double>(&fracChargedPatch)->default_value(0.5),
             "Surface area fraction of the patch (if N = 2, otherwise irrelevant).")
            ("externalPattern,E", value<string>(&externalPattern)->default_value("None"),
             "Filename of pattern to be imported (without the '.dat' extension).")
                // Physical Parameters for counterion condensation:
            ("counterionFlag,I", value<char>(&counterionFlag)->default_value('n'),
             "Flag specifying if explicit counterions should be present, 'y' for yes.")
            ("packingFraction,P", value<double>(&packing_fraction)->default_value(0.01),
             "NP packing fraction determining total box size.")
            ("counterionValency,Z", value<int>(&counterion_valency)->default_value(1),
             "Valency of the counterions.")
                // Virtual Parameters:
            ("totalTime,S", value<int>(&mdremote.steps)->default_value(250000),
             "Duration of the simulation (total timesteps).")
            ("timestep,d", value<double>(&mdremote.timestep)->default_value(0.001),
             "Time step used in the simulation.")
            ("initTemp,T", value<double>(&T)->default_value(0.001),
             "Target initial temperature (prior to annealing).")
            ("disc2,D", value<unsigned int>(&disc2)->default_value(8),
             "Discretization parameter 2")
            ("chain_length,C", value<unsigned int>(&chain_length)->default_value(5),
             "Number of thermostat particles in the Nose-Hoover chain (plus 1).")
            ("meshThermostatMass,Q", value<double>(&Q_mesh)->default_value(.01),
             "Mass of the mesh thermostat (higher means less strongly coupled).")
            ("ionsThermostatMass,Y", value<double>(&Q_ion)->default_value(.01),
             "Mass of the ions thermostat (higher means less strongly coupled).")
                 // Annealing Parameters:
            ("annealFlag,a", value<char>(&mdremote.anneal)->default_value('y'),
             "If annealing occurs or not ('y' for yes).")
            ("anneal_freq,f", value<int>(&mdremote.annealfreq)->default_value(50000),
             "Interval number of steps between which annealing commences.")
            ("anneal_dura,u", value<int>(&mdremote.annealDuration)->default_value(5000),
             "Number of steps over which to reduce (T = fT*T) & (Q_mesh = fQ*Q_mesh), same as before if 1")
            ("anneal_Tfac,e", value<double>(&mdremote.TAnnealFac)->default_value(1.25),
             "Fold reduction to temperature at initial annealing call.  Was 10 before (and constant throughout)")
                 // Runtime & Output Parameters:
            ("randomFlag,F", value<char>(&randomFlag)->default_value('y'),
             "If predefined shape used or not ('y' for yes).")
            ("offFlag,o", value<char>(&offFlag)->default_value('n'),
             "If predefined shape used or not ('y' for yes).")
            ("moviestart,m", value<int>(&mdremote.moviestart)->default_value(1),
             "The starting point of the movie")
            ("offfreq,O", value<int>(&mdremote.offfreq)->default_value(2500),
             "The frequency of making off files")
            ("moviefreq,M", value<int>(&mdremote.moviefreq)->default_value(1000),
             "The frequency of shooting the movie")
            ("povfreq,V", value<int>(&mdremote.povfreq)->default_value(100000),
             "The frequency of making povray files")
            ("writedata,W", value<int>(&mdremote.writedata)->default_value(500),
             "frequency of dumping thermo time series (after 1K steps)");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    // Compute the box radius (where NP radius is unit length so it is not found here):
    if(counterionFlag != 'y') packing_fraction = 1.0;
    box_radius = pow(1 / packing_fraction,1.0/3.0);
    // Reducing the counterion diameter:
    counterion_diameter = (counterion_diameter / unit_radius_sphere);

    // Compute the 2D Young's Modulus (in kB*T/nm^2) from the reduced stretching constant specified as input:
    youngsModulus = boundary.sconstant / (unit_radius_sphere * unit_radius_sphere);  // For information purposes only.

    // The dimensionless form of the tension factors are a function of radius:
    double dynePerCm_Scalefactor = pow(10,-14)*(pow(unit_radius_sphere, 2)/unitenergy);
    boundary.sigma_a = (dynePerCm_Scalefactor * boundary.sigma_a);    // Coefficient is 1 (dyne/cm) in (kB T_room/(R0 nm^2)).
    //double dynePerUmFifth_Scalefactor = pow(10,-4)*pow(10,-18)*(pow(unit_radius_sphere, 6)/unitenergy);
    //boundary.sigma_v = (dynePerUmFifth_Scalefactor * boundary.sigma_v);
    double atmPerNmThird_Scalefactor = (1.01325 * pow(10,5)) * pow(10,-27) * (pow(unit_radius_sphere, 6) / (unitenergy * pow(10, -7)));
    boundary.sigma_v = atmPerNmThird_Scalefactor * boundary.sigma_v;

    // The scale factor for electrostatic interactions is a function of radius:
    const double scalefactor = epsilon_water * lB_water / unit_radius_sphere;

    disc1 = 3;                                          //  Discretization parameters (h, k).
    alpha = 1.0;                                        //  Fractional charge occupancy.
    if (chain_length == 1) Q_mesh = 0;                  //  Reduced mass of the thermostat particle(s)
    if (chain_length == 1) Q_ion = 0;                   //  Reduced mass of the thermostat particle(s)

    mdremote.QAnnealFac =
            1.00 + (mdremote.TAnnealFac / 10.0);        // Maintain previous scaling, same as before (fQ = 2) if fT = 10.

    // Flags for rigid geometric constraint and form:
    constraintForm = 'L';                               // Enforce a linear constraint (not quadratic).

    // ### Simulation setup: ###
    // Unused quantities for antiquated constraints & energy functionals:
    lambda_v = 0;                        // hardwiring volume constraint. Unused.
    lambda_a = 0;                        // hardwiring area constraint. Unused.
    boundary.lambda_l = 0;            // forgot what was this?!  NB:  I think it's an edge line tension force. Unused.

    // Electrostatics setup (membrane, implicit salt ions):
    z_out = 1;                        // Valency of implicit (salt) ions.
    boundary.ein = ein;               // Permittivity inside the membrane.
    boundary.eout = eout;             // Permittivity outside the membrane.  (Both equal by default.)
    boundary.set_up(unit_radius_sphere);
    // Assigning the Debye length outside (in reduced units):
    if (conc_out == 0)
        boundary.inv_kappa_out = 100000000;  // Effectively infinite, no screening.
    else
        boundary.inv_kappa_out =
                (0.257 / sqrt(z_out * boundary.lB_out * unit_radius_sphere * conc_out)) / unit_radius_sphere;

    // Generate the membrane, assign charge to vertices:
    boundary.discretize(disc1, disc2);            // discretize the interface
    if (disc1 != 0 || disc2 != 0) {
        if (externalPattern == "None") 
			  boundary.assign_random_q_values(q_strength, alpha, numPatches, fracChargedPatch, randomFlag, functionFlag);
        else 
			  boundary.assign_external_q_values(q_strength, externalPattern);
	 }

    //  Assess total actual charge on the membrane for record & to place the correct number of counterions:
    double q_actual = 0;
    for (unsigned int i = 0; i < boundary.V.size(); i++) {
        q_actual += boundary.V[i].q;
    }

    // Populate the box with counterions, if requested:
    if(counterionFlag == 'y') {
        boundary.put_counterions(q_actual, unit_radius_sphere, counterion_diameter, box_radius, counterions, counterion_valency);
        cout << "Counterions have been placed inside the box." << endl;
    }

    boundary.dressup(lambda_a, lambda_v);            // Compute initial normals, areas, and volumes.

    boundary.ref_area = boundary.total_area; // NB changed from (4 * pi), discrete membrane differs.
    boundary.ref_volume = boundary.total_volume;
    //boundary.ref_Area_Vertices = boundary.total_Area_Vertices;    // (2017.09.19 NB added.)  Initial area by vertices.

    // LJ parameters (mesh-mesh almost never invoked, present for security reasons):
    double lj_length_mesh_mesh_coeff = 0.1;
    boundary.lj_length_mesh_mesh = lj_length_mesh_mesh_coeff * boundary.avg_edge_length;
    double meshPointRadius = 0.5 * boundary.avg_edge_length;
    boundary.lj_length_mesh_ions = ((2.0 * meshPointRadius + counterion_diameter) / 2.0);
    boundary.elj = 1.0;

    // Thermostat for MD of the interface:
    vector<THERMOSTAT> mesh_bath, ions_bath;        // thermostat for real system (interface + ions)
    if (chain_length == 1) {
        mesh_bath.push_back(THERMOSTAT(0, T, 3 * boundary.V.size(), 0.0, 0, 0));
        ions_bath.push_back(THERMOSTAT(0, T, 3 * counterions.size(), 0.0, 0, 0));
    }
    else {
        mesh_bath.push_back(THERMOSTAT(Q_mesh, T, 3 * boundary.V.size(), 0, 0, 0));
        ions_bath.push_back(THERMOSTAT(Q_ion, T, 3 * counterions.size(), 0, 0, 0));
        while (mesh_bath.size() != chain_length - 1)
            mesh_bath.push_back(THERMOSTAT(Q_mesh, T, 1, 0, 0, 0));
        while (ions_bath.size() != chain_length - 1)
            ions_bath.push_back(THERMOSTAT(Q_ion, T, 1, 0, 0, 0));
        mesh_bath.push_back(THERMOSTAT(0, T, 3 * boundary.V.size(), 0.0, 0, 0));    // dummy bath (has 0 mass)
        ions_bath.push_back(THERMOSTAT(0, T, 3 * counterions.size(), 0.0, 0, 0));
    }

    // Assign masses to vertices
    for (unsigned int i = 0; i < boundary.V.size(); i++)
        boundary.V[i].m = 1.0;      // NB:  Energy conservation fails invariably if not = 1.0 .

    // NB added following lines to load off-file (stops if the file is not found; specify correctly):
    if (offFlag == 'y') boundary.load_configuration("380000.off");

    int numOfNodes = world.size();
    if (world.rank() == 0) {
#pragma omp parallel default(shared)
        {
            if (omp_get_thread_num() == 0) {
                printf("Number of MPI processes used %d\n", numOfNodes);
                printf("Number of OpenMP threads per MPI process %d\n", omp_get_num_threads());
            }
        }
    }

    if (world.rank() == 0) {
        // Define an output stream to provide the same things as the (*.out) file on the HPC cluster (minus run outputs):
        // Useful in automated analysis when console output is not saved by default on local PC.
        ofstream list_out("outfiles/simulation.out", ios::app);

        // ### Output to the console (or HPC *.out file) & backup "simulation.out" file) the simulation parameters: ###

        cout << "\n===================\nMembrane\n===================\n";
        cout << "Unit radius of the sphere (in nm): " << unit_radius_sphere << endl;
        cout << "Number of vertices: " << boundary.V.size() << endl;
        cout << "Number of edges: " << boundary.E.size() << endl;
        cout << "Number of faces: " << boundary.F.size() << endl;
        cout << "Lambda_v: " << lambda_v << endl;
        cout << "Lambda_a: " << lambda_a << endl;
        cout << "Lambda_l: " << boundary.lambda_l << endl;
        cout << "Total intial volume: " << boundary.total_volume << "  Sphere volume: " << (4. / 3.) * pi
             << "  Ref volume: " << boundary.ref_volume << endl;
        cout << "Total intial face area: " << boundary.total_area << "  Sphere area: " << 4 * pi << "  Ref area: "
             << boundary.ref_area << endl;
        cout << "Mesh point radius: " << meshPointRadius << endl;
        cout << "Bending rigidity: " << boundary.bkappa << endl;
        cout << "Young's Modulus (2D): " << youngsModulus << endl;
        cout << "Stretching constant: " << boundary.sconstant << endl;
        cout << "Surface Tension constant: " << (boundary.sigma_a / dynePerCm_Scalefactor) << endl;
        cout << "Volume Tension constant: " << (boundary.sigma_v / atmPerNmThird_Scalefactor)<< endl;
        cout << "Unstretched edge length: " << boundary.avg_edge_length
             << endl; // NB uncommented & replaced the output value.
        cout << "LJ strength: " << boundary.elj << endl;
        cout << "LJ mesh-mesh distance cutoff: " << boundary.lj_length_mesh_mesh << endl;
        cout << "LJ mesh-ion distance cutoff: " << boundary.lj_length_mesh_ions << endl;
        cout << "Rigid geometric constraint: " << geomConstraint << endl;

        cout << "\n===================\nElectrostatics\n===================\n";
        cout << "Total charge on membrane (if homogeneously charged): " << q_strength << endl;
        cout << "Final actual charge on membrane: " << q_actual << endl;
        cout << "Number of simulated counterions: " << counterions.size() << endl;
        cout << "Fractional charge occupancy (alpha): " << alpha << endl;
        cout << "The number of patches is: " << numPatches << endl;
        cout << "The size of each patch is: " << fracChargedPatch << endl;
        cout << "The name of the external file used to specify charge patterning is: " << externalPattern << endl;
        cout << "Salt concentration: " << conc_out << endl;
        cout << "Debye length: " << boundary.inv_kappa_out << endl;
        cout << "Bjerrum length: " << boundary.lB_out << endl;

        cout << "\n===================\nThermodynamics\n===================\n";
        cout << "Initial Temperature (T): " << T << endl;
        cout << "Thermostat mass (Q_mesh): " << mesh_bath[0].Q << endl;
        cout << "Number of thermostats (C): " << mesh_bath.size() << endl;
        cout << "Timestep (dt): " << mdremote.timestep << "  &  Total steps: " << mdremote.steps << "." << endl;
        cout << "Annealing: " << mdremote.anneal << endl;
        cout << "Anneal frequency: " << mdremote.annealfreq << endl;
        cout << "Anneal duration: " << mdremote.annealDuration << endl;
        cout << "Temperature annealing decrement: " << mdremote.TAnnealFac << endl;
        cout << "Thermostat mass decrement: " << mdremote.QAnnealFac << endl;
        cout << "Randomness during run: " << randomFlag << endl;
        cout << "Off-file loaded: " << offFlag << endl;

        // Same thing in a file even if you're not on an HPC that can store the console output automatically on request:
        list_out << "===================\nMembrane\n===================\n";
        list_out << "Unit radius of the sphere (in nm): " << unit_radius_sphere
                 << endl; //NB added in observing size change.
        list_out << "Number of vertices: " << boundary.V.size() << endl;
        list_out << "Number of edges: " << boundary.E.size() << endl;
        list_out << "Number of faces: " << boundary.F.size() << endl;
        list_out << "Lambda_v: " << lambda_v << endl;
        list_out << "Lambda_a: " << lambda_a << endl;
        list_out << "Lambda_l: " << boundary.lambda_l << endl;
        list_out << "Total intial volume: " << boundary.total_volume << "  Sphere volume: " << (4. / 3.) * pi
                 << "  Ref volume: " << boundary.ref_volume << endl;
        list_out << "Total intial face area: " << boundary.total_area << "  Sphere area: " << 4 * pi << "  Ref area: "
                 << boundary.ref_area << endl;
        list_out << "Mesh point radius: " << meshPointRadius << endl;
        list_out << "Bending rigidity: " << boundary.bkappa << endl;
        list_out << "Young's Modulus (2D): " << youngsModulus << endl;
        list_out << "Stretching constant: " << boundary.sconstant << endl;
        list_out << "Surface Tension constant: " << (boundary.sigma_a / dynePerCm_Scalefactor) << endl;
        list_out << "Volume Tension constant: " << (boundary.sigma_v / atmPerNmThird_Scalefactor) << endl;
        list_out << "LJ strength: " << boundary.elj << endl;
        list_out << "LJ distance cutoff: " << boundary.lj_length_mesh_mesh << endl;
        list_out << "Rigid geometric constraint: " << geomConstraint << endl;

        list_out << "\n===================\nElectrostatics\n===================\n";
        list_out << "Total charge on membrane: " << q_strength << endl;
        list_out << "Final actual charge on membrane: " << q_actual << endl;
        list_out << "Number of simulated counterions: " << counterions.size() << endl;
        list_out << "Fractional charge occupancy (alpha): " << alpha << endl;
        list_out << "The number of patches is: " << numPatches << endl;
        list_out << "The size of each patch is: " << fracChargedPatch << endl;
        list_out << "The name of the external file used to specify charge patterning is: " << externalPattern << endl;
        list_out << "Salt concentration: " << conc_out << endl;
        list_out << "Debye length: " << boundary.inv_kappa_out << endl;
        list_out << "Bjerrum length: " << boundary.lB_out << endl;

        list_out << "\n===================\nThermodynamics\n===================\n";
        list_out << "Initial Temperature (T): " << T << endl;
        list_out << "Thermostat mass (Q_mesh): " << mesh_bath[0].Q << endl;
        list_out << "Number of thermostats (C): " << mesh_bath.size() << endl;
        list_out << "Timestep (dt): " << mdremote.timestep << "  &  Total steps: " << mdremote.steps << "." << endl;
        list_out << "Annealing: " << mdremote.anneal << endl;
        list_out << "Anneal frequency: " << mdremote.annealfreq << endl;
        list_out << "Anneal duration: " << mdremote.annealDuration << endl;
        list_out << "Temperature annealing decrement: " << mdremote.TAnnealFac << endl;
        list_out << "Thermostat mass decrement: " << mdremote.QAnnealFac << endl;
        list_out << "Randomness during run: " << randomFlag << endl;
        list_out << "Off-file loaded: " << offFlag << endl;

        // (NB added) Print the initial state in the movie file, prior to MD:
        interface_movie(0, boundary.V, counterions, box_radius);  // NB added.  Note, calls this "time step -1" in movie file.

    }
    // Initiate MD of the boundary/membrane:
    boundary.dressup(lambda_a, lambda_v);
    boundary.output_configuration();    // NB added.  More information on the initial state.
    /*boundary.compute_local_energies(scalefactor);  // Returns 'nan' at the moment as energies... should fix.
    boundary.compute_local_energies_by_component();
    if (world.rank() == 0) {
        rename("outfiles/local_electrostatic_E.off", "outfiles/local_electrostatic_E_initial.off");
        rename("outfiles/local_elastic_E.off", "outfiles/local_elastic_E_initial.off");
        rename("outfiles/local_stretching_E.off", "outfiles/local_stretching_E_initial.off");
        rename("outfiles/local_bending_E.off", "outfiles/local_bending_E_initial.off");
    }*/
    //MPI Boundary calculation for ions
    /*
    unsigned int rangeMesh = boundary.V.size() / world.size() + 1.5;
    lowerBoundMesh = world.rank() * rangeMesh;
    upperBoundMesh = (world.rank() + 1) * rangeMesh - 1;
    extraElementsMesh = world.size() * rangeMesh - boundary.V.size();
    sizFVecMesh = upperBoundMesh - lowerBoundMesh + 1;
    if (world.rank() == world.size() - 1) {
        upperBoundMesh = boundary.V.size() - 1;
        sizFVecMesh = upperBoundMesh - lowerBoundMesh + 1 + extraElementsMesh;
    }
    if (world.size() == 1) {
        lowerBoundMesh = 0;
        upperBoundMesh = boundary.V.size() - 1;
    }
*/
    // MPI Boundary calculations for the mesh:
    unsigned int rangeMesh = (boundary.V.size() + world.size() - 1) / (1.0 * world.size());
    lowerBoundMesh = world.rank() * rangeMesh;
    upperBoundMesh = (world.rank() + 1) * rangeMesh - 1;
    //extraElementsMesh = world.size() * rangeMesh - boundary.V.size();
    sizFVecMesh = rangeMesh;
    if (world.rank() == world.size() - 1) {
        upperBoundMesh = boundary.V.size() - 1;
    }

    // MPI Boundary calculations for the ions:
    int rangeIons = (counterions.size() + world.size() - 1) / (1.0 * world.size());
    lowerBoundIons = world.rank() * rangeIons;
    upperBoundIons = (world.rank() + 1) * rangeIons - 1;
    extraElementsIons = world.size() * rangeIons - counterions.size();
    sizFVecIons = rangeIons;
    if (world.rank() == world.size() - 1) {
        upperBoundIons = counterions.size() - 1; // With zero counterions, this is negative, preventing loop evaluation.
    }

    md_interface(boundary, counterions, mesh_bath, ions_bath, mdremote, geomConstraint, bucklingFlag, constraintForm, scalefactor, box_radius);
    boundary.compute_local_energies(scalefactor);
    boundary.compute_local_energies_by_component();

    return 0;
}
// End of main

