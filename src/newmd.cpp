#include "utility.h"
#include "interface.h"
#include "thermostat.h"
#include "particle.h"
#include "functions.h"
#include "newforces.h"
#include "newenergies.h"

void md_interface(INTERFACE &boundary, vector<PARTICLE> &counterions, vector<THERMOSTAT> &mesh_bath, vector<THERMOSTAT> &ions_bath, CONTROL &cpmdremote, char geomConstraint, char bucklingFlag,
                  char constraintForm, const double scalefactor, double box_radius) {

    double percentage = 0, percentagePre = -1;

    // Initialize the velocities to either the initial temperature or zero (choose which to comment out):
    //initialize_vertex_velocities(boundary.V, mesh_bath);
    initialize_velocities_to_zero(boundary.V, counterions);

    // ### Calculate the intial net force on all vertices: ###
    force_calculation_init(boundary, counterions, scalefactor, bucklingFlag);

    double expfac_mesh, expfac_ions; // Factors needed to update velocities for thermostatting.

    // Define output streams for data files:
    ofstream list_energy("outfiles/energy_nanomembrane.dat", ios::app);
    //ofstream list_variables ("outfiles/variables.dat", ios::app);
    ofstream list_area("outfiles/area.dat", ios::app);
    ofstream list_volume("outfiles/volume.dat", ios::app);
    ofstream list_temperature("outfiles/temperature.dat", ios::app);
    ofstream list_momentum("outfiles/total_momentum.dat", ios::app);
    ofstream list_netEnergyDrift("outfiles/net_Energy_Drift.dat",ios::app);

    // SHAKE & RATTLE prior to MD:
    if (geomConstraint == 'V') {
        if (world.rank() == 0) {
            if (constraintForm == 'Q') cout << "Constraint:  net volume, quadratically." << endl << endl;
            else cout << "Constraint:  net volume, linearly." << endl << endl;
        }
        SHAKE(boundary, cpmdremote, constraintForm);
        RATTLE(boundary, constraintForm);
    } else if (geomConstraint == 'A') {
        if (world.rank() == 0) {
            if (constraintForm == 'Q') cout << "Constraint:  net area, quadratically." << endl << endl;
            else cout << "Constraint:  net area, linearly." << endl << endl;
        }
        SHAKE_for_area(boundary, cpmdremote, constraintForm);
        RATTLE_for_area(boundary, constraintForm);
    }
    else if (world.rank() == 0)
        cout << endl << "No rigid geometric constraint will be enforced, soft constraints only." << endl << endl;

    // NB added to provide the initial quantities (pre-MD):
    long double vertex_ke = compute_kinetic_energy(boundary.V); // Store KE of vertices
    long double ions_ke = compute_kinetic_energy(counterions); // Store KE of ions
    boundary.compute_energy(0, counterions, scalefactor, bucklingFlag);
    double mesh_bath_ke = bath_kinetic_energy(mesh_bath);
    double mesh_bath_pe = bath_potential_energy(mesh_bath);
    double ions_bath_ke = bath_kinetic_energy(ions_bath);
    double ions_bath_pe = bath_potential_energy(ions_bath);
    VECTOR3D total_momentum = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < boundary.V.size(); i++)
        total_momentum += boundary.V[i].velvec;
    double initNetEnergy = (boundary.energy + bath_kinetic_energy(mesh_bath) + bath_potential_energy(mesh_bath));
    double globalNetEnergy = boundary.energy + mesh_bath_ke + mesh_bath_pe + ions_bath_ke + ions_bath_pe;
    if (world.rank() == 0) {
        list_energy << 0
        << setw(15) << boundary.netMeshKE + boundary.netIonKE // Net kinetic energy (mesh + ions)
        << setw(15) << boundary.penergy // Net potential energy (penergy of all)
        << setw(15) << boundary.energy  // Local net energy (netMeshKE + netIonKE + penergy of all)
        << setw(15) << globalNetEnergy // Global (conserved) net energy.
        << setw(15) << mesh_bath_ke + ions_bath_ke // Net bath KE
        << setw(15) << mesh_bath_pe + ions_bath_pe << endl; // Net bath PE

        // Compute and output the face-based net area & volumes:
        list_area << 0 << setw(15) << boundary.total_area << endl;
        list_volume << 0 << setw(15) << boundary.total_volume << endl;

        // Output the temperature (and optionally, thermostat parameters):
        list_temperature << 0 << setw(15) << 2 * boundary.netMeshKE / (mesh_bath[0].dof * kB) << endl;
    }

    // For annealing, compute per-step fractional decrement in {T, Q} required for desired net decrement per procedure.
    /* VJ suggested not changing annealing procedure from abrupt, but simply changing magnitude decremented per call.
     * In this case, just leave annealDuration = 1 in main. */
    long double fT = pow((1.0 / cpmdremote.TAnnealFac), (1.0 /
                                                         cpmdremote.annealDuration));  // Fractional decrement that, when f^(annealDuration), yields (T -> T/10).
    long double fQ = pow((1.0 / cpmdremote.QAnnealFac), (1.0 /
                                                         cpmdremote.annealDuration));  // Fractional decrement that, when f^(annealDuration), yields (Q -> Q/2)
    if (world.rank() == 0)
        cout << "Initial temperature annealing multiplier: " << fT << " Initial thermostat mass annealing multiplier: "
             << fQ << endl << endl;

    // A counter that will be incremented each write step after the system has drifted, to abort after providing info:
    unsigned int abortCounter = 0;

    // Run the MD:
    for (int num = 1; num <= cpmdremote.steps; num++) {
        // Reverse update of Nose-Hoover chain
        // xi
        for (int j = mesh_bath.size() - 1; j > -1; j--)
            update_chain_xi(j, mesh_bath, cpmdremote.timestep, vertex_ke);
        for (int j = ions_bath.size() - 1; j > -1; j--)
            update_chain_xi(j, ions_bath, cpmdremote.timestep, ions_ke);
        // eta
        for (unsigned int j = 0; j < mesh_bath.size(); j++)
            mesh_bath[j].update_eta(cpmdremote.timestep);
        for (unsigned int j = 0; j < ions_bath.size(); j++)
            ions_bath[j].update_eta(cpmdremote.timestep);
        // factor needed to update velocity
        expfac_mesh = exp(-0.5 * cpmdremote.timestep * mesh_bath[0].xi);
        expfac_ions = exp(-0.5 * cpmdremote.timestep * ions_bath[0].xi);

        // ### Velocity-Verlet Section ###
        // Propagate velocity (half time step):
        for (unsigned int i = 0; i < boundary.V.size(); i++)
            boundary.V[i].update_real_velocity(cpmdremote.timestep, mesh_bath[0], expfac_mesh);
        for (unsigned int i = 0; i < counterions.size(); i++)
            counterions[i].update_real_velocity(cpmdremote.timestep, ions_bath[0], expfac_ions);

        // Propagate position (full time step):
        for (unsigned int i = 0; i < boundary.V.size(); i++)
            boundary.V[i].update_real_position(cpmdremote.timestep);
        for (unsigned int i = 0; i < counterions.size(); i++)
            counterions[i].update_real_position(cpmdremote.timestep);

        // SHAKE to ensure constraint is satisfied for newly updated positions:
        if (geomConstraint == 'V') SHAKE(boundary, cpmdremote, constraintForm);
        else if (geomConstraint == 'A') SHAKE_for_area(boundary, cpmdremote, constraintForm);

        // After SHAKE, dress up interface again as vertex positions have changed:
        boundary.dressup(boundary.lambda_a, boundary.lambda_v);

        /*	// uncomment if flips desired
        for (int i=0; i<boundary.number_of_edges; i++)
          boundary.E[i].flipIfFavorable();
        for (int i=0; i<boundary.number_of_edges; i++)
          boundary.E[i].flipIfFavorable();
        for (int i=0; i<boundary.number_of_edges; i++)
          boundary.E[i].flipIfFavorable();
        boundary.dressup(boundary.lambda_a, boundary.lambda_v);
        */

        // ### Compute Forces ###
        force_calculation(boundary, counterions, scalefactor, bucklingFlag);

        // Propagate velocity (second half time step):
        for (unsigned int i = 0; i < boundary.V.size(); i++)
            boundary.V[i].update_real_velocity(cpmdremote.timestep, mesh_bath[0], expfac_mesh);
        for (unsigned int i = 0; i < counterions.size(); i++)
            counterions[i].update_real_velocity(cpmdremote.timestep, mesh_bath[0], expfac_ions);

        // RATTLE the system to enforce time-derivative of constraint (gradients update not required, positions unchanged):
        if (geomConstraint == 'V') RATTLE(boundary, constraintForm);
        else if (geomConstraint == 'A') RATTLE_for_area(boundary, constraintForm);

        // Compute KEs needed to set canonical ensemble:
        vertex_ke = compute_kinetic_energy(boundary.V);
        ions_ke = compute_kinetic_energy(counterions);

        // Invoke the thermostat (forward Nose-Hoover chain):
        // xi
        for (unsigned int j = 0; j < mesh_bath.size(); j++)
            mesh_bath[j].update_eta(cpmdremote.timestep);
        for (unsigned int j = 0; j < ions_bath.size(); j++)
            ions_bath[j].update_eta(cpmdremote.timestep);

        // eta
        for (unsigned int j = 0; j < mesh_bath.size(); j++)
            update_chain_xi(j, mesh_bath, cpmdremote.timestep, vertex_ke);
        for (unsigned int j = 0; j < ions_bath.size(); j++)
            update_chain_xi(j, ions_bath, cpmdremote.timestep, ions_ke);

        // ### Output file updates: ###
        // Output the series data at the specified interval:
        if (num % cpmdremote.writedata == 0 || num < 3) {
            // Compute & output the membrane-wide and global energies (energy_nanomembrane quantities):

            boundary.compute_energy(num, counterions, scalefactor, bucklingFlag);
            double mesh_bath_ke = bath_kinetic_energy(mesh_bath);
            double mesh_bath_pe = bath_potential_energy(mesh_bath);
            double ions_bath_ke = bath_kinetic_energy(ions_bath);
            double ions_bath_pe = bath_potential_energy(ions_bath);
            // Compute the global net energy to check for conservation:
            globalNetEnergy = boundary.energy + mesh_bath_ke + mesh_bath_pe + ions_bath_ke + ions_bath_pe;

            if (world.rank() == 0) {
                list_energy << num
                << setw(15) << boundary.netMeshKE + boundary.netIonKE // Net kinetic energy (mesh + ions)
                << setw(15) << boundary.penergy // Net potential energy (of all)
                << setw(15) << boundary.energy  // Local net energy (netMeshKE + netIonKE + penergy of all)
                << setw(15) << globalNetEnergy // Global (conserved) net energy.
                << setw(15) << mesh_bath_ke + ions_bath_ke // Mesh bath KE
                << setw(15) << mesh_bath_pe + ions_bath_pe << endl; // Mesh bath PE

                // Dump the energy drift and components explicitly (in file "net_Energy_Drift.dat"):
                list_netEnergyDrift << num << "\t" << (globalNetEnergy / initNetEnergy) << "\t" << globalNetEnergy << "\t" << initNetEnergy << "\t" << abortCounter << endl;

                // Compute and output the face-based net area & volumes:
                list_area << num << setw(15) << boundary.total_area << endl;
                list_volume << num << setw(15) << boundary.total_volume << endl;

                // Output the temperature (and optionally, thermostat parameters):
                list_temperature << num << setw(15) << 2 * boundary.netMeshKE / (mesh_bath[0].dof * kB) << endl;
                //list_variables << num << setw(15) << mesh_bath[0].xi << setw(15) << mesh_bath[0].eta << endl; // NOTE: commented out as not that important

                //  Abort the entire program if the global net energy has drifted upward by more than 5%:
                if (1.05 < (globalNetEnergy / initNetEnergy)) {
                    abortCounter++;
                    if (abortCounter == 10) {
                        cout << "Aborting due to excessive (>5%) net energy drift upwards." << endl;
                        abort();
                    }
                }
                //  Also abort if annealing hasn't yet begun, but the global energy has drifted downward by more than 5%:
                else if ((num < cpmdremote.annealfreq) && ((globalNetEnergy / initNetEnergy) < 0.95)) {
                    abortCounter++;
                    if (abortCounter == 10) {
                        cout << "Aborting due to excessive (>5%) net energy drift downwards prior to annealing." << endl;
                        abort();
                    }
                }
            }
            // Compute & output the net momentum (informative only, from velocities):
            total_momentum = VECTOR3D(0, 0, 0);
            for (unsigned int i = 0; i < boundary.V.size(); i++)
                total_momentum += boundary.V[i].velvec;
            if (world.rank() == 0)
                list_momentum << num << setw(15) << total_momentum.GetMagnitude() << endl;
        }

        // Output the dump (movie) file data at the specified interval:
        if (num % cpmdremote.moviefreq == 0) {
            for (unsigned int i = 0; i < boundary.F.size(); i++)
                boundary.F[i].compute_area_normal();
            for (unsigned int i = 0; i < boundary.V.size(); i++)
                boundary.V[i].compute_area_normal();
            if (num % cpmdremote.offfreq == 0)
                interface_off(num, boundary);
            if (world.rank() == 0)
                interface_movie(num, boundary.V, counterions, box_radius);
        }

        /*							// NOTE COMMENTING OUT THE POVRAY BLOCK AS IT IS NOT NECESSARY FOR INVESTIGATIONS
        if (num % cpmdremote.povfreq == 0)
        {
          for(unsigned int i=0; i<boundary.F.size(); i++)
            boundary.F[i].compute_area_normal();
          for(unsigned int i=0; i<boundary.V.size(); i++)
            boundary.V[i].compute_area_normal();
          interface_pov_smooth(num,boundary);
          interface_pov(num,boundary);
        }
        */

        // Conduct simulated annealing of the MD if above the annealing threshold (for the specified duration):
        /*  NB added a gradual annealing procedure that begins after the standard annealing interval ("freq").
              This proceeds for a specified number of steps "annealDuration".  Tested in Mathematica and here.
                VJ cautions against using gradual annealing, in which case just leave annealDuration = 1 in main. */
        if (cpmdremote.anneal == 'y' && num > cpmdremote.annealDuration &&
            ((num % cpmdremote.annealfreq == 0) || (num % cpmdremote.annealfreq < cpmdremote.annealDuration)) &&
            mesh_bath[0].T >= 0.000000000000001) {
            //  If annealing is commencing, update the annealing fractional decrements and print relevant information:
            if (cpmdremote.anneal == 'y' && num % cpmdremote.annealfreq == 0 &&
                mesh_bath[0].T >= 0.000000000000001) { // Decide which fractional reduction for {T, Q} to use:
                if (num == cpmdremote.annealfreq) // If in the first stage, decrement temperature by 1.25x.
                {
                    cpmdremote.TAnnealFac = 1.25;
                    cpmdremote.QAnnealFac = 1.0 + (cpmdremote.TAnnealFac / 10.0);
                } else if (num <= (3 * cpmdremote.annealfreq))  // If in the 2-3rd stages, decrement temperature by 2.0x.
                {
                    cpmdremote.TAnnealFac = 2;
                    cpmdremote.QAnnealFac = 1.0 + (cpmdremote.TAnnealFac / 10.0);
                } else if (num <= (6 * cpmdremote.annealfreq))  // If in the 4th-5th stages, decrement temperature by 4.0x.
                {
                    cpmdremote.TAnnealFac = 4;
                    cpmdremote.QAnnealFac = 1.0 + (cpmdremote.TAnnealFac / 10.0);
                }
                else if (num <= (7 * cpmdremote.annealfreq)) // If in 8th, decrement temperature by 8.0x.
                {
                    cpmdremote.TAnnealFac = 8;
                    cpmdremote.QAnnealFac = 1.0 + (cpmdremote.TAnnealFac / 10.0);
                } else  // Any further annealing calls will decrement the temperature 10.0x.
                {
                    cpmdremote.TAnnealFac = 10;
                    cpmdremote.QAnnealFac = 1.0 + (cpmdremote.TAnnealFac / 10.0);
                }

                // Fractional decrement per-step required to achieve net per-call decrement of (XAnnealFac):
                fT = pow((1.0 / cpmdremote.TAnnealFac),
                         (1.0 / cpmdremote.annealDuration));  // (fT != TAnnealFac if annealDuration != 1)
                fQ = pow((1.0 / cpmdremote.QAnnealFac), (1.0 / cpmdremote.annealDuration));
                if (world.rank() == 0)
                    cout << "Beginning annealing for " << cpmdremote.annealDuration
                         << " steps decrementing by a factor of fT = " << fT << " and fQ = " << fQ << "." << endl;
            }
            for (unsigned int b = 0; b < mesh_bath.size(); b++) {
                mesh_bath[b].T = fT * mesh_bath[b].T;
                mesh_bath[b].Q = fQ * mesh_bath[b].Q;
            }
        }

        // VJ's instantaneous annealing for constant decrement (fT, fQ):
        /*    if (cpmdremote.anneal == 'y' && num % cpmdremote.annealfreq == 0 && mesh_bath[0].T >= 0.000000000000001)
        {
          cout << "annealing ... " << endl;
          for (unsigned int b = 0; b < mesh_bath.size(); b++)
          {
            mesh_bath[b].T = mesh_bath[b].T / 10;
            mesh_bath[b].Q = mesh_bath[b].Q / 2;
          }
        }*/

        // Percentage calculation (for progress bar):
        if (world.rank() == 0) {
            percentage = roundf(num / (double) cpmdremote.steps * 100 * 10) / 10;
            //percentage output
            if (percentage != percentagePre) {
                double fraction_completed = percentage / 100;
                progressBar(fraction_completed);
                percentagePre = percentage;
            }
        }
    }
}

/* double sum trick
     for (unsigned int i = 0; i < boundary.V.size(); i++)
     {
       for (unsigned int j = 0; j < i; j++)
       {
 	VECTOR3D ljforce = lj_vertex_vertex(boundary.V[i], boundary.V[j], boundary.avg_edge_length);
 	boundary.V[i].forvec += ljforce;
 	boundary.V[j].forvec -= ljforce;
       }
     }*/