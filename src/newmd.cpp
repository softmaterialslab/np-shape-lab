#include "utility.h"
#include "interface.h"
#include "thermostat.h"
#include "functions.h"
#include "newforces.h"
#include "newenergies.h"

void md_interface(INTERFACE &boundary, vector<THERMOSTAT> &real_bath, CONTROL &cpmdremote, char geomConstraint,
                  char constraintForm, const double scalefactor) {

    double percentage = 0, percentagePre = -1;

    // Initialize the velocities to either the initial temperature or zero (choose which to comment out):
    //initialize_vertex_velocities(boundary.V, real_bath);
    initialize_vertex_velocities_to_zero(boundary.V);

    // ### Calculate the intial net force on all vertices: ###
    force_calculation_init(boundary, scalefactor);

    double expfac_real;

    // Define output streams for data files:
    ofstream list_energy("outfiles/energy_nanomembrane.dat", ios::app);
    //ofstream list_variables ("outfiles/variables.dat", ios::app);
    ofstream list_area("outfiles/area.dat", ios::app);
    ofstream list_volume("outfiles/volume.dat", ios::app);
    ofstream list_temperature("outfiles/temperature.dat", ios::app);
    ofstream list_momentum("outfiles/total_momentum.dat", ios::app);

    // (NB added.) Compute each & save pre-constraint, pre-MD {momentum, temp, net E}:
    vertex_kinetic_energy(boundary.V);
    boundary.compute_energy(0, scalefactor);
    VECTOR3D total_momentum = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < boundary.V.size(); i++)
        total_momentum += boundary.V[i].velvec;
    // Define the original temperature pre-MD SHAKE & RATTLE:
    double preConTemp = 2 * boundary.kenergy / (real_bath[0].dof * kB);
    double preConNetEnergy = (boundary.energy + bath_kinetic_energy(real_bath) + bath_potential_energy(real_bath));

    // SHAKE & RATTLE prior to MD:  (NB moved to below energy initialization to observe any temperature shift from before & after constraint function calls.)
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

    // Recompute the same quantities as above, assign post-constraint variables for comparison:
    long double vertex_ke = vertex_kinetic_energy(boundary.V);

    boundary.compute_energy(0, scalefactor);
    total_momentum = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < boundary.V.size(); i++)
        total_momentum += boundary.V[i].velvec;
    double initNetEnergy = (boundary.energy + bath_kinetic_energy(real_bath) + bath_potential_energy(real_bath));

    if (world.rank() == 0) {
        // Report the initial temperature and any shift due to constraints, for troubleshooting:
        cout << "Initial Temperature (prior to MD): " << 2 * boundary.kenergy / (real_bath[0].dof * kB)
             << "\t  &  constraint-induced shift: " << (2 * boundary.kenergy / (real_bath[0].dof * kB) - preConTemp)
             << "."
             << endl;
        cout << "Initial Net Energy (prior to MD):  " << initNetEnergy << "\t\t  &  constraint-induced shift: "
             << (initNetEnergy - preConNetEnergy) << "." << endl;
    }
    // NB added to provide the initial quantities (pre-MD) [those not dumped had scope issues]:
    //list_area << 0 << setw(15) << boundary.total_area << endl;
    //list_volume << 0 << setw(15) << boundary.total_volume << endl;
    //list_temperature << 0 << setw(15) << 2*boundary.kenergy/(real_bath[0].dof*kB) << endl;
    //list_momentum << 0 << setw(15) << total_momentum << endl;

    // For annealing, compute per-step fractional decrement in {T, Q} required for desired net decrement per procedure.
    /* VJ suggested not changing annealing procedure from abrupt, but simply changing magnitude decremented per call.
     * In this case, just leave annealDuration = 1 in main. */
    long double fT = pow((1.0 / cpmdremote.TAnnealFac), (1.0 /
                                                         cpmdremote.annealDuration));  // Fractional decrement that, when f^(annealDuration), yields (T -> T/10).
    long double fQ = pow((1.0 / cpmdremote.QAnnealFac), (1.0 /
                                                         cpmdremote.annealDuration));  // Fractional decrement that, when f^(annealDuration), yields (Q -> Q/2)
    if (world.rank() == 0)
        cout << "Initial temperature annealing multiplier: " << fT << " Initial thermostat mass anenealing multiplier: "
             << fQ << endl << endl;

    // Run the MD:
    for (int num = 1; num <= cpmdremote.steps; num++) {
        // Reverse update of Nose-Hoover chain
        // xi
        for (int j = real_bath.size() - 1; j > -1; j--)
            update_chain_xi(j, real_bath, cpmdremote.timestep, vertex_ke);
        // eta
        for (unsigned int j = 0; j < real_bath.size(); j++)
            real_bath[j].update_eta(cpmdremote.timestep);
        // factor needed to update velocity
        expfac_real = exp(-0.5 * cpmdremote.timestep * real_bath[0].xi);

        // ### Velocity-Verlet Section ###
        // Propagate velocity (half time step):
        for (unsigned int i = 0; i < boundary.V.size(); i++)
            boundary.V[i].update_real_velocity(cpmdremote.timestep, real_bath[0],
                                               expfac_real);    // update particle velocity half time step

        // Propagate position (full time step):
        for (unsigned int i = 0; i < boundary.V.size(); i++)
            boundary.V[i].update_real_position(cpmdremote.timestep);

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
        force_calculation(boundary, scalefactor);

        // Propagate velocity (second half time step):
        for (unsigned int i = 0; i < boundary.V.size(); i++)
            boundary.V[i].update_real_velocity(cpmdremote.timestep, real_bath[0], expfac_real);

        // RATTLE the system to enforce time-derivative of constraint (gradients update not required, positions unchanged):
        if (geomConstraint == 'V') RATTLE(boundary, constraintForm);
        else if (geomConstraint == 'A') RATTLE_for_area(boundary, constraintForm);

        // kinetic energies needed to set canonical ensemble
        vertex_ke = vertex_kinetic_energy(boundary.V);

        // Invoke the thermostat (forward Nose-Hoover chain):
        // xi
        for (unsigned int j = 0; j < real_bath.size(); j++)
            real_bath[j].update_eta(cpmdremote.timestep);
        // eta
        for (unsigned int j = 0; j < real_bath.size(); j++)
            update_chain_xi(j, real_bath, cpmdremote.timestep, vertex_ke);

        // ### Output file updates: ###
        // Output the series data at the specified interval:
        if (num % cpmdremote.writedata == 0 || num < 1000) {
            // Compute & output the membrane-wide and global energies (energy_nanomembrane quantities):

            boundary.compute_energy(num, scalefactor);
            double real_bath_ke = bath_kinetic_energy(real_bath);
            double real_bath_pe = bath_potential_energy(real_bath);

            if (world.rank() == 0) {
                list_energy << num << setw(15) << boundary.kenergy << setw(15) << boundary.penergy << setw(15)
                            << boundary.energy << setw(15) << boundary.energy
                                                              + real_bath_ke + real_bath_pe << setw(15) << real_bath_ke
                            << setw(15) << real_bath_pe << endl;

                // Print the initial and current net energy, along with the drift (their ratio):
/*                cout << "(" << num << ")" << "\tInitial Net Energy: " << initNetEnergy << "\t & \t Current Net Energy: "
                     << (boundary.energy + bath_kinetic_energy(real_bath) + bath_potential_energy(real_bath))
                     << "\t & \t Drift: "
                     << ((boundary.energy + bath_kinetic_energy(real_bath) + bath_potential_energy(real_bath)) /
                         initNetEnergy) << "." << endl;*/

                // Compute and output the face-based net area & volumes:
                list_area << num << setw(15) << boundary.total_area << endl;
                list_volume << num << setw(15) << boundary.total_volume << endl;

                // Output the temperature (and optionally, thermostat parameters):
                list_temperature << num << setw(15) << 2 * boundary.kenergy / (real_bath[0].dof * kB) << endl;
                //list_variables << num << setw(15) << real_bath[0].xi << setw(15) << real_bath[0].eta << endl; // NOTE: commented out as not that important
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
                boundary.F[i].computenormal();
            for (unsigned int i = 0; i < boundary.V.size(); i++)
                boundary.V[i].computenormal();
            if (num % cpmdremote.offfreq == 0)
                interface_off(num, boundary);
            if (world.rank() == 0)
                interface_movie(num, boundary.V, boundary);
        }

        /*							// NOTE COMMENTING OUT THE POVRAY BLOCK AS IT IS NOT NECESSARY FOR INVESTIGATIONS
        if (num % cpmdremote.povfreq == 0)
        {
          for(unsigned int i=0; i<boundary.F.size(); i++)
            boundary.F[i].computenormal();
          for(unsigned int i=0; i<boundary.V.size(); i++)
            boundary.V[i].computenormal();
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
            real_bath[0].T >= 0.000000000000001) {
            //  If annealing is commencing, update the annealing fractional decrements and print relevant information:
            if (cpmdremote.anneal == 'y' && num % cpmdremote.annealfreq == 0 &&
                real_bath[0].T >= 0.000000000000001) { // Decide which fractional reduction for {T, Q} to use:
                if (num == cpmdremote.annealfreq) // If in the first stage, decrement temperature by 1.25x.
                {
                    cpmdremote.TAnnealFac = 1.25;
                    cpmdremote.QAnnealFac = 1.0 + (cpmdremote.TAnnealFac / 10.0);
                } else if (num <=
                           (3 * cpmdremote.annealfreq))  // If in the 2-4th stages, decrement temperature by 2.0x.
                {
                    cpmdremote.TAnnealFac = 2;
                    cpmdremote.QAnnealFac = 1.0 + (cpmdremote.TAnnealFac / 10.0);
                } else if (num < 6 * cpmdremote.annealfreq)  // If in the 5-7th stages, decrement temperature by 5.0x.
                {
                    cpmdremote.TAnnealFac = 5;
                    cpmdremote.QAnnealFac = 1.0 + (cpmdremote.TAnnealFac / 10.0);
                } else if ((6 * cpmdremote.annealfreq) <= num <
                           (7 * cpmdremote.annealfreq)) // If in 8th, decrement temperature by 8.0x.
                { //  This phase will bring the (T = (10^-4)*T_0).
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
            for (unsigned int b = 0; b < real_bath.size(); b++) {
                real_bath[b].T = fT * real_bath[b].T;
                real_bath[b].Q = fQ * real_bath[b].Q;
            }
        }

        // VJ's instantaneous annealing for constant decrement (fT, fQ):
        /*    if (cpmdremote.anneal == 'y' && num % cpmdremote.annealfreq == 0 && real_bath[0].T >= 0.000000000000001)
        {
          cout << "annealing ... " << endl;
          for (unsigned int b = 0; b < real_bath.size(); b++)
          {
            real_bath[b].T = real_bath[b].T / 10;
            real_bath[b].Q = real_bath[b].Q / 2;
          }
        }*/

        //  2017.01.01 NB:  Abort the entire program if the global energy has drifted upward by more than 5%:
        if (1.05 <
            ((boundary.energy + bath_kinetic_energy(real_bath) + bath_potential_energy(real_bath)) / initNetEnergy))
            abort();
            //  Also abort if annealing hasn't yet begun, but the global energy has drifted downward by more than 5%:
        else if ((num < cpmdremote.annealfreq) &&
                 (((boundary.energy + bath_kinetic_energy(real_bath) + bath_potential_energy(real_bath)) /
                   initNetEnergy) < 0.95))
            abort();

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