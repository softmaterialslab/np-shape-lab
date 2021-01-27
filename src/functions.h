// This is a header file containing functions

#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include "utility.h"
#include "interface.h"
#include "vertex.h"
#include "control.h"
#include "thermostat.h"
#include "mpi_utility.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// ### General functions: ###

// Overloading << to print VECTOR3D:
ostream &operator<<(ostream &, VECTOR3D);

// post analysis : compute R
double compute_MD_trust_factor_R(int);

// post analysis : auto correlation function
void auto_correlation_function();

// display progress bar (code from the internet)
void progressBar(double);

// ### Functions useful in implementing constraints: ###

// SHAKE to ensure the volume constraint is true:  (2017.09.06 NB added conditional bypass if constraint is true already.)
inline void SHAKE(INTERFACE &boundary, CONTROL &mdremote, char constraintForm) {
    // (2019.09.17 NB added.)  Conditionally render the constraint quadratic, if desired (otherwise default to linear):
    if (constraintForm == 'Q')  // Quadratic slightly better for moderate parameterizations, blows up for extreme ones.
    {
        for (unsigned int k = 0; k < boundary.F.size(); k++) {
            for (unsigned int i = 0; i < boundary.F[k].itsV.size(); i++)
                boundary.F[k].itsV[i]->neg_GradConi = ((2 * (boundary.total_volume - boundary.ref_volume)) ^
                                                       boundary.F[k].itsV[i]->neg_GradVi);
        }
    } else // If not quadratic, then linear.  This is the defaut; any non-'Q' flag yields original linear constraint.
    {
        for (unsigned int k = 0; k < boundary.F.size(); k++) {
            for (unsigned int i = 0; i < boundary.F[k].itsV.size(); i++)
                boundary.F[k].itsV[i]->neg_GradConi = boundary.F[k].itsV[i]->neg_GradVi;
        }
    }

    double f0 = 0;
    double fa1 = 0;
    double fa2 = 0;
    double fa3 = 0;
    double fb1 = 0;
    double fb2 = 0;
    double fb3 = 0;
    double fc = 0;
    for (unsigned int i = 0; i < boundary.F.size(); i++) {
        f0 = f0 + boundary.F[i].itsV[0]->neg_GradConi *
                  (boundary.F[i].itsV[1]->neg_GradConi % boundary.F[i].itsV[2]->neg_GradConi);
        fa1 = fa1 + boundary.F[i].itsV[0]->neg_GradConi *
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->neg_GradConi);
        fa2 = fa2 + boundary.F[i].itsV[0]->neg_GradConi *
                    (boundary.F[i].itsV[1]->neg_GradConi % boundary.F[i].itsV[2]->posvec);
        fa3 = fa3 + boundary.F[i].itsV[0]->posvec *
                    (boundary.F[i].itsV[1]->neg_GradConi % boundary.F[i].itsV[2]->neg_GradConi);
        fb1 = fb1 +
              boundary.F[i].itsV[0]->neg_GradConi * (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->posvec);
        fb2 = fb2 +
              boundary.F[i].itsV[0]->posvec * (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->neg_GradConi);
        fb3 = fb3 +
              boundary.F[i].itsV[0]->posvec * (boundary.F[i].itsV[1]->neg_GradConi % boundary.F[i].itsV[2]->posvec);
        fc = fc + boundary.F[i].itsV[0]->posvec * (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->posvec);
    }

    double a, b, c;
    // Assuming mass = 1.0:
    a = (1.0 / mdremote.timestep) * (fa1 + fa2 + fa3) / f0;
    b = (1.0 / mdremote.timestep) * (1.0 / mdremote.timestep) * (fb1 + fb2 + fb3) / f0;
    c = (1.0 / mdremote.timestep) * (1.0 / mdremote.timestep) * (1.0 / mdremote.timestep) *
        (fc - 6 * boundary.ref_volume) / f0;

    // If the constraint is already exactly satisfied (c = 0), bypass the SHAKE:
    if ((c == 0) || (f0 ==
                     0)) { // If scaled volume discrepancy (c = 0), already satisfied (rare). If constraint gradient = 0, already satisfied.
        //Using in-principle-equivalent (boundary.total_volume - boundary.ref_volume) diff. as the if argument leads to less precision, best to use this.
        if (world.rank() == 0) {
            cout << "Volume_SHAKE:  Constraint already satisfied: " << c * f0 * pow(mdremote.timestep, 3) / 6
                 << ".  Bypassing SHAKE." << endl;
        }
        return;  // This is necessary in the quadratic case to avoid propagation of "nan" errors at t = 0.
    }

    // Equation is lambda^3  +  a lambda^2  +  b lambda + c = 0 , solving for lambda_V:

    // Declaring the 3 roots as required by GSL function (always puts a real in x0):
    double x0;
    double x1;
    double x2;
    int n; // number of real roots (one or three)

    n = gsl_poly_solve_cubic(a, b, c, &x0, &x1, &x2);
	 
	 if (n == 2)
		 cout << "INTERNAL WARNING: Unexpected number of roots in SHAKE" << endl;	// complex roots always come in pairs!

    //cout << "Volume_SHAKE:  Constraint NOT satisfied: " << c*pow((mdremote.timestep),3)*f0/6/*boundary.total_volume - boundary.ref_volume*/ << ". \tSHAKING, multiplier:  " << x0 << "." << endl;

    // Propagate positions and velocities with the additional constraint force (assuming mass = 1.0):
    for (unsigned int i = 0; i < boundary.V.size(); i++)
        boundary.V[i].velvec = boundary.V[i].velvec + (x0 ^ boundary.V[i].neg_GradConi);
    for (unsigned int i = 0; i < boundary.V.size(); i++) {
        boundary.V[i].posvec = boundary.V[i].posvec + ((x0 * mdremote.timestep) ^ boundary.V[i].neg_GradConi);
    }

    return;
}

// RATTLE to ensure time derivative of the volume constraint is true:  (2017.09.06 NB added conditional bypass, see comments below.)
inline void RATTLE(INTERFACE &boundary, char constraintForm) {
    // (2019.09.17 NB added.)  Conditional rendering constraint quadratic, if desired (otherwise default to linear):
    if (constraintForm == 'Q') {
        for (unsigned int k = 0; k < boundary.F.size(); k++) {
            for (unsigned int i = 0; i < boundary.F[k].itsV.size(); i++)
                boundary.F[k].itsV[i]->neg_GradConi = ((2.0 * (boundary.total_volume - boundary.ref_volume)) ^
                                                       boundary.F[k].itsV[i]->neg_GradVi);
        }
    } else // If not quadratic, then linear.  This is the defaut; any non-'Q' flag yields original linear constraint.
    {
        for (unsigned int k = 0; k < boundary.F.size(); k++) {
            for (unsigned int i = 0; i < boundary.F[k].itsV.size(); i++)
                boundary.F[k].itsV[i]->neg_GradConi = boundary.F[k].itsV[i]->neg_GradVi;
        }
    }

    double mu;
    double num, den;

    num = 0;
    den = 0;


    for (unsigned int i = 0; i < boundary.F.size(); i++) {
        // NB:  Notice the numerator is d/dt(V_tot):
        num = num + boundary.F[i].itsV[0]->velvec * (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->posvec)
              + boundary.F[i].itsV[0]->posvec * (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->velvec)
              + boundary.F[i].itsV[0]->posvec * (boundary.F[i].itsV[1]->velvec % boundary.F[i].itsV[2]->posvec);

        den = den +
              boundary.F[i].itsV[0]->neg_GradConi * (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->posvec)
              + boundary.F[i].itsV[0]->posvec * (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->neg_GradConi)
              + boundary.F[i].itsV[0]->posvec * (boundary.F[i].itsV[1]->neg_GradConi % boundary.F[i].itsV[2]->posvec);
    }

    // If the constraint is already exactly satisfied (only occurs initially for quadratic), bypass the SHAKE:
    if (num == 0 || den ==
                    0) {   // num = 0 if constraint is already satisfied (d/dt(V_tot) = 0) ; den = 0 if constraint gradient (g_i(t)) = 0.
        if (world.rank() == 0)
            cout << "Volume_RATTLE: Derivative constraint satisfied.  Bypassing RATTLE." << endl;
        return; // This is necessary in the quadratic case to avoid propagation of "nan" errors at t = 0.
    }

    mu = -(num / (den));

    //cout << "Volume_RATTLE: Derivative constraint is NOT satisfied. \tRATTLING, multiplier: " << mu << ". ";

    for (unsigned int i = 0; i < boundary.V.size(); i++) {
        boundary.V[i].velvec = boundary.V[i].velvec + (mu ^ boundary.V[i].neg_GradConi);
    }
    return;
}

// ### Functions useful in Nose-Hoover chain implementation: ###

// update bath xi value
inline void update_chain_xi(unsigned int j, vector<THERMOSTAT> &bath, double dt, long double ke) {
    if (bath[j].Q == 0)
        return;
    if (j != 0)
        bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j + 1].xi) + 0.5 * dt * (1.0 / bath[j].Q) *
                                                                    (bath[j - 1].Q * bath[j - 1].xi * bath[j - 1].xi -
                                                                     bath[j].dof * kB * bath[j].T) *
                                                                    exp(-0.25 * dt * bath[j + 1].xi);
    else
        bath[j].xi = bath[j].xi * exp(-0.5 * dt * bath[j + 1].xi) +
                     0.5 * dt * (1.0 / bath[j].Q) * (2 * ke - bath[j].dof * kB * bath[j].T) *
                     exp(-0.25 * dt * bath[j + 1].xi);
    return;
}

// ### New functions (as of ? prior to NB): ###

// interface movie
void interface_movie(int num, vector<VERTEX> &s, vector<PARTICLE> &counterions, double box_radius);

void create_input_coordinate(vector<VERTEX>& s, vector<PARTICLE>& counterions, double box_radius);

void interface_pov(int num, INTERFACE &dsphere);

void interface_pov_smooth(int num, INTERFACE &dsphere);

void interface_off(int num, INTERFACE &dsphere);

// initialize velocities of vertices
void initialize_vertex_velocities(vector<VERTEX> &V, vector<THERMOSTAT> &bath);

// initialize velocities of vertics to be zero
void initialize_velocities_to_zero(vector<VERTEX> &V, vector<PARTICLE> &counterions);

// ### Unused functions: ###

// computes gradient of green's fn
inline VECTOR3D Grad(VECTOR3D &vec1, VECTOR3D &vec2) {
    long double r = (vec1 - vec2).GetMagnitude();
    long double r3 = r * r * r;
    return (vec1 - vec2) ^ ((-1.0) / r3);
}

// computes gradient of normal dot gradient of 1/r
inline VECTOR3D GradndotGrad(VECTOR3D &vec1, VECTOR3D &vec2, VECTOR3D &normal) {
    long double r = (vec1 - vec2).GetMagnitude();
    long double r3 = r * r * r;
    long double r5 = r3 * r * r;
    return ((normal ^ (1.0 / r3)) - ((vec1 - vec2) ^ (3 * (normal * (vec1 - vec2)) / r5)));
}

// constraint equation
inline long double constraint(vector<VERTEX> &s, INTERFACE &dsphere) {
//   return dsphere.total_induced_charge(s) - dsphere.total_charge_inside(ion) * (1/dsphere.eout - 1/dsphere.ein);
    return 0;
}

// dot constraint equation
inline long double dotconstraint(vector<VERTEX> &s) {
    long double sigmadot = 0;
//   for (unsigned int k = 0; k < s.size(); k++)
//     sigmadot += s[k].vw * s[k].a;
    return sigmadot;
}

//  ### Works in progress: ###

// SHAKE_for_area to ensure the area constraint is true:  (Unfinished, vector issues.)
inline void SHAKE_for_area(INTERFACE &boundary, CONTROL &mdremote, char constraintForm) {
    double fa = 0;
    double fb = 0;
    double fc = 0;

    double x0 = 0;
    double x1 = 0;

    // If-statement bypassing of the SHAKE if the constraint is already satisfied:
    if (boundary.total_area == boundary.ref_area) {
        if (world.rank() == 0)
            cout << "SHAKE:  Constraint already satisfied, value: " << boundary.total_volume - boundary.ref_volume
                 << ".  Bypassing SHAKE, multiplier:  " << x0 << "." << endl;
    } else {   // Running SHAKE to regulate positions, velocities:
        for (unsigned int i = 0; i < boundary.F.size(); i++) {

            fa = fa + (
                    (boundary.F[i].itsV[1]->neg_GradAi % boundary.F[i].itsV[2]->neg_GradAi) -
                    (boundary.F[i].itsV[1]->neg_GradAi % boundary.F[i].itsV[0]->neg_GradAi) -
                    (boundary.F[i].itsV[0]->neg_GradAi % boundary.F[i].itsV[2]->neg_GradAi)).GetMagnitude();
            fb = fb + (
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->neg_GradAi) +
                    (boundary.F[i].itsV[1]->neg_GradAi % boundary.F[i].itsV[2]->posvec) -
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[0]->neg_GradAi) -
                    (boundary.F[i].itsV[1]->neg_GradAi % boundary.F[i].itsV[0]->posvec) -
                    (boundary.F[i].itsV[0]->posvec % boundary.F[i].itsV[2]->neg_GradAi) -
                    (boundary.F[i].itsV[0]->neg_GradAi % boundary.F[i].itsV[2]->posvec)).GetMagnitude();
            fc = fc + (
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->posvec) -
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[0]->posvec) -
                    (boundary.F[i].itsV[0]->posvec % boundary.F[i].itsV[2]->posvec)).GetMagnitude();
        }

        double a, b, c;

        // Definitions without scaling for constraint gradient:
        //a = fa;
        //b = (1.0/mdremote.timestep) * fb;
        //c = (1.0/mdremote.timestep) * (1.0/mdremote.timestep) * (fc - 2*boundary.ref_area);

        // Scaled definitions:
        a = fa / 4.;
        b = (1.0 / mdremote.timestep) * fb / 2.;
        c = (1.0 / mdremote.timestep) * (1.0 / mdremote.timestep) * (fc - 2 * boundary.ref_area);

        // Eqn (lambda^2  +  b lambda + c = 0) -> la_A:
        x0 = (-b - sqrt(pow(b, 2) - 4 * a * c)) /
             (2.0 * a); // (Init: {x0 = -25509.2}; typical {x0 <= this}) for qStr = 600.
        x1 = (-b + sqrt(pow(b, 2) - 4 * a * c)) /
             (2.0 * a); // (Init: {x1 = 0, typically x1 ~ -0.1 linear}) for qStr = 600.

        // Regulating positions and velocities with the additional constraint force.  (Mass assumed to be 1.)
        for (unsigned int i = 0; i < boundary.V.size(); i++)
            //boundary.V[i].velvec = boundary.V[i].velvec + (abs(x1) ^ boundary.V[i].dummy_aforce);
            boundary.V[i].velvec = boundary.V[i].velvec + (x1 ^ (0.5 ^ boundary.V[i].neg_GradAi));
        //boundary.V[i].velvec = boundary.V[i].velvec + (x1 * 0.5 ^ boundary.V[i].dummy_aforce);
        for (unsigned int i = 0; i < boundary.V.size(); i++) {
            //boundary.V[i].posvec = boundary.V[i].posvec + ((abs(x1) * mdremote.timestep) ^ boundary.V[i].dummy_aforce);
            boundary.V[i].posvec = boundary.V[i].posvec + ((x1 * mdremote.timestep) ^ (0.5 ^ boundary.V[i].neg_GradAi));
            //boundary.V[i].posvec = boundary.V[i].posvec + ((x1 * mdremote.timestep * 0.5) ^ boundary.V[i].dummy_aforce);  OLD
        }
        if (world.rank() == 0)
            cout << "SHAKE:  Constraint is NOT satisfied, value: " << boundary.total_area - boundary.ref_area
                 << ".  SHAKING, multiplier:  " << x1 << "." << endl;
    }

    return;
}

// RATTLE to ensure time derivative of the area constraint is true:  (Unfinished, vector issues.)
inline void RATTLE_for_area(INTERFACE &boundary, char constraintForm) {
    double mu, num, den;

    num = 0;
    den = 0;

    // RATTLE only needs bypassed if a quadratic constraint is used (latter condition), as linear does not render the denominator zero:
    if (boundary.total_area == boundary.ref_area) {
        if (world.rank() == 0)
            cout << "RATTLE: Constraint gradient zero.  Bypassing RATTLE to avoid infinite denominator." << endl;
    } else {
        for (unsigned int i = 0; i < boundary.F.size(); i++) {

            num = num + (
                    (boundary.F[i].itsV[1]->velvec % boundary.F[i].itsV[2]->posvec) +
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->velvec) -
                    (boundary.F[i].itsV[1]->velvec % boundary.F[i].itsV[0]->posvec) -
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[0]->velvec) -
                    (boundary.F[i].itsV[0]->velvec % boundary.F[i].itsV[2]->posvec) -
                    (boundary.F[i].itsV[0]->posvec % boundary.F[i].itsV[2]->velvec)).GetMagnitude();

            den = den + (
                    (boundary.F[i].itsV[1]->neg_GradAi % boundary.F[i].itsV[2]->posvec) +
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[2]->neg_GradAi) -
                    (boundary.F[i].itsV[1]->neg_GradAi % boundary.F[i].itsV[0]->posvec) -
                    (boundary.F[i].itsV[1]->posvec % boundary.F[i].itsV[0]->neg_GradAi) -
                    (boundary.F[i].itsV[0]->neg_GradAi % boundary.F[i].itsV[2]->posvec) -
                    (boundary.F[i].itsV[0]->posvec % boundary.F[i].itsV[2]->neg_GradAi)).GetMagnitude();
        }

        mu = -(num / (den / 2.));
        if (world.rank() == 0)
            cout << "RATTLE: Constraint is NOT satisfied, value: " << boundary.total_area - boundary.ref_area
                 << ".  RATTLING, multiplier:  " << mu << "." << endl;

        for (unsigned int i = 0; i < boundary.V.size(); i++)
            boundary.V[i].velvec = boundary.V[i].velvec + (mu ^ (0.5 ^ boundary.V[i].neg_GradAi));
        //boundary.V[i].velvec = boundary.V[i].velvec + (mu * 0.5 ^ boundary.V[i].dummy_aforce); OLD
    }
    return;
}

#endif

