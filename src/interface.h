// This is a header file for the INTERFACE class.  

#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "utility.h"
#include "vector3d.h"
#include "vertex.h"
#include "edge.h"
#include "face.h"
#include "particle.h"
#include <map>

class VERTEX;

class EDGE;

class FACE;

class INTERFACE {
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & V;
        ar & E;
        ar & F;
        ar & areaGradient;
        ar & ein;
        ar & eout;
        ar & em;
        ar & ed;
        ar & lambda_a;
        ar & lambda_v;
        ar & lambda_l;
        ar & bkappa;
        ar & ref_area;
        ar & ref_volume;
        ar & ref_Area_Vertices;
        ar & energy;
        ar & netMeshKE;
        ar & penergy;
        ar & lB_out;
        ar & inv_kappa_out;
        ar & number_of_vertices;
        ar & number_of_edges;
        ar & number_of_faces;
        ar & total_area;
        ar & total_volume;
        //ar & total_Area_Vertices;
        ar & rest_volume;
        ar & avg_edge_length;
        ar & elj;
        ar & lj_length_mesh_mesh;
        ar & sconstant;
        ar & sigma_a;
    }

public:

    vector<VERTEX> V;
    vector<VERTEX> Dual;

    vector<EDGE> E;
    vector<FACE> F;
    map<pair<VERTEX *, FACE *>, VECTOR3D> areaGradient;

    double ein;        // permittivity of inside medium
    double eout;        // permittivity of outside medium
    double em;          // permittivity at the interface, mean
    double ed;          // permittivity change at the inteface, difference scaled by 4 pi

    double lambda_a, lambda_v, lambda_l;
    double bkappa;
    double ref_area, ref_volume, ref_Area_Vertices;         // (2017.09.19 NB added ref_Area_Vertices.)
    double energy, netMeshKE, netIonKE, penergy;

    double lB_out;            // Bjerrum length outside
    double inv_kappa_out;            // debye length outside
//   double lB_in;			// Bjerrum length inside
//   double inv_kappa_in;		// debye length inside
//   double mean_sep_in;		// mean separation inside
//   double mean_sep_out;		// mean separation outside

    unsigned int number_of_vertices;
    unsigned int number_of_edges;
    unsigned int number_of_faces;
    double total_area, total_volume;
    //double total_Area_Vertices;   // (2017.09.19 NB added total_Area_Vertices.)
    double rest_volume;

    double avg_edge_length;
    double elj;            // strength of the lj potential
    double lj_length_mesh_mesh;
    double lj_length_mesh_ions;

    // stretching energy parameters
    double sconstant;
    double sigma_a;            // Strength of the surface tension penalty ((kB T) / L^2).

    // soft-volume constraint parameter:
    double sigma_v;

    EDGE *edge_between(VERTEX *, VERTEX *);

    // ### Geometric operations: ###
    void set_up(double);

    void load_configuration(string filename);

    void discretize(unsigned int, unsigned int);

    void dressup(double, double);

    void assign_dual_boundary_edges();

    void assign_dual_initial();

    void
    output_configuration();                                // Function to output information on the initial membrane.
    void assign_boundary_edges();                            // Not used.
    void update_velocity_interface();                        // Not used.
    void update_position_interface();                        // Not used.

    // ### Charge assignment functions: ###
    void assign_q_values(int num_divisions, double q_strength);

    void assign_q_spots(int num_divisions, int division_size, double q_strength);

    void assign_random_q_values(double q_strength, double alpha, int num_divisions, double fracChargedPatch, char randomFlag, char functionFlag);

    void assign_external_q_values(double q_strength, string externalPattern);

    void put_counterions(double q_actual, double unit_radius_sphere, double ion_diameter, double box_radius, vector<PARTICLE> &counterions, int counterion_valency);

    // ###  Energy computation operations: ###
    void compute_local_energies(const double scalefactor);  // Computes the local energetics profiles (creates local_*_E.off files).
    void compute_local_energies_by_component();             // Computes local elastic energy profiles (creates similar files to above).
    void compute_energy(int, vector<PARTICLE> &counterions, const double scalefactor, char bucklingFlag);     // Computes system energy.

    INTERFACE(double _ein = 1,
              double _eout = 1, double _lambda_a = 1, double _lambda_v = 1,
              double _ref_area = 1, double _ref_volume = 1) :
            ein(_ein), eout(_eout), lambda_a(_lambda_a),
            lambda_v(_lambda_v), ref_area(_ref_area), ref_volume(_ref_volume) {
    }

};

#endif
