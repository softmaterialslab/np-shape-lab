// This file contains the routines 

#include "functions.h"
#include <sstream>

// NB removed an extraneous factor of 2 from the denominator of p_sigma in the "initialize_vertex_velocities" command. With the volume constraint, the factor of 2 in the deminator led to starting at T = .5*T_target .

// overload out
ostream &operator<<(ostream &os, VECTOR3D vec) {
    if (world.rank() == 0)
        os << vec.x << setw(15) << vec.y << setw(15) << vec.z;
    return os;

}

// initialize velocities of vertices to start simulation
void initialize_vertex_velocities(vector<VERTEX> &V, vector<THERMOSTAT> &bath) {
    if (bath.size() == 1) {
        for (unsigned int i = 0; i < V.size(); i++)
            V[i].velvec = VECTOR3D(0, 0, 0);                    // initialized velocities
        if (world.rank() == 0)
            cout << "No thermostat for real system" << endl;
        return;
    }
    double p_sigma = sqrt(kB * bath[0].T / (/* 2 * */ V[0].m)); // Maxwell dist width; NB removed denom factor of 2.
    UTILITY ugsl;                                               // With that denom factor of 2, it starts at (T_target/2).
    for (unsigned int i = 0;
         i < V.size(); i++)                 // Starting at T_target seems to help stability at times.
        V[i].velvec = VECTOR3D(gsl_ran_gaussian(ugsl.r,p_sigma), gsl_ran_gaussian(ugsl.r,p_sigma), gsl_ran_gaussian(ugsl.r,p_sigma));	// initialized velocities
    VECTOR3D average_velocity_vector = VECTOR3D(0, 0, 0);
    for (unsigned int i = 0; i < V.size(); i++)
        average_velocity_vector = average_velocity_vector + V[i].velvec;
    average_velocity_vector = average_velocity_vector ^ (1.0 / V.size());
    for (unsigned int i = 0; i < V.size(); i++)
        V[i].velvec = V[i].velvec - average_velocity_vector;
    return;
}

// NB 2017.09.04:  initialize vertex/mesh velocities to zero
void initialize_velocities_to_zero(vector<VERTEX> &V, vector<PARTICLE> &counterions) {
    for (unsigned int i = 0; i < V.size(); i++) {
        V[i].velvec = VECTOR3D(0, 0, 0);
    }
    for (unsigned int i = 0; i < counterions.size(); i++) {
        counterions[i].velvec = VECTOR3D(0, 0, 0);
    }
}

// interface movie for VMD
void interface_movie(int num, vector<VERTEX> &V, vector<PARTICLE> &counterions, double box_halflength, double radius) {
    if (world.rank() == 0) {
        ofstream outdump("outfiles/p.lammpstrj", ios::app);

        // NB has formatted this to only output of decimal values (instead of Fortran format e.g. .00009 vs. 9e-5) at a certain precision using the following two lines, respectively.  The spacing in setw(~) also changed slightly.
        outdump << setprecision(12);
        outdump << std::fixed;

        outdump << "ITEM: TIMESTEP" << endl;
        outdump << num - 1 << endl;
        outdump << "ITEM: NUMBER OF ATOMS" << endl;
        outdump << V.size() + counterions.size() << endl;
        outdump << "ITEM: BOX BOUNDS" << endl;
        outdump << -box_halflength << "\t" << box_halflength << endl;
        outdump << -box_halflength << "\t" << box_halflength << endl;
        outdump << -box_halflength << "\t" << box_halflength << endl;
        outdump << "ITEM: ATOMS index type x y z Vq Va Vq/a" << endl;
        string type;
        for (unsigned int i = 0; i < V.size(); i++) {
            if (V[i].q > 0)
                type = "1";
            else if (V[i].q < 0)
                type = "-1";
            else
                type = "0";
            outdump << i << "\t" << type << "\t" << V[i].posvec.x << "\t" << V[i].posvec.y << "\t" << V[i].posvec.z
                    << "\t"
                    << V[i].q << "\t" << V[i].itsarea << "\t" << (V[i].q / (V[i].itsarea*radius*radius)) << endl;
        }
        for (unsigned int i = 0; i < counterions.size(); i++) {
            outdump << i + V.size() << "\t" << 2 << "\t" << counterions[i].posvec.x << "\t" << counterions[i].posvec.y
            << "\t" << counterions[i].posvec.z << "\t" << counterions[i].q << "\t" << 0 << "\t" << 0 << endl;
        }
        outdump.close();
    }
    return;
}


//input coordinate generated from orignal + dual meshes. 
void create_input_coordinate(vector<VERTEX>& V, vector<VERTEX>& Dual,vector<PARTICLE>& counterions, double box_halflength, double qLJ, double diameter) {
    if (world.rank() == 0) {
        ofstream outdump("outfiles/initCoordi.ShapeCondensation", ios::app);

        outdump << setprecision(12);
        outdump << std::fixed;

        outdump << "LAMMPS data file" << endl;
        outdump << V.size() + Dual.size() +counterions.size() << "\t" << "atoms" <<endl;
        outdump << "3 atom types" << endl;
        outdump << -box_halflength << "\t" << box_halflength << "\t" << "xlo" << "\t" << "xhi" << endl;
        outdump << -box_halflength << "\t" << box_halflength << "\t" << "ylo" << "\t" << "yhi" << endl;
        outdump << -box_halflength << "\t" << box_halflength << "\t" << "zlo" << "\t" << "zhi" << endl;
        outdump << endl;
        outdump << "Atoms" << endl;
        outdump << endl;
        //string type;
        for (unsigned int i = 0; i < V.size(); i++) {
            outdump << i+1 << "\t" << 1 << "\t" << V[i].posvec.x << "\t" << V[i].posvec.y << "\t" << V[i].posvec.z
                << "\t"
                << V[i].q * qLJ  << "\t" << diameter << "\t" << 1.0/((4.0/3.0)*pi* pow(diameter/2.0, 3)) << "\t" << endl;
        }
        for (unsigned int i = 0; i < Dual.size(); i++) {
            outdump << i + 1 + V.size() << "\t" << 1 << "\t" << Dual[i].posvec.x << "\t" << Dual[i].posvec.y << "\t" << Dual[i].posvec.z
                << "\t"
                << Dual[i].q * qLJ << "\t" << diameter << "\t" << 1.0 / ((4.0 / 3.0) * pi * pow(diameter / 2.0, 3))
                << "\t" << endl;
        }
        for (unsigned int i = 0; i < counterions.size(); i++) {
            string type;
            if (counterions[i].q > 0) {
                type = "2";
            }
            else {
                type = "3";
            }

            outdump << i + 1 + V.size()+ Dual.size() << "\t" << type << "\t" << counterions[i].posvec.x << "\t" << counterions[i].posvec.y
                << "\t" << counterions[i].posvec.z << "\t" << counterions[i].q * qLJ << "\t" << diameter << "\t" << 1.0 / ((4.0 / 3.0) * pi * pow(diameter / 2.0, 3)) << "\t" << endl;
        }
        outdump.close();
    }
    return;
}









// interface movie off files
void interface_off(int num, INTERFACE &dsphere) {
    stringstream ss;
    if (world.rank() == 0) {
        ss << "outfiles/" << num << ".off";
        ofstream outdump(ss.str().c_str(), ios::out);
        outdump << "OFF\n"
                << dsphere.number_of_vertices << " "
                << dsphere.number_of_faces << " "
                << dsphere.number_of_edges << "\n";
        for (unsigned int i = 0; i < dsphere.V.size(); i++) {
            outdump << dsphere.V[i].posvec.x << " "
                    << dsphere.V[i].posvec.y << " "
                    << dsphere.V[i].posvec.z << "\n";
        }
        for (unsigned int i = 0; i < dsphere.F.size(); i++) {
            outdump << "3 "
                    << dsphere.F[i].itsV[0]->index << " "
                    << dsphere.F[i].itsV[1]->index << " "
                    << dsphere.F[i].itsV[2]->index << " ";
            if (dsphere.F[i].itsV[0]->q > 0 &&
                dsphere.F[i].itsV[1]->q > 0 &&
                dsphere.F[i].itsV[2]->q > 0)
                outdump << "255 0 0 1\n";
            else if (dsphere.F[i].itsV[2]->q < 0 &&
                     dsphere.F[i].itsV[1]->q < 0 &&
                     dsphere.F[i].itsV[0]->q < 0)
                outdump << "0 0 255 1\n";
            else
                outdump << "0 0 0 0\n";
        }
        outdump.close();
    }
    return;
}

// ### Unused functions: ###

// compute MD trust factor R
double compute_MD_trust_factor_R(int hiteqm) {
    char filename[200];
    sprintf(filename, "outfiles/energy.dat");
    ifstream in(filename, ios::in);
    if (!in) {
        if (world.rank() == 0)
            cout << "File could not be opened" << endl;
        return 0;
    }

    int col1;
    double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11;
    vector<double> ext, ke, pe, fake, real;
    while (in >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11) {
        ext.push_back(col2);
        ke.push_back(col3);
        pe.push_back(col4);
        fake.push_back(col6);
        real.push_back(col5);
    }

    double ext_mean = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_mean += ext[i];
    ext_mean = ext_mean / ext.size();
    double ke_mean = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
        ke_mean += ke[i];
    ke_mean = ke_mean / ke.size();

    double ext_sd = 0;
    for (unsigned int i = 0; i < ext.size(); i++)
        ext_sd += (ext[i] - ext_mean) * (ext[i] - ext_mean);
    ext_sd = ext_sd / ext.size();
    ext_sd = sqrt(ext_sd);

    double ke_sd = 0;
    for (unsigned int i = 0; i < ke.size(); i++)
        ke_sd += (ke[i] - ke_mean) * (ke[i] - ke_mean);
    ke_sd = ke_sd / ke.size();
    ke_sd = sqrt(ke_sd);

    double R = ext_sd / ke_sd;
    if (world.rank() == 0) {
        ofstream out("outfiles/R.dat");
        out << "Sample size " << ext.size() << endl;
        out << "Sd: ext, kinetic energy and R" << endl;
        out << ext_sd << setw(15) << ke_sd << setw(15) << R << endl;
    }
    return R;
}

// auto correlation function
void auto_correlation_function() {
    char filename[200];
    sprintf(filename, "outfiles/for_auto_corr.dat");
    ifstream in(filename, ios::in);
    if (!in) {
        cout << "File could not be opened" << endl;
        return;
    }

    double col1, col2;
    vector<double> n, autocorr;
    while (in >> col1 >> col2)
        n.push_back(col2);

    double avg = 0;
    for (unsigned int j = 0; j < n.size(); j++)
        avg = avg + n[j];
    avg = avg / n.size();

    int ntau = 5000;            // time to which the auto correlation function is computed

    for (int i = 0; i < ntau; i++) {
        double A = 0;
        for (unsigned int j = 0; j < n.size(); j++)
            A = A + n[j + i] * n[j];
        A = A / n.size();
        autocorr.push_back(A - avg * avg);
    }
    if (world.rank() == 0) {
        ofstream out("outfiles/auto_correlation.dat");
        for (int i = 0; i < ntau; i++)
            out << i << setw(15) << autocorr[i] / autocorr[0] << endl;

        cout << "Auto correlation function generated" << endl;
    }
    return;
}

// interface movie povray
void interface_pov(int num, INTERFACE &dsphere) {

    if (world.rank() == 0) {
        stringstream ss;
        ss << "outfiles/" << num << ".pov";
        ofstream outdump(ss.str().c_str(), ios::out);
        outdump << "background{ rgb<1,1,1> }\n"
                << "camera {\n"
                << "  location <3, 0, 0>\n"
                << "  look_at <0, 0, 0>\n"
                << "}\n"
                << "light_source { <50, 50, -50> color rgb<1, 1, 1> }\n"
                << "#declare Plus = texture {\n"
                << "  pigment { color rgb<1.0, 0.5, 0.5> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare Minus = texture {\n"
                << "  pigment { color rgb<0.5, 0.5, 1.0> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare BoundaryColor = texture {\n"
                << "  pigment { color rgb<0.9, 0.9, 0.9> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare VertexPlusColor = texture {\n"
                << "  pigment { color rgb<1.0, 0.0, 0.0> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare VertexMinusColor = texture {\n"
                << "  pigment { color rgb<0.0, 0.0, 1.0> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare VertexNeutralColor = texture {\n"
                << "  pigment { color rgb<0.3, 0.3, 0.3> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare EdgeColor = texture {\n"
                << "  pigment { color rgb<0.0, 0.0, 0.0> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare CylinderRadius = " << 0.1 * dsphere.avg_edge_length
                << ";\n"
                << "#declare SphereRadius   = " << 0.5 * dsphere.avg_edge_length
                << ";\n";

        outdump << "mesh {\n";
        for (unsigned int i = 0; i < dsphere.F.size(); i++) {
            outdump << "  triangle {\n"
                    << "    <"
                    << dsphere.F[i].itsV[0]->posvec.x << ", "
                    << dsphere.F[i].itsV[0]->posvec.y << ", "
                    << dsphere.F[i].itsV[0]->posvec.z << ">,\n    <"
                    << dsphere.F[i].itsV[1]->posvec.x << ", "
                    << dsphere.F[i].itsV[1]->posvec.y << ", "
                    << dsphere.F[i].itsV[1]->posvec.z << ">,\n    <"
                    << dsphere.F[i].itsV[2]->posvec.x << ", "
                    << dsphere.F[i].itsV[2]->posvec.y << ", "
                    << dsphere.F[i].itsV[2]->posvec.z << ">\n";
            if (dsphere.F[i].itsV[0]->q > 0 &&
                dsphere.F[i].itsV[1]->q > 0 &&
                dsphere.F[i].itsV[2]->q > 0)
                outdump << "    texture { Plus }\n";
            else if (dsphere.F[i].itsV[2]->q < 0 &&
                     dsphere.F[i].itsV[1]->q < 0 &&
                     dsphere.F[i].itsV[0]->q < 0)
                outdump << "    texture { Minus }\n";
            else
                outdump << "    texture { BoundaryColor }\n";
            outdump << "  }\n";
        }
        outdump << "}\n";
        for (unsigned int i = 0; i < dsphere.E.size(); i++) {
            outdump << "cylinder {\n"
                    << "  <"
                    << dsphere.E[i].itsV[0]->posvec.x << ", "
                    << dsphere.E[i].itsV[0]->posvec.y << ", "
                    << dsphere.E[i].itsV[0]->posvec.z << ">,\n  <"
                    << dsphere.E[i].itsV[1]->posvec.x << ", "
                    << dsphere.E[i].itsV[1]->posvec.y << ", "
                    << dsphere.E[i].itsV[1]->posvec.z << ">\n"
                    << "  CylinderRadius\n  texture { EdgeColor }\n}\n";
        }
        for (unsigned int i = 0; i < dsphere.V.size(); i++) {
            outdump << "sphere {\n"
                    << "  <"
                    << dsphere.V[i].posvec.x << ", "
                    << dsphere.V[i].posvec.y << ", "
                    << dsphere.V[i].posvec.z << ">\n";
            if (dsphere.V[i].q > 0)
                outdump << "  SphereRadius\n  texture { VertexPlusColor }\n}\n";
            else if (dsphere.V[i].q < 0)
                outdump << "  SphereRadius\n  texture { VertexMinusColor }\n}\n";
            else
                outdump << "  SphereRadius\n  texture { VertexNeutralColor }\n}\n";

        }
    }

}

// interface movie povray smooth
void interface_pov_smooth(int num, INTERFACE &dsphere) {
    if (world.rank() == 0) {
        stringstream ss;
        ss << "outfiles/" << num << "_smooth.pov";
        ofstream outdump(ss.str().c_str(), ios::out);
        outdump << "background{ rgb<1,1,1> }\n"
                << "camera {\n"
                << "  location <3, 0, 0>\n"
                << "  look_at <0, 0, 0>\n"
                << "}\n"
                << "light_source { <50, 50, -50> color rgb<1, 1, 1> }\n"
                << "#declare Plus = texture {\n"
                << "  pigment { color rgb<1.0, 0.5, 0.5> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare Minus = texture {\n"
                << "  pigment { color rgb<0.5, 0.5, 1.0> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare BoundaryColor = texture {\n"
                << "  pigment { color rgb<0.9, 0.9, 0.9> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare VertexPlusColor = texture {\n"
                << "  pigment { color rgb<1.0, 0.0, 0.0> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare VertexMinusColor = texture {\n"
                << "  pigment { color rgb<0.0, 0.0, 1.0> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare VertexNeutralColor = texture {\n"
                << "  pigment { color rgb<0.3, 0.3, 0.3> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare EdgeColor = texture {\n"
                << "  pigment { color rgb<0.0, 0.0, 0.0> }\n"
                << "  finish { ambient 0.2 diffuse 0.5 }\n"
                << "}\n"
                << "#declare CylinderRadius = " << 0 << ";\n"
                << "#declare SphereRadius   = " << 0 << ";\n";

        outdump << "mesh {\n";
        for (unsigned int i = 0; i < dsphere.F.size(); i++) {
            outdump << "  smooth_triangle {\n"
                    << "    <"
                    << dsphere.F[i].itsV[0]->posvec.x << ", "
                    << dsphere.F[i].itsV[0]->posvec.y << ", "
                    << dsphere.F[i].itsV[0]->posvec.z << ">,\n    <"
                    << dsphere.F[i].itsV[0]->itsnormal.x << ", "
                    << dsphere.F[i].itsV[0]->itsnormal.y << ", "
                    << dsphere.F[i].itsV[0]->itsnormal.z << ">,\n    <"
                    << dsphere.F[i].itsV[1]->posvec.x << ", "
                    << dsphere.F[i].itsV[1]->posvec.y << ", "
                    << dsphere.F[i].itsV[1]->posvec.z << ">,\n    <"
                    << dsphere.F[i].itsV[1]->itsnormal.x << ", "
                    << dsphere.F[i].itsV[1]->itsnormal.y << ", "
                    << dsphere.F[i].itsV[1]->itsnormal.z << ">,\n    <"
                    << dsphere.F[i].itsV[2]->posvec.x << ", "
                    << dsphere.F[i].itsV[2]->posvec.y << ", "
                    << dsphere.F[i].itsV[2]->posvec.z << ">,\n    <"
                    << dsphere.F[i].itsV[2]->itsnormal.x << ", "
                    << dsphere.F[i].itsV[2]->itsnormal.y << ", "
                    << dsphere.F[i].itsV[2]->itsnormal.z << ">\n";
            if (dsphere.F[i].itsV[0]->q > 0 &&
                dsphere.F[i].itsV[1]->q > 0 &&
                dsphere.F[i].itsV[2]->q > 0)
                outdump << "    texture { Plus }\n";
            else if (dsphere.F[i].itsV[2]->q < 0 &&
                     dsphere.F[i].itsV[1]->q < 0 &&
                     dsphere.F[i].itsV[0]->q < 0)
                outdump << "    texture { Minus }\n";
            else
                outdump << "    texture { BoundaryColor }\n";
            outdump << "  }\n";
        }
        outdump << "}\n";
        for (unsigned int i = 0; i < dsphere.E.size(); i++) {
            outdump << "cylinder {\n"
                    << "  <"
                    << dsphere.E[i].itsV[0]->posvec.x << ", "
                    << dsphere.E[i].itsV[0]->posvec.y << ", "
                    << dsphere.E[i].itsV[0]->posvec.z << ">,\n  <"
                    << dsphere.E[i].itsV[1]->posvec.x << ", "
                    << dsphere.E[i].itsV[1]->posvec.y << ", "
                    << dsphere.E[i].itsV[1]->posvec.z << ">\n"
                    << "  CylinderRadius\n  texture { EdgeColor }\n}\n";
        }
        for (unsigned int i = 0; i < dsphere.V.size(); i++) {
            outdump << "sphere {\n"
                    << "  <"
                    << dsphere.V[i].posvec.x << ", "
                    << dsphere.V[i].posvec.y << ", "
                    << dsphere.V[i].posvec.z << ">\n";
            if (dsphere.V[i].q > 0)
                outdump << "  SphereRadius\n  texture { VertexPlusColor }\n}\n";
            else if (dsphere.V[i].q < 0)
                outdump << "  SphereRadius\n  texture { VertexMinusColor }\n}\n";
            else
                outdump << "  SphereRadius\n  texture { VertexNeutralColor }\n}\n";
        }
    }

}

void progressBar(double fraction_completed) {

    if (world.rank() == 0) {
        int val = (int) (fraction_completed * 100);
        int lpad = (int) (fraction_completed * PBWIDTH);
        int rpad = PBWIDTH - lpad;
        printf("\r%3d%% |%.*s%*s|", val, lpad, PBSTR, rpad, "");
        fflush(stdout);
    }
}