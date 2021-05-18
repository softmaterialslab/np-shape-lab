// This file contains member functions for interface class

#include "interface.h"
#include "functions.h"
#include "newenergies.h"
#include "particle.h"
#include <algorithm>
#include <list>
#include <sstream>
#include <cmath>

void colormap(double x, double *r, double *g, double *b, double *a) {
    (*r) = 0;
    (*g) = x;
    (*b) = 1 - x;
    (*a) = 1;
	 return ;
}

void INTERFACE::set_up(double unit_radius_sphere) {
    // useful combinations of different dielectric constants (inside and outside)
    em = 0.5 * (ein + eout);
    ed = (eout - ein) / (4 * pi);

    // useful length scales signifying competition between electrostatics and entropy
    lB_out = (lB_water * epsilon_water / eout) / unit_radius_sphere;

    return;
}

// %%% NP Mesh Setup Functionality: %%%
// Dress up the mesh with normals areas and volume elements:
void INTERFACE::dressup(double _lambda_a, double _lambda_v) {
    for (unsigned int i = 0; i < F.size(); i++) {
        F[i].compute_area_normal();
        F[i].computevolume();
    }

    for (unsigned int i = 0; i < V.size(); i++) {
        V[i].compute_area_normal();
    }

    total_area = 0;
    for (unsigned int i = 0; i < F.size(); i++)
        total_area += F[i].itsarea;

    total_volume = 0;
    for (unsigned int i = 0; i < F.size(); i++)
        total_volume += F[i].subtends_volume;

    /*total_Area_Vertices = 0;                          // (2017.09.17 NB added.)  Computes current net area by vertices.
    for (unsigned int i = 0; i < V.size(); i++)
        total_Area_Vertices += V[i].itsarea;*/

    //   ofstream filenormals("outfiles/normal.dat", ios::out);
    //   ofstream fileareas("outfiles/area.dat", ios::out);
    //   ofstream filevolumes("outfiles/volume.dat", ios::out);
    //   for (unsigned int i = 0; i < F.size(); i++)
    //   {
    //     filenormals << F[i].itsnormal << endl;
    //     fileareas << F[i].itsarea << endl;
    //     filevolumes << F[i].subtends_volume << endl;
    //   }

    lambda_a = _lambda_a;
    lambda_v = _lambda_v;
}

// Discretize the mesh:
void INTERFACE::discretize(unsigned int disc1, unsigned int disc2) {
    char filename[200];
    if (disc1 == 0 && disc2 == 0)
        sprintf(filename, "infiles/final_configuration.dat");
    else
        sprintf(filename, "infiles/%d_%d.dat", disc1, disc2);
    ifstream in(filename, ios::in);
    if (!in) {
        if (world.rank() == 0)
            cout << "File could not be opened" << endl;
        exit(1);
    }

    // Vertices
    unsigned int vcol1;
    double vcol2, vcol3, vcol4;
    double cur_q;
    string dummy;
    in >> dummy >> dummy >> dummy >> dummy >> vcol1 >> dummy;
    number_of_vertices = vcol1;
    for (unsigned int i = 0; i < number_of_vertices; i++) {
        in >> vcol1 >> vcol2 >> vcol3 >> vcol4 >> cur_q;
        V.push_back(VERTEX(VECTOR3D(vcol2, vcol3, vcol4)));
        V[vcol1].index = vcol1;
        V[vcol1].q = cur_q;
    }
    assert(V.size() == number_of_vertices);

    // Edges
    unsigned int ecol1, ecol2, ecol3;
    in >> dummy >> dummy >> dummy >> dummy >> ecol1 >> dummy;
    number_of_edges = ecol1;
    E.resize(number_of_edges);
    for (unsigned int i = 0; i < number_of_edges; i++) {
        in >> ecol1 >> ecol2 >> ecol3 >> dummy;
        E[ecol1].index = ecol1;
        E[ecol1].itsV.push_back(&V[ecol2]);
        E[ecol1].itsV.push_back(&V[ecol3]);
        V[ecol2].itsE.push_back(&E[ecol1]);
        V[ecol3].itsE.push_back(&E[ecol1]);
    }

    // Faces
    unsigned int fcol1, fcol2, fcol3, fcol4;
    int sign;
    in >> dummy >> dummy >> dummy >> dummy >> fcol1 >> dummy;
    number_of_faces = fcol1;
    F.resize(number_of_faces);
    for (unsigned int i = 0; i < number_of_faces; i++) {
        in >> fcol1 >> fcol2 >> fcol3 >> fcol4 >> sign >> dummy;
        F[fcol1].index = fcol1;
        if (sign == 1) {
            F[fcol1].itsV.push_back(&V[fcol2]);
            F[fcol1].itsV.push_back(&V[fcol3]);
            F[fcol1].itsV.push_back(&V[fcol4]);
        } else {
            F[fcol1].itsV.push_back(&V[fcol4]);
            F[fcol1].itsV.push_back(&V[fcol3]);
            F[fcol1].itsV.push_back(&V[fcol2]);
        }
        V[fcol2].itsF.push_back(&F[fcol1]);
        V[fcol3].itsF.push_back(&F[fcol1]);
        V[fcol4].itsF.push_back(&F[fcol1]);
        EDGE *currentedge;
        currentedge = edge_between(F[fcol1].itsV[1], F[fcol1].itsV[2]);
        F[fcol1].itsE.push_back(currentedge);
        currentedge->itsF.push_back(&F[fcol1]);
        currentedge = edge_between(F[fcol1].itsV[2], F[fcol1].itsV[0]);
        F[fcol1].itsE.push_back(currentedge);
        currentedge->itsF.push_back(&F[fcol1]);
        currentedge = edge_between(F[fcol1].itsV[0], F[fcol1].itsV[1]);
        F[fcol1].itsE.push_back(currentedge);
        currentedge->itsF.push_back(&F[fcol1]);
    }

    for (unsigned int i = 0; i < V.size(); i++)
        assert(V[i].itsE.size() == V[i].itsF.size());        // check if e(v) = f(v)

    // compute the average edge length
    avg_edge_length = 0;
    for (unsigned int i = 0; i < E.size(); i++) {
        E[i].l0 = E[i].length();
        avg_edge_length += E[i].l0;
    }

    avg_edge_length /= E.size();

    //   for (unsigned int i = 0; i < E.size(); i++)
    //    E[i].l0 = avg_edge_length;


    // for local ref area
    for (unsigned int i = 0; i < F.size(); i++) {
        F[i].computearea();
    }

    for (unsigned int i = 0; i < V.size(); i++) {
        V[i].computearea();
        V[i].itsrefarea = V[i].itsarea;
    }

    total_area = 0;
    for (unsigned int i = 0; i < F.size(); i++)
        total_area += F[i].itsarea;

    if (world.rank() == 0) {
        ofstream listvertices("outfiles/initial_membrane.xyz", ios::out);
        listvertices << number_of_vertices << endl;
        listvertices << "interface" << endl;
        for (unsigned int k = 0; k < V.size(); k++)
            listvertices << "I" << setw(15) << V[k].posvec << endl;
        listvertices.close();

        ofstream filevertices("outfiles/initial_info_vertices.dat", ios::out);
        for (unsigned int i = 0; i < V.size(); i++) {
            filevertices << V[i].index << setw(15) << V[i].posvec << V[i].q << endl;
            for (unsigned int k = 0; k < V[i].itsE.size(); k++) {
                filevertices << setw(10) << (V[i].itsE[k])->index;
            }
            filevertices << endl;
            for (unsigned int k = 0; k < V[i].itsF.size(); k++) {
                filevertices << setw(10) << (V[i].itsF[k])->index;
            }
            filevertices << endl;
        }

        ofstream fileedges("outfiles/initial_info_edges.dat", ios::out);
        for (unsigned int i = 0; i < E.size(); i++) {
            fileedges << E[i].index << endl;
            for (unsigned int k = 0; k < E[i].itsV.size(); k++) {
                fileedges << setw(10) << (E[i].itsV[k])->index;
            }
            fileedges << endl;
            for (unsigned int k = 0; k < E[i].itsF.size(); k++) {
                fileedges << setw(10) << (E[i].itsF[k])->index;
            }
            fileedges << endl;
        }
        ofstream filefaces("outfiles/initial_info_faces.dat", ios::out);
        for (unsigned int i = 0; i < F.size(); i++) {
            filefaces << F[i].index << endl;
            for (unsigned int k = 0; k < F[i].itsV.size(); k++) {
                filefaces << setw(10) << (F[i].itsV[k])->index;
            }
            filefaces << endl;
            for (unsigned int k = 0; k < F[i].itsE.size(); k++) {
                filefaces << setw(10) << (F[i].itsE[k])->index;
            }
            filefaces << endl;
        }

        filevertices.close();
        fileedges.close();
        filefaces.close();
    }
    /*  int output_num = 1;
    for (unsigned int i=0; i<E.size(); i++)
    {
    //cout << "edge " << i << ": " << E[i].itsV[0]->index << " " << E[i].itsV[1]->index << "\n";
    //E[i].flipIfFavorable();
    //cout << "edge " << i << ": " << E[i].itsV[0]->index << " " << E[i].itsV[1]->index << "\n";
    stringstream ss;
    ss << "outfiles/edge" << output_num << ".dat";
    output_num++;
    ofstream fileedges(ss.str().c_str(), ios::out);
    fileedges << "\n\n--------------------------------------------------------------------------------\n\n";
    for (unsigned int j = 0; j < E.size(); j++)
    {
    fileedges << E[j].index << endl;
    for (unsigned int k = 0; k < E[j].itsV.size(); k++)
    {
    fileedges << setw(10) << (E[j].itsV[k])->index;
    }
    fileedges << endl;
    for (unsigned int k = 0; k < E[j].itsF.size(); k++)
    {
    fileedges << setw(10) << (E[j].itsF[k])->index;
    }
    fileedges << endl;
    }

    //cout << sqrt(E[i].crossingLengthSquared()) / sqrt(E[i].lengthSquared()) << "\n";
    }
    //exit(1);
    */
    return;
}

//FS added 2021/1/20
//Compute V's duals from initial input data file;
void INTERFACE::assign_dual_initial() {
    for (unsigned int i = 0; i < number_of_faces; i++) {
        Dual.push_back(VERTEX(VECTOR3D((F[i].itsV[0]->posvec.x + F[i].itsV[1]->posvec.x + F[i].itsV[2]->posvec.x) / 3.0,
            (F[i].itsV[0]->posvec.y + F[i].itsV[1]->posvec.y + F[i].itsV[2]->posvec.y) / 3.0,
            (F[i].itsV[0]->posvec.z + F[i].itsV[1]->posvec.z + F[i].itsV[2]->posvec.z) / 3.0)));
        Dual[i].index = i;
        Dual[i].q = (F[i].itsV[0]->q + F[i].itsV[1]->q + F[i].itsV[2]->q) / 3.0;
    }
}

//reassign charges based on the ratio of total charges (q_original + q_dual) and previous charges q_original. 
void INTERFACE::reassign_charges() {
    double q_original = 0.0;
    double q_dual = 0.0;
    //Compute total charges (q_dual + q_orignal)
    for (unsigned int i = 0; i < number_of_vertices; i++) {
        q_original += V[i].q;
    }
    for (unsigned int i = 0; i < number_of_faces; i++) {
        q_dual += Dual[i].q;
    }

    for (unsigned int i = 0; i < number_of_vertices; i++) {
        V[i].q *= (q_original / (q_original + q_dual));
    }
    for (unsigned int i = 0; i < number_of_faces; i++) {
        Dual[i].q *= (q_original / (q_original + q_dual));
    }
}

void INTERFACE::reassign_pm_charges() {
    double q_original_positive = 0.0;
    double q_original_negative = 0.0;
    double q_dual_positive = 0.0;
    double q_dual_negative = 0.0;

    //Compute total charges (q_dual + q_orignal)
    for (unsigned int i = 0; i < number_of_vertices; i++) {
        if (V[i].q > 0)
            q_original_positive += V[i].q;
        if (V[i].q < 0)
            q_original_negative += V[i].q;
    }

    for (unsigned int i = 0; i < number_of_faces; i++) {
        if (Dual[i].q > 0)
            q_dual_positive += Dual[i].q;
        if (Dual[i].q < 0)
            q_dual_negative += Dual[i].q;
    }


    for (unsigned int i = 0; i < number_of_vertices; i++) {
        if (V[i].q > 0)
            V[i].q *= (q_original_positive / (q_original_positive + q_dual_positive));
        if (V[i].q < 0)
            V[i].q *= (q_original_negative / (q_original_negative + q_dual_negative));
    }
    for (unsigned int i = 0; i < number_of_faces; i++) {
        if (Dual[i].q > 0)
            Dual[i].q *= (q_original_positive / (q_original_positive + q_dual_positive));
        if (Dual[i].q < 0)
            Dual[i].q *= (q_original_negative / (q_original_negative + q_dual_negative));
    }

    double q_original_positive_ref = 0.0;
    double q_original_negative_ref = 0.0;
    double q_dual_positive_ref = 0.0;
    double q_dual_negative_ref = 0.0;
    //testing
    for (unsigned int i = 0; i < number_of_vertices; i++) {
        if (V[i].q > 0)
            q_original_positive_ref += V[i].q;
        if (V[i].q < 0)
            q_original_negative_ref += V[i].q;
    }

    for (unsigned int i = 0; i < number_of_faces; i++) {
        if (Dual[i].q > 0)
            q_dual_positive_ref += Dual[i].q;
        if (Dual[i].q < 0)
            q_dual_negative_ref += Dual[i].q;
    }

    cout << "positive total charges(include dual): "<< q_original_positive_ref + q_dual_positive_ref << endl;
    cout  << "negative total charges(include dual): "<< q_original_negative_ref + q_dual_negative_ref << endl;

}




//  Load off-file for restarting pre-existing simulations (only works if there was no edge-flipping, same connectivity):
void INTERFACE::load_configuration(string filename) {
    ifstream in(filename.c_str(), ios::in);
    if (!in) {
        if (world.rank() == 0)
            cout << "File could not be opened" << endl;
        exit(1);
    }
    string dummy;
    unsigned int _nv, _ne, _nf;
    if (world.rank() == 0) {
        in >> dummy >> _nv >> _nf >> _ne;
        cout << V.size() << " " << E.size() << " " << F.size() << " ";
        cout << _nv << " " << _ne << " " << _nf << "\n";
    }
    assert(_nv == V.size() && _ne == E.size() && _nf == F.size());
    if (world.rank() == 0)
        for (unsigned int i = 0; i < V.size(); i++)
            in >> V[i].posvec.x >> V[i].posvec.y >> V[i].posvec.z;
}

EDGE *INTERFACE::edge_between(VERTEX *v1, VERTEX *v2) {
    for (unsigned int i = 0; i < v1->itsE.size(); i++) {
        if ((v1->itsE[i])->opposite(v1) == v2)
            return v1->itsE[i];
    }
    return NULL;
    //     if((*itsE[i]).opposite(v1)==v2)
}

void INTERFACE::assign_dual_boundary_edges() {
    for (unsigned int i = 0; i < number_of_edges; i++)
        E[i].isOnBoundary = 0;
    for (unsigned int i = 0; i < number_of_faces; i++) {
        if ((F[i].itsV[0]->q != F[i].itsV[1]->q) ||
            (F[i].itsV[0]->q != F[i].itsV[2]->q)) {
            // use the fact that edges and vertices are stored cyclically
            if (F[i].itsV[1]->q != F[i].itsV[2]->q)
                F[i].itsE[0]->isOnBoundary = 1;
            if (F[i].itsV[2]->q != F[i].itsV[0]->q)
                F[i].itsE[1]->isOnBoundary = 1;
            if (F[i].itsV[0]->q != F[i].itsV[1]->q)
                F[i].itsE[2]->isOnBoundary = 1;
        }
    }
}

//  NB added function for pH (to distribute charge randomly); uniform for a = 1.0, equivalent to legacy but shuffled.
void INTERFACE::assign_random_q_values(double q_strength, double alpha, int num_divisions, double fracChargedPatch, char randomFlag, char functionFlag) {
    vector<pair<double, int> > permutations;
    for (unsigned int i = 0; i < number_of_vertices; i++)
        permutations.push_back(pair<double, int>(V[i].posvec.z, i));
    sort(permutations.begin(), permutations.end());
    assert(num_divisions >= 1 && num_divisions <= 17); // Verify the prescribed number of divisions is supported.
    if (number_of_vertices % num_divisions != 0)
        if (world.rank() == 0)
            cout << "Warning:  The number of vertices modulo the number of divisions is not zero; the closest approximation will be used." << endl;
		  
        // Define a vector containing the to-be-randomized charge occupancy information:
    vector<int> chargeStateList;
    // Begin the list with values representing charged vertices:
    for (unsigned int i = 1; i <= alpha * number_of_vertices; i++) {
        chargeStateList.push_back(1);
    }
    // Finish the list with values representing uncharged vertices:
    for (unsigned int i = chargeStateList.size(); i < number_of_vertices; i++) {
        chargeStateList.push_back(0);
    }

    // Previously, permutations & mismatched indices were used for pseudorandom scrambling.  This enforces a more random one.
    vector<double> randomAreaList;
    for (unsigned int i = 0; i < number_of_vertices; i++) {
        randomAreaList.push_back(V[i].itsarea); // Create the to-be-randomized vertex area list:
    }

    // Randomize both the charge state list (for pH studies) and vertex area list (for normalization of all methods):
    if (randomFlag == 'y'){
        srand(123456);
        random_shuffle(chargeStateList.begin(), chargeStateList.end());
        random_shuffle(randomAreaList.begin(), randomAreaList.end());
    }

    unsigned int nVertPerPatch = (fracChargedPatch) * V.size();


    if (num_divisions == 1)              // one section: uniformly charged(default); cube formation if functionFlag = "c"
    {
	if (functionFlag == 'c'){
	    for (unsigned int i = 0; i < number_of_vertices; i++){
	      if((V[permutations[i].second].posvec.x * V[permutations[i].second].posvec.x + V[permutations[i].second].posvec.y * V[permutations[i].second].posvec.y <= 0.25) || (V[permutations[i].second].posvec.y * V[permutations[i].second].posvec.y + V[permutations[i].second].posvec.z * V[permutations[i].second].posvec.z <= 0.25) || (V[permutations[i].second].posvec.x * V[permutations[i].second].posvec.x + V[permutations[i].second].posvec.z * V[permutations[i].second].posvec.z <= 0.25)){
		    V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;	
	       }
	       else{
	           V[permutations[i].second].q = 0;
	       }
	    }
	}
	else {
	    for (unsigned int i = 0; i < number_of_vertices; i++)
		{
                     if (chargeStateList[i] == 1){
             		V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;}
                     else{
               		V[permutations[i].second].q = 0;
		  }
	      }
	}
    }
    if (num_divisions == 2) {           //  two patch Janus, specifiable fractional coverage
		 unsigned int i;
        for (i = 0; i < nVertPerPatch; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = 0;
	if (functionFlag == 'y'){
		for (i = 0; i < number_of_vertices; i++) {
		if (V[permutations[i].second].posvec.z <= 0 && (V[permutations[i].second].posvec.x-0.5) *( V[permutations[i].second].posvec.x-0.5)+V[permutations[i].second].posvec.z*V[permutations[i].second].posvec.z <= 0.25 ){
			V[permutations[i].second].q = 0;

		}
		if (V[permutations[i].second].posvec.z >= 0 && (V[permutations[i].second].posvec.x+0.5) *( V[permutations[i].second].posvec.x+0.5)+V[permutations[i].second].posvec.z*V[permutations[i].second].posvec.z <= 0.25){
			V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
		}
	}
	}
	//code for plus and minus
    /*
	    for (i = 0; i < nVertPerPatch; i++)
            V[permutations[i].second].q = 0.12 * randomAreaList[i];
        for (; i < (number_of_vertices - nVertPerPatch ); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = - 0.12 * randomAreaList[i];
    */

    }
    if (num_divisions == 3) {           //  n patch striped, approximately equal areas each about z-axis
		 unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
    }
    if (num_divisions == 4) {           //  n patch striped, approximately equal areas each about z-axis
		 unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = 0;

    }
    if (num_divisions == 5) {           //  n patch striped, approximately equal areas each about z-axis
		 unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
    }
    if (num_divisions == 6) {           //  n patch striped, approximately equal areas each about z-axis
		 unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = 0;
    }
    if (num_divisions == 7) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
    }
    if (num_divisions == 8) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = 0;
    }
    if (num_divisions == 9) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
    }
    if (num_divisions == 10) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 9*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = 0;
    }
    if (num_divisions == 11) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 9*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 10*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
    }
    if (num_divisions == 12) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 9*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 10*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 11*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = 0;
    }
    if (num_divisions == 13) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 9*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 10*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 11*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 12*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
    }
    if (num_divisions == 14) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 9*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 10*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 11*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 12*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 13*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = 0;
    }
    if (num_divisions == 15) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 9*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 10*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 11*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 12*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 13*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 14*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
    }
    if (num_divisions == 16) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 9*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 10*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 11*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 12*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 13*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 14*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 15*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = 0;
    }
    if (num_divisions == 17) {           //  n patch striped, approximately equal areas each about z-axis
        unsigned int i;
        for (i = 0; i < number_of_vertices * (alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 2*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 3*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 4*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 5*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 6*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 7*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 8*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 9*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 10*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 11*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 12*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 13*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 14*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices * 15*(alpha / num_divisions); i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
        for (; i < number_of_vertices * 16*(alpha / num_divisions); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
    }
    assign_dual_boundary_edges();

    //  Scale the vertices' charges to achieve the target net charge exactly:
    if (q_strength == 0){
        for (unsigned int i = 0; i < V.size(); i++) {
            V[i].q = 0;
        }
    }
    else {
        //  Assess total actual charge on the membrane (which may be less than the desired amount due to shuffling):
        double q_actual = 0.;
        for (unsigned int i = 0; i < V.size(); i++) {
            q_actual += V[i].q;
        }

        //  Assess the target total charge, to scale charged vertices by to achieve it:
        double q_target;
        if (num_divisions == 1) q_target = round(q_strength); // Round ensures electroneutrality with counterions.
        else if (num_divisions == 2) q_target = round(fracChargedPatch * q_strength);
        else if (num_divisions >= 3) q_target = round(ceil(0.5 * num_divisions) * fracChargedPatch * q_strength);
            // 'Ceil' necessary as number of charged patches is always >= number of uncharged patches by design above.

        //  Scale the vertices' charges to achieve the target net charge exactly:
        for (unsigned int i = 0; i < V.size(); i++) {
            V[i].q = V[i].q * (q_target / q_actual);
        }
    }

/*	if (0)
	{
		for (i = 0; i<number_of_vertices; i++)
			cout << i << " " << V[i].q << "\n";
		for (i = 0; i<number_of_edges; i++)
			cout << setw(2) << E[i].itsV[0]->index << " "
			<< setw(2) << E[i].itsV[1]->index << " "
			<< setw(2) << E[i].isOnBoundary << "\n";
		for (i = 0; i<number_of_faces; i++)
		{
			cout << setw(2) << F[i].itsV[0]->index << " "
				<< setw(2) << F[i].itsV[1]->index << " "
				<< setw(2) << F[i].itsV[2]->index << "\n";
		}
	}*/
}

void INTERFACE::assign_random_plusminus_values(double sigma, double radius, int num_divisions, double fracChargedPatch, char randomFlag, char functionFlag) {
    vector<pair<double, int> > permutations;
    for (unsigned int i = 0; i < number_of_vertices; i++)
        permutations.push_back(pair<double, int>(V[i].posvec.z, i));
    sort(permutations.begin(), permutations.end());
    //assert(num_divisions == 2); // Verify the prescribed number of divisions is supported.
    if (number_of_vertices % num_divisions != 0)
        if (world.rank() == 0)
            cout << "Warning:  The number of vertices modulo the number of divisions is not zero; the closest approximation will be used." << endl;

    // Previously, permutations & mismatched indices were used for pseudorandom scrambling.  This enforces a more random one.
    vector<double> randomAreaList;
    for (unsigned int i = 0; i < number_of_vertices; i++) {
        randomAreaList.push_back(V[i].itsarea); // Create the to-be-randomized vertex area list:
    }


    unsigned int nVertPerPatch = (fracChargedPatch)*V.size();
    double side_charges = sigma * (4 * pi * radius * radius) * fracChargedPatch;

    if (num_divisions == 2) {           //  two patch Janus, specifiable fractional coverage
        unsigned int i;
        //code for plus and minus
        for (i = 0; i < nVertPerPatch; i++)
            V[permutations[i].second].q = sigma * randomAreaList[i];
        for (; i < (number_of_vertices - nVertPerPatch); i++)
            V[permutations[i].second].q = 0;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = -1.0 * sigma * randomAreaList[i];
            //V[permutations[i].second].q = 0;
    }

    if (functionFlag == 'c') {
        for (unsigned int i = 0; i < number_of_vertices; i++) {
            if ((V[permutations[i].second].posvec.x * V[permutations[i].second].posvec.x + V[permutations[i].second].posvec.y * V[permutations[i].second].posvec.y <= 0.25) || (V[permutations[i].second].posvec.y * V[permutations[i].second].posvec.y + V[permutations[i].second].posvec.z * V[permutations[i].second].posvec.z <= 0.25) || (V[permutations[i].second].posvec.x * V[permutations[i].second].posvec.x + V[permutations[i].second].posvec.z * V[permutations[i].second].posvec.z <= 0.25)) {
                V[permutations[i].second].q = sigma * randomAreaList[i] * radius * radius;
            }
            else {
                V[permutations[i].second].q = 0;
            }
        }
    }
    else {
        //  Scale the vertices' charges to achieve the target net charge exactly:
        if (sigma == 0) {
            for (unsigned int i = 0; i < V.size(); i++) {
                V[i].q = 0;
            }
        }
        else {
            //  Assess total actual charge on the membrane (which may be less than the desired amount due to shuffling):
            double q_positive_actual = 0.;
            double q_negative_actual = 0.;
            for (unsigned int i = 0; i < nVertPerPatch; i++) {
                q_positive_actual += V[permutations[i].second].q;
            }

            for (unsigned int i = number_of_vertices - nVertPerPatch; i < number_of_vertices; i++) {
                q_negative_actual += V[permutations[i].second].q;
            }


            //  Assess the target total charge, to scale charged vertices by to achieve it:
            double q_target;
            q_target = round(side_charges); // Round ensures electroneutrality with counterions.

        //  Scale the vertices' charges to achieve the target net charge exactly:
            for (unsigned int i = 0; i < nVertPerPatch; i++) {
                V[permutations[i].second].q = V[permutations[i].second].q * (q_target / q_positive_actual);
            }

            for (unsigned int i = number_of_vertices - nVertPerPatch; i < number_of_vertices; i++) {
                V[permutations[i].second].q = V[permutations[i].second].q * (-1.0 * q_target / q_negative_actual);
            }
        }
    }
 
}




//  NB added function to import charge values/patterns from a file:
void INTERFACE::assign_external_q_values(double q_strength, string externalPattern) {
    // Read data from a single specified file:
    stringstream fileNameStream;
    fileNameStream << "infiles/Charge_Patterns/" << externalPattern << ".dat";
    string fileName = fileNameStream.str();
    int col1;
    double col2;
    ifstream inStream(fileName, ios::in);
    if (!inStream)    // Verify the file could be opened.
    {
        cout << "Charge assignment file could not be opened.  Mesh will be uncharged." << endl;
        for(unsigned int i = 0; i < V.size(); i++) {
            V[i].q = 0;
        }
    } else
    {
        cout << "Charge assignment file opened successfully.  Charging the mesh." << endl;
        while(inStream >> col1 >> col2) {
            V[col1].q = (q_strength / total_area)*col2;
        }
    }

    //  Scale the vertices' charges to achieve an integer charge (for if ions are present):
    if (q_strength == 0){
        for (unsigned int i = 0; i < V.size(); i++) {
            V[i].q = 0;
        }
    }
    else {
        //  Assess total actual charge on the membrane (which may be less than the desired amount due to shuffling):
        double q_actual = 0.;
        for (unsigned int i = 0; i < V.size(); i++) {
            q_actual += V[i].q;
        }

        //  Assess the target total charge, to scale charged vertices by to achieve it:
        double q_target = round(q_actual);
        // 'Ceil' necessary as number of charged patches is always >= number of uncharged patches by design above.

        //  Scale the vertices' charges to achieve the target net charge exactly:
        for (unsigned int i = 0; i < V.size(); i++) {
            V[i].q = V[i].q * (q_target / q_actual);
        }
    }
}

//  Output information on the initial state of the NP mesh:
void INTERFACE::output_configuration() {

    ofstream list_config("outfiles/final_configuration.dat", ios::out);
    if (world.rank() == 0) {
        list_config << "# Number of Vertices= " << this->number_of_vertices << "\n";
        list_config << "\nvertices index pos pos pos charge\n\n";
        for (unsigned int i = 0; i < number_of_vertices; i++) {
            list_config << V[i].index << "\t"
                        << V[i].posvec.x << "\t"
                        << V[i].posvec.y << "\t"
                        << V[i].posvec.z << "\t"
                        << V[i].q << "\n";
        }
        list_config << "\n\n\n# Number of Edges= " << this->number_of_edges << "\n";
        list_config << "\nedges: index vertex_index vertex_index length\n\n";
    }
    for (unsigned int i = 0; i < number_of_edges; i++) {
        assert(E[i].itsV.size() == 2);
        if (world.rank() == 0)
            list_config << E[i].index << "\t"
                        << E[i].itsV[0]->index << "\t"
                        << E[i].itsV[1]->index << "\t"
                        << E[i].length() << "\t"
                        << "1\n";
    }
    if (world.rank() == 0) {
        list_config << "\n\n\n# Number of Faces= " << this->number_of_faces << "\n";
        list_config << "\nfaces\n\n";
    }
    for (unsigned int i = 0; i < number_of_faces; i++) {
        assert(F[i].itsV.size() == 3);
        if (world.rank() == 0)
            list_config << F[i].index << "\t"
                        << F[i].itsV[0]->index << "\t"
                        << F[i].itsV[1]->index << "\t"
                        << F[i].itsV[2]->index << "\t"
                        << "1\t1\n";
    }
}

// %%% Counterion Setup Functionality: %%%

void INTERFACE::put_counterions(double q_actual, double unit_radius_sphere, double counterion_diameter, double box_halflength, vector<PARTICLE> &counterions, int counterion_valency, bool counterion_flag) {

// % Populate the list of counterions:
// Verify the requested valency counterions can actually enforce electroneutrality (within precision):
assert(int(round(q_actual)) % counterion_valency == 0);
// Compute the total number needed for ideal integer valency:
unsigned int total_counterions = round(abs(q_actual / counterion_valency));
if (total_counterions % 2 != 0) total_counterions += 1;
// NP-Counterion touching distance:
double r0 = (2.5 + 0.5 * counterion_diameter);      //Assuming the biggest radius is based on the longest axis of the deformed shape which lambda = 2.5.
// Box-Counterion touching distance (spherical box):
double r0_box = box_halflength - 0.5 * counterion_diameter;
double counterion_count = 0.0;

UTILITY ugsl;

// Generate counterions in the box:
while (counterions.size() != total_counterions) {
    double x = gsl_rng_uniform(ugsl.r);
    x = (1 - x) * (-r0_box) + x * (r0_box);
    double y = gsl_rng_uniform(ugsl.r);
    y = (1 - y) * (-r0_box) + y * (r0_box);
    double z = gsl_rng_uniform(ugsl.r);
    z = (1 - z) * (-r0_box) + z * (r0_box);
    VECTOR3D posvec = VECTOR3D(x, y, z);
    if (posvec.GetMagnitude() < r0 + counterion_diameter) // Avoid placing counterions within the NP.
        continue;
    if (posvec.GetMagnitude() >= box_halflength - counterion_diameter) // Avoid placing counterions outside of the box.
        continue;
    bool continuewhile = false;
    for (unsigned int i = 0; i < counterions.size() && !continuewhile; i++)
        if ((posvec - counterions[i].posvec).GetMagnitude() <= counterion_diameter)
            continuewhile = true;        // Avoid overlapping with existing ions.
        if (continuewhile)
            continue;
    if (counterion_flag){
        if (counterion_count < (total_counterions / 2.0)) {
            counterions.push_back(PARTICLE(int(counterions.size()) + 1, counterion_diameter, counterion_valency, counterion_valency * -1.0, 1.0, posvec));
            counterion_count += 1.0;
        }
        else {
            counterions.push_back(PARTICLE(int(counterions.size()) + 1, counterion_diameter, counterion_valency, counterion_valency * 1.0, 1.0, posvec));
        }
    }
    else {
        counterions.push_back(PARTICLE(int(counterions.size()) + 1, counterion_diameter, counterion_valency, counterion_valency * -1.0, 1.0, posvec));
    }
}
if (world.rank() == 0) {
ofstream listcounterions("outfiles/counterions.xyz", ios::out);
listcounterions << counterions.size() << endl;
listcounterions << "counterions" << endl;

for (unsigned int i = 0; i < counterions.size(); i++)
listcounterions << "C" << setw(15) << counterions[i].posvec.x<< setw(15) << counterions[i].posvec.y
<< setw(15) << counterions[i].posvec.z << setw(15)
<< counterions[i].posvec.GetMagnitude() << endl;
listcounterions.close();
}

// Verify electroneutrality to within machine and type precision:
if (!counterion_flag) {
    assert(round(q_actual) == (counterion_valency * counterions.size()));
}

return;
}

// %%% Run-time & Post-processing Functionality: %%%

// Compute membrane-wide energies at a given step (num), component-wise {totKE, BE, SE, TE, totLJ, totES, VE}:
// This produces the "energy_in_parts" output file and calculates the quantities used in "energy_nanomembrane".
void INTERFACE::compute_energy(int num, vector<PARTICLE> &counterions, const double scalefactor, char bucklingFlag) {

    // Initialize & compute net KE for both the mesh and counterions:
    netMeshKE = 0;
    for (unsigned int i = 0; i < V.size(); i++)
        netMeshKE += V[i].ke;
    netIonKE = 0;
    for (unsigned int i = 0; i < counterions.size(); i++)
        netIonKE += counterions[i].ke;
    // Compute the total for output:
    double totalKE = netMeshKE + netIonKE;

/*	// linear area
	penergy = lambda_a * (total_area);
	//   penergy += lambda_v * (total_volume - ref_volume) * (total_volume - ref_volume);
	penergy += lambda_v * (total_volume);*/

    // Initiate the potential energy:
    penergy = 0;

    // Initiate the output stream:
    ofstream output("outfiles/energy_in_parts_kE_bE_sE_tE_ljE_esE.dat", ios::app);
    if (world.rank() == 0) {
        output << num << " " << totalKE << " ";
    }
    // Initialize, compute, & output the net bending energy (mesh-only):
    double benergy = 0;
    for (unsigned int i = 0; i < E.size(); i++) {
        benergy += bkappa * (1 - E[i].itsS
                                 / (4 * E[i].itsF[0]->itsarea * E[i].itsF[1]->itsarea));
    }
    penergy += benergy;
    if (world.rank() == 0)
        output << benergy << " ";

    // Initialize, compute, and output the net stretching energy (mesh-only):
    double senergy = 0;
    if(bucklingFlag == 'n'){
        for (unsigned int i = 0; i < E.size(); i++) {
            double stretched = (E[i].length() - E[i].l0);
            senergy += 0.5 * sconstant * stretched * stretched / (E[i].l0 * E[i].l0);
        }
        senergy = senergy * avg_edge_length * avg_edge_length;
    }
    else {
        for (unsigned int i = 0; i < E.size(); i++) {
            //double stretched = (E[i].length() - E[i].l0);
            double stretched = (E[i].length() - avg_edge_length);
            senergy += 0.5 * sconstant * stretched * stretched;
        }
    }
    penergy += senergy;
    if (world.rank() == 0)
        output << senergy << " ";

    // Initialize, compute, and output the net tension energy (mesh-only):
    double TEnergy = 0;
    //TEnergy = (sigma_a * (((total_Area_Vertices - ref_Area_Vertices) * (total_Area_Vertices - ref_Area_Vertices)) / ref_Area_Vertices)); // Quadratic in area difference from sphere.
    //TEnergy = (sigma_a * total_Area_Vertices * total_Area_Vertices / ref_Area_Vertices);		// Quadratic in absolute area.
    TEnergy = (sigma_a * total_area);    // Linear in absolute area.
    penergy += TEnergy;
    if (world.rank() == 0)
        output << TEnergy << " ";

    // Initialize, compute, and (later) output the net volume tension energy (mesh-only):
    double VolTEnergy = 0;
    VolTEnergy = (sigma_v * ((total_volume - ref_volume) * (total_volume - ref_volume))); // Quadratic in volume difference from sphere.
    //VolTEnergy = (sigma_v * (total_volume  * total_volume));		// Quadratic in absolute volume.
    //VolTEnergy = (sigma_v * total_volume);    // Linear in absolute volume.
    penergy += VolTEnergy; // Output below to occur as last column.

    // Initialize, compute, and output the net {LJ, ES} energies (mesh & counterions):
    //Common MPI Message objects
    vector<double> netMeshLJ(sizFVecMesh, 0.0);
    vector<double> netMeshES(sizFVecMesh, 0.0);
    vector<double> netIonsLJ(sizFVecIons,0.0);
    vector<double> netIonsES(sizFVecIons,0.0);
    double lj_total = 0,lj_totalT = 0;
    double es_total = 0,es_totalT = 0;
    int i=0, j=0;

    // The mesh-mesh and mesh-ion components:
    #pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = lowerBoundMesh; i <= upperBoundMesh; i++) {
        double mesh_mesh_net_LJ_temp = 0;
        double mesh_mesh_net_ES_temp = 0;
        double mesh_ions_net_LJ_temp = 0;
        double mesh_ions_net_ES_temp = 0;
        for (j = 0; j < V.size(); j++) {
            if (i == j) continue;
            mesh_mesh_net_LJ_temp += energy_lj_pair(V[i].posvec, V[j].posvec, lj_length_mesh_mesh, elj);
            mesh_mesh_net_ES_temp += energy_es_pair(V[i], V[j], em, inv_kappa_out, scalefactor);
        }
        mesh_mesh_net_LJ_temp = 0.5 * mesh_mesh_net_LJ_temp; // Cancel out double-counting in the mesh-mesh.
        mesh_mesh_net_ES_temp = 0.5 * mesh_mesh_net_ES_temp;
        for (j = 0; j < counterions.size(); j++) {
            mesh_ions_net_LJ_temp += energy_lj_pair(V[i].posvec, counterions[j].posvec, lj_length_mesh_ions, elj);
            mesh_ions_net_ES_temp += energy_es_pair(V[i], counterions[j], em, inv_kappa_out, scalefactor);
        }
        netMeshLJ[i - lowerBoundMesh] = mesh_mesh_net_LJ_temp + mesh_ions_net_LJ_temp;
        netMeshES[i - lowerBoundMesh] = mesh_mesh_net_ES_temp + mesh_ions_net_ES_temp;
    }

    for (i = lowerBoundMesh; i <= upperBoundMesh; i++) {
        lj_total += netMeshLJ[i - lowerBoundMesh];
        es_total += netMeshES[i - lowerBoundMesh];
    }

    // The ion-ion component:
    #pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        double ion_ion_net_LJ_temp = 0;
        double ion_ion_net_ES_temp = 0;
        for (j = 0; j < counterions.size(); j++) {
            if (i == j) continue;
            ion_ion_net_LJ_temp += energy_lj_pair(counterions[i].posvec, counterions[j].posvec, counterions[0].diameter, elj);
            ion_ion_net_ES_temp += energy_es_pair(counterions[i], counterions[j], em, inv_kappa_out, scalefactor);
        }
        netIonsLJ[i - lowerBoundIons] = 0.5 * ion_ion_net_LJ_temp ;
        netIonsES[i - lowerBoundIons] = 0.5 * ion_ion_net_ES_temp;
    }

    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        lj_total += netIonsLJ[i - lowerBoundIons];
        es_total += netIonsES[i - lowerBoundIons];
    }

    //MPI Operations
    if (world.size() > 1) {
        all_reduce(world, lj_total, lj_totalT,std::plus<double>());
        all_reduce(world, es_total, es_totalT,std::plus<double>());
    }else{
        lj_totalT=lj_total;
        es_totalT=es_total;
    }
    penergy += lj_totalT + es_totalT;

    /*
    if (world.rank() == 0) {
        if (num == 0)
            cout << "Initial ES energy: " << es_total << endl
                 << endl;  // boundary.PE(0) = ES(0) + LJ(0) = ES(0) as LJ_0 = 0.
        else
            cout << "ES energy: " << es_total << endl << endl;
    }
     */

    //penergy += lj_total + es_total;
    if (world.rank() == 0)
        output << lj_totalT << " " << es_totalT << " " << VolTEnergy << endl;

/*	// line tension energy

	//  double line_tension_total = 0;
	//  for (unsigned int i = 0; i < E.size(); i++)
	//    if (E[i].isOnBoundary)
	//      line_tension_total += lambda_l * E[i].length();

	//  penergy += line_tension_total;
	//  output << line_tension_total << "\n";*/

/*	double line_tension_total = 0;
	for (unsigned int i = 0; i < F.size(); i++)
		F[i].computecenter();
	for (unsigned int i = 0; i < E.size(); i++)
		if (E[i].isOnBoundary)
			line_tension_total += lambda_l * (E[i].itsF[0]->itscenter - E[i].itsF[1]->itscenter).GetMagnitude();
	penergy += line_tension_total;
	//  output << line_tension_total << "\n";*/

    // Record (for output) the local net energy (excluding thermostats):
    energy = netMeshKE + netIonKE + penergy;
    return;
}

// Compute the spatial energetics profiles on the membrane (elastic, electrostatic), say for visualization:
void INTERFACE::compute_local_energies(const double scalefactor) {
    ofstream es_output("outfiles/local_electrostatic_E.off", ios::out);
    if (world.rank() == 0)
        es_output << "OFF\n" << V.size() << " "
                  << F.size() << " " << E.size() << "\n";
    vector<double> es_energy(F.size(), 0);
    for (unsigned int i = 0; i < V.size(); i++)
        for (unsigned int j = i + 1; j < V.size(); j++) {
            double cur_E = energy_es_pair(V[i], V[j], em, inv_kappa_out, scalefactor);
            for (unsigned int k = 0; k < V[i].itsF.size(); k++)
                es_energy[V[i].itsF[k]->index] += cur_E;
            for (unsigned int k = 0; k < V[j].itsF.size(); k++)
                es_energy[V[j].itsF[k]->index] += cur_E;
        }

    // Hone in on the minimum and maximum values of per-face electrostatic energies:
    double min_E = 1e100;
    double max_E = -1e100;
    for (unsigned int i = 0; i < F.size(); i++) {
        if (es_energy[i] < min_E)
            min_E = es_energy[i];
        if (es_energy[i] > max_E)
            max_E = es_energy[i];
    }
    // Output the minimum and maximum values for the scale bar (NB commented out normalization, see below):
    if (world.rank() == 0)
        cout << "electrostatic range" << "\t" << min_E << "\t" << max_E << endl;

    // Dump the vertex coordinates alone (for all vertices):
    if (world.rank() == 0)
        for (unsigned int i = 0; i < V.size(); i++)
            es_output << V[i].posvec.x << " "
                      << V[i].posvec.y << " "
                      << V[i].posvec.z << "\n";

    // Dump the per-face vertex indices & energy per-face with given scaling:
    for (unsigned int i = 0; i < F.size(); i++) {
        double r, g, b, a;
        //colormap((es_energy[F[i].index] - min_E) / (max_E - min_E), &r, &g, &b, &a); // This is Vikram's original, shifted & normalized.
        colormap(es_energy[F[i].index], &r, &g, &b, &a); // NB added to compute non-relative quantities with no offset.
        if (world.rank() == 0)
            es_output << "3 "
                      << F[i].itsV[0]->index << " "
                      << F[i].itsV[1]->index << " "
                      << F[i].itsV[2]->index << " "
                      << r << " "
                      << g << " "
                      << b << " "
                      << a << "\n";
    }

    ofstream el_output("outfiles/local_elastic_E.off", ios::out);
    if (world.rank() == 0)
        el_output << "OFF\n" << V.size() << " "
                  << F.size() << " " << E.size() << "\n";

    vector<double> elastic_energy(E.size(), 0);
    for (unsigned int i = 0; i < E.size(); i++) {
        double stretched = (E[i].length() - E[i].l0);
        double senergy = 0.5 * sconstant * stretched * stretched / (E[i].l0 * E[i].l0);
        senergy *= avg_edge_length * avg_edge_length;
        for (unsigned int k = 0; k < E[i].itsF.size(); k++)
            elastic_energy[E[i].itsF[k]->index] += senergy;

        double benergy = bkappa * (1 - E[i].itsS
                                       / (4 * E[i].itsF[0]->itsarea * E[i].itsF[1]->itsarea));
        for (unsigned int k = 0; k < E[i].itsF.size(); k++)
            elastic_energy[E[i].itsF[k]->index] += benergy;
    }
    if (world.rank() == 0)
        for (unsigned int i = 0; i < V.size(); i++)
            el_output << V[i].posvec.x << " "
                      << V[i].posvec.y << " "
                      << V[i].posvec.z << "\n";
    min_E = 1e100;
    max_E = -1e100;
    for (unsigned int i = 0; i < F.size(); i++) {
        if (elastic_energy[i] < min_E)
            min_E = elastic_energy[i];
        if (elastic_energy[i] > max_E)
            max_E = elastic_energy[i];
    }

    // output min and max values for scale bar
    // elastic
    if (world.rank() == 0)
        cout << "elastic range" << "\t" << min_E << "\t" << max_E << endl;

    for (unsigned int i = 0; i < F.size(); i++) {
        double r, g, b, a;
        //colormap((elastic_energy[F[i].index] - min_E) / (max_E - min_E), &r, &g, &b, &a); // This is VJ's original version, shifted & normalized.
        colormap(elastic_energy[F[i].index], &r, &g, &b,
                 &a);   // NB absolute value version.
        if (world.rank() == 0)
            el_output << "3 "
                      << F[i].itsV[0]->index << " "
                      << F[i].itsV[1]->index << " "
                      << F[i].itsV[2]->index << " "
                      << r << " "
                      << g << " "
                      << b << " "
                      << a << "\n";
    }

}

// Compute the spatial elastic energetics profiles per each component, separately:
void INTERFACE::compute_local_energies_by_component() {
    // Open the streams for both components:
    ofstream el_output_bending("outfiles/local_bending_E.off", ios::out);
    if (world.rank() == 0)
        el_output_bending << "OFF\n" << V.size() << " "
                          << F.size() << " " << E.size() << "\n";
    ofstream el_output_stretching("outfiles/local_stretching_E.off", ios::out);
    if (world.rank() == 0)
        el_output_stretching << "OFF\n" << V.size() << " "
                             << F.size() << " " << E.size() << "\n";

    // Iterate over edges to determine both components:
    vector<double> stretching_energy(E.size(), 0);
    vector<double> bending_energy(E.size(), 0);
    for (unsigned int i = 0; i < E.size(); i++) {
        double stretched = (E[i].length() - E[i].l0);
        double senergy = 0.5 * sconstant * stretched * stretched / (E[i].l0 * E[i].l0);
        senergy *= avg_edge_length * avg_edge_length;
        for (unsigned int k = 0; k < E[i].itsF.size(); k++)
            stretching_energy[E[i].itsF[k]->index] += senergy;

        double benergy = bkappa * (1 - E[i].itsS
                                       / (4 * E[i].itsF[0]->itsarea * E[i].itsF[1]->itsarea));
        for (unsigned int k = 0; k < E[i].itsF.size(); k++)
            bending_energy[E[i].itsF[k]->index] += benergy;
    }
    if (world.rank() == 0)
        for (unsigned int i = 0; i < V.size(); i++)
            el_output_stretching << V[i].posvec.x << " "
                                 << V[i].posvec.y << " "
                                 << V[i].posvec.z << "\n";
    if (world.rank() == 0)
        for (unsigned int i = 0; i < V.size(); i++)
            el_output_bending << V[i].posvec.x << " "
                              << V[i].posvec.y << " "
                              << V[i].posvec.z << "\n";
    // Determine the maximum and minimum value for both components:
    double min_E_stretching = 1e100;
    double max_E_stretching = -1e100;
    for (unsigned int i = 0; i < F.size(); i++) {
        if (stretching_energy[i] < min_E_stretching)
            min_E_stretching = stretching_energy[i];
        if (stretching_energy[i] > max_E_stretching)
            max_E_stretching = stretching_energy[i];
    }
    double min_E_bending = 1e100;
    double max_E_bending = -1e100;
    for (unsigned int i = 0; i < F.size(); i++) {
        if (bending_energy[i] < min_E_bending)
            min_E_bending = bending_energy[i];
        if (bending_energy[i] > max_E_bending)
            max_E_bending = bending_energy[i];
    }

    // Output min and max values for scale bars for each component:
    if (world.rank() == 0)
        cout << "stretching range" << "\t" << min_E_stretching << "\t" << max_E_stretching << endl;
    // stretching
    if (world.rank() == 0)
        cout << "bending range" << "\t" << min_E_bending << "\t" << max_E_bending << endl;

    // Output the final local profiles for each file:
    for (unsigned int i = 0; i < F.size(); i++) {
        double r, g, b, a;
        //colormap((stretching_energy[F[i].index] - min_E_stretching) / (max_E_stretching - min_E_stretching), &r, &g, &b, &a); // This is VJ's original version, shifted & normalized.
        colormap(stretching_energy[F[i].index], &r, &g, &b, &a);   // NB absolute value version.
        if (world.rank() == 0)
            el_output_stretching << "3 "
                                 << F[i].itsV[0]->index << " "
                                 << F[i].itsV[1]->index << " "
                                 << F[i].itsV[2]->index << " "
                                 << r << " "
                                 << g << " "
                                 << b << " "
                                 << a << "\n";
    }
    for (unsigned int i = 0; i < F.size(); i++) {
        double r, g, b, a;
        //colormap((bending_energy[F[i].index] - min_E_bending) / (max_E_bending - min_E_bending), &r, &g, &b, &a); // This is VJ's original version, shifted & normalized.
        colormap(bending_energy[F[i].index], &r, &g, &b, &a);   // NB absolute value version.
        if (world.rank() == 0)
            el_output_bending << "3 "
                              << F[i].itsV[0]->index << " "
                              << F[i].itsV[1]->index << " "
                              << F[i].itsV[2]->index << " "
                              << r << " "
                              << g << " "
                              << b << " "
                              << a << "\n";
    }
}

// ### Unused/Legacy Functions:  ###

void INTERFACE::assign_boundary_edges() {
    for (unsigned int i = 0; i < number_of_edges; i++)
        E[i].isOnBoundary = 0;
    for (unsigned int i = 0; i < number_of_faces; i++) {
        if ((F[i].itsV[0]->q != F[i].itsV[1]->q) ||
            (F[i].itsV[0]->q != F[i].itsV[2]->q)) {
            // use the fact that edges and vertices are stored cyclically
            if (F[i].itsV[1]->q == F[i].itsV[2]->q)
                F[i].itsE[0]->isOnBoundary = 1;
            if (F[i].itsV[2]->q == F[i].itsV[0]->q)
                F[i].itsE[1]->isOnBoundary = 1;
            if (F[i].itsV[0]->q == F[i].itsV[1]->q)
                F[i].itsE[2]->isOnBoundary = 1;
        }
    }
}

//void INTERFACE::assign_q_spots(int num_spots, int spot_size, double q_strength)
// {
// 	int    NChargedV = num_spots * spot_size;
//
//     assert(NChargedV <= V.size());                       // Verify it isn't > possible.
//
//     double qEach = q_strength / V.size();                // Assign charge in consistent fashion with pH.
// 	cout << "Number of charged vertices:  " << NChargedV << "\n"; // Note:  q_str = charge were all vertices charged.
//                                                                   // Net charge here = qEach*NChargedV .
// 	for (unsigned int i = 0; i<number_of_vertices; i++)
// 		V[i].q = 0;                                               // At first, assign all vertices zero charge.
//
//
//
//     int size = 0;
//     list<int> vertexQueue;                                    // The single entity "vector" to be charged.
//     vertexQueue.push_back(0);                                 // Adding an entry to the "vector".
//     while (!vertexQueue.empty() && size < spot_size)
//     {
//         int cur = vertexQueue.front();                        // Reference to the first element in the "vector".
//         vertexQueue.pop_front();                              // Deletes the first element.
//         if (V[cur].q == 0)                                    // Checks to see if vertex (q) has already been updated.
//         {
//             size++;
//             V[cur].q = qplus;
//             for (unsigned int j = 0; j<V[cur].itsE.size(); j++)
//                 vertexQueue.push_back(V[cur].itsE[j]->opposite(&V[cur])->index);
//         }
//     }
// 	//assign_boundary_edges();
// 	assign_dual_boundary_edges();
// }

//  NB has not changed this functional at all; it is exactly as originally received in May 2017.
void INTERFACE::assign_q_values(int num_divisions, double q_strength) {
    unsigned int i;
    vector<pair<double, int> > permutations;
    for (i = 0; i < number_of_vertices; i++)
        permutations.push_back(pair<double, int>(V[i].posvec.z, i));
    sort(permutations.begin(), permutations.end());
    assert(num_divisions >= 1 && num_divisions <= 4);
    double q = q_strength / number_of_vertices;
    if (num_divisions == 1)              // only one section: uniformly charged
    {
//     for (i=0; i<number_of_vertices; i++)
//       V[permutations[i].second].q = q;
        for (i = 0; i < number_of_vertices; i++)
            //V[permutations[i].second].q = q_strength * V[i].itsarea / total_area;  // Original with mismatched.
            V[permutations[i].second].q = q_strength * V[permutations[i].second].itsarea / total_area;  // NB matched.
    }
    if (num_divisions == 2) {
        assert(number_of_vertices % 2 == 0);
        for (i = 0; i < number_of_vertices / 2; i++)
            V[permutations[i].second].q = q;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = -q;
    }
    if (num_divisions == 3) {
        assert(number_of_vertices % 4 == 0);
        for (i = 0; i < number_of_vertices / 4; i++)
            V[permutations[i].second].q = q;
        for (; i < number_of_vertices * 3 / 4; i++)
            V[permutations[i].second].q = -q;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = q;
    }
    if (num_divisions == 4) {
        assert(number_of_vertices % 4 == 0);
        for (i = 0; i < number_of_vertices / 4; i++)
            V[permutations[i].second].q = q;
        for (; i < number_of_vertices / 2; i++)
            V[permutations[i].second].q = -q;
        for (; i < number_of_vertices * 3 / 4; i++)
            V[permutations[i].second].q = q;
        for (; i < number_of_vertices; i++)
            V[permutations[i].second].q = -q;
    }

    //assign_boundary_edges();
    assign_dual_boundary_edges();

    if (0) {
        for (i = 0; i < number_of_vertices; i++)
            cout << i << " " << V[i].q << "\n";
        for (i = 0; i < number_of_edges; i++)
            cout << setw(2) << E[i].itsV[0]->index << " "
                 << setw(2) << E[i].itsV[1]->index << " "
                 << setw(2) << E[i].isOnBoundary << "\n";
        for (i = 0; i < number_of_faces; i++) {
            cout << setw(2) << F[i].itsV[0]->index << " "
                 << setw(2) << F[i].itsV[1]->index << " "
                 << setw(2) << F[i].itsV[2]->index << "\n";
        }
    }
}
