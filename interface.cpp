// This file contains member functions for interface class

#include "interface.h"
#include "functions.h"
#include "newenergies.h"
#include <algorithm>
#include <list>
#include <sstream>

//extern double alpha;

double colormap(double x, double* r, double* g, double* b, double* a)
{
	(*r) = 0;
	(*g) = x;
	(*b) = 1 - x;
	(*a) = 1;
}

void INTERFACE::set_up()
{
	// useful combinations of different dielectric constants (inside and outside)
	em = 0.5 * (ein + eout);
	ed = (eout - ein) / (4 * pi);

	// useful length scales signifying competition between electrostatics and entropy
	lB_out = (lB_water * epsilon_water / eout) / unit_radius_sphere;

	return;
}

// dress up the interface with normals areas and volume elements
void INTERFACE::dressup(double _lambda_a, double _lambda_v)
{
	for (unsigned int i = 0; i < F.size(); i++)
	{
		F[i].computenormal();
		F[i].computearea();
		F[i].computevolume();
	}

	for (unsigned int i = 0; i < V.size(); i++)
	{
		V[i].computenormal();
		V[i].computearea();
	}

	total_area = 0;
	for (unsigned int i = 0; i < F.size(); i++)
		total_area += F[i].itsarea;

	total_volume = 0;
	for (unsigned int i = 0; i < F.size(); i++)
		total_volume += F[i].subtends_volume;

	total_Area_Vertices = 0;												// (2017.09.17 NB added.)  Computes current net area by vertices.
	for (unsigned int i = 0; i < V.size(); i++)
		total_Area_Vertices += V[i].itsarea;

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

/*
* This function reads an off file and moves all the vertices to the points
* specified in the off file.  This will only work if there was no edge
* flipping, as the connectivity is assumed to be identical to the initial
* connectivity (for simplicity).
*/
void INTERFACE::load_configuration(string filename)
{
	ifstream in(filename.c_str(), ios::in);
	if (!in)
	{
		cout << "File could not be opened" << endl;
		exit(1);
	}
	string dummy;
	int _nv, _ne, _nf;
	in >> dummy >> _nv >> _nf >> _ne;
	cout << V.size() << " " << E.size() << " " << F.size() << " ";
	cout << _nv << " " << _ne << " " << _nf << "\n";

	assert(_nv == V.size() && _ne == E.size() && _nf == F.size());
	for (int i = 0; i<V.size(); i++)
		in >> V[i].posvec.x >> V[i].posvec.y >> V[i].posvec.z;
}

// Discretize the interface:
void INTERFACE::discretize(unsigned int disc1, unsigned int disc2)
{
	char filename[200];
	if (disc1 == 0 && disc2 == 0)
		sprintf(filename, "infiles/final_configuration.dat");
	else
		sprintf(filename, "infiles/%d_%d.dat", disc1, disc2);
	ifstream in(filename, ios::in);
	if (!in)
	{
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
	for (unsigned int i = 0; i < number_of_vertices; i++)
	{
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
	for (unsigned int i = 0; i < number_of_edges; i++)
	{
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
	for (unsigned int i = 0; i < number_of_faces; i++)
	{
		in >> fcol1 >> fcol2 >> fcol3 >> fcol4 >> sign >> dummy;
		F[fcol1].index = fcol1;
		if (sign == 1)
		{
			F[fcol1].itsV.push_back(&V[fcol2]);
			F[fcol1].itsV.push_back(&V[fcol3]);
			F[fcol1].itsV.push_back(&V[fcol4]);
		}
		else
		{
			F[fcol1].itsV.push_back(&V[fcol4]);
			F[fcol1].itsV.push_back(&V[fcol3]);
			F[fcol1].itsV.push_back(&V[fcol2]);
		}
		V[fcol2].itsF.push_back(&F[fcol1]);
		V[fcol3].itsF.push_back(&F[fcol1]);
		V[fcol4].itsF.push_back(&F[fcol1]);
		EDGE* currentedge;
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
		assert(V[i].itsE.size() == V[i].itsF.size());		// check if e(v) = f(v)

															// compute the average edge length 
	avg_edge_length = 0;
	for (unsigned int i = 0; i < E.size(); i++)
	{
		E[i].l0 = E[i].length();
		avg_edge_length += E[i].l0;
	}

	avg_edge_length /= E.size();

	//   for (unsigned int i = 0; i < E.size(); i++)
	//    E[i].l0 = avg_edge_length;


	// for local ref area
	for (unsigned int i = 0; i < F.size(); i++)
	{
		F[i].computearea();
	}

	for (unsigned int i = 0; i < V.size(); i++)
	{
		V[i].computearea();
		V[i].itsrefarea = V[i].itsarea;
	}

	total_area = 0;
	for (unsigned int i = 0; i < F.size(); i++)
		total_area += F[i].itsarea;

	ofstream listvertices("outfiles/initial_membrane.xyz", ios::out);
	listvertices << number_of_vertices << endl;
	listvertices << "interface" << endl;
	for (unsigned int k = 0; k < V.size(); k++)
		listvertices << "I" << setw(15) << V[k].posvec << endl;
	listvertices.close();

	ofstream filevertices("outfiles/initial_info_vertices.dat", ios::out);
	for (unsigned int i = 0; i < V.size(); i++)
	{
		filevertices << V[i].index << setw(15) << V[i].posvec << V[i].q << endl;
		for (unsigned int k = 0; k < V[i].itsE.size(); k++)
		{
			filevertices << setw(10) << (V[i].itsE[k])->index;
		}
		filevertices << endl;
		for (unsigned int k = 0; k < V[i].itsF.size(); k++)
		{
			filevertices << setw(10) << (V[i].itsF[k])->index;
		}
		filevertices << endl;
	}

	ofstream fileedges("outfiles/initial_info_edges.dat", ios::out);
	for (unsigned int i = 0; i < E.size(); i++)
	{
		fileedges << E[i].index << endl;
		for (unsigned int k = 0; k < E[i].itsV.size(); k++)
		{
			fileedges << setw(10) << (E[i].itsV[k])->index;
		}
		fileedges << endl;
		for (unsigned int k = 0; k < E[i].itsF.size(); k++)
		{
			fileedges << setw(10) << (E[i].itsF[k])->index;
		}
		fileedges << endl;
	}
	ofstream filefaces("outfiles/initial_info_faces.dat", ios::out);
	for (unsigned int i = 0; i < F.size(); i++)
	{
		filefaces << F[i].index << endl;
		for (unsigned int k = 0; k < F[i].itsV.size(); k++)
		{
			filefaces << setw(10) << (F[i].itsV[k])->index;
		}
		filefaces << endl;
		for (unsigned int k = 0; k < F[i].itsE.size(); k++)
		{
			filefaces << setw(10) << (F[i].itsE[k])->index;
		}
		filefaces << endl;
	}

	filevertices.close();
	fileedges.close();
	filefaces.close();

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

EDGE* INTERFACE::edge_between(VERTEX* v1, VERTEX* v2)
{
	for (unsigned int i = 0; i < v1->itsE.size(); i++)
	{
		if ((v1->itsE[i])->opposite(v1) == v2)
			return v1->itsE[i];
	}
	return NULL;
	//     if((*itsE[i]).opposite(v1)==v2)
}

// void INTERFACE::assign_q_spots(int num_spots, int spot_size, double q_strength)
// {
// 	int    NChargedV = num_spots * spot_size;
//     
//     assert(NChargedV <= boundary.V.size());                       // Verify it isn't > possible.
//     
//     double qEach = q_strength / boundary.V.size();                // Assign charge in consistent fashion with pH.
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

void INTERFACE::assign_boundary_edges()
{
	for (unsigned int i = 0; i<number_of_edges; i++)
		E[i].isOnBoundary = 0;
	for (unsigned int i = 0; i<number_of_faces; i++)
	{
		if ((F[i].itsV[0]->q != F[i].itsV[1]->q) ||
			(F[i].itsV[0]->q != F[i].itsV[2]->q))
		{
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

void INTERFACE::assign_dual_boundary_edges()
{
	for (unsigned int i = 0; i<number_of_edges; i++)
		E[i].isOnBoundary = 0;
	for (unsigned int i = 0; i<number_of_faces; i++)
	{
		if ((F[i].itsV[0]->q != F[i].itsV[1]->q) ||
			(F[i].itsV[0]->q != F[i].itsV[2]->q))
		{
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

// NB has not changed this functional at all; it is exactly as originally received in May 2017.
void INTERFACE::assign_q_values(int num_divisions, double q_strength)
{
	unsigned int i;
	vector< pair<double,int> > permutations;
	for (i=0; i<number_of_vertices; i++)
		permutations.push_back(pair<double,int>(V[i].posvec.z, i));
	sort(permutations.begin(),permutations.end());
	assert(num_divisions >= 1 && num_divisions <=4);
	double q = q_strength / number_of_vertices;
	if (num_divisions == 1)              // only one section: uniformly charged
	{
//     for (i=0; i<number_of_vertices; i++)
//       V[permutations[i].second].q = q;
		for (i=0; i<number_of_vertices; i++)
			//V[permutations[i].second].q = q_strength * V[i].itsarea / total_area;  // Original with mismatched.
		  V[permutations[i].second].q = q_strength * V[permutations[i].second].itsarea / total_area;  // NB matched.
	}
	if (num_divisions == 2)
	{
		assert(number_of_vertices % 2 == 0);
		for (i=0; i<number_of_vertices/2; i++)
			V[permutations[i].second].q = q;
		for (; i<number_of_vertices; i++)
			V[permutations[i].second].q = -q;
	}
	if (num_divisions == 3)
	{
		assert(number_of_vertices % 4 == 0);
		for (i=0; i<number_of_vertices/4; i++)
			V[permutations[i].second].q = q;
		for (; i<number_of_vertices*3/4; i++)
			V[permutations[i].second].q = -q;
		for (; i<number_of_vertices; i++)
			V[permutations[i].second].q = q;
	}
	if (num_divisions == 4)
	{
		assert(number_of_vertices % 4 == 0);
		for (i=0; i<number_of_vertices/4; i++)
			V[permutations[i].second].q = q;
		for (; i<number_of_vertices/2; i++)
			V[permutations[i].second].q = -q;
		for (; i<number_of_vertices*3/4; i++)
			V[permutations[i].second].q = q;
		for (; i<number_of_vertices; i++)
			V[permutations[i].second].q = -q;
	}

	//assign_boundary_edges();
	assign_dual_boundary_edges();

	if (0)
	{
		for (i=0; i<number_of_vertices; i++)
			cout << i << " " << V[i].q << "\n";
		for (i=0; i<number_of_edges; i++)
			cout << setw(2) << E[i].itsV[0]->index << " "
					 << setw(2) << E[i].itsV[1]->index << " "
					 << setw(2) << E[i].isOnBoundary << "\n";
		for (i=0; i<number_of_faces; i++)
		{
			cout << setw(2) << F[i].itsV[0]->index << " "
					 << setw(2) << F[i].itsV[1]->index << " "
					 << setw(2) << F[i].itsV[2]->index << "\n";
		}
	}
}

//  NB added function for pH (to distribute charge randomly); uniform for a = 1.0, equivalent to above but shuffled.
void INTERFACE::assign_random_q_values(int num_divisions, double q_strength, double alpha)
{
	unsigned int i;
	vector< pair<double, int> > permutations;
	for (i = 0; i<number_of_vertices; i++)
		permutations.push_back(pair<double, int>(V[i].posvec.z, i));
	sort(permutations.begin(), permutations.end());
	assert(num_divisions >= 1 && num_divisions <= 3); // Verify the prescribed number of divisions is supported.
    if(number_of_vertices % num_divisions == 0) cout << "Warning:  The number of vertices modulo the number of divisions is not zero; the closest approximation will be used." << endl;
	double q = q_strength / number_of_vertices;

    // Define a vector containing the to-be-randomized charge occupancy information:
	vector<int> chargeStateList;
	// Begin the list with values representing charged vertices:
	for (i = 1; i <= alpha*number_of_vertices; i++)
	{
		chargeStateList.push_back(1);
	}
	// Finish the list with values representing uncharged vertices:
	for (i = chargeStateList.size(); i < number_of_vertices; i++)
	{
		chargeStateList.push_back(0);
	}
	// Shuffle the list so it is randomly distributed:
        srand(time(0));
	random_shuffle(chargeStateList.begin(), chargeStateList.end());
    
    // Previously, permutations & mismatched indices were used for pseudorandom scrambling.  This enforces a different and more random one each time.
    vector<double> randomAreaList;
    for(int i = 0; i < number_of_vertices; i++)
    {
        randomAreaList.push_back(V[i].itsarea); // Create the to-be-randomized vertex area list:
    }
          
    // Once generated, randomize the vertex area list (uses the same random seed based on time above):
    random_shuffle(randomAreaList.begin(), randomAreaList.end());
    
	if (num_divisions == 1)              // one section: uniformly charged
	{
		for (i = 0; i < number_of_vertices; i++)
			if(chargeStateList[i]==1)
				V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
			else
				V[permutations[i].second].q = 0;
	}
	if (num_divisions == 2)
	{
		for (i = 0; i<number_of_vertices / 2; i++)
			V[permutations[i].second].q = q_strength * V[i].itsarea / total_area;
		for (; i<number_of_vertices; i++)
			V[permutations[i].second].q = 0;
	}
	if (num_divisions == 3)
	{
		for (i = 0; i<number_of_vertices*(alpha/(num_divisions-1)); i++)
			V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
		for (; i<number_of_vertices * (1-alpha/(num_divisions-1)); i++)
			V[permutations[i].second].q = 0;
		for (; i<number_of_vertices; i++)
			V[permutations[i].second].q = q_strength * randomAreaList[i] / total_area;
	}

	//assign_boundary_edges();
	assign_dual_boundary_edges();
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

// Compute the spatial energetics profiles on the membrane (elastic, electrostatic), say for visualization:
void INTERFACE::compute_local_energies()
{
	ofstream es_output("outfiles/local_electrostatic_E.off", ios::out);
	es_output << "OFF\n" << V.size() << " "
		<< F.size() << " " << E.size() << "\n";
	vector<double> es_energy(F.size(), 0);
	for (unsigned int i = 0; i < V.size(); i++)
		for (unsigned int j = i + 1; j < V.size(); j++)
		{
			double cur_E = energy_es_vertex_vertex(V[i], V[j], em, inv_kappa_out);
			for (int k = 0; k<V[i].itsF.size(); k++)
				es_energy[V[i].itsF[k]->index] += cur_E;
			for (int k = 0; k<V[j].itsF.size(); k++)
				es_energy[V[j].itsF[k]->index] += cur_E;
		}
		
    // Hone in on the minimum and maximum values of per-face electrostatic energies:
	double min_E = 1e100;
	double max_E = -1e100;
	for (unsigned int i = 0; i < F.size(); i++)
	{
		if (es_energy[i] < min_E)
			min_E = es_energy[i];
		if (es_energy[i] > max_E)
			max_E = es_energy[i];
	}
	// Output the minimum and maximum values for the scale bar (NB commented out normalization, see below):
	cout << "electrostatic range" << "\t" << min_E << "\t" << max_E << endl;

    // Dump the vertex coordinates alone (for all vertices):
	for (unsigned int i = 0; i<V.size(); i++)
		es_output << V[i].posvec.x << " "
		<< V[i].posvec.y << " "
		<< V[i].posvec.z << "\n";

    // Dump the per-face vertex indices & energy per-face with given scaling:
	for (unsigned int i = 0; i<F.size(); i++)
	{
		double r, g, b, a;
		//colormap((es_energy[F[i].index] - min_E) / (max_E - min_E), &r, &g, &b, &a); // This is Vikram's original, shifted & normalized.
		colormap(es_energy[F[i].index], &r, &g, &b, &a); // NB added to compute non-relative quantities with no offset.
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
	el_output << "OFF\n" << V.size() << " "
		<< F.size() << " " << E.size() << "\n";
	vector<double> elastic_energy(E.size(), 0);
	for (unsigned int i = 0; i < E.size(); i++)
	{
		double stretched = (E[i].length() - E[i].l0);
		double senergy = 0.5*sconstant*stretched*stretched / (E[i].l0*E[i].l0);
		senergy *= avg_edge_length * avg_edge_length;
		for (int k = 0; k<E[i].itsF.size(); k++)
			elastic_energy[E[i].itsF[k]->index] += senergy;

		double benergy = bkappa * (1 - E[i].itsS
			/ (4 * E[i].itsF[0]->itsarea * E[i].itsF[1]->itsarea));
		for (int k = 0; k<E[i].itsF.size(); k++)
			elastic_energy[E[i].itsF[k]->index] += benergy;
	}

	for (unsigned int i = 0; i<V.size(); i++)
		el_output << V[i].posvec.x << " "
		<< V[i].posvec.y << " "
		<< V[i].posvec.z << "\n";
	min_E = 1e100;
	max_E = -1e100;
	for (unsigned int i = 0; i < F.size(); i++)
	{
		if (elastic_energy[i] < min_E)
			min_E = elastic_energy[i];
		if (elastic_energy[i] > max_E)
			max_E = elastic_energy[i];
	}

	// output min and max values for scale bar
	// elastic
	cout << "elastic range" << "\t" << min_E << "\t" << max_E << endl;

	for (unsigned int i = 0; i<F.size(); i++)
	{
		double r, g, b, a;
		//colormap((elastic_energy[F[i].index] - min_E) / (max_E - min_E), &r, &g, &b, &a); // This is VJ's original version, shifted & normalized.
		colormap((elastic_energy[F[i].index] - min_E) / (max_E - min_E), &r, &g, &b, &a);   // NB absolute value version.
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

// Compute the membrane-wide energies at a given step (num), component-wise {kinetic, BE, SE, TE, LJ, ES} respectively:
		// This produces the "energy_in_parts" output file and calculates the quantities used in "energy_nanomembrane".
void INTERFACE::compute_energy(int num)
{
	// Initialize & compute net kinetic energy:
	kenergy = 0;
	for (unsigned int i = 0; i < V.size(); i++)
		kenergy += V[i].Realke;

/*	// linear area
	penergy = lambda_a * (total_area);
	//   penergy += lambda_v * (total_volume - ref_volume) * (total_volume - ref_volume);
	penergy += lambda_v * (total_volume);*/

	// Initiate the potential energy:
	penergy = 0;

	// Initiate the output stream:
	ofstream output("outfiles/energy_in_parts_kE_bE_sE_tE_ljE_esE.dat", ios::app);
	output << num << " ";
	output << kenergy << " ";
	// NOTE: constraints are not explicit like below, so commented out the energies; also for teaching purpose
	//  output << lambda_a * (total_area) << " ";
	//  output <<  lambda_v * (total_volume - ref_volume) * (total_volume - ref_volume) << " ";

	// Initialize, compute, & output the net bending energy:
	double benergy = 0;
	for (unsigned int i = 0; i < E.size(); i++)
	{
		benergy += bkappa * (1 - E[i].itsS
														 / (4 * E[i].itsF[0]->itsarea * E[i].itsF[1]->itsarea));
	}
	penergy += benergy;
	output << benergy << " ";

	// Initialize, compute, and output the net stretching energy:
	double senergy = 0;
	for (unsigned int i = 0; i < E.size(); i++)
	{
		double stretched = (E[i].length() - E[i].l0);
		//     senergy += 0.5 * sconstant * stretched * stretched;
		senergy += 0.5 * sconstant * stretched * stretched / (E[i].l0 * E[i].l0);
	}
	senergy = senergy * avg_edge_length * avg_edge_length;
	penergy += senergy;
	output << senergy << " ";

	double TEnergy = 0;
	//TEnergy = (sigma_a * (((total_Area_Vertices - ref_Area_Vertices) * (total_Area_Vertices - ref_Area_Vertices)) / ref_Area_Vertices)); // Quadratic in area difference from sphere.
	//TEnergy = (sigma_a * total_Area_Vertices * total_Area_Vertices / ref_Area_Vertices);		// Quadratic in absolute area.
	TEnergy = (sigma_a * total_Area_Vertices);	// Linear in absolute area.
	penergy += TEnergy;

	output << TEnergy << " ";

	// Initialize, compute, and output the net {LJ, ES} energies:
	double lj_total = 0;
	double es_total = 0;
	for (unsigned int i = 0; i < V.size(); i++)
	{
		for (unsigned int j = i + 1; j < V.size(); j++)
		{
			lj_total += energy_lj_vertex_vertex(V[i], V[j], lj_length, elj);
			es_total += energy_es_vertex_vertex(V[i], V[j], em, inv_kappa_out);
    }
  }

  cout << "Initial ES energy: " << es_total << endl << endl;  // boundary.PE(0) = ES(0) + LJ(0) = ES(0) as LJ_0 = 0.

	penergy += lj_total + es_total;
	output << lj_total << " " << es_total << endl;

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

	energy = kenergy + penergy;
}

void INTERFACE::output_configuration()
{
	ofstream list_config("outfiles/final_configuration.dat", ios::out);
	list_config << "# Number of Vertices= " << this->number_of_vertices << "\n";
	list_config << "\nvertices index pos pos pos charge\n\n";
	for (unsigned int i = 0; i<number_of_vertices; i++)
	{
		list_config << V[i].index << "\t"
			<< V[i].posvec.x << "\t"
			<< V[i].posvec.y << "\t"
			<< V[i].posvec.z << "\t"
			<< V[i].q << "\n";
	}
	list_config << "\n\n\n# Number of Edges= " << this->number_of_edges << "\n";
	list_config << "\nedges: index vertex_index vertex_index length\n\n";
	for (unsigned int i = 0; i<number_of_edges; i++)
	{
		assert(E[i].itsV.size() == 2);
		list_config << E[i].index << "\t"
			<< E[i].itsV[0]->index << "\t"
			<< E[i].itsV[1]->index << "\t"
            << E[i].length() << "\t"
			<< "1\n";
	}
	list_config << "\n\n\n# Number of Faces= " << this->number_of_faces << "\n";
	list_config << "\nfaces\n\n";
	for (unsigned int i = 0; i<number_of_faces; i++)
	{
		assert(F[i].itsV.size() == 3);
		list_config << F[i].index << "\t"
			<< F[i].itsV[0]->index << "\t"
			<< F[i].itsV[1]->index << "\t"
			<< F[i].itsV[2]->index << "\t"
			<< "1\t1\n";
	}
}
