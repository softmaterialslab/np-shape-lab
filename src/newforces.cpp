#include "newforces.h"
#include "functions.h"

void force_calculation_init(INTERFACE &boundary, vector<PARTICLE> &counterions, const double scalefactor, char bucklingFlag) {
    unsigned int i=0, j=0;

// Compute constraint gradients: (called 'forces' in prev comment & in SHAKE, RATTLE.)
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].compute_neg_VGrads(&boundary/*, constraintForm*/);

    // Initialize forces on all vertices, ions:
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].forvec = VECTOR3D(0, 0, 0);
    for (i = 0; i < counterions.size(); i++)
        counterions[i].forvec = VECTOR3D(0, 0, 0);

    // Contribute LJ & electrostatics forces (initialization implicit):
    #pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = 0; i < boundary.V.size(); i++) {
        for (j = 0; j < boundary.V.size(); j++) { // The mesh-mesh component:
            if (i == j) continue;
            VECTOR3D ljForce;
            VECTOR3D esForce;
            VECTOR3D r_vec = boundary.V[i].posvec - boundary.V[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = boundary.lj_length_mesh_mesh * boundary.lj_length_mesh_mesh;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljForce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esForce = ((-scalefactor * boundary.V[i].q * boundary.V[j].q / boundary.em) *
                       (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1 / r2)) ^ ((-1.0 / r) - (1.0 / boundary.inv_kappa_out)));

            boundary.V[i].forvec += ljForce + esForce;
        }

        for (j = 0; j < counterions.size(); j++) { // The mesh-ion component (for mesh):
            VECTOR3D ljForce;
            VECTOR3D esForce;
            VECTOR3D r_vec = boundary.V[i].posvec - counterions[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = boundary.lj_length_mesh_ions * boundary.lj_length_mesh_ions;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljForce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esForce = ((-scalefactor * boundary.V[i].q * counterions[j].q / boundary.em) *
                       (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1 / r2)) ^ ((-1.0 / r) - (1.0 / boundary.inv_kappa_out)));

            boundary.V[i].forvec += ljForce + esForce;
        }
    }
    #pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = 0; i < counterions.size(); i++) {
        for (j = 0; j < counterions.size(); j++) { // The ion-ion component:
            if (i == j) continue;
            VECTOR3D ljForce;
            VECTOR3D esForce;
            VECTOR3D r_vec = counterions[i].posvec - counterions[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = counterions[0].diameter * counterions[0].diameter;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljForce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esForce = ((-scalefactor * counterions[i].q * counterions[j].q / boundary.em) *
                       (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1 / r2)) ^ ((-1.0 / r) - (1.0 / boundary.inv_kappa_out)));

            counterions[i].forvec += ljForce + esForce;
        }
        for (j = 0; j < boundary.V.size(); j++) { // The mesh-ion component:
            VECTOR3D ljForce;
            VECTOR3D esForce;
            VECTOR3D r_vec = counterions[i].posvec - boundary.V[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = boundary.lj_length_mesh_ions * boundary.lj_length_mesh_ions;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljForce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esForce = ((-scalefactor * counterions[i].q * boundary.V[j].q / boundary.em) *
                       (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1 / r2)) ^ ((-1.0 / r) - (1.0 / boundary.inv_kappa_out)));

            counterions[i].forvec += ljForce + esForce;
        }
    }

    // Initialize & contribute initial bending forces:
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].bforce = VECTOR3D(0, 0, 0);
    for (i = 0; i < boundary.E.size(); i++)
        boundary.E[i].bending_forces(&boundary);
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].forvec += boundary.V[i].bforce;
    // Contribute initial stretching forces (2017.09.17 NB added, initialization implicit):
    for (i = 0; i < boundary.V.size(); i++) {
        boundary.V[i].stretching_forces(&boundary, bucklingFlag);
        boundary.V[i].forvec += boundary.V[i].sforce;
    }
    // (2017.09.17 NB added.) Contribute initial tension forces:
    for (i = 0; i < boundary.V.size(); i++) { // Re-calculation of geometric gradients not required as positions haven't changed.
        boundary.V[i].tension_forces(&boundary);
        boundary.V[i].forvec += boundary.V[i].TForce;
        boundary.V[i].volume_tension_forces(&boundary);
        boundary.V[i].forvec += boundary.V[i].VolTForce;
    }
}

void force_calculation(INTERFACE &boundary, vector<PARTICLE> &counterions, const double scalefactor, char bucklingFlag) {

    //Common MPI Message objects
    vector<VECTOR3D> forVecMesh(sizFVecMesh, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> forVecMeshGather(boundary.V.size(), VECTOR3D(0, 0, 0));
    vector<VECTOR3D> forVecIons(sizFVecIons, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> forVecIonsGather(counterions.size(), VECTOR3D(0, 0, 0));
    int i=0, j=0;
    
// ### Compute Forces ###
    //boundary.areaGradient.clear();
    // Recalculate negated geometric gradients:
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].compute_neg_VGrads(&boundary/*, constraintForm*/);

    // Update LJ & ES forces:
    #pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = lowerBoundMesh; i <= upperBoundMesh; i++) {
        // Mesh-mesh component:
        for (j = 0; j < boundary.V.size(); j++) {
            if (i == j) continue;
            VECTOR3D ljForce;
            VECTOR3D esForce;
            VECTOR3D r_vec = boundary.V[i].posvec - boundary.V[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = boundary.lj_length_mesh_mesh * boundary.lj_length_mesh_mesh;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljForce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esForce = ((-scalefactor * boundary.V[i].q * boundary.V[j].q / boundary.em) * (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1/r2)) ^ ((-1.0/r)-(1.0/boundary.inv_kappa_out)));

            forVecMesh[i - lowerBoundMesh] += ljForce + esForce;
        }
        // Mesh-ion component (for mesh):
        for (j = 0; j < counterions.size(); j++) {
            VECTOR3D ljForce;
            VECTOR3D esForce;
            VECTOR3D r_vec = boundary.V[i].posvec - counterions[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = boundary.lj_length_mesh_ions * boundary.lj_length_mesh_ions;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljForce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esForce = ((-scalefactor * boundary.V[i].q * counterions[j].q / boundary.em) *
                       (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1 / r2)) ^ ((-1.0 / r) - (1.0 / boundary.inv_kappa_out)));

            forVecMesh[i - lowerBoundMesh] += ljForce + esForce;
        }
    }
    // Ion-ion component:
    #pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = lowerBoundIons; i <= upperBoundIons; i++) {
        // Ion-ion component:
        for (j = 0; j < counterions.size(); j++) {
            if (i == j) continue;
            VECTOR3D ljForce;
            VECTOR3D esForce;
            VECTOR3D r_vec = counterions[i].posvec - counterions[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = counterions[0].diameter * counterions[0].diameter;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljForce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esForce = ((-scalefactor * counterions[i].q * counterions[j].q / boundary.em) * (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1/r2)) ^ ((-1.0/r)-(1.0/boundary.inv_kappa_out)));

            forVecIons[i - lowerBoundIons] += ljForce + esForce;
        }
        // Mesh-ion component (for ions):
        for (j = 0; j < boundary.V.size(); j++) {
            VECTOR3D ljForce;
            VECTOR3D esForce;
            VECTOR3D r_vec = counterions[i].posvec - boundary.V[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = boundary.lj_length_mesh_ions * boundary.lj_length_mesh_ions;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljForce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esForce = ((-scalefactor * counterions[i].q * boundary.V[j].q / boundary.em) * (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1/r2)) ^ ((-1.0/r)-(1.0/boundary.inv_kappa_out)));

            forVecIons[i - lowerBoundIons] += ljForce + esForce;
        }
    }

    // %%% Additional mesh-only updates: %%%
    // Initialize and update bending forces:
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].bforce = VECTOR3D(0,0,0);
    for (i = 0; i < boundary.E.size(); i++)
        boundary.E[i].bending_forces(&boundary);
    // Update stretching forces:
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].stretching_forces(&boundary, bucklingFlag);
    // Surface & volume tension forces (2017.09.17, 2019.08.16 NB added):
    for (i = 0; i < boundary.V.size(); i++) { // Recalculation of negated gradients already handled above.
        boundary.V[i].tension_forces(&boundary);
        boundary.V[i].volume_tension_forces(&boundary);
    }

    // Collect MPI forVec computations via broadcasting using (all_gather = gather + broadcast):
    if (world.size() > 1) {
        all_gather(world, &forVecMesh[0], forVecMesh.size(), forVecMeshGather);
        all_gather(world, &forVecIons[0], forVecIons.size(), forVecIonsGather);
    }
    else {
        for (i = lowerBoundMesh; i <= upperBoundMesh; i++)
            forVecMeshGather[i] = forVecMesh[i - lowerBoundMesh];
        for (i = lowerBoundIons; i <= upperBoundIons; i++)
            forVecIonsGather[i] = forVecIons[i - lowerBoundMesh];
    }

    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].forvec = forVecMeshGather[i] + boundary.V[i].bforce + boundary.V[i].sforce + boundary.V[i].TForce + boundary.V[i].VolTForce;

    for (i = 0; i < counterions.size(); i++)
        counterions[i].forvec = forVecIonsGather[i];

    forVecMesh.clear();
    forVecIons.clear();
    forVecMeshGather.clear();
    forVecIonsGather.clear();
    }