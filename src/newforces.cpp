#include "newforces.h"
#include "functions.h"

void force_calculation_init(INTERFACE &boundary, const double scalefactor, char bucklingFlag) {
    unsigned int i=0, j=0;

// Compute constraint gradients: (called 'forces' in prev comment & in SHAKE, RATTLE.)
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].compute_neg_VGrads(&boundary/*, constraintForm*/);

    // Initialize forces on all vertices (NB added): (Verified they autoset to zero before.)
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].forvec = VECTOR3D(0, 0, 0);

    // Contribute LJ & electrostatics forces (initialization implicit):
    #pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = 0; i < boundary.V.size(); i++) {
        for (j = 0; j < boundary.V.size(); j++) {
            if (i == j) continue;
            VECTOR3D ljforce;
            VECTOR3D esforce;
            VECTOR3D r_vec = boundary.V[i].posvec - boundary.V[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = boundary.lj_length * boundary.lj_length;

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljforce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esforce = ((-scalefactor * boundary.V[i].q * boundary.V[j].q / boundary.em) * (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                       (r_vec ^ (1/r2)) ^ ((-1.0/r)-(1.0/boundary.inv_kappa_out)));

            boundary.V[i].forvec += ljforce + esforce;
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
    for (i = 0; i < boundary.V.size(); i++) {   // Re-calculation of the numerical geometric gradients not required as positions haven't changed from above.
        boundary.V[i].tension_forces(&boundary);
        boundary.V[i].forvec += boundary.V[i].TForce;
        boundary.V[i].volume_tension_forces(&boundary);
        boundary.V[i].forvec += boundary.V[i].VolTForce;
    }
}

void force_calculation(INTERFACE &boundary, const double scalefactor, char bucklingFlag) {

    //Common MPI Message objects
    vector<VECTOR3D> forvec(sizFVec, VECTOR3D(0, 0, 0));
    vector<VECTOR3D> forvecGather(boundary.V.size() + extraElements, VECTOR3D(0, 0, 0));
    unsigned int i=0, j=0;
    
// ### Compute Forces ###
    //boundary.areaGradient.clear();
    // Recalculate negated geometric gradients:
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].compute_neg_VGrads(&boundary/*, constraintForm*/);

    // Update vertex-vertex LJ & electrostatic forces:
    //i = 0; i < boundary.V.size(); i++
    //i = lowerBound; i <= upperBound; i++
    #pragma omp parallel for schedule(dynamic) default(shared) private(i, j)
    for (i = lowerBound; i <= upperBound; i++) {
        for (j = 0;j < boundary.V.size(); j++) {
        //for (j = i + 1; j < boundary.V.size(); j++) {
            // These do not need re-initialized because they are being set each time.
            if (i == j) continue;
            VECTOR3D ljforce;
            VECTOR3D esforce;
            VECTOR3D r_vec = boundary.V[i].posvec - boundary.V[j].posvec;
            double r = r_vec.GetMagnitude();
            double r2 = r_vec.GetMagnitudeSquared();
            double d2 = boundary.lj_length * boundary.lj_length;

            //VECTOR3D grad_vec = r_vec ^ ((-1.0) / r3);

            if (r2 < dcut2 * d2) {
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double d6 = d2 * d2 * d2;
                double d12 = d6 * d6;
                ljforce = (r_vec ^ (48 * boundary.elj * ((d12 / r12) - 0.5 * (d6 / r6)) * (1 / r2)));
            }

            esforce = ((-scalefactor * boundary.V[i].q * boundary.V[j].q / boundary.em) * (exp(-(1.0 / boundary.inv_kappa_out) * r)) ^
                        (r_vec ^ (1/r2)) ^ ((-1.0/r)-(1.0/boundary.inv_kappa_out)));

            forvec[i-lowerBound] += ljforce + esforce;
            //boundary.V[i].forvec += ljforce + esforce;
            //boundary.V[j].forvec -= ljforce + esforce;
        }
    }

    // Update bending forces:
    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].bforce = VECTOR3D(0,0,0);
    for (i = 0; i < boundary.E.size(); i++)
        boundary.E[i].bending_forces(&boundary);
    //for (i = 0; i < boundary.V.size(); i++)
    //    boundary.V[i].forvec += boundary.V[i].bforce;
    // Update stretching forces:
    for (i = 0; i < boundary.V.size(); i++)
    {
        boundary.V[i].stretching_forces(&boundary, bucklingFlag);
        //boundary.V[i].forvec += boundary.V[i].sforce;
    }
    // Surface & volume tension forces (2017.09.17, 2019.08.16 NB added):
    for (i = 0; i < boundary.V.size(); i++)
    {   // Recalculation of negated gradients already handled above.
        boundary.V[i].tension_forces(&boundary);
        boundary.V[i].volume_tension_forces(&boundary);
    }

    //forvec broadcasting using all gather = gather + broadcast
    if (world.size() > 1) {
        all_gather(world, &forvec[0], forvec.size(), forvecGather);
    } else {
        for (i = lowerBound; i <= upperBound; i++)
            forvecGather[i] = forvec[i - lowerBound];
    }

    for (i = 0; i < boundary.V.size(); i++)
        boundary.V[i].forvec = forvecGather[i] + boundary.V[i].bforce + boundary.V[i].sforce + boundary.V[i].TForce + boundary.V[i].VolTForce;

    forvec.clear();
    forvecGather.clear();



}