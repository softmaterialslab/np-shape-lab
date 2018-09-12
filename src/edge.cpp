// Edge cpp file

#include "edge.h"

VERTEX *EDGE::opposite(VERTEX *currentv) {
    assert(itsV.size() == 2);
    if (itsV[0] == currentv)
        return itsV[1];
    if (itsV[1] == currentv)
        return itsV[0];
    return NULL;
}


FACE *EDGE::opposite(FACE *currentf) {
    assert(itsF.size() == 2);
    if (itsF[0] == currentf)
        return itsF[1];
    if (itsF[1] == currentf)
        return itsF[0];
    return NULL;
}

double EDGE::lengthSquared() {
    VECTOR3D difference = itsV[1]->posvec - itsV[0]->posvec;
    return difference * difference;
}

double EDGE::length() {
    return sqrt(lengthSquared());
}

double EDGE::compute_S() {
    itsS = (itsF[0]->itsdirection) * (itsF[1]->itsdirection);
    return itsS;
}

void EDGE::bending_forces(INTERFACE *boundary) {
    compute_S();
    vector<VERTEX *> v;
    v.resize(4);
    v[0] = itsV[0];
    v[1] = itsV[1];
    v[2] = itsF[0]->across(this);
    v[3] = itsF[1]->across(this);


    for (int i = 0; i < 4; i++) {
        VECTOR3D grad0 = (boundary->areaGradient[make_pair(v[i], itsF[0])]);
        VECTOR3D grad1 = (boundary->areaGradient[make_pair(v[i], itsF[1])]);
        v[i]->bforce += boundary->bkappa
                        * (0.25 / (itsF[0]->itsarea * itsF[1]->itsarea))
                        ^ (compute_gradS(v[i])
                           - ((grad0 ^ (itsS / itsF[0]->itsarea))
                              + (grad1 ^ (itsS / itsF[1]->itsarea))));
    }

}

VECTOR3D EDGE::compute_gradS(VERTEX *wrt) {
    VERTEX *v0;
    VERTEX *v1;
    VERTEX *v2;
    VERTEX *v3;
    VECTOR3D result;

    if (wrt == itsV[0] || wrt == itsV[1]) {
        v0 = wrt;
        v2 = opposite(wrt);
        v1 = itsF[0]->across(this);
        v3 = itsF[1]->across(this);

        result.x =
                (v2->posvec.y - v3->posvec.y) *
                (-v2->posvec.x * v1->posvec.y + v0->posvec.x * v1->posvec.y + v1->posvec.x * v2->posvec.y -
                 v0->posvec.x * v2->posvec.y - v1->posvec.x * v0->posvec.y + v2->posvec.x * v0->posvec.y)
                + (v1->posvec.y - v2->posvec.y) *
                  (-v3->posvec.x * v2->posvec.y + v0->posvec.x * v2->posvec.y + v2->posvec.x * v3->posvec.y -
                   v0->posvec.x * v3->posvec.y - v2->posvec.x * v0->posvec.y + v3->posvec.x * v0->posvec.y)
                + (-v2->posvec.z + v3->posvec.z) *
                  (v2->posvec.x * v1->posvec.z - v0->posvec.x * v1->posvec.z - v1->posvec.x * v2->posvec.z +
                   v0->posvec.x * v2->posvec.z + v1->posvec.x * v0->posvec.z - v2->posvec.x * v0->posvec.z)
                + (-v1->posvec.z + v2->posvec.z) *
                  (v3->posvec.x * v2->posvec.z - v0->posvec.x * v2->posvec.z - v2->posvec.x * v3->posvec.z +
                   v0->posvec.x * v3->posvec.z + v2->posvec.x * v0->posvec.z - v3->posvec.x * v0->posvec.z);
        result.y =
                (-v2->posvec.x + v3->posvec.x) *
                (-v2->posvec.x * v1->posvec.y + v0->posvec.x * v1->posvec.y + v1->posvec.x * v2->posvec.y -
                 v0->posvec.x * v2->posvec.y - v1->posvec.x * v0->posvec.y + v2->posvec.x * v0->posvec.y)
                + (-v1->posvec.x + v2->posvec.x) *
                  (-v3->posvec.x * v2->posvec.y + v0->posvec.x * v2->posvec.y + v2->posvec.x * v3->posvec.y -
                   v0->posvec.x * v3->posvec.y - v2->posvec.x * v0->posvec.y + v3->posvec.x * v0->posvec.y)
                + (v2->posvec.z - v3->posvec.z) *
                  (-v2->posvec.y * v1->posvec.z + v0->posvec.y * v1->posvec.z + v1->posvec.y * v2->posvec.z -
                   v0->posvec.y * v2->posvec.z - v1->posvec.y * v0->posvec.z + v2->posvec.y * v0->posvec.z)
                + (v1->posvec.z - v2->posvec.z) *
                  (-v3->posvec.y * v2->posvec.z + v0->posvec.y * v2->posvec.z + v2->posvec.y * v3->posvec.z -
                   v0->posvec.y * v3->posvec.z - v2->posvec.y * v0->posvec.z + v3->posvec.y * v0->posvec.z);
        result.z =
                (v2->posvec.x - v3->posvec.x) *
                (v2->posvec.x * v1->posvec.z - v0->posvec.x * v1->posvec.z - v1->posvec.x * v2->posvec.z +
                 v0->posvec.x * v2->posvec.z + v1->posvec.x * v0->posvec.z - v2->posvec.x * v0->posvec.z)
                + (v1->posvec.x - v2->posvec.x) *
                  (v3->posvec.x * v2->posvec.z - v0->posvec.x * v2->posvec.z - v2->posvec.x * v3->posvec.z +
                   v0->posvec.x * v3->posvec.z + v2->posvec.x * v0->posvec.z - v3->posvec.x * v0->posvec.z)
                + (-v2->posvec.y + v3->posvec.y) *
                  (-v2->posvec.y * v1->posvec.z + v0->posvec.y * v1->posvec.z + v1->posvec.y * v2->posvec.z -
                   v0->posvec.y * v2->posvec.z - v1->posvec.y * v0->posvec.z + v2->posvec.y * v0->posvec.z)
                + (-v1->posvec.y + v2->posvec.y) *
                  (-v3->posvec.y * v2->posvec.z + v0->posvec.y * v2->posvec.z + v2->posvec.y * v3->posvec.z -
                   v0->posvec.y * v3->posvec.z - v2->posvec.y * v0->posvec.z + v3->posvec.y * v0->posvec.z);

    } else {
        v1 = wrt;
        for (unsigned int i = 0; i < itsF[0]->itsV.size(); i++)
            if (itsF[0]->itsV[i] != itsV[0] && itsF[0]->itsV[i] != itsV[1]
                && itsF[0]->itsV[i] != wrt)
                v3 = itsF[0]->itsV[i];
        for (unsigned int i = 0; i < itsF[1]->itsV.size(); i++)
            if (itsF[1]->itsV[i] != itsV[0] && itsF[1]->itsV[i] != itsV[1]
                && itsF[1]->itsV[i] != wrt)
                v3 = itsF[1]->itsV[i];
        v0 = itsV[0];
        v2 = itsV[1];

        result.x =
                (-v0->posvec.y + v2->posvec.y) *
                (-v2->posvec.x * v0->posvec.y + v3->posvec.x * v0->posvec.y + v0->posvec.x * v2->posvec.y -
                 v3->posvec.x * v2->posvec.y - v0->posvec.x * v3->posvec.y + v2->posvec.x * v3->posvec.y)
                + (v0->posvec.z - v2->posvec.z) *
                  (v2->posvec.x * v0->posvec.z - v3->posvec.x * v0->posvec.z - v0->posvec.x * v2->posvec.z +
                   v3->posvec.x * v2->posvec.z + v0->posvec.x * v3->posvec.z - v2->posvec.x * v3->posvec.z);
        result.y =
                (v0->posvec.x - v2->posvec.x) *
                (-v2->posvec.x * v0->posvec.y + v3->posvec.x * v0->posvec.y + v0->posvec.x * v2->posvec.y -
                 v3->posvec.x * v2->posvec.y - v0->posvec.x * v3->posvec.y + v2->posvec.x * v3->posvec.y)
                + (-v0->posvec.z + v2->posvec.z) *
                  (-v2->posvec.y * v0->posvec.z + v3->posvec.y * v0->posvec.z + v0->posvec.y * v2->posvec.z -
                   v3->posvec.y * v2->posvec.z - v0->posvec.y * v3->posvec.z + v2->posvec.y * v3->posvec.z);
        result.z =
                (-v0->posvec.x + v2->posvec.x) *
                (v2->posvec.x * v0->posvec.z - v3->posvec.x * v0->posvec.z - v0->posvec.x * v2->posvec.z +
                 v3->posvec.x * v2->posvec.z + v0->posvec.x * v3->posvec.z - v2->posvec.x * v3->posvec.z)
                + (v0->posvec.y - v2->posvec.y) *
                  (-v2->posvec.y * v0->posvec.z + v3->posvec.y * v0->posvec.z + v0->posvec.y * v2->posvec.z -
                   v3->posvec.y * v2->posvec.z - v0->posvec.y * v3->posvec.z + v2->posvec.y * v3->posvec.z);

    }
    return result;
}

double EDGE::crossingLengthSquared() {
    VERTEX *v1;
    VERTEX *v2;
    for (unsigned int i = 0; i < itsF[0]->itsV.size(); i++)
        if (itsF[0]->itsV[i] != itsV[0] && itsF[0]->itsV[i] != itsV[1])
            v1 = itsF[0]->itsV[i];
    for (unsigned int i = 0; i < itsF[1]->itsV.size(); i++) //NOTE replaced F[0] by F[1], technically true.
        if (itsF[1]->itsV[i] != itsV[0] && itsF[1]->itsV[i] != itsV[1])
            v2 = itsF[1]->itsV[i];
    return ((v2->posvec) - (v1->posvec)).GetMagnitudeSquared();
}

void EDGE::flipIfFavorable() {
    if (((itsF[0]->itsnormal) * (itsF[1]->itsnormal)) < 0) {
        flip();
        itsF[0]->computenormal();
        itsF[1]->computenormal();
        return;
    }
    if (crossingLengthSquared() < lengthSquared()) {
        flip();
        itsF[0]->computenormal();
        itsF[1]->computenormal();
        if (((itsF[0]->itsnormal) * (itsF[1]->itsnormal)) < 0) {
            flip();
            itsF[0]->computenormal();
            itsF[1]->computenormal();
            return;
        }
        return;
    }
    return;
}

//
//  v0   e1   v2        v0   e1   v2
//    *------*            *------*
//    |f0   /|            |\   f1|
//  e2|  e0  |e3   -->  e2|  e0  |e3
//    |/   f1|            |f0   \|
//    *------*            *------*
//  v1   e4   v3        v1   e4   v3
//
// do we do clockwise or counterclockwise on our faces?  does it matter?
// convention: v1 comes 1 before v2 on f0; 1 after v2 on f1 (both cyclically)
void EDGE::flip() {
    FACE *f0 = itsF[0];
    FACE *f1 = itsF[1];
    EDGE *e0 = this;
    unsigned int i;
    for (i = 0; i < f0->itsE.size() && f0->itsE[i] != e0; i++);
    assert(f0->itsE[i] == e0);
    EDGE *e1 = f0->itsE[(i + 1) % f0->itsE.size()];
    EDGE *e2 = f0->itsE[(i + 2) % f0->itsE.size()];
    VERTEX *v0 = f0->itsV[(i + 0) % f0->itsV.size()];
    VERTEX *v1 = f0->itsV[(i + 1) % f0->itsV.size()];
    VERTEX *v2 = f0->itsV[(i + 2) % f0->itsV.size()];
    assert(v0 != e0->itsV[0] && v0 != e0->itsV[1]);
    assert(v0 == e1->itsV[0] || v0 == e1->itsV[1]);
    assert(v0 == e2->itsV[0] || v0 == e2->itsV[1]);

    for (i = 0; i < f1->itsE.size() && f1->itsE[i] != e0; i++);
    assert(f1->itsE[i] == e0);
    EDGE *e4 = f1->itsE[(i + 1) % f1->itsE.size()];
    EDGE *e3 = f1->itsE[(i + 2) % f1->itsE.size()];
    VERTEX *v3 = f1->itsV[(i + 0) % f1->itsV.size()];
    assert(v2 == f1->itsV[(i + 1) % f1->itsV.size()]);
    assert(v1 == f1->itsV[(i + 2) % f1->itsV.size()]);
    assert(v3 != e0->itsV[0] && v3 != e0->itsV[1]);
    assert(v3 == e4->itsV[0] || v3 == e4->itsV[1]);
    assert(v3 == e3->itsV[0] || v3 == e3->itsV[1]);

    if (v2->itsF.size() <= 3 && v1->itsF.size() <= 3)
        return;
    if (e0->index == e1->index ||
        e0->index == e2->index ||
        e0->index == e3->index ||
        e0->index == e4->index ||
        e1->index == e2->index ||
        e1->index == e3->index ||
        e1->index == e4->index ||
        e2->index == e3->index ||
        e2->index == e4->index ||
        e3->index == e4->index) {
        //cout << "duplicate edges on flip, refusing to flip\n";
        return;
    }
    if (v0->index == v1->index ||
        v0->index == v2->index ||
        v0->index == v3->index ||
        v1->index == v2->index ||
        v1->index == v3->index ||
        v2->index == v3->index) {
        //cout << "duplicate vertices on flip, refusing to flip\n";
        return;
    }
    if (f0->index == f1->index) {
        //cout << "duplicate faces on flip, refusing to flip\n";
        return;
    }

//  cout << e0->index << "\t" 
//       << e1->index << "\t" 
//       << e2->index << "\t" 
//       << e3->index << "\t" 
//       << e4->index << "\t" 
//       << v0->index << "\t" 
//       << v1->index << "\t" 
//       << v2->index << "\t" 
//       << v3->index << "\t" 
//       << f0->index << "\t" 
//       << f1->index << "\n";

    // start rewiring everything
    assert(e0->itsV.size() == 2);

    // rewire the crossing edge, e0
    e0->itsV[0] = v3;
    e0->itsV[1] = v0;
    // the faces shouldn't change for e0 (both will sill be f0 and f1)

    assert(f0->itsE.size() == 3);
    assert(f0->itsV.size() == 3);
    assert(f1->itsE.size() == 3);
    assert(f1->itsV.size() == 3);

    // the two faces are wired correctly here
    f0->itsE[0] = e0;
    f0->itsE[1] = e2;
    f0->itsE[2] = e4;
    f0->itsV[0] = v1;
    f0->itsV[1] = v3;
    f0->itsV[2] = v0;

    f1->itsE[0] = e0;
    f1->itsE[1] = e3;
    f1->itsE[2] = e1;
    f1->itsV[0] = v2;
    f1->itsV[1] = v0;
    f1->itsV[2] = v3;

    // replace f0 and f1 on the two edges where something changes
    int count = 0;
    for (i = 0; i < e1->itsF.size(); i++)
        if (e1->itsF[i] == f0) {
            e1->itsF[i] = f1;
            count++;
        }
    assert (count == 1);
    count = 0;
    for (i = 0; i < e4->itsF.size(); i++)
        if (e4->itsF[i] == f1) {
            e4->itsF[i] = f0;
            count++;
        }
    assert (count == 1);

    // fix the vertices
    vector<FACE *> fnew;
    vector<EDGE *> enew;
    // v0
    for (i = 0; i < v0->itsF.size(); i++) {
        fnew.push_back(v0->itsF[i]);
        if (v0->itsF[i] == f0)
            fnew.push_back(f1);
    }
    v0->itsF.swap(fnew);
    fnew.clear();
    for (i = 0; i < v0->itsE.size(); i++) {
        enew.push_back(v0->itsE[i]);
        if (v0->itsE[i] == e2)
            enew.push_back(e0);
    }
    v0->itsE.swap(enew);
    enew.clear();
    // v1
    for (i = 0; i < v1->itsF.size(); i++) {
        if (v1->itsF[i] != f1)
            fnew.push_back(v1->itsF[i]);
    }
    v1->itsF.swap(fnew);
    fnew.clear();
    for (i = 0; i < v1->itsE.size(); i++) {
        if (v1->itsE[i] != e0)
            enew.push_back(v1->itsE[i]);
    }
    v1->itsE.swap(enew);
    enew.clear();
    // v2
    for (i = 0; i < v2->itsF.size(); i++) {
        if (v2->itsF[i] != f0)
            fnew.push_back(v2->itsF[i]);
    }
    v2->itsF.swap(fnew);
    fnew.clear();
    for (i = 0; i < v2->itsE.size(); i++) {
        if (v2->itsE[i] != e0)
            enew.push_back(v2->itsE[i]);
    }
    v2->itsE.swap(enew);
    enew.clear();
    // v3
    for (i = 0; i < v3->itsF.size(); i++) {
        fnew.push_back(v3->itsF[i]);
        if (v3->itsF[i] == f1)
            fnew.push_back(f0);
    }
    v3->itsF.swap(fnew);
    fnew.clear();
    for (i = 0; i < v3->itsE.size(); i++) {
        enew.push_back(v3->itsE[i]);
        if (v3->itsE[i] == e3)
            enew.push_back(e0);
    }
    v3->itsE.swap(enew);
    enew.clear();

}
