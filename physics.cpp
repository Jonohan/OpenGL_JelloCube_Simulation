/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

    Haohan Wang
*/

#include "jello.h"
#include "physics.h"
#include "input.h"
#include <math.h>
#include <stdlib.h>


// 1*1*1/7
const double restLength = 0.14285714;


/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world* jello, struct point a[8][8][8])
{
    struct point totalForce = { 0.0, 0.0, 0.0 };
    struct point tempForce = { 0.0, 0.0, 0.0 };
    struct point testExternalForce = { -0.1, 0.0, 0.0 }; // Hardcode test

    double mass = jello->mass;

    double forceAdjustFactor = 1.0; // Manual adjust, if move too fast

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
            {
                a[i][j][k].x = 0;
                a[i][j][k].y = 0;
                a[i][j][k].z = 0;

                StructuralSpringsForce(jello, i, j, k, tempForce);
                pSUM(a[i][j][k], tempForce, a[i][j][k]);

                ShearSpringsForce(jello, i, j, k, tempForce);
                pSUM(a[i][j][k], tempForce, a[i][j][k]);

                BendSpringsForce(jello, i, j, k, tempForce);
                pSUM(a[i][j][k], tempForce, a[i][j][k]);

                CollsionForce(jello, i, j, k, tempForce);
                pSUM(a[i][j][k], tempForce, a[i][j][k]);

                ForceField(jello, i, j, k, tempForce);
                pSUM(a[i][j][k], tempForce, a[i][j][k]);

                pSUM(a[i][j][k], mouseForce, a[i][j][k]);

                //pSUM(a[i][j][k], testExternalForce, a[i][j][k]);

                pMULTIPLY(a[i][j][k], forceAdjustFactor, a[i][j][k]);

                pMULTIPLY(a[i][j][k], (1 / mass), a[i][j][k]);

            }
}


// Elastic force exerted on a , neighbor is b
// Hook's law: F = -k(|L|-R)(L/|L|)
void HookForce(struct point a, struct point b, double kHook, double restLength, struct point &hF)
{
    hF = { 0.0,0.0,0.0 };
    struct point L;
    struct point unitL;
    double len;   
    double length;// for pNORMALIZE()


    // L be the vector pointing from B to A
    // L = vector a - vecotor b
    pDIFFERENCE(a, b, L); 

    // |L|    
    pLENGTH(L, len); 

    //-k(|L|-R)
    double forceCal = -1.0 * kHook * (len - restLength);

    // L/|L|
    pCPY(L, unitL);
    pNORMALIZE(unitL);

    //result
    pMULTIPLY(unitL, forceCal, hF);
}

// Damping force exerted on A, neighbor is B
// F=-k((velA-velB).L)/|L|*L/|L|
void DampForce(struct point a, struct point b, struct point velA , struct point velB, double kDamp, struct point &dF)
{
    dF = { 0.0,0.0,0.0 };
    struct point velRel;
    struct point L;
    struct point unitL;

    double len;
    double length;
    double dotP;

    //relative velocity velA-velB
    pDIFFERENCE(velA, velB, velRel);

    // L be the vector pointing from B to A
    // L = vector a - vecotor b
    pDIFFERENCE(a, b, L);

    //(velA-velB).L
    DOTPRODUCT(velRel, L, dotP);

    // |L|
    pLENGTH(L, len);

    // -k((velA-velB).L)/|L|
    double forceCal = -1.0 * kDamp * dotP / len;

    // L/|L|
    pCPY(L, unitL);
    pNORMALIZE(unitL);

    // result
    pMULTIPLY(unitL, forceCal, dF);
}

/* # of Structural spring :
    central: 6
    surface: 5
    edge: 4
    diagnal: 3 
*/ 

void StructuralSpringsForce(struct world *jello, int i, int j, int k, struct point& totalsSF)
{
    struct point hookF;
    struct point dampF;
    struct point sumF;

    totalsSF = { 0.0,0.0,0.0 };

    // front
    if (j + 1 <= 7 && j >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j + 1][k], jello->kElastic, restLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j + 1][k], jello->v[i][j][k], jello->v[i][j + 1][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalsSF, sumF, totalsSF);
    }

    // back
    if (j - 1 >= 0 && j <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j - 1][k], jello->kElastic, restLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j - 1][k], jello->v[i][j][k], jello->v[i][j - 1][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalsSF, sumF, totalsSF);
    }


    // right
    if (i - 1 >= 0 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j][k], jello->kElastic, restLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j][k], jello->v[i][j][k], jello->v[i - 1][j][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalsSF, sumF, totalsSF);
    }

    // left
    if (i + 1 <=7 && i >=0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j][k], jello->kElastic, restLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j][k], jello->v[i][j][k], jello->v[i + 1][j][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalsSF, sumF, totalsSF);
    }

    // up
    if (k + 1 <= 7 && k >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j][k + 1], jello->kElastic, restLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j][k + 1], jello->v[i][j][k], jello->v[i][j][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalsSF, sumF, totalsSF);
    }

    // down
    if (k - 1 >= 0 && k <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j][k - 1], jello->kElastic, restLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j][k - 1], jello->v[i][j][k], jello->v[i][j][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalsSF, sumF, totalsSF);
    }
}

/* # of shear springs
* 27-6-1=20
*/
void ShearSpringsForce(struct world* jello, int i, int j, int k, struct point& totalssF)
{
    struct point hookF;
    struct point dampF;
    struct point sumF;

    totalssF = { 0.0,0.0,0.0 };

    // sqrt 2 * length
    double shrestLength = 1.414213 * restLength;

    // Assume when looking from above on -k direction, +i point to south, +j point to east
    // 
    // Upper layer: 8
    // NW
    if (k + 1 <= 7 && j - 1 >= 0 && i - 1 >= 0 && k >= 0 && j <= 7 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    // NE
    if (k + 1 <= 7 && j + 1 <= 7 && i - 1 >= 0 && k >= 0 && j >= 0 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //SE
    if (k + 1 <= 7 && j + 1 <= 7 && i + 1 <= 7 && k >= 0 && j >= 0 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //SW
    if (k + 1 <= 7 && j - 1 >= 0 && i + 1 <= 7 && k >= 0 && j <= 7 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //N
    if (k + 1 <= 7 && i - 1 >= 0 && k >= 0 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j][k + 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j][k + 1], jello->v[i][j][k], jello->v[i - 1][j][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //S
    if (k + 1 <= 7 && i + 1 <= 7 && k >= 0 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->v[i][j][k], jello->v[i + 1][j][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //W
    if (k + 1 <= 7 && j - 1 >= 0 && k >= 0 && j <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->v[i][j][k], jello->v[i][j - 1][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //E
    if (k + 1 <= 7 && j + 1 <= 7 && k >= 0 && j >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->v[i][j][k], jello->v[i][j + 1][k + 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    // Lower layer: 8
    // NW
    if (k - 1 >= 0 && j - 1 >= 0 && i - 1 >= 0 && k <= 7 && j <= 7 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    // NE
    if (k - 1 >= 0 && j + 1 <= 7 && i - 1 >= 0 && k <= 7 && j >= 0 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //SE
    if (k - 1 >= 0 && j + 1 <= 7 && i + 1 <= 7 && k <= 7 && j >= 0 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //SW
    if (k - 1 >= 0 && j - 1 >= 0 && i + 1 <= 7 && k <= 7 && j <= 7 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //N
    if (k - 1 >= 0 && i - 1 >= 0 && k <= 7 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j][k - 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j][k - 1], jello->v[i][j][k], jello->v[i - 1][j][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //S
    if (k - 1 >= 0 && i + 1 <= 7 && k <= 7 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->v[i][j][k], jello->v[i + 1][j][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //W
    if (k - 1 >= 0 && j - 1 >= 0 && k <= 7 && j <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j - 1][k - 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j - 1][k - 1], jello->v[i][j][k], jello->v[i][j - 1][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //E
    if (k - 1 >= 0 && j + 1 <= 7 && k <= 7 && j >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->v[i][j][k], jello->v[i][j + 1][k - 1], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    // middle layer: 4
    // NW
    if (j - 1 >= 0 && i - 1 >= 0 && j <= 7 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->v[i][j][k], jello->v[i - 1][j - 1][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    // NE
    if (j + 1 <= 7 && i - 1 >= 0 && j >= 0 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->v[i][j][k], jello->v[i - 1][j + 1][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //SE
    if (j + 1 <= 7 && i + 1 <= 7 && j >= 0 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->v[i][j][k], jello->v[i + 1][j + 1][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }

    //SW
    if (j - 1 >= 0 && i + 1 <= 7 && j <= 7 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->kElastic, shrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->v[i][j][k], jello->v[i + 1][j - 1][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalssF, sumF, totalssF);
    }
}

/* # of bend springs:
* 6
*/
void BendSpringsForce(struct world* jello, int i, int j, int k, struct point& totalbSF)
{
    struct point hookF;
    struct point dampF;
    struct point sumF;

    double bendAdj = 0.6;
    double bendrestLength = bendAdj * 2.0 * restLength; // I adjest it a little bit to confront external pressure

    totalbSF = { 0.0,0.0,0.0 };

    // front
    if (j + 2 <= 7 && j >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j + 2][k], jello->kElastic, bendrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j + 2][k], jello->v[i][j][k], jello->v[i][j + 2][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalbSF, sumF, totalbSF);
    }

    // back
    if (j - 2 >= 0 && j <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j - 2][k], jello->kElastic, bendrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j - 2][k], jello->v[i][j][k], jello->v[i][j - 2][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalbSF, sumF, totalbSF);
    }

    // right
    if (i - 2 >= 0 && i <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i - 2][j][k], jello->kElastic, bendrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i - 2][j][k], jello->v[i][j][k], jello->v[i - 2][j][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalbSF, sumF, totalbSF);
    }

    // left
    if (i + 2 <= 7 && i >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i + 2][j][k], jello->kElastic, bendrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i + 2][j][k], jello->v[i][j][k], jello->v[i + 2][j][k], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalbSF, sumF, totalbSF);
    }

    // up
    if (k + 2 <= 7 && k >= 0)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j][k + 2], jello->kElastic, bendrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j][k + 2], jello->v[i][j][k], jello->v[i][j][k + 2], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalbSF, sumF, totalbSF);
    }

    // down
    if (k - 2 >= 0 && k <= 7)
    {
        HookForce(jello->p[i][j][k], jello->p[i][j][k - 2], jello->kElastic, bendrestLength, hookF);
        DampForce(jello->p[i][j][k], jello->p[i][j][k - 2], jello->v[i][j][k], jello->v[i][j][k - 2], jello->dElastic, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalbSF, sumF, totalbSF);
    }
}

/*
* Collison with boundary box (4*4)
* from -2 to 2
* # of surface: 6
*/
void CollsionForce(struct world* jello, int i, int j, int k, struct point& totalcF)
{
    struct point hookF;
    struct point dampF;
    struct point sumF;

    totalcF = { 0.0,0.0,0.0 };

    struct point velSur = { 0.0, 0.0, 0.0 }; // velocity of surface is 0

    double adjRestLen = -0.0; // mannually adjust, length ¡ý deformation ¡ü, linear

    // with left surface
    if (jello->p[i][j][k].x > 2.0)
    {
        struct point SurfP = { 2.0, jello->p[i][j][k].y , jello->p[i][j][k].z };
        HookForce(jello->p[i][j][k], SurfP, jello->kCollision, adjRestLen, hookF);
        DampForce(jello->p[i][j][k], SurfP, jello->v[i][j][k], velSur, jello->dCollision, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalcF, sumF, totalcF);
    }

    // with right surface
    if (jello->p[i][j][k].x < -2.0)
    {
        struct point SurfP = { -2.0, jello->p[i][j][k].y , jello->p[i][j][k].z };
        HookForce(jello->p[i][j][k], SurfP, jello->kCollision, adjRestLen, hookF);
        DampForce(jello->p[i][j][k], SurfP, jello->v[i][j][k], velSur, jello->dCollision, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalcF, sumF, totalcF);
    }

    // with front surface
    if (jello->p[i][j][k].y > 2.0)
    {
        struct point SurfP = { jello->p[i][j][k].x, 2.0 , jello->p[i][j][k].z };
        HookForce(jello->p[i][j][k], SurfP, jello->kCollision, adjRestLen, hookF);
        DampForce(jello->p[i][j][k], SurfP, jello->v[i][j][k], velSur, jello->dCollision, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalcF, sumF, totalcF);
    }

    // with front surface
    if (jello->p[i][j][k].y < -2.0)
    {
        struct point SurfP = { jello->p[i][j][k].x, -2.0 , jello->p[i][j][k].z };
        HookForce(jello->p[i][j][k], SurfP, jello->kCollision, adjRestLen, hookF);
        DampForce(jello->p[i][j][k], SurfP, jello->v[i][j][k], velSur, jello->dCollision, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalcF, sumF, totalcF);
    }

    // with top surface
    if (jello->p[i][j][k].z > 2.0)
    {
        struct point SurfP = { jello->p[i][j][k].x, jello->p[i][j][k].y , 2.0 };
        HookForce(jello->p[i][j][k], SurfP, jello->kCollision, adjRestLen, hookF);
        DampForce(jello->p[i][j][k], SurfP, jello->v[i][j][k], velSur, jello->dCollision, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalcF, sumF, totalcF);
    }

    // with bottom surface
    if (jello->p[i][j][k].z < -2.0)
    {
        struct point SurfP = { jello->p[i][j][k].x, jello->p[i][j][k].y , -2.0 };
        HookForce(jello->p[i][j][k], SurfP, jello->kCollision, adjRestLen, hookF);
        DampForce(jello->p[i][j][k], SurfP, jello->v[i][j][k], velSur, jello->dCollision, dampF);
        pSUM(hookF, dampF, sumF);
        pSUM(totalcF, sumF, totalcF);
    }

    // Computed inclined plane collision force
    if (jello->incPlanePresent)
    {
        
        double posPoint;
        double length;
        // a=-1;b = 1;c = 1;d = 2;
        // in rotate.w it's a=1;b = -1;c = -1; d = 4;
        struct point planeNormalVector = { jello->a, jello->b, jello->c };
        // F(i,j,k) = ai+bj+ck+d
        // if F(i,j,k) > 0, on the positive side of vector, if < 0 on the negative side (inside the corne)
        posPoint = jello->a * jello->p[i][j][k].x + jello->b * jello->p[i][j][k].y + jello->c * jello->p[i][j][k].z + jello->d;

        // distance of point to plane = |ax+by+cz+d|/sqrt(a^2+b^2+c^2)
        double normalVectorDist;
        double pointDist;
        pLENGTH(planeNormalVector, normalVectorDist);
        pointDist = posPoint / normalVectorDist;  // distance with direction (can be negative)

        // when cube get out off boundray, collision force
        struct point posProjOnPlane;
        struct point planePointVector;
        struct point planeUnitNormalVector;
        if (posPoint < 0)
        {
            //printf("Vector: (%f, %f, %f), d = %f\n", planeNormalVector.x, planeNormalVector.y, planeNormalVector.z , jello->d);
            pCPY(planeNormalVector, planeUnitNormalVector);
            pNORMALIZE(planeUnitNormalVector);
            // scale plane unit normal vector to plane vector point to point
            pMULTIPLY(planeUnitNormalVector, pointDist, planePointVector);
            // position of point - vector(from plane to point) = position of projection point(vector begin)
            pDIFFERENCE(jello->p[i][j][k], planePointVector, posProjOnPlane);

            HookForce(jello->p[i][j][k], posProjOnPlane, jello->kCollision, adjRestLen, hookF);
            DampForce(jello->p[i][j][k], posProjOnPlane, jello->v[i][j][k], velSur, jello->dCollision, dampF);
            pSUM(hookF, dampF, sumF);
            pSUM(totalcF, sumF, totalcF);
        }
    }
}



void ForceField(struct world* jello, int i, int j, int k, struct point& totalfF)
{
    totalfF = { 0.0,0.0,0.0 };

    // length of Force Field each edge = total length / # of edges
    double fFLength = 4.0 / (jello->resolution-1.0);

    // find the approximation point of p[i][j][k] on the bottom corner on -i,-j,-k direction in the jello.forceField
    // e.g. p[i][j][k] = {-1.9,-1.9,-1.9}, then approximation point = [0][0][0] in jello.forceField ({-2.0,-2.0,-2.0} in world coordinate)   
    int ffcoorX;
    int ffcoorY;
    int ffcoorZ;

    if (jello->p[i][j][k].x <= -2.0)
    {
        ffcoorX = 0;
    }
    else if (jello->p[i][j][k].x >= 2.0)
    {
        ffcoorX = (jello->resolution - 1)-1; //since one of 8 points coordinate can be ffPos+1
    }
    else
    {
        ffcoorX = (int)((jello->p[i][j][k].x + 2.0) / fFLength);
    }

    if (jello->p[i][j][k].y <= -2.0)
    {
        ffcoorY = 0;
    }
    else if (jello->p[i][j][k].y >= 2.0)
    {
        ffcoorY = (jello->resolution - 1) - 1;
    }
    else
    {
        ffcoorY = (int)((jello->p[i][j][k].y + 2.0) / fFLength);
    }

    if (jello->p[i][j][k].z <= -2.0)
    {
        ffcoorZ = 0;
    }
    else if (jello->p[i][j][k].z >= 2.0)
    {
        ffcoorZ = (jello->resolution - 1) - 1;
    }
    else
    {
        ffcoorZ = (int)((jello->p[i][j][k].z + 2.0) / fFLength);
    }

    //Assume ffPoint on -i,-j,-k direction is A000, other can represent by +1 on corresponding axis
    //A000 represent force on A000 point
    //e.g. A100 is p[i+1][j][k]
    //# of ffpoint: 8
    //e.g. jello.forceField[i * jello.resolution * jello.resolution +j * jello.resolution + k].x
    struct point A000 = jello->forceField[ffcoorX * jello->resolution * jello->resolution + ffcoorY * jello->resolution + ffcoorZ];
    struct point A100 = jello->forceField[(ffcoorX + 1) * jello->resolution * jello->resolution + ffcoorY * jello->resolution + ffcoorZ];
    struct point A110 = jello->forceField[(ffcoorX + 1) * jello->resolution * jello->resolution + (ffcoorY + 1) * jello->resolution + ffcoorZ];
    struct point A111 = jello->forceField[(ffcoorX + 1) * jello->resolution * jello->resolution + (ffcoorY + 1) * jello->resolution + (ffcoorZ + 1)];
    struct point A011 = jello->forceField[ffcoorX * jello->resolution * jello->resolution + (ffcoorY + 1) * jello->resolution + (ffcoorZ + 1)];
    struct point A001 = jello->forceField[ffcoorX * jello->resolution * jello->resolution + ffcoorY * jello->resolution + (ffcoorZ + 1)];
    struct point A101 = jello->forceField[(ffcoorX + 1) * jello->resolution * jello->resolution + ffcoorY * jello->resolution + (ffcoorZ + 1)];
    struct point A010 = jello->forceField[ffcoorX * jello->resolution * jello->resolution + (ffcoorY + 1) * jello->resolution + ffcoorZ];

    // relative coordinate compute = p[i][j][k] - position of approximation point
    double iRel = (jello->p[i][j][k].x - (ffcoorX * fFLength - 2.0))/ fFLength;
    double jRel = (jello->p[i][j][k].y - (ffcoorY * fFLength - 2.0))/ fFLength;
    double kRel = (jello->p[i][j][k].z - (ffcoorZ * fFLength - 2.0))/ fFLength;

    //printf("Relative coor: %f", iRel);

    //totalfF = (1-iRel)(1-jRel)(1-kRel)A000 + ...
    struct point tempForce = { 0.0,0.0,0.0 };

    pMULTIPLY(A000, (1 - iRel) * (1 - jRel) * (1 - kRel), tempForce);
    pSUM(totalfF, tempForce, totalfF);

    pMULTIPLY(A100, iRel * (1 - jRel) * (1 - kRel), tempForce);
    pSUM(totalfF, tempForce, totalfF);

    pMULTIPLY(A110, iRel * jRel * (1 - kRel), tempForce);
    pSUM(totalfF, tempForce, totalfF);

    pMULTIPLY(A111, iRel * jRel * kRel, tempForce);
    pSUM(totalfF, tempForce, totalfF);

    pMULTIPLY(A011, (1 - iRel) * jRel * kRel, tempForce);
    pSUM(totalfF, tempForce, totalfF);

    pMULTIPLY(A001, (1 - iRel) * (1 - jRel) * kRel, tempForce);
    pSUM(totalfF, tempForce, totalfF);

    pMULTIPLY(A101, iRel * (1 - jRel) * kRel, tempForce);
    pSUM(totalfF, tempForce, totalfF);

    pMULTIPLY(A010, (1 - iRel) * jRel * (1 - kRel), tempForce);
    pSUM(totalfF, tempForce, totalfF);
}


/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;
      }
}



/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
