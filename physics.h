/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);

void HookForce(struct point a, struct point b, double kHook, double restLength, struct point& hF);
void DampForce(struct point a, struct point b, struct point velA, struct point velB, double kDamp, struct point& dF);

void StructuralSpringsForce(struct world* jello, int i, int j, int k, struct point& totalsSF);
void ShearSpringsForce(struct world* jello, int i, int j, int k, struct point& totalssF);
void BendSpringsForce(struct world* jello, int i, int j, int k, struct point& totalbSF);

void CollsionForce(struct world* jello, int i, int j, int k, struct point& totalcF);

void ForceField(struct world* jello, int i, int j, int k, struct point& totalfF);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);

#endif

