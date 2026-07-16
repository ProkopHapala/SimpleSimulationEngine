#ifndef KinematicSolver_h
#define KinematicSolver_h

/// @file KinematicSolver.h
/// @brief Newton-Raphson constraint solver for assemblies of rigid bodies with anchor and slider joints.
///
/// Given N rigid bodies (6 DOF each: 3 position + 3 rotation) and M constraints, solves for
/// body poses satisfying all constraints. Constraint types: Anchor (coincident points on two
/// bodies) and Slider (point constrained to a B-spline path). Uses Levenberg-Marquardt with
/// numerical Jacobian and line search. Slider parameters (t along path) are free unknowns,
/// adding to the state vector. sweepSolve() iterates slider parameters through a full cycle.

#include <math.h>
#include <stdio.h>
#include <vector>
#include <functional>

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

// ====================================================================
//   Kinematic Constraint Solver
//   Given N rigid bodies (6 DOF each) and M constraints (anchors,
//   sliders), solves for body poses satisfying all constraints
//   using Newton-Raphson with numerical Jacobian.
//
//   Body pose: {pos (Vec3d), rot (Mat3d)} — 6 DOF
//   Perturbation: {dp (Vec3d), domega (Vec3d)} — tangent space
//   New pose: pos' = pos + dp, rot' = rot * R(domega)
//
//   Constraint types:
//   - Anchor: P_world(bodyA, lposA) = P_world(bodyB, lposB)  [3 eqs]
//   - Slider: P_world(bodyA, lposA) = pathPoint(t)           [3 eqs]
//     where pathPoint(t) uses cubic B-spline over control points
//
//   Slider parameters t are free unknowns (added to the state vector).
//   Total unknowns: 6*N_bodies + N_sliders
//   Total equations: 3*N_constraints
// ====================================================================

// --- B-spline interpolation (ported from JS SplineCubic.js) ---
// ps: vertex indices into verts array
// closed: wrap indices for closed path
// Returns interpolated position at parameter t in [0,1]
inline Vec3d bsplineInterpolate( double t, const std::vector<int>& ps, bool closed, const std::vector<Vec3d>& verts ){
    int n = ps.size();
    if(n < 2) return Vec3d{0,0,0};
    int totalSegments = closed ? n : n - 1;
    if(totalSegments <= 0) return Vec3d{0,0,0};

    double scaledT = t * totalSegments;
    int i = (int)floor(scaledT);
    double f = scaledT - i;

    if(i >= totalSegments){
        if(closed){ i = i % n; }
        else      { i = totalSegments - 1; f = 1.0; }
    }

    auto getPt = [&](int idx) -> Vec3d {
        if(closed){ idx = ((idx % n) + n) % n; }
        else      { idx = std::max(0, std::min(n - 1, idx)); }
        return verts[ps[idx]];
    };

    Vec3d p0 = getPt(i-1), p1 = getPt(i), p2 = getPt(i+1), p3 = getPt(i+2);
    double f2 = f*f, f3 = f2*f;
    double b0 = (1 - 3*f + 3*f2 - f3)   / 6.0;
    double b1 = (4 - 6*f2 + 3*f3)       / 6.0;
    double b2 = (1 + 3*f + 3*f2 - 3*f3) / 6.0;
    double b3 = f3                      / 6.0;
    return p0*b0 + p1*b1 + p2*b2 + p3*b3;
}

// --- Body pose: 6 DOF (3 pos + 3 rot as rotation matrix) ---
struct KinematicPose {
    Vec3d pos = {0,0,0};
    Mat3d rot = Mat3dIdentity;
    bool bFixed = false;  // if true, body DOF are not part of unknowns
};

// Apply perturbation {dp, domega} to pose
inline KinematicPose applyPerturbation( const KinematicPose& pose, const double* dq ){
    KinematicPose out;
    out.pos = pose.pos + Vec3d{dq[0], dq[1], dq[2]};
    Vec3d omega{dq[3], dq[4], dq[5]};
    Mat3d dR = Mat3dIdentity;
    double r2 = omega.norm2();
    if(r2 > 1e-20){
        double r = sqrt(r2);
        double ca = cos(r), sa = sin(r);
        dR.rotate_csa(ca, sa, omega * (1.0/r));
    }
    out.rot.set_mmul(pose.rot, dR);  // right-multiply: local-frame perturbation
    return out;
}

// Transform local point to world using pose
inline Vec3d localToWorld( const KinematicPose& pose, const Vec3d& lp ){
    Vec3d wp;
    pose.rot.dot_to_T(lp, wp);
    wp.add(pose.pos);
    return wp;
}

// --- Constraint types ---
enum class ConstraintType { Anchor, Slider };

struct KinematicConstraint {
    ConstraintType type;
    // Body indices
    int bodyA = -1;   // -1 means world-fixed
    int bodyB = -1;   // -1 means world-fixed (for slider: the "path body" is fixed)
    // Local positions on each body
    Vec3d lposA = {0,0,0};
    Vec3d lposB = {0,0,0};
    // For Slider: path control points (vertex indices into external vert array)
    // The slider parameter t is a free unknown
    std::vector<int> pathIndices;
    bool pathClosed = false;
    const std::vector<Vec3d>* pathVerts = nullptr;  // pointer to vertex array
    int sliderParamIndex = -1;  // index into slider parameter vector
};

// --- Solver ---
struct KinematicSystem {
    std::vector<KinematicPose> poses;
    std::vector<KinematicConstraint> constraints;
    std::vector<double> sliderParams;  // free parameters, one per slider constraint
    std::vector<bool> sliderFixed;     // if true, slider param is not an unknown (set externally)

    int nBodies()    const { return poses.size(); }
    int nSliders()   const { return sliderParams.size(); }
    int nFreeBodies() const { int n=0; for(const auto& p: poses) if(!p.bFixed) n++; return n; }
    int nFreeSliders() const { int n=0; for(int i=0; i<(int)sliderFixed.size(); i++) if(!sliderFixed[i]) n++; return n; }
    int nUnknowns()  const { return 6 * nFreeBodies() + nFreeSliders(); }
    int nEqPerCon()  const { return 3; }
    int nEquations() const { return 3 * constraints.size(); }

    // Map: unknown index -> (bodyIndex, dofIndex) or (slider, -1)
    std::vector<std::pair<int,int>> unknownMap;  // (ib, k) for body DOF, (-1, si) for slider
    void buildUnknownMap(){
        unknownMap.clear();
        for(int ib=0; ib<nBodies(); ib++){
            if(poses[ib].bFixed) continue;
            for(int k=0; k<6; k++) unknownMap.push_back({ib, k});
        }
        for(int si=0; si<nSliders(); si++){
            if(si < (int)sliderFixed.size() && sliderFixed[si]) continue;
            unknownMap.push_back({-1, si});
        }
    }

    // Get world position of constraint point on body A
    Vec3d getPointA( int ic ) const {
        const auto& c = constraints[ic];
        if(c.bodyA < 0) return c.lposA;  // world-fixed
        return localToWorld(poses[c.bodyA], c.lposA);
    }

    // Get world position of constraint point on body B (anchor) or path point (slider)
    Vec3d getPointB( int ic ) const {
        const auto& c = constraints[ic];
        if(c.type == ConstraintType::Anchor){
            if(c.bodyB < 0) return c.lposB;  // world-fixed
            return localToWorld(poses[c.bodyB], c.lposB);
        } else {  // Slider
            double t = sliderParams[c.sliderParamIndex];
            return bsplineInterpolate(t, c.pathIndices, c.pathClosed, *c.pathVerts);
        }
    }

    // Evaluate all constraints: residual[i*3..i*3+2] = pointA - pointB for constraint i
    void evalResiduals( double* residuals ) const {
        for(int ic = 0; ic < (int)constraints.size(); ic++){
            Vec3d pa = getPointA(ic);
            Vec3d pb = getPointB(ic);
            Vec3d d  = pa - pb;
            residuals[ic*3+0] = d.x;
            residuals[ic*3+1] = d.y;
            residuals[ic*3+2] = d.z;
        }
    }

    // Solve: find poses + slider params such that all residuals = 0
    // Uses Levenberg-Marquardt with line search
    // Returns number of iterations, or -1 if not converged
    int solve( int maxIters = 50, double tol = 1e-8, bool bVerbose = false ){
        buildUnknownMap();
        int n = nUnknowns();
        int m = nEquations();
        if(bVerbose) printf("KinematicSolver: %d unknowns, %d equations (mobility = %d)\n", n, m, n - m);

        std::vector<double> residuals(m, 0.0);
        std::vector<double> Jacobian(m * n, 0.0);
        std::vector<double> dq(n, 0.0);
        double eps = 1e-6;
        double lambda = 1e-3;  // LM damping, adapted per iteration

        evalResiduals(residuals.data());
        double curRes = 0;
        for(int i = 0; i < m; i++) curRes += residuals[i]*residuals[i];

        for(int iter = 0; iter < maxIters; iter++){
            double maxRes = 0;
            for(int i = 0; i < m; i++) maxRes = std::max(maxRes, fabs(residuals[i]));
            if(bVerbose) printf("  iter %d: maxRes = %g, lambda = %g\n", iter, maxRes, lambda);
            if(maxRes < tol){
                if(bVerbose) printf("KinematicSolver: converged in %d iterations (maxRes=%g)\n", iter, maxRes);
                return iter;
            }

            // 1. Compute Jacobian by finite differences
            for(int j = 0; j < n; j++){
                auto [ib, k] = unknownMap[j];
                std::vector<KinematicPose> savedPoses = poses;
                std::vector<double> savedSliders = sliderParams;

                double eps_j = eps;
                if(ib >= 0){
                    std::vector<double> perturb(6, 0.0);
                    perturb[k] = eps_j;
                    poses[ib] = applyPerturbation(savedPoses[ib], perturb.data());
                } else {
                    sliderParams[k] += eps_j;
                }

                std::vector<double> resPert(m, 0.0);
                evalResiduals(resPert.data());
                for(int i = 0; i < m; i++) Jacobian[i * n + j] = (resPert[i] - residuals[i]) / eps_j;

                poses = savedPoses;
                sliderParams = savedSliders;
            }

            // 2. Solve for dq
            // For over-determined (m > n): normal equations (J^T J) dq = -J^T r + lambda*diag(J^T J)
            // For square (m == n): direct solve J dq = -r with Tikhonov
            if(m >= n){
                // Normal equations with LM damping
                std::vector<double> JTJ(n * n, 0.0);
                std::vector<double> JTr(n, 0.0);
                for(int i = 0; i < n; i++){
                    for(int j = 0; j < n; j++){
                        double s = 0;
                        for(int kk = 0; kk < m; kk++) s += Jacobian[kk*n+i] * Jacobian[kk*n+j];
                        JTJ[i*n+j] = s;
                    }
                    double s = 0;
                    for(int kk = 0; kk < m; kk++) s += Jacobian[kk*n+i] * (-residuals[kk]);
                    JTr[i] = s;
                }
                // LM damping
                for(int i = 0; i < n; i++) JTJ[i*n+i] += lambda * (JTJ[i*n+i] + 1e-12);

                std::vector<double> aug(n * (n+1), 0.0);
                for(int i = 0; i < n; i++){
                    for(int j = 0; j < n; j++) aug[i*(n+1)+j] = JTJ[i*n+j];
                    aug[i*(n+1)+n] = JTr[i];
                }
                bool singular = false;
                for(int k = 0; k < n; k++){
                    int piv = k;
                    for(int i = k+1; i < n; i++){
                        if(fabs(aug[i*(n+1)+k]) > fabs(aug[piv*(n+1)+k])) piv = i;
                    }
                    if(piv != k){ for(int j = 0; j <= n; j++) std::swap(aug[k*(n+1)+j], aug[piv*(n+1)+j]); }
                    double diag = aug[k*(n+1)+k];
                    if(fabs(diag) < 1e-15){ singular = true; break; }
                    for(int i = k+1; i < n; i++){
                        double f = aug[i*(n+1)+k] / diag;
                        for(int j = k; j <= n; j++) aug[i*(n+1)+j] -= f * aug[k*(n+1)+j];
                    }
                }
                if(singular){ lambda *= 10; continue; }
                for(int i = n-1; i >= 0; i--){
                    double s = aug[i*(n+1)+n];
                    for(int j = i+1; j < n; j++) s -= aug[i*(n+1)+j] * dq[j];
                    dq[i] = s / aug[i*(n+1)+i];
                }
            } else {
                // Under-determined: augmented system [J; sqrt(lambda)*I] * dq = [-r; 0]
                int nrows = m + n;
                std::vector<double> aug(nrows * (n + 1), 0.0);
                for(int i = 0; i < m; i++){
                    for(int j = 0; j < n; j++) aug[i*(n+1)+j] = Jacobian[i*n+j];
                    aug[i*(n+1)+n] = -residuals[i];
                }
                double sqrtLambda = sqrt(lambda);
                for(int i = 0; i < n; i++){
                    aug[(m+i)*(n+1)+i] = sqrtLambda;
                    aug[(m+i)*(n+1)+n] = 0;
                }
                bool singular = false;
                for(int k = 0; k < n; k++){
                    int piv = k;
                    for(int i = k+1; i < nrows; i++){
                        if(fabs(aug[i*(n+1)+k]) > fabs(aug[piv*(n+1)+k])) piv = i;
                    }
                    if(piv != k){ for(int j = 0; j <= n; j++) std::swap(aug[k*(n+1)+j], aug[piv*(n+1)+j]); }
                    double diag = aug[k*(n+1)+k];
                    if(fabs(diag) < 1e-15){ singular = true; break; }
                    for(int i = k+1; i < nrows; i++){
                        double f = aug[i*(n+1)+k] / diag;
                        if(f == 0) continue;
                        for(int j = k; j <= n; j++) aug[i*(n+1)+j] -= f * aug[k*(n+1)+j];
                    }
                }
                if(singular){ lambda *= 10; continue; }
                for(int i = n-1; i >= 0; i--){
                    double s = aug[i*(n+1)+n];
                    for(int j = i+1; j < n; j++) s -= aug[i*(n+1)+j] * dq[j];
                    dq[i] = s / aug[i*(n+1)+i];
                }
            }

            // 3. Line search: try step, accept if residual decreases
            std::vector<KinematicPose> savedPoses = poses;
            std::vector<double> savedSliders = sliderParams;

            double alpha = 1.0;
            bool accepted = false;
            for(int ls = 0; ls < 10; ls++){
                // Apply correction with step alpha
                poses = savedPoses;
                sliderParams = savedSliders;
                std::vector<std::vector<double>> bodyDQ(nBodies(), std::vector<double>(6, 0.0));
                std::vector<bool> bodyTouched(nBodies(), false);
                for(int j = 0; j < n; j++){
                    auto [ib, kk] = unknownMap[j];
                    if(ib >= 0){ bodyDQ[ib][kk] = alpha * dq[j]; bodyTouched[ib] = true; }
                    else       { sliderParams[kk] += alpha * dq[j]; }
                }
                for(int ib = 0; ib < nBodies(); ib++){
                    if(bodyTouched[ib]) poses[ib] = applyPerturbation(poses[ib], bodyDQ[ib].data());
                }

                std::vector<double> newRes(m, 0.0);
                evalResiduals(newRes.data());
                double newResSq = 0;
                for(int i = 0; i < m; i++) newResSq += newRes[i]*newRes[i];

                if(newResSq < curRes){
                    curRes = newResSq;
                    residuals = newRes;
                    accepted = true;
                    lambda = std::max(lambda * 0.5, 1e-12);
                    break;
                }
                alpha *= 0.5;
            }

            if(!accepted){
                // Line search failed, increase damping
                poses = savedPoses;
                sliderParams = savedSliders;
                lambda *= 10;
                if(lambda > 1e10){
                    if(bVerbose) printf("KinematicSolver: lambda exploded at iter %d, giving up\n", iter);
                    return -1;
                }
            }
        }

        evalResiduals(residuals.data());
        double maxRes = 0;
        for(int i = 0; i < m; i++) maxRes = std::max(maxRes, fabs(residuals[i]));
        if(bVerbose) printf("KinematicSolver: max residual after %d iters = %g\n", maxIters, maxRes);
        if(maxRes < tol * 100) return maxIters;
        return -1;
    }

    // Sweep slider parameters and output poses at each step
    // Each slider gets offset + t (mod 1), so sliders maintain their relative spacing
    // Sliders are fixed during solve — only body poses are solved for
    void sweepSolve( int nSteps, std::function<void(int step, const std::vector<KinematicPose>&)> callback,
                     int maxIters = 50, double tol = 1e-8, bool bVerbose = false ){
        std::vector<double> offsets = sliderParams;  // save initial params as offsets
        sliderFixed.assign(nSliders(), true);  // fix all sliders during sweep
        for(int step = 0; step <= nSteps; step++){
            double t = (double)step / nSteps;
            for(int si = 0; si < nSliders(); si++){
                sliderParams[si] = offsets[si] + t;
                if(sliderParams[si] >= 1.0) sliderParams[si] -= 1.0;
            }
            int iters = solve(maxIters, tol, bVerbose);
            if(bVerbose) printf("sweepSolve step %d (t=%.3f): iters=%d\n", step, t, iters);
            callback(step, poses);
        }
    }

    // Sweep from tStart to tEnd (no wrapping, for open paths)
    void sweepSolveRange( double tStart, double tEnd, int nSteps, std::function<void(int step, const std::vector<KinematicPose>&)> callback,
                          int maxIters = 50, double tol = 1e-8, bool bVerbose = false ){
        sliderFixed.assign(nSliders(), true);
        for(int step = 0; step <= nSteps; step++){
            double t = tStart + (tEnd - tStart) * (double)step / nSteps;
            for(int si = 0; si < nSliders(); si++) sliderParams[si] = t;
            int iters = solve(maxIters, tol, bVerbose);
            if(bVerbose) printf("sweepSolveRange step %d (t=%.3f): iters=%d\n", step, t, iters);
            callback(step, poses);
        }
    }
};

#endif // KinematicSolver_h
