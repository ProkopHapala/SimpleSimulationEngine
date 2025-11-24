#include <stdio.h>
#include <stdlib.h>

#include "globals.h"

#include "IO_utils.h"
#include "argparse.h"

#include "TrussDynamics_d.h"

static void print_usage(const char* prog){
    printf("Usage: %s -i <truss.txt> -o <traj.xyz> [-n steps] [-dt dt] [-psave k] [-v verbosity]\n", prog);
}

static void write_xyz_frame(FILE* f, const TrussDynamics_d& sim, int istep, double time){
    int np = sim.nPoint;
    fprintf(f, "%d\n", np);
    fprintf(f, "step %d time %g\n", istep, time);
    for(int i=0; i<np; i++){
        const Quat4d& p = sim.points[i];
        fprintf(f, "% .16g % .16g % .16g\n", p.f.x, p.f.y, p.f.z);
    }
}

int main(int argc, char** argv){
    setbuf(stdout, NULL);

    const char* inTruss  = NULL;
    const char* outTraj  = "traj.xyz";
    int   nSteps  = 1000;
    double dt     = 5.0e-4;   // default consistent with TrussDynamics_d dt scale
    int   perSave = 10;
    int   method  = 0;        // reserved for future (different solvers)

    LambdaDict opts;
    opts["-i"]    = {1, [&](const char** ss){ inTruss  = ss[0]; if(verbosity>0) printf("ARG -i %s\n", inTruss);  }};
    opts["-o"]    = {1, [&](const char** ss){ outTraj  = ss[0]; if(verbosity>0) printf("ARG -o %s\n", outTraj);  }};
    opts["-n"]    = {1, [&](const char** ss){ int v=0; sscanf(ss[0], "%d", &v); nSteps  = v; if(verbosity>0) printf("ARG -n %d\n", nSteps); }};
    opts["-dt"]   = {1, [&](const char** ss){ double v=0; sscanf(ss[0], "%lf", &v); dt = v; if(verbosity>0) printf("ARG -dt %g\n", dt); }};
    opts["-psave"]={1, [&](const char** ss){ int v=0; sscanf(ss[0], "%d", &v); perSave = v; if(verbosity>0) printf("ARG -psave %d\n", perSave); }};
    opts["-m"]    = {1, [&](const char** ss){ int v=0; sscanf(ss[0], "%d", &v); method  = v; if(verbosity>0) printf("ARG -m %d\n", method); }};
    opts["-v"]    = {1, [&](const char** ss){ int v=1; sscanf(ss[0], "%d", &v); verbosity = v; printf("ARG -v %d\n", verbosity); }};

    process_args(argc, argv, opts);

    if(!inTruss){
        print_usage(argv[0]);
        return 1;
    }

    if(verbosity>0){
        printf("trussSimBatch: in='%s' out='%s' nSteps=%d dt=%g perSave=%d method=%d verbosity=%d\n",
               inTruss, outTraj, nSteps, dt, perSave, method, verbosity);
    }

    TrussDynamics_d sim;
    loadSimFromFile(inTruss, sim);
    if(sim.nPoint <= 0 || sim.nBonds <= 0){
        printf("ERROR: trussSimBatch: loaded sim has nPoint=%d nBonds=%d\n", sim.nPoint, sim.nBonds);
        return 2;
    }

    // basic integrator selection hook (for now, just use run)
    double damp = sim.damping;

    FILE* fout = fopen(outTraj, "w");
    if(!fout){
        printf("ERROR: trussSimBatch: cannot open trajectory '%s' for writing\n", outTraj);
        return 3;
    }

    // initial frame at step 0
    write_xyz_frame(fout, sim, 0, 0.0);

    double time = 0.0;
    for(int step=1; step<=nSteps; step++){
        // One time step of dynamics; method switch reserved for future
        switch(method){
            default:
                sim.run(1, dt, damp);
                break;
        }
        time += dt;

        if(perSave > 0 && (step % perSave) == 0){
            write_xyz_frame(fout, sim, step, time);
        }
    }

    fclose(fout);

    if(verbosity>0){
        printf("trussSimBatch: DONE, wrote trajectory '%s'\n", outTraj);
    }

    return 0;
}
