/*
 *  A minimal C interface to Pythia 6.4
 *
 *  I. Volobouev
 *  April 2008
 */

#ifndef PYTHIA6_H_
#define PYTHIA6_H_

#ifdef __cplusplus
extern "C" {
#endif

/* The PYJETS common block -- main event record in PYTHIA
 *
 * n            The number of particles in the event
 *
 * npad         Unused
 *
 * k[0][i]      Gives information on whether or not a parton or particle
 *                number i has fragmented or decayed
 * k[1][i]      Gives the particle code
 * k[2][i]      Gives the particle origin
 * k[3,4][i]    These give the positions of fragmentation/decay products
 *
 * p[0,1,2][i]  These give the x, y, and z components of the momentum (GeV/c)
 * p[3][i]      Particle energy (GeV)
 * p[4][i]      Particle mass (GeV/c^2)
 *
 * v[0,1,2][i]  Production vertex location (mm)
 * v[3][i]      Production vertex time (units mm/c, approximately 3.33e-12 s)
 * v[4][i]      Particle invariant lifetime (c*tau, distance from
 *                production to decay in mm)
 */
extern struct {
    int n, npad;
    int k[5][4000];
    double p[5][4000], v[5][4000];
} pyjets_;

extern struct {
    int idbmup[2];
    double ebmup[2];
    int pdfgup[2], pdfsup[2], idwtup, nprup;
    double xsecup[100], xerrup[100], xmaxup[100];
    int lprup[100];
} heprup_;

extern struct {
    int nup, idprup;
    double xwgtup, scalup, aqedup, aqcdup;
    int idup[500], istup[500], mothup[500][2], icolup[500][2];
    double pup[500][5], vtimup[500], spinup[500];
} hepeup_;

// Functions with C-style semantics equivalent to
// their corresponding Fortran analogs
void pyinit(const char*, const char*, const char*, double);
void pygive(const char*);
void pyupin(void);
void pyupev(void);
void pyevnt(void);
void pyevnw(void);
void pylist(int);
void pystat(int);
void pyedit(int);
void py1ent(int, int, double e, double theta, double phi);
void py2ent(int, int, int, double ecms);
void pyname(int, char name[16]);
double pyr(void);
double pymass(int);

/* The following function returns three times the charge of a particle */
int pychge(int);


#ifdef __cplusplus
}
#endif

#endif /* PYTHIA6_H_ */
