#include <string.h>
#include "pythia6.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void pyupin_(void);
extern void pyupev_(void);
extern void pylist_(int*);
extern void pystat_(int*);
extern void pyedit_(int*);
extern void py1ent_(int*, int*, double*, double*, double*);
extern void py2ent_(int*, int*, int*, double*);
extern void pygive_(const char*, int);
extern void pyinit_(const char*, const char*, const char*, double*,
                    int, int, int);
extern void pyname_(int*, char*);
extern double pyr_(int*);
extern double pymass_(int*);
extern int pychge_(int*);
extern void pyevnt_(void);
extern void pyevnw_(void);

void pyinit(const char* s1, const char* s2, const char* s3, double d)
{
    pyinit_(s1, s2, s3, &d, strlen(s1), strlen(s2), strlen(s3));
}

void pyedit(int i)
{
    pyedit_(&i);
}

void pygive(const char* p)
{
    pygive_(p, strlen(p));
}

void py2ent(int cmd, int p1code, int p2code, double ecms)
{
    py2ent_(&cmd, &p1code, &p2code, &ecms);
}

void pylist(int i)
{
    pylist_(&i);
}
    
void pystat(int i)
{
    pystat_(&i);
}

void pyupev(void)
{
    pyupev_();
}

void pyupin(void)
{
    pyupin_();
}

void pyevnt(void)
{
    pyevnt_();
}

void pyevnw(void)
{
    pyevnw_();
}

void pyname(int i, char name[16])
{
    pyname_(&i, &name[0]);
    name[15] = '\0';
}

double pyr(void)
{
    int i = 0;
    return pyr_(&i);
}

double pymass(int i)
{
    return pymass_(&i);
}

int pychge(int i)
{
    return pychge_(&i);
}

void py1ent(int i, int j, double e, double theta, double phi)
{
    py1ent_(&i, &j, &e, &theta, &phi);
}




#ifdef __cplusplus
}
#endif
