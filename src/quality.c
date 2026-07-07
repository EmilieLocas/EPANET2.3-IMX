#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "mempool.h"
#include "types.h"
#include "funcs.h"

const double Q_STAGNANT = 0.005 / GPMperCFS;

double  findsourcequal(Project *, int, double, long);
/* BAM start */
int     findcrossjuncs(Project *pr);
/* BAM end */
/* IMX start */
double  imxadjoutconc(Project *pr, Cjunc *cj);
double  imxoppoutconc(Project *pr, Cjunc *cj);
void    assigncontamination(Project *pr);
double  getlinkangle(Project *pr, int lnk, int node);
void    assigncontaminationnode(Project *pr, int n);
/* IMX end */

extern char    setreactflag(Project *);
extern double  getucf(double);
extern void    ratecoeffs(Project *);
extern void    initsegs(Project *);
extern void    reversesegs(Project *, int);
extern int     sortnodes(Project *);
extern void    transport(Project *, long);

static double  sourcequal(Project *, Psource);
static void    evalmassbalance(Project *);
static double  findstoredmass(Project *);
static int     flowdirchanged(Project *);
/* BAM start */
static double angle(double x1, double y1, double x2, double y2);
static void getlinkcoords(Project *pr, int lnk, int node, double *x, double *y);
static int nonadjlink(Project *pr, int tolnk, int lnk2, int lnk3, int lnk4, int node);
/* BAM end */

int openqual(Project *pr)
{
    Network  *net = &pr->network;
    Quality *qual = &pr->quality;

    int errcode = 0;
    int n;

    if (qual->Qualflag == NONE) return errcode;

    if (net->Adjlist == NULL)
    {
        if (net->Nnodes < 2) return 223;
        if (net->Ntanks == 0) return 224;
    
        errcode = buildadjlists(net);
        if (errcode ) return errcode;

        if (errcode = unlinked(pr)) return errcode;
    }

    qual->OutOfMemory = FALSE;
    qual->SegPool = mempool_create();
    if (qual->SegPool == NULL) errcode = 101;

    n = net->Nlinks + 1;
    qual->FlowDir = (FlowDirection *)calloc(n, sizeof(FlowDirection));
    qual->PipeRateCoeff = (double *)calloc(n, sizeof(double));

    n = net->Nlinks + net->Ntanks + 1;
    qual->FirstSeg = (Pseg *)calloc(n, sizeof(Pseg));
    qual->LastSeg = (Pseg *)calloc(n, sizeof(Pseg));

    qual->SortedNodes = (int *)calloc(n, sizeof(int));

    /* BAM start */
    // Allocate memory for cross junctions
    qual->Crossjuncs = (Cjunc *)calloc(net->Njuncs + 1, sizeof(Cjunc));
    /* BAM end */

    ERRCODE(MEMCHECK(qual->FlowDir));
    ERRCODE(MEMCHECK(qual->PipeRateCoeff));
    ERRCODE(MEMCHECK(qual->FirstSeg));
    ERRCODE(MEMCHECK(qual->LastSeg));
    ERRCODE(MEMCHECK(qual->SortedNodes));
    /* BAM start */
    ERRCODE(MEMCHECK(qual->Crossjuncs));
    /* BAM end */
    return errcode;
}

int initqual(Project *pr)
{
    Network *net = &pr->network;
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;
    Times   *time = &pr->times;

    int i;
    int errcode = 0;

    if (!hyd->OpenHflag)
    {
        fseek(pr->outfile.HydFile, pr->outfile.HydOffset, SEEK_SET);
    }

    time->Qtime = 0;
    time->Htime = 0;
    time->Rtime = time->Rstart;
    pr->report.Nperiods = 0;

    for (i = 1; i <= net->Nnodes; i++)
    {
        if (qual->Qualflag == TRACE) qual->NodeQual[i] = 0.0;
        else                         qual->NodeQual[i] = net->Node[i].C0;
        if (net->Node[i].S != NULL) net->Node[i].S->Smass = 0.0;
    }
    if (qual->Qualflag == NONE) return errcode;

    for (i = 1; i <= net->Ntanks; i++)
    {
        net->Tank[i].C = qual->NodeQual[net->Tank[i].Node];
    }

    if (qual->Qualflag == TRACE) qual->NodeQual[qual->TraceNode] = 100.0;

    if (qual->Diffus > 0.0) qual->Sc = hyd->Viscos / qual->Diffus;
    else                    qual->Sc = 0.0;

    qual->Bucf = getucf(qual->BulkOrder);
    qual->Tucf = getucf(qual->TankOrder);

    qual->Reactflag = setreactflag(pr);

    qual->FreeSeg = NULL;
    mempool_reset(qual->SegPool);

    initsegs(pr);

    for (i = 1; i <= net->Nlinks; i++) qual->FlowDir[i] = ZERO_FLOW;

    qual->Wbulk = 0.0;
    qual->Wwall = 0.0;
    qual->Wtank = 0.0;
    qual->Wsource = 0.0;

    qual->MassBalance.initial = findstoredmass(pr);
    qual->MassBalance.inflow = 0.0;
    qual->MassBalance.outflow = 0.0;
    qual->MassBalance.reacted = 0.0;
    qual->MassBalance.final = 0.0;
    qual->MassBalance.ratio = 0.0;
    qual->MassBalance.segCount = 0;
    return errcode;
}

int runqual(Project *pr, long *t)
{
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;
    Times   *time = &pr->times;

    long hydtime = 0; 
    long hydstep = 0;
    int errcode = 0;

    *t = time->Qtime;

    if (time->Qtime == time->Htime)
    {
        if (!hyd->OpenHflag)
        {
            if (!readhyd(pr, &hydtime)) return 307;
            if (!readhydstep(pr, &hydstep)) return 307;
            time->Htime = hydtime;
        }

        if (time->Htime >= time->Rtime)
        {
            if (pr->outfile.Saveflag)
            {
                errcode = saveoutput(pr);
                pr->report.Nperiods++;
            }
            time->Rtime += time->Rstep;
        }
        if (errcode) return errcode;

        if (qual->Qualflag != NONE && time->Qtime < time->Dur)
        {
            if (qual->Reactflag && qual->Qualflag != AGE) ratecoeffs(pr);

            if (flowdirchanged(pr) == TRUE)
            {
                errcode = sortnodes(pr);
            }
        }
        if (!hyd->OpenHflag) time->Htime = hydtime + hydstep;
    }
    return errcode;
}

int nextqual(Project *pr, long *tstep)
{
    Quality *qual = &pr->quality;
    Times   *time = &pr->times;

    long hydstep;
    long dt, qtime;
    int errcode = 0;

    *tstep = 0;
    hydstep = 0;
    if (time->Htime <= time->Dur) hydstep = time->Htime - time->Qtime;

    if (qual->Qualflag != NONE && hydstep > 0)
    {
        qtime = 0;
        while (!qual->OutOfMemory && qtime < hydstep)
        {
            dt = MIN(time->Qstep, hydstep - qtime);
            qtime += dt;
            transport(pr, dt);
        }
        if (qual->OutOfMemory) errcode = 101;
    }

    evalmassbalance(pr);

    if (!errcode) *tstep = hydstep;
    time->Qtime += hydstep;

    if (!errcode && *tstep == 0)
    {
        if (qual->Qualflag != NONE && pr->report.Statflag)
        {
            writemassbalance(pr);
        }

        if (pr->outfile.Saveflag) errcode = savefinaloutput(pr);
    }
    return errcode;
}

int stepqual(Project *pr, long *tleft)
{
    Quality *qual = &pr->quality;
    Times   *time = &pr->times;

    long dt, hstep, t, tstep;
    int errcode = 0;

    tstep = time->Qstep;
    do
    {
        dt = tstep;

        hstep = time->Htime - time->Qtime;

        if (hstep < dt)
        {
            dt = hstep;

            if (qual->Qualflag != NONE) transport(pr, dt);
            time->Qtime += dt;

            if (pr->hydraul.OpenHflag) break;

            errcode = runqual(pr, &t);
            time->Qtime = t;
        }

        else
        {
            if (qual->Qualflag != NONE) transport(pr, dt);
            time->Qtime += dt;
        }

        tstep -= dt;
        if (qual->OutOfMemory) errcode = 101;

    } while (!errcode && tstep > 0);

    evalmassbalance(pr);

    *tleft = time->Dur - time->Qtime;

    if (!errcode && *tleft == 0)
    {
        if (qual->Qualflag != NONE && pr->report.Statflag)
        {
            writemassbalance(pr);
        }

        if (pr->outfile.Saveflag) errcode = savefinaloutput(pr);
    }
    return errcode;
}

int closequal(Project *pr)
{
    Quality *qual = &pr->quality;
    int errcode = 0;

    if (qual->Qualflag != NONE)
    {
        if (qual->SegPool) mempool_delete(qual->SegPool);
        FREE(qual->FirstSeg);
        FREE(qual->LastSeg);
        FREE(qual->PipeRateCoeff);
        FREE(qual->FlowDir);
        FREE(qual->SortedNodes);
        /* BAM start */
        FREE(qual->Crossjuncs);
        /* BAM end */
    }
    freeadjlists(&pr->network);
    return errcode;
}

double avgqual(Project *pr, int k)
{
    Network  *net = &pr->network;
    Quality  *qual = &pr->quality;

    double vsum = 0.0, msum = 0.0;
    Pseg seg;

    if (qual->Qualflag == NONE) return 0.0;

    if (qual->FirstSeg != NULL)
    {
        seg = qual->FirstSeg[k];
        while (seg != NULL)
        {
            vsum += seg->v;
            msum += (seg->c) * (seg->v);
            seg = seg->prev;
        }
    }

    if (vsum > 0.0) return (msum / vsum);

    else
    {
        return ((qual->NodeQual[net->Link[k].N1] +
            qual->NodeQual[net->Link[k].N2]) / 2.);
    }
}

double findsourcequal(Project *pr, int n, double volout, long tstep)
{
    Network *net = &pr->network;
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;
    Times   *time = &pr->times;

    double massadded = 0.0, c;
    Psource source;

    if (qual->Qualflag != CHEM) return 0.0;

    source = net->Node[n].S;
    if (source == NULL)    return 0.0;
    if (source->C0 == 0.0) return 0.0;
    if (volout / tstep <= Q_STAGNANT) return 0.0;

    c = sourcequal(pr, source);
    switch (source->Type)
    {
        case CONCEN:
        if (net->Node[n].Type == JUNCTION)
        {
            if (hyd->NodeDemand[n] < 0.0)
            {
                c = -c * hyd->NodeDemand[n] * tstep / volout;
            }
            else c = 0.0;
        }
        break;

        case MASS:
            c = c * tstep / volout;
            break;

        case SETPOINT:
            c = MAX(c - qual->NodeQual[n], 0.0);
            break;

        case FLOWPACED:
            break;
    }

    massadded = c * volout;

    source->Smass += massadded;

    if (time->Htime >= time->Rstart)
    {
        qual->Wsource += massadded;
    }
    return c;
}

double sourcequal(Project *pr, Psource source)
{
    Network *net = &pr->network;
    Times   *time = &pr->times;

    int i;
    long k;
    double c;

    c = source->C0;

    if (source->Type == MASS) c /= 60.0;
    else                      c /= pr->Ucf[QUALITY];

    i = source->Pat;
    if (i == 0)  return c;
    k = ((time->Qtime + time->Pstart) / time->Pstep) %
        (long)net->Pattern[i].Length;
    return (c * net->Pattern[i].F[k]);
}

void  evalmassbalance(Project *pr)
{
    Quality *qual = &pr->quality;

    double massin;
    double massout;
    double massreacted;

    if (qual->Qualflag == NONE) qual->MassBalance.ratio = 1.0;
    else
    {
        qual->MassBalance.final = findstoredmass(pr);
        massin = qual->MassBalance.initial + qual->MassBalance.inflow;
        massout = qual->MassBalance.outflow + qual->MassBalance.final;
        massreacted = qual->MassBalance.reacted;
        if (massreacted > 0.0) massout += massreacted;
        else                   massin -= massreacted;
        if (massin == 0.0) qual->MassBalance.ratio = 1.0;
        else               qual->MassBalance.ratio = massout / massin;
    }
}

double  findstoredmass(Project *pr)
{
    Network  *net = &pr->network;
    Quality  *qual = &pr->quality;

    int    i, k;
    double totalmass = 0.0;
    Pseg   seg;

    for (k = 1; k <= net->Nlinks; k++)
    {
        seg = qual->FirstSeg[k];
        while (seg != NULL)
        {
            totalmass += (seg->c) * (seg->v);
            seg = seg->prev;
        }
    }

    for (i = 1; i <= net->Ntanks; i++)
    {
        if (net->Tank[i].A == 0.0) continue;

        else
        {
            k = net->Nlinks + i;
            seg = qual->FirstSeg[k];
            while (seg != NULL)
            {
                totalmass += seg->c * seg->v;
                seg = seg->prev;
            }
        }
    }
    return totalmass;
}

int flowdirchanged(Project *pr)
{
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;

    int k;
    int result = FALSE;
    int newdir;
    int olddir;
    double q;

    for (k = 1; k <= pr->network.Nlinks; k++)
    {
        olddir = qual->FlowDir[k];
        q = (hyd->LinkStatus[k] <= CLOSED) ? 0.0 : hyd->LinkFlow[k];
        newdir = SGN(q);

        if (fabs(q) < Q_STAGNANT) newdir = 0;

        if (newdir * olddir < 0) reversesegs(pr, k);

        if (newdir != olddir) result = TRUE;

        qual->FlowDir[k] = newdir;
    }
    return result;
}

/* BAM start */
int findcrossjuncs(Project *pr)
{
    Network *net = &pr->network;
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;

    int j, k;
    int Nupnode, Ndownnode;
    int inlink1, inlink2, outlink1, outlink2;
    int opplink1, opplink2;

    for (j = 1; j <= net->Njuncs; j++)
    {
        if (net->Node[j].X == MISSING || net->Node[j].Y == MISSING)
            return 256;
    }

    int prevNcross = qual->Ncrossjuncs;
    qual->Ncrossjuncs = 0;

    for (j = 1; j <= net->Njuncs; j++)
    {
        Nupnode = 0;
        Ndownnode = 0;

        for (k = 1; k <= net->Nlinks; k++)
        {
            if (net->Link[k].Len == 0.0) continue;

            if (hyd->LinkFlow[k] > 0.0 && net->Link[k].N1 == j)
            {
                Nupnode++;
                if (Nupnode == 1) outlink1 = k;
                else if (Nupnode == 2) outlink2 = k;
            }
            else if (hyd->LinkFlow[k] > 0.0 && net->Link[k].N2 == j)
            {
                Ndownnode++;
                if (Ndownnode == 1) inlink1 = k;
                else if (Ndownnode == 2) inlink2 = k;
            }
            else if (hyd->LinkFlow[k] < 0.0 && net->Link[k].N2 == j)
            {
                Nupnode++;
                if (Nupnode == 1) outlink1 = k;
                else if (Nupnode == 2) outlink2 = k;
            }
            else if (hyd->LinkFlow[k] < 0.0 && net->Link[k].N1 == j)
            {
                Ndownnode++;
                if (Ndownnode == 1) inlink1 = k;
                else if (Ndownnode == 2) inlink2 = k;
            }
        }

        qual->Crossjuncs[j].iscrossjunc = 0;

        if (Nupnode == 2 && Ndownnode == 2)
        {
            opplink1 = nonadjlink(pr, inlink1, inlink2, outlink1, outlink2, j);

            if (opplink1 != inlink2)
            {
                opplink2 = nonadjlink(pr, inlink2, inlink1, outlink1, outlink2, j);

                double a_in1 = getlinkangle(pr, inlink1, j);
                double a_in2 = getlinkangle(pr, inlink2, j);
                double crossAngle = fabs(a_in2 - a_in1);
                if (crossAngle > M_PI) crossAngle = 2.0 * M_PI - crossAngle;

                double TOL = 5.0 * M_PI / 180.0;
                if (fabs(crossAngle - M_PI / 2.0) > TOL) continue;

                int cl = inlink1, pl = inlink2;
                if ((fabs(hyd->LinkFlow[inlink2]) + fabs(hyd->LinkFlow[opplink2])) >
                    (fabs(hyd->LinkFlow[inlink1]) + fabs(hyd->LinkFlow[opplink1])))
                {
                    cl = inlink2;
                    pl = inlink1;
                }

                double a_cl      = getlinkangle(pr, cl, j);
                double diff_out1 = fabs(fmod(fabs(getlinkangle(pr, outlink1, j) - a_cl), 2.0*M_PI) - M_PI);
                double diff_out2 = fabs(fmod(fabs(getlinkangle(pr, outlink2, j) - a_cl), 2.0*M_PI) - M_PI);

                int opp_out = (diff_out1 < diff_out2) ? outlink1 : outlink2;
                int adj_out = (diff_out1 < diff_out2) ? outlink2 : outlink1;

                qual->Ncrossjuncs++;
                qual->Crossjuncs[j].iscrossjunc  = 1;
                qual->Crossjuncs[j].nodeindex    = j;
                qual->Crossjuncs[j].contaminlink = cl;
                qual->Crossjuncs[j].purelink     = pl;
                qual->Crossjuncs[j].oppoutlink   = opp_out;
                qual->Crossjuncs[j].adjoutlink   = adj_out;
            }
        }
    }
    if (qual->Ncrossjuncs != prevNcross)
    {
        printf("Nombre de cross-junctions detectees: %d\n", qual->Ncrossjuncs);
        for (j = 1; j <= net->Njuncs; j++)
        {
            if (!qual->Crossjuncs[j].iscrossjunc) continue;
            printf("  Noeud %d | contaminlink=%d purelink=%d | adjoutlink=%d oppoutlink=%d\n",
                j, qual->Crossjuncs[j].contaminlink, qual->Crossjuncs[j].purelink,
                qual->Crossjuncs[j].adjoutlink, qual->Crossjuncs[j].oppoutlink);
        }
    }
    return 0;
}

/* IMX start */
double imxadjoutconc(Project *pr, Cjunc *cj)
{
    Hydraul *hyd = &pr->hydraul;

    double Q1 = fabs(hyd->LinkFlow[cj->contaminlink]);
    double Q2 = fabs(hyd->LinkFlow[cj->purelink]);
    double Q3 = fabs(hyd->LinkFlow[cj->adjoutlink]);
    double Ccont = cj->contaminconc;
    double Cpur  = cj->pureconc;
    double ratio = Q3 / Q1;

    if (fabs(Ccont - Cpur) < 1e-10) return Ccont;

    double C3_star;
    if (ratio <= 0.85)
        C3_star = 1.0;
    else
        C3_star = 0.22 * log(ratio) + 0.91 * pow(Q3/Q2, -0.79);

    /* IMX : borner C3* entre 0 et 1 pour éviter concentrations
       hors bornes physiques */
    if (C3_star < 0.0) C3_star = 0.0;
    if (C3_star > 1.0) C3_star = 1.0;

    double result = C3_star * (Ccont - Cpur) + Cpur;

    /* Clamp final entre Cpur et Ccont */
    if (result < Cpur)  result = Cpur;
    if (result > Ccont) result = Ccont;

    return result;
}

double imxoppoutconc(Project *pr, Cjunc *cj)
{
        Hydraul *hyd = &pr->hydraul;

    double Q1 = fabs(hyd->LinkFlow[cj->contaminlink]);
    double Q2 = fabs(hyd->LinkFlow[cj->purelink]);
    double Q3 = fabs(hyd->LinkFlow[cj->adjoutlink]);
    double Q4 = fabs(hyd->LinkFlow[cj->oppoutlink]);
    double Ccont = cj->contaminconc;
    double Cpur  = cj->pureconc;

    if (fabs(Ccont - Cpur) < 1e-10) return Ccont;
    if (Q4 < 1e-10) return Cpur;

    double C3 = imxadjoutconc(pr, cj);

    /* Bilan de masse : C4 = (C1*Q1 + C2*Q2 - C3*Q3) / Q4 */
    double result = (Ccont * Q1 + Cpur * Q2 - C3 * Q3) / Q4;

    /* Clamp physique : la concentration ne peut pas sortir
       de l'intervalle [Cpur, Ccont] */
    if (result < Cpur)  result = Cpur;
    if (result > Ccont) result = Ccont;

    return result;
}
/* IMX end */

double angle(double x1, double y1, double x2, double y2)
{
    double relptx = x1 - x2;
    double relpty = y1 - y2;
    double arctanpt = atan(relpty / relptx);

    if ((relptx >= 0.0) && (relpty >= 0.0)) return arctanpt; // 1st quadrant
    else if ((relptx <= 0.0) && (relpty >= 0.0)) return arctanpt + M_PI; // 2nd quadrant
    else if ((relptx < 0.0) && (relpty <= 0.0)) return arctanpt + M_PI; // 3rd quadrant
    else if ((relptx >= 0.0) && (relpty <= 0.0)) return arctanpt + 2 * M_PI; // 4th quadrant
    else return 0.0; // should never happen
}

void getlinkcoords(Project *pr, int lnk, int node, double *x, double *y)
{
    Network *net = &pr->network;
    Slink *link = &net->Link[lnk];
    Pvertices verts = link->Vertices;
    
    if (verts != NULL && verts->Npts > 0)
    {
        if (link->N1 == node) 
        {
            *x = verts->X[0]; 
            *y = verts->Y[0]; 
        }
        else
        { 
            *x = verts->X[verts->Npts-1]; 
            *y = verts->Y[verts->Npts-1]; 
        }
    }
    else
    {
        int other = (link->N1 == node) ? link->N2 : link->N1;
        *x = net->Node[other].X;
        *y = net->Node[other].Y;
    }
}

int nonadjlink(Project *pr, int tolnk, int lnk2, int lnk3, int lnk4, int node)
{
    /* IMX start */
    double a1 = getlinkangle(pr, tolnk, node);
    double a2 = getlinkangle(pr, lnk2,  node);
    double a3 = getlinkangle(pr, lnk3,  node);
    double a4 = getlinkangle(pr, lnk4,  node);
    /* IMX end */

    a2 -= a1;
    a3 -= a1;
    a4 -= a1;

    if (a2 < 0.0) a2 += (2.0 * M_PI);
    if (a3 < 0.0) a3 += (2.0 * M_PI);
    if (a4 < 0.0) a4 += (2.0 * M_PI);

    if (((a2 >= a3) && (a2 <= a4)) || ((a2 <= a3) && (a2 >= a4))) return lnk2;
    if (((a3 >= a2) && (a3 <= a4)) || ((a3 <= a2) && (a3 >= a4))) return lnk3;
    if (((a4 >= a2) && (a4 <= a3)) || ((a4 <= a2) && (a4 >= a3))) return lnk4;
    return lnk2;
}
/* BAM end */

/* IMX start */
double getlinkangle(Project *pr, int lnk, int node)
{
    double x, y;
    getlinkcoords(pr, lnk, node, &x, &y);
    return angle(x, y, pr->network.Node[node].X, pr->network.Node[node].Y);
}

void assigncontaminationnode(Project *pr, int n)
{
    Quality *qual = &pr->quality;
    Cjunc *cj = &qual->Crossjuncs[n];

    if (!cj->iscrossjunc) return;

    if (cj->pureconc > cj->contaminconc)
    {
        int tmplink      = cj->contaminlink;
        cj->contaminlink = cj->purelink;
        cj->purelink     = tmplink;
        double tmpc      = cj->contaminconc;
        cj->contaminconc = cj->pureconc;
        cj->pureconc     = tmpc;

        int oldadj = cj->adjoutlink;
        int oldopp = cj->oppoutlink;
        double a_contam = getlinkangle(pr, cj->contaminlink, n);
        double a_adj    = getlinkangle(pr, oldadj, n);
        double a_opp    = getlinkangle(pr, oldopp, n);
        double diff_adj = fabs(fmod(fabs(a_adj - a_contam), 2*M_PI) - M_PI);
        double diff_opp = fabs(fmod(fabs(a_opp - a_contam), 2*M_PI) - M_PI);

        if (diff_adj < diff_opp)
        {
            cj->oppoutlink = oldadj;
            cj->adjoutlink = oldopp;
        }
    }
}
/* IMX end */
