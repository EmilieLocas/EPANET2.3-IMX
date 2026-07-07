/*
******************************************************************************
Project:      OWA EPANET
Version:      2.3
Module:       qualroute.c
Description:  computes water quality transport over a single time step
Authors:      see AUTHORS
Copyright:    see AUTHORS
License:      see LICENSE
Last Updated: 02/14/2025
******************************************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "mempool.h"
#include "types.h"

// Macro to compute the volume of a link
#define LINKVOL(k) (0.785398 * net->Link[(k)].Len * SQR(net->Link[(k)].Diam))

// Macro to get link flow compatible with flow saved to hydraulics file
// Structure hydraul
#define LINKFLOW(k) ((hyd->LinkStatus[k] <= CLOSED) ? 0.0 : hyd->LinkFlow[k])

// Exported functions
int     sortnodes(Project *);
void    transport(Project *, long);
void    initsegs(Project *);
void    reversesegs(Project *, int);
void    addseg(Project *, int, double, double);

// Imported functions
extern double  findsourcequal(Project *, int, double, long);
extern void    reactpipes(Project *, long);
extern void    reacttanks(Project *, long);
extern double  mixtank(Project *, int, double, double, double);
/* BAM start */
extern int     findcrossjuncs(Project *pr);
extern double getlinkangle(Project *pr, int lnk, int node);
/* BAM end */
/* IMX start */
extern double  imxadjoutconc(Project *pr, Cjunc *cj);
extern double  imxoppoutconc(Project *pr, Cjunc *cj);
extern void    assigncontamination(Project *pr);
extern void    assigncontaminationnode(Project *pr, int n);
/* IMX end */

// Local functions
static void    evalnodeinflow(Project *, int, long, double *, double *);
static void    evalnodeoutflow(Project *, int, double, long);
static double  findnodequal(Project *, int, double, double, double, long);
static double  noflowqual(Project *, int);
static void    updatemassbalance(Project *, int, double, double, long);
static int     selectnonstacknode(Project *, int, int *);


void transport(Project *pr, long tstep)
{
    Network *net = &pr->network;
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;

    int j, k, m, n;
    double volin, massin, volout, nodequal;
    Padjlist alink;

    if (qual->Reactflag)
    {
        reactpipes(pr, tstep);
        reacttanks(pr, tstep);
    }

    if (findcrossjuncs(pr) > 0) return;

    for (j = 1; j <= net->Nnodes; j++)
    {
        n = qual->SortedNodes[j];
        volin  = 0.0;
        massin = 0.0;
        volout = 0.0;

        for (alink = net->Adjlist[n]; alink != NULL; alink = alink->next)
        {
            k = alink->link;
            m = (qual->FlowDir[k] < 0) ? net->Link[k].N1 : net->Link[k].N2;
            if (m == n)
            {
                evalnodeinflow(pr, k, tstep, &volin, &massin);
            }
            else volout += fabs(LINKFLOW(k));
        }

        if (net->Node[n].Type == JUNCTION)
            volout += MAX(0.0, hyd->NodeDemand[n]);
        volout *= tstep;

        /* IMX : capturer les concentrations des liens entrants APRES evalnodeinflow,
           en lisant NodeQual du nœud amont (déjà calculé en ordre topologique) */
        if (n <= net->Njuncs && qual->Crossjuncs[n].iscrossjunc == 1)
        {
            Cjunc *cj = &qual->Crossjuncs[n];
            int inlinks[2] = { cj->contaminlink, cj->purelink };
            double concs[2];

            for (int ii = 0; ii < 2; ii++)
            {
                k = inlinks[ii];
                /* LastSeg contient la concentration IMX poussée par evalnodeoutflow
                du nœud amont — c'est la valeur correcte après mélange incomplet. */
                if (qual->LastSeg[k] != NULL)
                    concs[ii] = qual->LastSeg[k]->c;
                else
                {
                    int upnode = (qual->FlowDir[k] >= 0)
                                ? net->Link[k].N1
                                : net->Link[k].N2;
                    concs[ii] = qual->NodeQual[upnode];
                }
            }
            cj->contaminconc = concs[0];
            cj->pureconc     = concs[1];

            assigncontaminationnode(pr, n);
        }

        nodequal = findnodequal(pr, n, volin, massin, volout, tstep);

        for (alink = net->Adjlist[n]; alink != NULL; alink = alink->next)
        {
            k = alink->link;
            m = (qual->FlowDir[k] < 0) ? net->Link[k].N2 : net->Link[k].N1;
            if (m == n)
            {
                double outconc = nodequal;

                if (n <= net->Njuncs && qual->Crossjuncs[n].iscrossjunc == 1)
                {
                    Cjunc *cj = &qual->Crossjuncs[n];
                    if (k == cj->adjoutlink)
                        outconc = imxadjoutconc(pr, cj);
                    else if (k == cj->oppoutlink)
                        outconc = imxoppoutconc(pr, cj);
                }
                if (n <= net->Njuncs && qual->Crossjuncs[n].iscrossjunc == 1)
                {
                    Cjunc *cj = &qual->Crossjuncs[n];
                    double ucf = pr->Ucf[QUALITY];
                    printf("Cross-junc n=%d | contaminlink=%d purelink=%d | contaminconc=%.4f pureconc=%.4f mg/L\n",
                        n, cj->contaminlink, cj->purelink,
                        cj->contaminconc / ucf, cj->pureconc / ucf);
                    printf("  adjoutlink=%d oppoutlink=%d\n", cj->adjoutlink, cj->oppoutlink);
                    printf("  IMX_adj=%.4f IMX_opp=%.4f mg/L\n",
                        imxadjoutconc(pr, cj) / ucf, imxoppoutconc(pr, cj) / ucf);
                    printf("  NodeQual[upnode_contam]=%.4f NodeQual[upnode_pur]=%.4f mg/L\n",
                        qual->NodeQual[net->Link[cj->contaminlink].N1] / ucf,
                        qual->NodeQual[net->Link[cj->purelink].N1] / ucf);
                }
                evalnodeoutflow(pr, k, outconc, tstep);
            }
        }
        updatemassbalance(pr, n, massin, volout, tstep);
    }
}

void  evalnodeinflow(Project *pr, int k, long tstep, double *volin,
                     double *massin)
// Cumule masse et volume arrivant au noeud/jonction
/*
**--------------------------------------------------------------
**   Input:   k = link index -> unique, integer-based identifier (from 1) 
**   assigned to each link (pipe, pump, or valve)
**            tstep = quality routing time step (time interval where software
**                     computes/updates results)
**   Output:  volin = flow volume entering a node/junction
**            Flow volume : volume of a fluid that passes through a given surface 
**            per unit of time ==>  Flowrate in article (l/s) **Careful : units
**            massin = constituent mass entering a node/junction
**            Dimensionless concentration ? ==> *Concentration = massin / volin ?
**            **Pressure does not have a significant impact on mixing phenomenon
**            **Output ==> dimensionless concentration or volin and massin ?
**   Purpose: adds the contribution of a link's (pipe's) outflow volume
**            and constituent mass to the total inflow into its
**            downstream node/junction over a time step.
**--------------------------------------------------------------
*/
{
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;

    double q, v, vseg;
    Pseg seg;

    // Get flow rate (q) and flow volume (v) through link
    q = LINKFLOW(k);
    v = fabs(q) * tstep;

    // Transport flow volume v from link's leading segments into downstream
    // node, removing segments once their full volume is consumed
    while (v > 0.0)
    {
        seg = qual->FirstSeg[k];
        if (!seg) break;

        // ... volume transported from first segment is smaller of
        //     remaining flow volume & segment volume
        vseg = seg->v;
        vseg = MIN(vseg, v);

        // ... update total volume & mass entering downstream node
        *volin += vseg;
        *massin += vseg * seg->c;

        // ... reduce remaining flow volume by amount transported
        v -= vseg;

        // ... if all of segment's volume was transferred
        if (v >= 0.0 && vseg >= seg->v)
        {
            // ... replace this leading segment with the one behind it
            qual->FirstSeg[k] = seg->prev;
            if (qual->FirstSeg[k] == NULL) qual->LastSeg[k] = NULL;

            // ... recycle the used up segment
            seg->prev = qual->FreeSeg;
            qual->FreeSeg = seg;
            qual->MassBalance.segCount--;                                         
        }

        // ... otherwise just reduce this segment's volume
        else seg->v -= vseg;
    }
}


double findnodequal(Project *pr, int n, double volin,
                    double massin, double volout, long tstep)
{
    Network *net = &pr->network;
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;

    // Junction : mélange complet
    if (net->Node[n].Type == JUNCTION)
    {
        volin -= MIN(0.0, hyd->NodeDemand[n]) * tstep;
        if (volin > 0.0) qual->NodeQual[n] = massin / volin;
        else if (qual->Reactflag) qual->NodeQual[n] = noflowqual(pr, n);

        /* IMX : cross-junction — NodeQual[n] garde le mélange complet
           (bilan de masse), evalnodeoutflow appliquera IMX par lien */
        if (n <= net->Njuncs && qual->Crossjuncs[n].iscrossjunc == 1)
            return qual->NodeQual[n];
    }

    // Tank
    else if (net->Node[n].Type == TANK)
    {
        qual->NodeQual[n] = mixtank(pr, n, volin, massin, volout);
    }

    // Initialiser SourceQual avant toute source externe
    qual->SourceQual = 0.0;

    // Traçage
    if (qual->Qualflag == TRACE)
    {
        if (n == qual->TraceNode)
        {
            if (net->Node[n].Type == RESERVOIR) qual->SourceQual = 100.0;
            else qual->SourceQual = MAX(100.0 - qual->NodeQual[n], 0.0);
            qual->NodeQual[n] = 100.0;
        }
        return qual->NodeQual[n];
    }

    // Source chimique externe
    qual->SourceQual = findsourcequal(pr, n, volout, tstep);
    if (qual->SourceQual == 0.0) return qual->NodeQual[n];

    switch (net->Node[n].Type)
    {
        case JUNCTION:
            qual->NodeQual[n] += qual->SourceQual;
            return qual->NodeQual[n];
        case TANK:
            return qual->NodeQual[n] + qual->SourceQual;
        case RESERVOIR:
            qual->NodeQual[n] = qual->SourceQual;
            return qual->SourceQual;
    }
    return qual->NodeQual[n];
}

double  noflowqual(Project *pr, int n)
/*
**--------------------------------------------------------------
**   Input:   n = node index
**   Output:  quality for node n
**   Purpose: sets the quality for a junction node that has no
**            inflow to the average of the quality in its
**            adjoining link segments.
**   Note:    this function is only used for reactive substances.
**--------------------------------------------------------------
*/
{
    Network *net = &pr->network;
    Quality *qual = &pr->quality;

    int k, inflow, kount = 0;
    double c = 0.0;
    FlowDirection dir;
    Padjlist  alink;

    // Examine each link incident on the node
    for (alink = net->Adjlist[n]; alink != NULL; alink = alink->next)
    {
        // ... index of an incident link
        k = alink->link;
        dir = qual->FlowDir[k];

        // Node n is link's downstream node - add quality
        // of link's first segment to average
        if (net->Link[k].N2 == n && dir >= 0) inflow = TRUE;
        else if (net->Link[k].N1 == n && dir < 0)  inflow = TRUE;
        else inflow = FALSE;
        if (inflow == TRUE && qual->FirstSeg[k] != NULL)
        {
            c += qual->FirstSeg[k]->c;
            kount++;
        }

        // Node n is link's upstream node - add quality
        // of link's last segment to average
        else if (inflow == FALSE && qual->LastSeg[k] != NULL)
        {
            c += qual->LastSeg[k]->c;
            kount++;
        }
    }
    if (kount > 0) c = c / (double)kount;
    return c;
}


void evalnodeoutflow(Project *pr, int k, double c, long tstep)
// Renvoie la concentration calculée dans les tuyaux sortant (output)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**            c = quality from upstream node
**            tstep = time step
**   Output:  none
**   Purpose: releases flow volume and mass from the upstream
**            node of a link over a time step.
**--------------------------------------------------------------
*/
{
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;

    double v;
    Pseg seg;

    // Find flow volume (v) released over time step
    v = fabs(LINKFLOW(k)) * tstep;
    if (v == 0.0) return;

    // Release flow and mass into upstream end of the link

    // ... case where link has a last (most upstream) segment
    seg = qual->LastSeg[k];
    if (seg)
    {
        // ... if node quality close to segment quality then mix
        //     the nodal outflow volume with the segment's volume
        if (fabs(seg->c - c) < qual->Ctol)
        {
            seg->c = (seg->c*seg->v + c*v) / (seg->v + v);
            seg->v += v;
        }

        // ... otherwise add a new segment at upstream end of link
        else addseg(pr, k, v, c);
    }

    // ... link has no segments so add one
    else addseg(pr, k, v, c);
}


void updatemassbalance(Project *pr, int n, double massin,
                       double volout, long tstep)
/*
**--------------------------------------------------------------
**   Input:   n = node index
**            massin = mass inflow to node
**            volout = outflow volume from node
**   Output:  none
**   Purpose: Adds a node's external mass inflow and outflow
**            over the current time step to the network's
**            overall mass balance.
**--------------------------------------------------------------
*/
{
    Network *net = &pr->network;
    Hydraul *hyd = &pr->hydraul;
    Quality *qual = &pr->quality;

    double masslost = 0.0,
           massadded = 0.0;

    switch (net->Node[n].Type)
    {
        // Junctions lose mass from outflow demand & gain it from source inflow
    case JUNCTION:
        masslost = MAX(0.0, hyd->NodeDemand[n]) * tstep * qual->NodeQual[n];
        massadded = qual->SourceQual * volout;
        break;

        // Reservoirs add mass from quality source if specified or from a fixed
        // initial quality
    case RESERVOIR:
        masslost = massin;
        if (qual->SourceQual > 0.0) massadded = qual->SourceQual * volout;
        else                        massadded = qual->NodeQual[n] * volout;
        break;

        // Tanks add mass only from external source inflow
    case TANK:
        massadded = qual->SourceQual * volout;
        break;
    }
    qual->MassBalance.outflow += masslost;
    qual->MassBalance.inflow += massadded;
}


int sortnodes(Project *pr)
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  returns an error code
**   Purpose: topologically sorts nodes from upstream to downstream.
**   Note:    links with negligible flow are ignored since they can
**            create spurious cycles that cause the sort to fail.
**--------------------------------------------------------------
*/
{
    Network *net = &pr->network;
    Quality *qual = &pr->quality;

    int i, j, k, n;
    int *indegree = NULL;
    int *stack = NULL;
    int stacksize = 0;
    int numsorted = 0;
    int errcode = 0;
    FlowDirection dir;
    Padjlist  alink;

    // Allocate an array to count # links with inflow to each node
    // and for a stack to hold nodes waiting to be processed
    indegree = (int *)calloc(net->Nnodes + 1, sizeof(int));
    stack = (int *)calloc(net->Nnodes + 1, sizeof(int));
    if (indegree && stack)
    {
        // Count links with "non-negligible" inflow to each node
        for (k = 1; k <= net->Nlinks; k++)
        {
            dir = qual->FlowDir[k];
            if (dir == POSITIVE) n = net->Link[k].N2;
            else if (dir == NEGATIVE) n = net->Link[k].N1;
            else continue;
            indegree[n]++;
        }

        // Place nodes with no inflow onto a stack
        for (i = 1; i <= net->Nnodes; i++)
        {
            if (indegree[i] == 0)
            {
                stacksize++;
                stack[stacksize] = i;
            }
        }

        // Examine each node on the stack until none are left
        while (numsorted < net->Nnodes)
        {
            // ... if stack is empty then a cycle exists
            if (stacksize == 0)
            {
                //  ... add a non-sorted node connected to a sorted one to stack
                j = selectnonstacknode(pr, numsorted, indegree);
                if (j == 0) break;  // This shouldn't happen.
                indegree[j] = 0;
                stacksize++;
                stack[stacksize] = j;
            }

            // ... make the last node added to the stack the next
            //     in sorted order & remove it from the stack
            i = stack[stacksize];
            stacksize--;
            numsorted++;
            qual->SortedNodes[numsorted] = i;

            // ... for each outflow link from this node reduce the in-degree
            //     of its downstream node
            for (alink = net->Adjlist[i]; alink != NULL; alink = alink->next)
            {
                // ... k is the index of the next link incident on node i
                k = alink->link;

                // ... skip link if flow is negligible
                if (qual->FlowDir[k] == 0) continue;

                // ... link has flow out of node (downstream node n not equal to i)
                n = net->Link[k].N2;
                if (qual->FlowDir[k] < 0) n = net->Link[k].N1;

                // ... reduce degree of node n
                if (n != i && indegree[n] > 0)
                {
                    indegree[n]--;

                    // ... no more degree left so add node n to stack
                    if (indegree[n] == 0)
                    {
                        stacksize++;
                        stack[stacksize] = n;
                    }
                }
            }
        }
    }
    else errcode = 101;
    if (numsorted < net->Nnodes) errcode = 120;
    FREE(indegree);
    FREE(stack);
    return errcode;
}


int selectnonstacknode(Project *pr, int numsorted, int *indegree)
/*
**--------------------------------------------------------------
**   Input:   numsorted = number of nodes that have been sorted
**            indegree = number of inflow links to each node
**   Output:  returns a node index
**   Purpose: selects a next node for sorting when a cycle exists.
**--------------------------------------------------------------
*/
{
    Network  *net = &pr->network;
    Quality  *qual = &pr->quality;

    int i, m, n;
    Padjlist  alink;

    // Examine each sorted node in last in - first out order
    for (i = numsorted; i > 0; i--)
    {
        // For each link connected to the sorted node
        m = qual->SortedNodes[i];
        for (alink = net->Adjlist[m]; alink != NULL; alink = alink->next)
        {
            // ... n is the node of link k opposite to node m
            n = alink->node;

            // ... select node n if it still has inflow links
            if (indegree[n] > 0) return n;
        }
    }

    // If no node was selected by the above process then return the
    // first node that still has inflow links remaining
    for (i = 1; i <= net->Nnodes; i++)
    {
        if (indegree[i] > 0) return i;
    }

    // If all else fails return 0 indicating that no node was selected
    return 0;
}


void initsegs(Project *pr)
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: initializes water quality volume segments in each
**            pipe and tank.
**--------------------------------------------------------------
*/
{
    Network *net = &pr->network;
    Quality *qual = &pr->quality;

    int j, k;
    double c, v, v1;

    // Add one segment with assigned downstream node quality to each pipe
    for (k = 1; k <= net->Nlinks; k++)
    {
        qual->FirstSeg[k] = NULL;
        qual->LastSeg[k] = NULL;
        if (net->Link[k].Type == PIPE)
        {
            v = LINKVOL(k);
            j = net->Link[k].N2;
            c = qual->NodeQual[j];
            addseg(pr, k, v, c);
        }
    }

    // Initialize segments in tanks
    for (j = 1; j <= net->Ntanks; j++)
    {
        // Skip reservoirs
        if (net->Tank[j].A == 0.0) continue;

        // Establish initial tank quality & volume
        k = net->Tank[j].Node;
        c = net->Node[k].C0;
        v = net->Tank[j].V0;

        // Create one volume segment for entire tank
        k = net->Nlinks + j;
        qual->FirstSeg[k] = NULL;
        qual->LastSeg[k] = NULL;
        addseg(pr, k, v, c);

        // Create a 2nd segment for the 2-compartment tank model
        if (!qual->OutOfMemory && net->Tank[j].MixModel == MIX2)
        {
            // ... mixing zone segment
            v1 = MAX(0, v - net->Tank[j].V1frac * net->Tank[j].Vmax);
            qual->FirstSeg[k]->v = v1;

            // ... stagnant zone segment
            v = v - v1;
            addseg(pr, k, v, c);
        }
    }
}


void reversesegs(Project *pr, int k)
/*
**--------------------------------------------------------------
**   Input:   k = link index
**   Output:  none
**   Purpose: re-orients a link's segments when flow reverses.
**--------------------------------------------------------------
*/
{
    Quality *qual = &pr->quality;
    Pseg  seg, nseg, pseg;

    seg = qual->FirstSeg[k];
    qual->FirstSeg[k] = qual->LastSeg[k];
    qual->LastSeg[k] = seg;
    pseg = NULL;
    while (seg != NULL)
    {
        nseg = seg->prev;
        seg->prev = pseg;
        pseg = seg;
        seg = nseg;
    }
}


void addseg(Project *pr, int k, double v, double c)
/*
**-------------------------------------------------------------
**   Input:   k = segment chain index
**            v = segment volume
**            c = segment quality
**   Output:  none
**   Purpose: adds a segment to the start of a link
**            upstream of its current last segment.
**-------------------------------------------------------------
*/
{
    Quality *qual = &pr->quality;
    Pseg seg;

    // Grab the next free segment from the segment pool if available
    if (qual->FreeSeg != NULL)
    {
        seg = qual->FreeSeg;
        qual->FreeSeg = seg->prev;
    }

    // Otherwise allocate a new segment
    else
    {
        seg = (struct Sseg *) mempool_alloc(qual->SegPool, sizeof(struct Sseg));
        if (seg == NULL)
        {
            qual->OutOfMemory = TRUE;
            return;
        }
    }

    // Assign volume and quality to the segment
    seg->v = v;
    seg->c = c;

    // Add the new segment to the end of the segment chain
    seg->prev = NULL;
    if (qual->FirstSeg[k] == NULL) qual->FirstSeg[k] = seg;
    if (qual->LastSeg[k] != NULL)  qual->LastSeg[k]->prev = seg;
    qual->LastSeg[k] = seg;
    qual->MassBalance.segCount++;                                     
}


