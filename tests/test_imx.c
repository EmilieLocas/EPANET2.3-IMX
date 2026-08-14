#include <stdio.h>
#include "epanet2_2.h"

int main()
{
    EN_Project p;
    int errcode = 0;
    char errmsg[256];

    // Créer et ouvrir le projet
    errcode = EN_createproject(&p);
    if (errcode) { printf("Erreur création projet: %d\n", errcode); return 1; }

    errcode = EN_open(p, "imx_warning_test.inp", "test.rpt", "test.bin");
    if (errcode)
    {
        EN_geterror(errcode, errmsg, 256);
        printf("Erreur ouverture: %s\n", errmsg);
        return 1;
    }

    // Rouler la simulation hydraulique et qualité
    errcode = EN_solveH(p);
    if (errcode) { printf("Erreur hydraulique: %d\n", errcode); return 1; }

    errcode = EN_solveQ(p);
    if (errcode) { printf("Erreur qualité: %d\n", errcode); return 1; }

    // Lire les concentrations aux noeuds
    int nodeCount;
    EN_getcount(p, EN_NODECOUNT, &nodeCount);

    printf("\nConcentrations aux noeuds:\n");
    for (int i = 1; i <= nodeCount; i++)
    {
        char id[256];
        double qual;
        EN_getnodeid(p, i, id);
        EN_getnodevalue(p, i, EN_QUALITY, &qual);
        printf("  Noeud %s : %.4f mg/L\n", id, qual);
    }

    // Fermer et détruire le projet
    EN_close(p);
    EN_deleteproject(p);

    printf("\nSimulation terminée sans erreur.\n");
    return 0;
}