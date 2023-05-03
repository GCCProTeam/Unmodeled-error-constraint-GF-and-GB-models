/*------------------------------------------------------------------------------
*
*         Unmodeled error constraint GF model
*
*-----------------------------------------------------------------------------*/
#include <direct.h>
#include "NavYPPP.h"

/* Code and phase geometry-free measurement -------------------------------------------*/
extern double gfmeas(const obsd_t* obs, const nav_t* nav)
{
    double freq1;
    freq1 = sat2freq(obs->sat, obs->code[0], nav);
    if (freq1 == 0.0 || obs->L[0] == 0.0 || obs->P[0] == 0.0) return 0.0;
    return (obs->L[0] / freq1) * CLIGHT - obs->P[0];
}

/* detect cycle slip by geometry-free measurement or Unmodeled error constraint GF model-----------------------------*/
extern void detslp_gf(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
    double g0, g1, el, thres, GFResult=0.0;
    int i, j;
    rtk->opt.cs_gf = 0.8486;

    for (i = 0; i < n && i < MAXOBS; i++)
    {
        if ((g1 = gfmeas(obs + i, nav)) == 0.0)
        {
            continue;
        }
        g0 = rtk->ssat[obs[i].sat - 1].gf[0];
        if (g0 == 0.0) continue;
        rtk->ssat[obs[i].sat - 1].gf[0] = g1;

        el = rtk->ssat[obs[i].sat - 1].azel[1] * R2D;
        if (el < rtk->opt.elmin * R2D) continue;
        else  thres = rtk->opt.cs_gf;

        rtk->ssat[obs[i].sat - 1].delta_gf[0] = fabs(g1 - g0);
        rtk->ssat[obs[i].sat - 1].delta_gf[1] = thres;

        if (g0 != 0.0 && fabs(g1 - g0) >= thres)
        {
            GFResult = (g1 - rtk->ssat[obs[i].sat - 1].GFPrediction);
        }
        else
        {
            GFResult = g1 - g0;
        }
    }
    return;
}