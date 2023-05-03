
/*------------------------------------------------------------------------------
*
*         Unmodeled error constraint GB model
*
*-----------------------------------------------------------------------------*/
#include <direct.h>
#include "NavYPPP.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    50               /* max number of iterations */
#define MAX_STD_FIX 0.15             /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4               /* min satellite number for solution */
#define THRES_REJECT 4.0             /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS     SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.06               /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.007             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */
#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */


/* number and index of states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#define NC(opt)     (NSYS)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))

#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))

#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

static int GBFlag = 0;
static int GPSBDSNvNum[2];
static int REExclude;
static int CycleSlipRESat[MAXSAT];
#define K0YUZHI 1.5
#define K1YUZHI 3.0
#define REYUZHI 0.05
#define CycleREYuzhi         1.50
#define CycleEstimationYuzhi 0.75
static double GPSClock = 0.0;
static double BDSClock = 0.0;
static int CSSat = 0;
static int CSSatNum[10];
static int CSSatCount = 0;

/* PR */
static int StandPosRes(double* H, double* v, double* R, int nv, rtk_t* rtk)
{
    /* Tq  */
    int i, j, k, n, q, flag = 0, count = 0, STFlag = 0, m, iter;
    double Tq, Max = 0.0, Index = 0, u = 1.96, W;
    double* V1, * V2, * QL;
    /* w */
    double* w;
    int wi;
    w = mat(nv, 1);
    int exflag = 0, jishu = 0;
    double fact1 = 1.0, fact2, fact3 = 100000000, fact4;
    n = nv;
    q = nv - (5);
    V1 = mat(n, 1); V2 = mat(1, n); QL = mat(n, n);
    matcpy(V1, v, n, 1);
    for (i = 0; i < n; i++)
    {
        V1[i] = -V1[i];
    }
    matcpy(QL, R, n, n);
    matinv(QL, n);
    matmul("TN", 1, n, n, 1.0, V1, QL, 0.0, V2);
    matmul("NN", 1, 1, n, 1.0, V2, V1, 0.0, &Tq);
    Tq = Tq / q;
    /* Identification */
    if (Tq > chisqr[q])
    {
        for (iter = 0; ; iter++)
        {
        }
    }
    free(V1); free(V2); free(QL); free(w);
    return 1000;
}
/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t* rtk, const obsd_t* obs, int n, int stat)
{
    const prcopt_t* opt = &rtk->opt;
    int i, j;
    /* test # of valid satellites */
    rtk->sol.ns = 0;
    for (i = 0; i < n && i < MAXOBS; i++)
    {
        for (j = 0; j < opt->nf; j++)
        {

            if (!rtk->ssat[obs[i].sat - 1].vsat[j]) continue;

            rtk->ssat[obs[i].sat - 1].lock[j]++;
            rtk->ssat[obs[i].sat - 1].outc[j] = 0;
            if (j == 0) rtk->sol.ns++;
        }
    }
    rtk->sol.stat = rtk->sol.ns < MIN_NSAT_SOL ? SOLQ_NONE : stat;
    rtk->sol.PostSigma = PostSigma_2;
    rtk->sol.DeltrGPS = rtk->x[3];
    rtk->sol.DeltrBDS = rtk->x[4];

    if (rtk->sol.stat == SOLQ_FIX)
    {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->xa[i];
            rtk->sol.qr[i] = (float)rtk->Pa[i + i * rtk->na];
        }
        rtk->sol.qr[3] = (float)rtk->Pa[1];
        rtk->sol.qr[4] = (float)rtk->Pa[1 + 2 * rtk->na];
        rtk->sol.qr[5] = (float)rtk->Pa[2];
    }
    else
    {
        for (i = 0; i < 3; i++)
        {
            rtk->sol.rr[i] = rtk->x[i];
            rtk->sol.qr[i] = (float)rtk->P[i + i * rtk->nx];
        }
        rtk->sol.qr[3] = (float)rtk->P[1];
        rtk->sol.qr[4] = (float)rtk->P[2 + rtk->nx];
        rtk->sol.qr[5] = (float)rtk->P[2];
    }
    for (i = 0; i < MAXSAT; i++) for (j = 0; j < opt->nf; j++)
    {
        if (rtk->ssat[i].slip[j] & 3) rtk->ssat[i].slipc[j]++;
        if (rtk->ssat[i].fix[j] == 2 && stat != SOLQ_FIX) rtk->ssat[i].fix[j] = 1;
    }
}
/* least square estimation -----------------------------------------------------
* least square estimation by solving normal equation (x=(A*A')^-1*A*y)
* args   : double *A        I   transpose of (weighted) design matrix (n x m)
*          double *y        I   (weighted) measurements (m x 1)
*          int    n,m       I   number of parameters and measurements (n<=m)
*          double *x        O   estmated parameters (n x 1)
*          double *Q        O   esimated parameters covariance matrix (n x n)
* return : status (0:ok,0>:error)
* notes  : for weighted least square, replace A and y by A*w and w*y (w=W^(1/2))
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
static int PPPlsq_(const double* A, const double* y, int n, int m, double* x, double* Q)
{
    double* Ay;
    int info;

    if (m < n) return -1;
    Ay = mat(n, 1);
    matmul("NN", n, 1, m, 1.0, A, y, 0.0, Ay); /* Ay=A*y */
    matmul("NT", n, n, m, 1.0, A, A, 0.0, Q);  /* Q=A*A' */
    if (!(info = matinv(Q, n))) matmul("NN", n, 1, n, 1.0, Q, Ay, 0.0, x); /* x=Q^-1*Ay */
    free(Ay);
    return info;
}
static int MyPPPlsq(double* x, double* P, const double* H, const double* v, const double* R, int n, int m, double* QVV, int flag)
{
    double* x_, * xp_, * P_, * Pp_, * H_, * v_;
    int i, j, k, info, * ix;
    double sig;
    int ii, jj;
    double* B1, * NBB, * QLL, * QLL_, * LSQH, * QV, * QV_, * QVV_;

    //ix = imat(n, 1); for (i = k = 0; i < n; i++) if (x[i] != 0.0 && P[i + i * n] > 0.0) ix[k++] = i;
    ix = imat(n, 1); for (i = k = 0; i < n; i++) if (P[i + i * n] > 0.0) ix[k++] = i;
    x_ = mat(k, 1);
    xp_ = mat(k, 1);
    P_ = mat(k, k);
    Pp_ = mat(k, k);
    H_ = mat(k, m);
    v_ = mat(m, 1);

    for (i = 0; i < k; i++)
    {
        x_[i] = x[ix[i]];
        //for (j = 0; j < k; j++) P_[i + j * k] = P[ix[i] + ix[j] * n];
        for (j = 0; j < m; j++) H_[i + j * k] = H[ix[i] + j * n];
    }
    QLL = mat(m, m);
    matcpy(QLL, R, m, m);
    LSQH = mat(k, m);
    matcpy(LSQH, H_, k, m);
    QLL_ = mat(m, m);
    for (i = 0; i < m; i++)
    {
        v_[i] = v[i];
    }
    for (i = 0; i < m; i++)
    {
        sig = sqrt(R[i + i * m]);
        v_[i] /= sig;
        for (j = 0; j < k; j++) H_[j + i * k] /= sig;
    }
    info = PPPlsq_(H_, v_, k, m, xp_, Pp_);
    for (i = 0; i < k; i++)
    {
        x[ix[i]] = x[ix[i]] + xp_[i];
        for (j = 0; j < k; j++) P[ix[i] + ix[j] * n] = Pp_[i + j * k];
    }
    free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
    return info;
}
static int PPPlsq(double* x, double* P, const double* H, const double* v, const double* R, int n, int m, double* QVV, int flag)
{
    double* x_, * xp_, * P_, * Pp_, * H_, * v_;
    int i, j, k, info, * ix;
    double sig;
    int ii, jj;
    double* B1, * NBB, * QLL, * QLL_, * LSQH, * QV, * QV_, * QVV_;

    //ix = imat(n, 1); for (i = k = 0; i < n; i++) if (x[i] != 0.0 && P[i + i * n] > 0.0) ix[k++] = i;
    ix = imat(n, 1); for (i = k = 0; i < n; i++) if (P[i + i * n] > 0.0) ix[k++] = i;
    x_ = mat(k, 1);
    xp_ = mat(k, 1);
    P_ = mat(k, k);
    Pp_ = mat(k, k);
    H_ = mat(k, m);
    v_ = mat(m, 1);
    for (i = 0; i < k; i++)
    {
        x_[i] = x[ix[i]];
        for (j = 0; j < m; j++) H_[i + j * k] = H[ix[i] + j * n];
    }
    QLL = mat(m, m);
    matcpy(QLL, R, m, m);
    LSQH = mat(k, m);
    matcpy(LSQH, H_, k, m);
    QLL_ = mat(m, m);
    for (i = 0; i < m; i++)
    {
        v_[i] = v[i];
    }
    for (i = 0; i < m; i++)
    {
        sig = sqrt(R[i + i * m]);
        v_[i] /= sig;
        for (j = 0; j < k; j++) H_[j + i * k] /= sig;
    }
    info = PPPlsq_(H_, v_, k, m, xp_, Pp_);
    if (flag)
    {
        for (ii = 0; ii < m; ii++)
        {
            for (jj = 0; jj < m; jj++)
            {
                QLL_[jj + ii * m] = QLL[jj + ii * m] / (0.003 * 0.003);
            }
        }
        QV = mat(m, k);
        QV_ = mat(m, m);
        QVV_ = mat(m, m);
        B1 = mat(k, m);
        NBB = mat(k, k);
        matinv(QLL_, m);
        matmul("NN", k, m, m, 1.0, LSQH, QLL, 0.0, B1);
        matmul("NT", k, k, m, 1.0, B1, LSQH, 0.0, NBB);
        matinv(NBB, k);
        matmul("TN", m, k, k, 1.0, LSQH, NBB, 0.0, QV);
        matmul("NN", m, m, k, 1.0, QV, LSQH, 0.0, QV_);
        for (ii = 0; ii < m; ii++)
        {
            for (jj = 0; jj < m; jj++)
            {
                QVV_[jj + ii * m] = QLL[jj + ii * m] - QV_[jj + ii * m];
            }
        }
        matcpy(QVV, QVV_, m, m);
        free(B1); free(NBB); free(QLL); free(QLL_); free(LSQH); free(QV); free(QV_); free(QVV_);
    }
    for (i = 0; i < k; i++)
    {
        x[ix[i]] = x[ix[i]] + xp_[i];
        for (j = 0; j < k; j++) P[ix[i] + ix[j] * n] = Pp_[i + j * k];
    }
    free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
    return info;
}
/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t* rtk, int i)
{
    if (rtk->sol.stat == SOLQ_FIX) return SQRT(rtk->Pa[i + i * rtk->nx]);
    return SQRT(rtk->P[i + i * rtk->nx]);
}
/* write solution status for PPP ---------------------------------------------*/
extern int pppoutstat(rtk_t* rtk, char* buff)
{
    ssat_t* ssat;
    double tow, pos[3], vel[3], acc[3], * x;
    int i, j, week;
    char id[32], * p = buff;

    if (!rtk->sol.stat) return 0;
    tow = time2gpst(rtk->sol.time, &week);
    x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;
    /* receiver position */
    p += sprintf(p, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,rtk->sol.stat, x[0], x[1], x[2], STD(rtk, 0), STD(rtk, 1), STD(rtk, 2));
    /* receiver clocks */
    i = IC(0, &rtk->opt);
    p += sprintf(p, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",week, tow, rtk->sol.stat, 1, x[i] * 1E9 / CLIGHT, x[i + 1] * 1E9 / CLIGHT,STD(rtk, i) * 1E9 / CLIGHT, STD(rtk, i + 1) * 1E9 / CLIGHT);
    return (int)(p - buff);
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta) < 1E-12 && fabs(mu) < 1E-12) return PI;
    return atan2(-tan(beta), sin(mu)) + PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char* type, int opt, double beta, double mu,double* yaw)
{
    *yaw = yaw_nominal(beta, mu);
    return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char* type, int opt,const double* rs, double* exs, double* eys)
{
    double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
    double yaw, cosy, siny, erpv[5] = { 0 };
    int i;

    sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);
    /* beta and orbit angle */
    matcpy(ri, rs, 6, 1);
    ri[3] -= OMGE * ri[1];
    ri[4] += OMGE * ri[0];
    cross3(ri, ri + 3, n);
    cross3(rsun, n, p);
    if (!normv3(rs, es) || !normv3(rsun, esun) || !normv3(n, en) ||
        !normv3(p, ep)) return 0;
    beta = PI / 2.0 - acos(dot(esun, en, 3));
    E = acos(dot(es, ep, 3));
    mu = PI / 2.0 + (dot(es, esun, 3) <= 0 ? -E : E);
    if (mu < -PI / 2.0) mu += 2.0 * PI;
    else if (mu >= PI / 2.0) mu -= 2.0 * PI;

    /* yaw-angle of satellite */
    if (!yaw_angle(sat, type, opt, beta, mu, &yaw)) return 0;
    /* satellite fixed x,y-vector */
    cross3(en, es, ex);
    cosy = cos(yaw);
    siny = sin(yaw);
    for (i = 0; i < 3; i++) {
        exs[i] = -siny * en[i] + cosy * ex[i];
        eys[i] = -cosy * en[i] - siny * ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char* type, int opt,const double* rs, const double* rr, double* phw)
{
    double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
    double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
    int i;

    if (opt <= 0) return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!sat_yaw(time, sat, type, opt, rs, exs, eys)) return 0;
    /* unit vector satellite to receiver */
    for (i = 0; i < 3; i++) r[i] = rr[i] - rs[i];
    if (!normv3(r, ek)) return 0;
    /* unit vectors of receiver antenna */
    ecef2pos(rr, pos);
    xyz2enu(pos, E);
    exr[0] = E[1]; exr[1] = E[4]; exr[2] = E[7]; /* x = north */
    eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

    /* phase windup effect */
    cross3(ek, eys, eks);
    cross3(ek, eyr, ekr);
    for (i = 0; i < 3; i++) {
        ds[i] = exs[i] - ek[i] * dot(ek, exs, 3) - eks[i];
        dr[i] = exr[i] - ek[i] * dot(ek, exr, 3) + ekr[i];
    }
    cosp = dot(ds, dr, 3) / norm(ds, 3) / norm(dr, 3);
    if (cosp < -1.0) cosp = -1.0;
    else if (cosp > 1.0) cosp = 1.0;
    ph = acos(cosp) / 2.0 / PI;
    cross3(ds, dr, drs);
    if (dot(ek, drs, 3) < 0.0) ph = -ph;

    *phw = ph + floor(*phw - ph + 0.5); /* in cycle */
    return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int idx, int type, const prcopt_t* opt)
{
    double fact = 1.0, sinel = sin(el);

    if (type == 1) fact *= opt->eratio[idx == 0 ? 0 : 1];

    fact *= sys == SYS_GLO ? EFACT_GLO : (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);

    return SQR(fact * opt->err[1]) + SQR(fact * opt->err[2] / sinel);
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t* rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i] = xi;
    for (j = 0; j < rtk->nx; j++)
    {
        rtk->P[i + j * rtk->nx] = rtk->P[j + i * rtk->nx] = i == j ? var : 0.0;
    }
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t* obs, const nav_t* nav, const double* azel, const prcopt_t* opt, const double* dantr, const double* dants, double phw, double* L, double* P, double* Lc, double* Pc)
{
    double freq[NFREQ] = { 0 }, C1, C2;
    int i, sys = satsys(obs->sat, NULL);

    for (i = 0; i < 1; i++)
    {
        L[i] = P[i] = 0.0;
        freq[i] = sat2freq(obs->sat, obs->code[i], nav);
        if (freq[i] == 0.0 || obs->L[i] == 0.0 || obs->P[i] == 0.0) continue;
        if (testsnr(0, 0, azel[1], obs->SNR[i] * SNR_UNIT, &opt->snrmask)) continue;
        /* antenna phase center and phase windup correction */
        L[i] = obs->L[i] * CLIGHT / freq[i] - dants[i] - phw * CLIGHT / freq[i];
        P[i] = obs->P[i] - dants[i];
    }
}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t* rtk)
{
    double* x, * xp, pos[3];
    int i, j, * ix, nx;

    /* fixed mode */
    if (rtk->opt.mode == PMODE_PPP_FIXED)
    {
        for (i = 0; i < 3; i++) initx(rtk, rtk->opt.ru[i], 1E-8, i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x, 3) <= 0.0)
    {
        for (i = 0; i < 3; i++) initx(rtk, rtk->sol.rr[i], VAR_POS, i);
    }
    /* static ppp mode */
    if (rtk->opt.mode == PMODE_PPP_STATIC)
    {
        for (i = 0; i < 3; i++)
        {
            rtk->P[i * (1 + rtk->nx)] += SQR(rtk->opt.prn[5]) * fabs(rtk->tt);
        }
        return;
    }
    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics)
    {
        for (i = 0; i < 3; i++)
        {
            initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        }
        return;
    }
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t* rtk)
{
    int i;
    double DeltaGPSClock, DeltaBDSClock;
    DeltaGPSClock = 0.0;
    DeltaBDSClock = 0.0;
    initx(rtk, DeltaGPSClock, VAR_CLK, IC(0, &rtk->opt));
    initx(rtk, DeltaBDSClock, VAR_CLK, IC(1, &rtk->opt));
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav, int epoch)
{
    double L[NFREQ], P[NFREQ], Lc, Pc, bias[MAXOBS], offset = 0.0, pos[3] = { 0 };
    double freq1, freq2, ion, dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
    int i, j, k, f, sat, slip[MAXOBS] = { 0 }, clk_jump = 0;

    for (i = 0; i < MAXSAT; i++) for (j = 0; j < rtk->opt.nf; j++)
    {
        rtk->ssat[i].slip[j] = 0;
    }
    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(rtk, obs, n, nav);
    for (f = 0; f < NF(&rtk->opt); f++)
    {
        /* reset phase-bias if expire obs outage counter */
        for (i = 0; i < MAXSAT; i++)
        {
            if (++rtk->ssat[i].outc[f] > (uint32_t)rtk->opt.maxout || rtk->opt.modear == ARMODE_INST || clk_jump)
            {
                if (epoch > 2)
                {
                    rtk->ssat[i].SDReset = 1;
                }
            }
        }
        for (i = 0; i < n && i < MAXOBS; i++)
        {
            if (rtk->ssat[obs[i].sat - 1].slip[0] == 1)
            {
                sat = obs[i].sat;
                j = IB(sat, f, &rtk->opt);
                // reinitialize phase-bias if detecting cycle slip
                initx(rtk, 0, VAR_BIAS, IB(sat, f, &rtk->opt));
            }
        }
    }
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav, int epoch)
{
    /* temporal update of position */
    udpos_ppp(rtk);
    /* temporal update of clock */
    udclk_ppp(rtk);
    /* temporal update of phase-bias */
    udbias_ppp(rtk, obs, n, nav, epoch);
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double* rs, const double* rr, const pcv_t* pcv, double* dant)
{
    double ru[3], rz[3], eu[3], ez[3], nadir, cosa;
    int i;
    for (i = 0; i < 3; i++) {
        ru[i] = rr[i] - rs[i];
        rz[i] = -rs[i];
    }
    if (!normv3(ru, eu) || !normv3(rz, ez)) return;

    cosa = dot(eu, ez, 3);
    cosa = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
    nadir = acos(cosa);
    antmodel_s(pcv, nadir, dant);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double* pos, const double* azel, const prcopt_t* opt, const double* x, double* dtdx, const nav_t* nav, double* dtrp, double* var)
{
    double trp[3] = { 0 };
    if (opt->tropopt == TROPOPT_SAAS)
    {
        *dtrp = tropmodel(time, pos, azel, REL_HUMI);
        *var = SQR(ERR_SAAS);
        return 1;
    }
    return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double* pos, const double* azel, const prcopt_t* opt, int sat, const double* x, const nav_t* nav, double* dion, double* var)
{
    if (opt->ionoopt == IONOOPT_BRDC)
    {
        *dion = ionmodel(time, nav->ion_gps, pos, azel);
        *var = SQR(*dion * ERR_BRDCI);
        return 1;
    }
    return 0;
}
static int Myppp_res_3(int post, const obsd_t* obs, int n, const double* rs, const double* dts, const double* var_rs, const int* svh, const double* dr, int* exc, const nav_t* nav, const double* x, rtk_t* rtk, double* v, double* H, double* R, double* azel)
{
    prcopt_t* opt = &rtk->opt;
    double y, r, GPScdtr = 0.0, BDScdtr = 0.0, bias = 0.0, C = 0.0, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc, Pc;
    double var[MAXOBS * 2], dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, dcb, freq;
    double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
    double ve[MAXOBS * 2 * NFREQ] = { 0 }, vmax = 0;
    char str[32];
    int ne = 0, obsi[MAXOBS * 2 * NFREQ] = { 0 }, frqi[MAXOBS * 2 * NFREQ], maxobs, maxfrq, rej;
    int i, j, k, sat, sys, nv = 0, nx = rtk->nx, stat = 1;
    double Lastp = 0.0, Laste[3];
    double Fixedp = 0.0, Fixedpos[3];
    double Fixedrr[3], Fixede[3];
    double* Fixedazel;
    double Fixeddtrp = 0.0, Fixeddion = 0.0, Fixedvart = 0.0, Fixedvari = 0.0, FixedC = 0.0;
    double DionCorrected = 0.0;
    int week, ionsat1, ionsat2;
    double tow, TimeEpoch;
    tow = time2gpst(obs[0].time, &week);
    Fixedazel = zeros(2, n);

    time2str(obs[0].time, str, 2);
    for (i = 0; i < MAXSAT; i++) for (j = 0; j < opt->nf; j++) rtk->ssat[i].vsat[j] = 0;
    for (i = 0; i < 3; i++) rr[i] = x[i] + dr[i];
    ecef2pos(rr, pos);
    for (i = 0; i < n && i < MAXOBS; i++)
    {
        sat = obs[i].sat;
        if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 || satazel(pos, e, azel + i * 2) < opt->elmin)
        {
            exc[i] = 1; continue;
        }
        if (!(sys = satsys(sat, NULL)) || !rtk->ssat[sat - 1].vs || satexclude(obs[i].sat, var_rs[i], svh[i], opt) || exc[i])
        {
            exc[i] = 1; continue;
        }
        if (!model_trop(obs[i].time, pos, azel + i * 2, opt, x, dtdx, nav, &dtrp, &vart) || !model_iono(obs[i].time, pos, azel + i * 2, opt, sat, x, nav, &dion, &vari))
        {
            continue;
        }
        //satellite and receiver antenna model
        if (opt->posopt[0]) satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants);
        //phase windup model
        if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type, opt->posopt[2] ? 2 : 0, rs + i * 6, rr, &rtk->ssat[sat - 1].phw))
        {
            continue;
        }
        //corrected phase and code measurements
        corr_meas(obs + i, nav, azel + i * 2, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, L, P, NULL, NULL);
        for (j = 0; j < 1; j++)
        {
            if ((y = j % 2 == 0 ? L[j / 2] : P[j / 2]) == 0.0) continue;
            if ((freq = sat2freq(sat, obs[i].code[j / 2], nav)) == 0.0) continue;
            C = SQR(FREQ1 / freq) * ionmapf(pos, azel + i * 2) * (j % 2 == 0 ? -1.0 : 1.0);
            for (k = 0; k < nx; k++) H[k + nx * nv] = k < 3 ? -e[k] : 0.0;
            GPScdtr = x[3];
            BDScdtr = x[4];
            switch (sys)
            {
                case SYS_CMP:
                {                
                    H[IC(0, opt) + nx * nv] = 1.0;
                    H[IC(1, opt) + nx * nv] = 1.0;
                    break;
                }
                default:
                {            
                    H[IC(0, opt) + nx * nv] = 1.0;
                    break;
                }
            }
            //Unmodeled error constraint
            DionCorrected = rtk->ssat[sat - 1].GBPrediction;
            if (rtk->ssat[sat - 1].slip[0] == 1)
            {
                bias = x[IB(sat, 0, opt)];
                H[IB(sat, 0, opt) + nx * nv] = CLIGHT / freq;
            }
            if ((rtk->ssat[sat - 1].ph[0][0] == 0.0) || rtk->ssat[sat - 1].SDReset == 1)
            {
                if (GBFlag)
                {
                    Fixedp = geodist(rs + i * 6, Fixedrr, Fixede);
                    if (obs[i].L[0] == 0.0) continue;
                    ecef2pos(Fixedrr, Fixedpos);
                    satazel(Fixedpos, Fixede, Fixedazel + i * 2);
                    if (!model_trop(obs[i].time, Fixedpos, Fixedazel + i * 2, opt, x, dtdx, nav, &Fixeddtrp, &Fixedvart) || !model_iono(obs[i].time, Fixedpos, Fixedazel + i * 2, opt, sat, x, nav, &Fixeddion, &Fixedvari))
                    {
                        continue;
                    }
                    //satellite and receiver antenna model
                    if (opt->posopt[0]) satantpcv(rs + i * 6, Fixedrr, nav->pcvs + sat - 1, dants);
                    //phase windup model
                    if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type, opt->posopt[2] ? 2 : 0, rs + i * 6, Fixedrr, &rtk->ssat[sat - 1].phw))
                    {
                        continue;
                    }
                    //corrected phase and code measurements
                    corr_meas(obs + i, nav, Fixedazel + i * 2, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, L, P, NULL, NULL);
                    if ((y = j % 2 == 0 ? L[j / 2] : P[j / 2]) == 0.0) continue;
                    if ((freq = sat2freq(sat, obs[i].code[j / 2], nav)) == 0.0) continue;
                    FixedC = SQR(FREQ1 / freq) * ionmapf(Fixedpos, Fixedazel + i * 2) * (-1.0);
                    rtk->ssat[obs[i].sat - 1].pt[0][0] = obs[i].time;
                    rtk->ssat[obs[i].sat - 1].ph[0][0] = (y - (Fixedp - CLIGHT * dts[i * 2] + Fixeddtrp + FixedC * Fixeddion));
                    rtk->ssat[sat - 1].SDReset = 0;
                    rtk->ssat[sat - 1].outc[j] = 0;
                }
                exc[i] = 1;
                continue;
            }
            switch (sys)
            {
                case SYS_CMP:
                {
                    v[nv] = (y - (r + GPScdtr + BDScdtr - CLIGHT * dts[i * 2] + dtrp + C * dion + bias * (CLIGHT / freq) + DionCorrected)) - rtk->ssat[obs[i].sat - 1].ph[0][0];
                    break;
                }
                default:
                {
                    v[nv] = (y - (r + GPScdtr - CLIGHT * dts[i * 2] + dtrp + C * dion + bias * (CLIGHT / freq) + DionCorrected)) - rtk->ssat[obs[i].sat - 1].ph[0][0];
                    break;
                }
            }
            bias = 0.0;
            DionCorrected = 0.0;
            if (post)
            {
                if (j % 2 == 0) rtk->ssat[sat - 1].resc[j / 2] = v[nv];
            }
            /* variance */
            var[nv] = varerr(obs[i].sat, sys, azel[1 + i * 2], j / 2, j % 2, opt);
            var[nv] = var[nv] * 2;
            if (j % 2 == 0) rtk->ssat[sat - 1].vs = 1;
            if (j % 2 == 0) rtk->ssat[sat - 1].vsat[j / 2] = 1;
            nv++;
            if (GBFlag && post)
            {
                Fixedp = geodist(rs + i * 6, Fixedrr, Fixede);
                if (obs[i].L[0] == 0.0) continue;
                ecef2pos(Fixedrr, Fixedpos);
                satazel(Fixedpos, Fixede, Fixedazel + i * 2);
                if (!model_trop(obs[i].time, Fixedpos, Fixedazel + i * 2, opt, x, dtdx, nav, &Fixeddtrp, &Fixedvart) || !model_iono(obs[i].time, Fixedpos, Fixedazel + i * 2, opt, sat, x, nav, &Fixeddion, &Fixedvari))
                {
                    continue;
                }
                //satellite and receiver antenna model
                if (opt->posopt[0]) satantpcv(rs + i * 6, Fixedrr, nav->pcvs + sat - 1, dants);
                //phase windup model
                if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type, opt->posopt[2] ? 2 : 0, rs + i * 6, Fixedrr, &rtk->ssat[sat - 1].phw))
                {
                    continue;
                }
                //corrected phase and code measurements
                corr_meas(obs + i, nav, Fixedazel + i * 2, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, L, P, NULL, NULL);
                if ((y = j % 2 == 0 ? L[j / 2] : P[j / 2]) == 0.0)          continue;
                if ((freq = sat2freq(sat, obs[i].code[j / 2], nav)) == 0.0) continue;
                FixedC = SQR(FREQ1 / freq) * ionmapf(Fixedpos, Fixedazel + i * 2) * (-1.0);
                rtk->ssat[obs[i].sat - 1].pt[0][0] = obs[i].time;
                rtk->ssat[obs[i].sat - 1].ph[0][0] = (y - (Fixedp - CLIGHT * dts[i * 2] + Fixeddtrp + FixedC * Fixeddion));
            }
        }
    }
    for (i = 0; i < nv; i++) for (j = 0; j < nv; j++)
    {
        R[i + j * nv] = i == j ? var[i] : 0.0;
    }
    return nv;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t* opt)
{
    return NX(opt);
}
/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav, int epoch)
{
    const prcopt_t* opt = &rtk->opt;
    double* rs, * dts, * var, * v, * H, * R, * azel, * xp, * Pp, dr[3] = { 0 }, std[3];
    char str[32];
    int i, j, nv, info, svh[MAXOBS], exc[MAXOBS] = { 0 }, stat = SOLQ_SINGLE;

    time2str(obs[0].time, str, 2);
    rs = mat(6, n); dts = mat(2, n); var = mat(1, n); azel = zeros(2, n);
    /* temporal update of ekf states */
    udstate_ppp(rtk, obs, n, nav, epoch);
    /* satellite positions and clocks */
    satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);
    /* earth tides correction */
    if (opt->tidecorr)
    {
        tidedisp(gpst2utc(obs[0].time), rtk->x, opt->tidecorr == 1 ? 1 : 7, &nav->erp, opt->odisp[0], dr);
    }
    nv = n * rtk->opt.nf * 2 + MAXSAT + 3;
    xp = mat(rtk->nx, 1);
    Pp = zeros(rtk->nx, rtk->nx);
    v = mat(nv, 1);
    H = mat(rtk->nx, nv);
    R = mat(nv, nv);

    double* Qvv = mat(nv, nv);
    double* LastX, * LastPp, DeltaX[3] = { 1000, 1000, 1000 };
    LastX = mat(rtk->nx, 1);
    LastPp = zeros(rtk->nx, rtk->nx);
    GBFlag = 1;
    for (i = 0; i < 1; i++)
    {
        if (i == 0)
        {
            matcpy(xp, rtk->x, rtk->nx, 1);
            matcpy(Pp, rtk->P, rtk->nx, rtk->nx);
        }
        //prefit residuals
        if (!(nv = Myppp_res_3(0, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel)))
        {
            break;
        }
        if ((info = MyPPPlsq(xp, Pp, H, v, R, rtk->nx, nv, Qvv, 1)))
        {
            break;
        }
        if (i > 0)
        {
            for (j = 0; j < 3; j++)
            {
                DeltaX[j] = xp[j] - LastX[j];
            }
        }
        if (norm(DeltaX, 3) < REYUZHI)
        {
            break;
        }
        matcpy(LastX, xp, rtk->nx, 1);
        matcpy(LastPp, Pp, rtk->nx, rtk->nx);
    }
    if (i > 0)
    {
        //postfit residuals
        if (!(nv = Myppp_res_3(1, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel)))
        {
            trace(1, "%s ppp (%d) no valid obs data\n", str, i + 1);
        }
        CSSat = StandPosRes(H, v, R, nv, rtk);
        for (ii = 0; ii < 10; ii++)
        {
            CSSatNum[ii] = 0;
        }
        CSSatCount = 0;
        matcpy(rtk->x, xp, rtk->nx, 1);
        matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
        stat = SOLQ_PPP;
        /* update solution status */
        update_stat(rtk, obs, n, stat);
    }
    free(rs); free(dts); free(var); free(azel);
    free(xp); free(Pp); free(v); free(H); free(R);
    return;
}