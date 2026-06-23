#include "Solvers.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

#define gg1 (ga - 1.0)
#define g2 (ga + 1.0)
#define gg2 (ga + 1.0)
#define gp ((g2/ga)/2.0)
#define gm ((g1/ga)/2.0)
#define gga ga
#define ER_S std::cout << "\n---------------------\nStandart error in file: Solvers.cpp\n" << endl
#define watch(x) cout << (#x) << " is " << (x) << endl
#define M(x) cout << (#x)  << endl

using namespace std;


Solvers::Solvers()
{
    //ctor
}

Solvers::~Solvers()
{
    //dtor
}

void Solvers::perpendicular(const vector<double>& a, vector<double>& b, vector<double>& c, bool t)
{
    if (t == false)
    {
       double A = a[0]*a[0] + a[1]*a[1];
       if( A > 0.01*(A + a[2]*a[2]))
       {
           double B = sqrt(A);
           b[0] = -a[1]/B;
           b[1] = a[0]/B;
           b[2] = 0.0;
           double C = sqrt(A* (A + a[2]*a[2]));
           c[0] = -a[0]*a[2]/C;
           c[1] = -a[1]*a[2]/C;
           c[2] = A/C;
           return;
       }
       A = a[0]*a[0] + a[2]*a[2];
       if( A > 0.01*(A + a[1]*a[1]))
       {
           double B = sqrt(A);
           b[0] = -a[2]/B;
           b[1] = 0.0;
           b[2] = a[0]/B;
           double C = sqrt( A* (A + a[1]*a[1]) );
           c[0] = a[0]*a[1]/C;
           c[1] = -A/C;
           c[2] = a[1]*a[2]/C;
           return;
       }
    }
    else
    {
       double A = a[0]*a[0] + a[1]*a[1];
       if( A > 0.01)
       {
           double B = sqrt(A);
           b[0] = -a[1]/B;
           b[1] = a[0]/B;
           b[2] = 0.0;;
           c[0] = -a[0]*a[2]/B;
           c[1] = -a[1]*a[2]/B;
           c[2] = A/B;
           return;
       }
       A = a[0]*a[0] + a[2]*a[2];
       if( A > 0.01)
       {
           double B = sqrt(A);
           b[0] = -a[2]/B;
           b[1] = 0.0;
           b[2] = a[0]/B;

           c[0] = a[0]*a[1]/B;
           c[1] = -A/B;
           c[2] = a[1]*a[2]/B;
           return;
       }
    }

}

double Solvers::HLLD_Solver_Alexashov(const std::vector<double>& qqq1, const std::vector<double>& qqq2,//
                                  const std::vector<double>& n, std::vector<double>& qqq,//
                double& dsl, double& dsp, double& dsc, double w, const double& DIST)
{
    int id_bn = 1;
    int n_state = 1;
    double FR[8],FL[8];
    double FW[8],UL[8],UZ[8],UR[8];
    double UZL[8],UZR[8];
    double UZZL[8],UZZR[8];

    double vL[3],vR[3],bL[3],bR[3];
    double vzL[3],vzR[3],bzL[3],bzR[3];
    double vzzL[3],vzzR[3],bzzL[3],bzzR[3];
    double qv[3],qb[3];
    double aco[3][3];

    double wv = w;
    double r1 = qqq1[0];
    double u1 = qqq1[2];
    double v1 = qqq1[3];
    double w1 = qqq1[4];
    double p1 = qqq1[1];
    double bx1 = qqq1[5]/spi4;
    double by1 = qqq1[6]/spi4;
    double bz1 = qqq1[7]/spi4;


    double r2 = qqq2[0];
    double u2 = qqq2[2];
    double v2 = qqq2[3];
    double w2 = qqq2[4];
    double p2 = qqq2[1];
    double bx2 = qqq2[5]/spi4;
    double by2 = qqq2[6]/spi4;
    double bz2 = qqq2[7]/spi4;

    double ro = (r2+r1)/2.0;
    double au = (u2+u1)/2.0;
    double av = (v2+v1)/2.0;
    double aw = (w2+w1)/2.0;
    double ap = (p2+p1)/2.0;
    double abx = (bx2+bx1)/2.0;
    double aby = (by2+by1)/2.0;
    double abz = (bz2+bz1)/2.0;


    double bk = abx*n[0] + aby*n[1]  + abz*n[2];
    double b2 = kv(abx) + kv(aby) + kv(abz) ;

    double d = b2 - kv(bk);
    aco[0][0] = n[0];
    aco[1][0] = n[1];
    aco[2][0] = n[2];
    if (d > eps)
    {
        d = sqrt(d);
        aco[0][1] = (abx - bk * n[0] )/d;
        aco[1][1] = (aby - bk * n[1] )/d;
        aco[2][1] = (abz - bk * n[2] )/d;
        aco[0][2] = (aby * n[2] - abz * n[1] )/d;
        aco[1][2] = (abz * n[0] - abx * n[2] )/d;
        aco[2][2] = (abx * n[1] - aby * n[0] )/d;
    }
    else
    {
        double aix, aiy, aiz;
        if ( (fabs( n[0] ) < fabs( n[1] )) && ( fabs( n[0] ) < fabs( n[2] ) ) )
        {
            aix =  1.0;
            aiy =  0.0;
            aiz =  0.0;
        }
        else if( fabs( n[1] ) < fabs( n[2] ))
        {
        aix =  0.0;
        aiy =  1.0;
        aiz =  0.0;
        }
        else
        {
        aix =  0.0;
        aiy =  0.0;
        aiz =  1.0;
        }

        double aik = aix* n[0] + aiy* n[1] + aiz* n[2];
        d = sqrt( 1.0 - kv(aik) );
        aco[0][1] = (aix - aik * n[0] )/d;
        aco[1][1] = (aiy - aik * n[1] )/d;
        aco[2][1] = (aiz - aik * n[2] )/d;
        aco[0][2] = (aiy * n[2] - aiz * n[1] )/d;
        aco[1][2] = (aiz * n[0] - aix * n[2] )/d;
        aco[2][2] = (aix * n[1] - aiy * n[0] )/d;
    }

    for (int i = 0; i < 3; i++)
    {
        vL[i] = aco[0][i] * u1 + aco[1][i] * v1 + aco[2][i] * w1;
        vR[i] = aco[0][i] * u2 + aco[1][i] * v2 + aco[2][i] * w2;
        bL[i] = aco[0][i] * bx1 + aco[1][i] * by1 + aco[2][i] * bz1;
        bR[i] = aco[0][i] * bx2 + aco[1][i] * by2 + aco[2][i] * bz2;
    }

    double aaL =  bL[0]/sqrt(r1);
    double b2L =  kv(bL[0]) + kv(bL[1]) + kv(bL[2]);
    double b21 = b2L/r1;
    double cL  = sqrt(ga*p1/r1);
    double qp = sqrt(b21 + cL*(cL + 2.0*aaL));
    double qm = sqrt(b21 + cL*(cL - 2.0*aaL));
    double cfL = (qp + qm)/2.0;
    double ptL = p1 + b2L/2.0;

    double aaR =  bR[0]/sqrt(r2);
    double b2R =  kv(bR[0]) + kv(bR[1]) + kv(bR[2]);
    double b22 = b2R/r2;
    double cR  = sqrt(ga*p2/r2);
    qp = sqrt(b22 + cR*(cR + 2.0*aaR));
    qm = sqrt(b22 + cR*(cR - 2.0*aaR));
    double cfR = (qp + qm)/2.0;
    double ptR = p2 + b2R/2.0;

    double aC  =  (aaL + aaR)/2.0;
    double b2o = (b22 + b21)/2.0;
    double cC  = sqrt(ga*ap/ro);
    qp = sqrt(b2o + cC*(cC + 2.0*aC));
    qm = sqrt(b2o + cC*(cC - 2.0*aC));
    double cfC = (qp + qm)/2.0;
    double vC1 = (vL[0] + vR[0])/2.0;

    double SL = min( (vL[0] - cfL),(vR[0] - cfR) );
    double SR = max( (vL[0] + cfL),(vR[0] + cfR) );

    double suR = SR - vR[0];
    double suL = SL - vL[0];
    double SM = (suR*r2*vR[0] - ptR+ptL - suL*r1*vL[0])/(suR*r2 - suL*r1);

    if(SR <= SL)
    {
        ER_S;
        cout << 231 << endl;
        exit(-1);
    }

    double SM00 = SM;
    double SR00 = SR;
    double SL00 = SL;
    double SM01, SR01, SL01;
    if ( (SM00 >=  SR00)||(SM00 <=  SL00) )
    {
        SL = min( (vL[0] - cfL),(vR[0] - cfR) );
        SR = max( (vL[0] + cfL),(vR[0] + cfR) );
        suR = SR - vR[0];
        suL = SL - vL[0];
        SM = (suR*r2*vR[0] - ptR + ptL - suL*r1*vL[0])/(suR*r2 - suL*r1);
        SM01 = SM;
        SR01 = SR;
        SL01 = SL;
        if ( (SM01 >=  SR01)||(SM01 <=  SL01) )
        {
            ER_S;
            cout << 251 << endl;
            exit(-1);
        }
    }
    dsl = SL;
    dsc = SM;
    dsp = SR;
     
    double UU = max( fabs(dsl), fabs(dsp) );
    double time = kurant * DIST/UU;

    double upt1 = ( kv(u1) + kv(v1) + kv(w1) )/2.0;
    double sbv1 = u1*bx1 + v1*by1 + w1*bz1;

    double upt2 = ( kv(u2) + kv(v2) + kv(w2) )/2.0;
    double sbv2 = u2*bx2 + v2*by2 + w2*bz2;

    double e1 = p1/g1 + r1*upt1 + b2L/2.0;
    double e2 = p2/g1 + r2*upt2 + b2R/2.0;

    FL[0] = r1*vL[0];
    FL[1] = r1*vL[0]*vL[0] + ptL - kv(bL[0]);
    FL[2] = r1*vL[0]*vL[1]     - bL[0]*bL[1];
    FL[3] = r1*vL[0]*vL[2]     - bL[0]*bL[2];
    FL[4] = (e1 + ptL)*vL[0]     - bL[0]*sbv1;
    FL[5] = 0.0;
    FL[6] = vL[0]*bL[1] - vL[1]*bL[0];
    FL[7] = vL[0]*bL[2] - vL[2]*bL[0];

    FR[0] = r2*vR[0];
    FR[1] = r2*vR[0]*vR[0] + ptR - kv(bR[0]);
    FR[2] = r2*vR[0]*vR[1] - bR[0]*bR[1];
    FR[3] = r2*vR[0]*vR[2] - bR[0]*bR[2];
    FR[4] = (e2 + ptR)*vR[0] - bR[0]*sbv2;
    FR[5] = 0.0;
    FR[6] = vR[0]*bR[1] - vR[1]*bR[0];
    FR[7] = vR[0]*bR[2] - vR[2]*bR[0];

    UL[0] = r1;
    UL[4] = e1;
    UR[0] = r2;
    UR[4] = e2;


    for (int ik = 0; ik < 3; ik++)
    {
        UL[ik+1] = r1*vL[ik];
        UL[ik+5] =   bL[ik];
        UR[ik+1] = r2*vR[ik];
        UR[ik+5] =   bR[ik];
    }

    for (int ik = 0; ik < 8; ik++)
    {
        UZ[ik] = (SR*UR[ik]-SL*UL[ik]+FL[ik]-FR[ik])/(SR-SL);
    }


    // if(id_bn == 1)
    // {
    //     UZ[5] = 0.0;
    // }
    if(n_state == 1)  /// HLL
    {
        double dq[8];
        for(int ik = 0; ik < 8; ik++)
        {
            dq[ik] = UR[ik] - UL[ik];
        }

        double TL = SL;
        double TR = SR;
        if(SL > wv)
        {
            TL = 0.0;
            for(int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv*UL[ik];
            }
        }
        if ( (SL <= wv) && (wv <= SR) )
        {
            for(int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv*UZ[ik];
            }
        }
        if(SR < wv)
        {
            TR = 0.0;
            for(int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv*UR[ik];
            }
        }


        double a = TR*TL;
        double b = TR - TL;

        qqq[0] = (TR*FL[0] - TL*FR[0] + a*dq[0])/b - FW[0];
        qqq[4] = (TR*FL[4] - TL*FR[4] + a*dq[4])/b - FW[4];
        
        for(int ik = 1; ik < 4; ik++)
        {
            qv[ik - 1] = (TR*FL[ik] - TL*FR[ik] + a*dq[ik])/b - FW[ik];
        }
        for(int ik = 5; ik < 8; ik++)
        {
            qb[ik - 5] = (TR*FL[ik] - TL*FR[ik] + a*dq[ik])/b - FW[ik];
        }

        for(int i = 0; i < 3; i++)
        {
            qqq[i + 1] = aco[i][0]*qv[0] + aco[i][1]*qv[1] + aco[i][2]*qv[2];
            qqq[i + 5] = aco[i][0]*qb[0] + aco[i][1]*qb[1] + aco[i][2]*qb[2];
            qqq[i + 5] = spi4*qqq[i + 5];
        }

        double SWAP = qqq[4];
        qqq[4] = qqq[5];
        qqq[5] = qqq[6];
        qqq[6] = qqq[7];
        qqq[7] = SWAP;

        return time;
    }
    if(n_state == 3)
    {
        double ptz = (suR*r2*ptL - suL*r1*ptR + r1*r2*suR*suL*(vR[0] - vL[0]))/(suR*r2 - suL*r1);

        vzL[0] = SM;
        vzR[0] = SM;
        vzzL[0] = SM;
        vzzR[0] = SM;
        double ptzL = ptz;
        double ptzR = ptz;
        double ptzzL = ptz;
        double ptzzR = ptz;

        double suRm = suR/(SR - SM);
        double suLm = suL/(SL - SM);
        double rzR = r2*suRm;
        double rzL = r1*suLm;

        double bn = UZ[5];
        double bn2 = bn*bn;
        bzL[0] = bn;
        bzR[0] = bn;
        bzzL[0] = bn;
        bzzR[0] = bn;
        
        double ttR = r2*suR*(SR - SM) - bn2;
        double tvR, tbR;
        
        if( fabs(ttR) <= eps)
        {
            tvR = 0.0;
            tbR = 0.0;
        }
        else
        {
            tvR = (SM - vR[0])/ttR;
            tbR = (r2*suR*suR - bn2)/ttR;
        }

        double ttL = r1*suL*(SL - SM) - bn2;
        double tvL, tbL;

        if(fabs(ttL) <= eps)
        {
            tvL = 0.0;
            tbL = 0.0;
        }
        else
        {
            tvL = (SM - vL[0])/ttL;
            tbL = (r1*suL*suL - bn2)/ttL;
        }

        vzL[1] = vL[1] - bn*bL[1]*tvL;
        vzL[2] = vL[2] - bn*bL[2]*tvL;
        vzR[1] = vR[1] - bn*bR[1]*tvR;
        vzR[2] = vR[2] - bn*bR[2]*tvR;

        bzL[1] = bL[1]*tbL;
        bzL[2] = bL[2]*tbL;
        bzR[1] = bR[1]*tbR;
        bzR[2] = bR[2]*tbR;

        double sbvL = bzL[0]*vzL[0] + bzL[1]*vzL[1] + bzL[2]*vzL[2];
        double sbvR = bzR[0]*vzR[0] + bzR[1]*vzR[1] + bzR[2]*vzR[2];

        double ezR = e2*suRm + (ptz*SM - ptR*vR[0] + bn*(sbv2 - sbvR))/(SR - SM);
        double ezL = e1*suLm + (ptz*SM - ptL*vL[0] + bn*(sbv1 - sbvL))/(SL - SM);

        double rzzR = rzR;
        double rzzL = rzL;
        double rzRs = sqrt(rzR);
        double rzLs = sqrt(rzL);
        double rzss = rzRs + rzLs;
        double rzps = rzRs*rzLs;

        double SZL = SM - fabs(bn)/rzLs;
        double SZR = SM + fabs(bn)/rzRs;

        int ibn = 0;
        double sbn;
        if(fabs(bn) > epsb)
        {
            sbn = fabs(bn)/bn;
            ibn = 1;
        }
        else
        {
            sbn = 0.0;
            ibn = 0;
            SZL = SM;
            SZR = SM;
        }

        vzzL[1] = (rzLs*vzL[1] + rzRs*vzR[1] + sbn*(bzR[1] - bzL[1]) )/rzss;
        vzzL[2] = (rzLs*vzL[2] + rzRs*vzR[2] + sbn*(bzR[2] - bzL[2]) )/rzss;
        vzzR[1] = vzzL[1];
        vzzR[2] = vzzL[2];

        bzzL[1] = (rzLs*bzR[1] + rzRs*bzL[1] + sbn*rzps*(vzR[1] - vzL[1]) )/rzss;
        bzzL[2] = (rzLs*bzR[2] + rzRs*bzL[2] + sbn*rzps*(vzR[2] - vzL[2]) )/rzss;
        bzzR[1] = bzzL[1];
        bzzR[2] = bzzL[2];

        double sbzz = bzzL[0]*vzzL[0] + bzzL[1]*vzzL[1] + bzzL[2]*vzzL[2];

        double ezzR = ezR + rzRs*sbn*(sbvR - sbzz);
        double ezzL = ezL - rzLs*sbn*(sbvL - sbzz);

        UZL[0] = rzL;
        UZL[4] = ezL;
        UZR[0] = rzR;
        UZR[4] = ezR;
        
        for (int ik = 0; ik < 3; ik++)
        {
            UZL[ik + 1] = vzL[ik]*rzL;
            UZL[ik + 5] = bzL[ik];
            UZR[ik + 1] = vzR[ik]*rzR;
            UZR[ik + 5] = bzR[ik];
        }

        UZZL[0] = rzzL;
        UZZL[4] = ezzL;
        UZZR[0] = rzzR;
        UZZR[4] = ezzR;

        for (int ik = 0; ik < 3; ik++)
        {
            UZZL[ik + 1] = vzzL[ik]*rzzL;
            UZZL[ik + 5] = bzzL[ik];
            UZZR[ik + 1] = vzzR[ik]*rzzR;
            UZZR[ik + 5] = bzzR[ik];
        }

        int j_ccs =  - 1;

        if(SL > wv)
        {
            qqq[0] = FL[0] - wv*UL[0];
            qqq[4] = FL[4] - wv*UL[4];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] - wv*UL[ik];
            }
            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] - wv*UL[ik];
            }
            j_ccs =  1;
        }

        if ( (SL <= wv)&&(SZL >= wv) )
        {
            int ik = 0;     
            qqq[ik] = FL[ik] + SL*(UZL[ik] - UL[ik]) - wv*UZL[ik];
            ik = 4;     
            qqq[ik] = FL[ik] + SL*(UZL[ik] - UL[ik]) - wv*UZL[ik];
            for (ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] + SL*(UZL[ik] - UL[ik]) - wv*UZL[ik];
            }
            for (ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] + SL*(UZL[ik] - UL[ik]) - wv*UZL[ik];
            }
            j_ccs =  2;
        }
        
       if(ibn == 1)
       {
           if ( (SZL <= wv)&&(SM >= wv) )
           {
                int ik = 0;
                qqq[ik] = FL[ik] + SZL*(UZZL[ik] - UZL[ik]) +  SL*( UZL[ik] -  UL[ik]) - wv*UZZL[ik];
                ik = 4;
                qqq[ik] = FL[ik] + SZL*(UZZL[ik] - UZL[ik]) +  SL*( UZL[ik] -  UL[ik]) - wv*UZZL[ik];

                for (ik = 1; ik < 4; ik++)
                {
                    qv[ik - 1] = FL[ik] + SZL*(UZZL[ik] - UZL[ik]) +  SL*( UZL[ik] -  UL[ik]) - wv*UZZL[ik];
                }
                for (ik = 5; ik < 8; ik++)
                {
                    qb[ik - 5] = FL[ik] + SZL*(UZZL[ik] - UZL[ik]) +  SL*( UZL[ik] -  UL[ik]) - wv*UZZL[ik];
                }
                j_ccs =  3;
           }

           if ((SM <= wv)&&(SZR >= wv))
           {
                int ik = 1;
                qqq[ik] = FR[ik] + SZR*(UZZR[ik] - UZR[ik])  +  SR*( UZR[ik] -  UR[ik]) - wv*UZZR[ik];
                ik = 5;
                qqq[ik] = FR[ik] + SZR*(UZZR[ik] - UZR[ik]) +  SR*( UZR[ik] -  UR[ik]) - wv*UZZR[ik];

                for (ik = 1; ik < 4; ik++)
                {
                    qv[ik - 1] = FR[ik] + SZR*(UZZR[ik] - UZR[ik]) +  SR*( UZR[ik] -  UR[ik]) - wv*UZZR[ik];
                }
                for (ik = 5; ik < 8; ik++)
                {
                    qb[ik - 5] = FR[ik] + SZR*(UZZR[ik] - UZR[ik]) +  SR*( UZR[ik] -  UR[ik]) - wv*UZZR[ik];
                }
                j_ccs =  4;
           }
       }
       
        if ( (SZR <= wv)&&(SR >= wv) )
        {
            int ik = 1;
            qqq[ik] = FR[ik] + SR*(UZR[ik] - UR[ik]) - wv*UZR[ik];
            ik = 5;
            qqq[ik] = FR[ik] + SR*(UZR[ik] - UR[ik]) - wv*UZR[ik];
            for (ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] + SR*(UZR[ik] - UR[ik]) - wv*UZR[ik];
            }
            for (ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] + SR*(UZR[ik] - UR[ik]) - wv*UZR[ik];
            }
            j_ccs =  5;
        }

        if(SR < wv)
        {
            qqq[0] = FR[0] - wv*UR[0];
            qqq[4] = FR[4] - wv*UR[4];
            for (int ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] - wv*UR[ik];
            }
            for (int ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] - wv*UR[ik];
            }
            j_ccs =  6;
        }

        if(j_ccs ==  - 1)
        {
            ER_S;
            cout << 559 << endl;
            watch(SL);
            watch(r1);
            watch(p1);
            watch(v1);
            watch(u1);
            watch(w1);
            watch(bx1);
            watch(by1);
            watch(bz1);
            watch(r2);
            watch(p2);
            watch(v2);
            watch(u2);
            watch(w2);
            watch(bx2);
            watch(by2);
            watch(bz2);
            exit(-1);
        }


        double SN  =  - max(fabs(SL),fabs(SR))/2.0;

        double wbn = 0.0;
        if(wv >= SR)
        {
            wbn = wv*bR[0];
        }
        else if (wv <= SL)
        {
            wbn = wv*bL[0];
        }
        else
        {
            wbn = wv*(bL[0] + bR[0])/2.0;
        }

        qb[0] =  SN*(bR[0] - bL[0]) - wbn;



        for (int i = 0; i < 3; i++)
        {
            qqq[i + 1] = aco[i][0]*qv[0] + aco[i][1]*qv[1] + aco[i][2]*qv[2];
            qqq[i + 5] = aco[i][0]*qb[0] + aco[i][1]*qb[1] + aco[i][2]*qb[2];
            qqq[i + 5] = spi4*qqq[i + 5];
        }

        double SWAP = qqq[4];
        qqq[4] = qqq[5];
        qqq[5] = qqq[6];
        qqq[6] = qqq[7];
        qqq[7] = SWAP;

        // watch(qqq[0]);
        // watch(qqq[1]);
        // watch(qqq[2]);
        // watch(qqq[3]);
        // watch(qqq[4]);
        // watch(qqq[5]);
        // watch(qqq[6]);
        // watch(qqq[7]);
        // watch(qqq[8]);

        return time;
    }
    return time;
}

double Solvers::Godunov_Solver_Alexashov(const std::vector<double>& qqq1, const std::vector<double>& qqq2,//
                                  const std::vector<double>& n, std::vector<double>& qqq,//
                double& dsl, double& dsp, double& dsc, double w, const double& DIST)
{
    double al = n[0];
    double be = n[1];
    double ge = n[2];
    double time = 0.0;

    vector<double> tau_1(3);
    vector<double> tau_2(3);

    //perpendicular(n, tau_1, tau_2, true);  // Находим базис,  Считаем нормаль единичной

    /*double al2 = tau_1[0];
    double be2 = tau_1[1];
    double ge2 = tau_1[2];
    double al3 = tau_2[0];
    double be3 = tau_2[1];
    double ge3 = tau_2[2];*/

    double al2 = tau_1[0] = be;
    double be2 = tau_1[1]  = -al;
    double ge2 = tau_1[2] =  0.0;
    double al3 = tau_2[0] = 0.0;
    double be3 = tau_2[1] = 0.0;
    double ge3 = tau_2[2] = 1.0;

    double enI  = al * qqq1[2] + be * qqq1[3] + ge * qqq1[4];
    // if (enI > 110500)
    // {
    //     ER_S;
    //     M("Problem enI");
    //     watch(al);
    //     watch(be);
    //     watch(ge);
    //     watch(qqq1[0]);
    //     watch(qqq1[1]);
    //     watch(qqq1[2]);
    //     watch(qqq1[3]);
    //     watch(qqq1[4]);
    //     watch(n[0]);
    //     watch(n[1]);
    //     watch(n[2]);
    // }
    double teI2 = al2* qqq1[2] + be2* qqq1[3] + ge2* qqq1[4];
    double teI3 = al3* qqq1[2] + be3* qqq1[3] + ge3* qqq1[4];
    double enII = al * qqq2[2] + be * qqq2[3] + ge * qqq2[4];
    double teII2= al2* qqq2[2] + be2* qqq2[3] + ge2* qqq2[4];
    double teII3= al3* qqq2[2] + be3* qqq2[3] + ge3* qqq2[4];

    double pI   = qqq1[1];
    double pII  = qqq2[1];
    double rI   = qqq1[0];
    double rII  = qqq2[0];

    int ipiz = 0;
    if(pI > pII)   // Смена местами величин
    {
        double eno2  = enII;;
        double teo22  = teII2;
        double teo23 = teII3;
        double p2    = pII;
        double r2    = rII;

        double eno1  = enI;
        double teo12 = teI2;
        double teo13 = teI3;
        double p1    = pI;
        double r1    = rI;

        enI  = -eno2;
        teI2 = teo22;
        teI3 = teo23;
        pI   = p2;
        rI   = r2;

        enII = -eno1;
        teII2= teo12;
        teII3= teo13;
        pII  = p1;
        rII  = r1;
        w    = -w;
        ipiz = 1;                                                                // ???? Он точно здесь должен быть?
    }

    double cI = 0.0;
    double cII = 0.0;
    if(pI != 0.0)
    {
        cI = sqrt(ga*pI /rI );
    }
    if(pII != 0.0)
    {
        cII = sqrt(ga*pII /rII );
    }

    // cout << "!! C2 = " << cII << "     = " << ga*pII /rII << endl;
    // cout << ga << " " << pII << " " << rII << endl;

    double a   = sqrt( rI * (  g2 * pII + g1 * pI ) / 2.0 );
    double Uud = ( pII - pI ) / a;
    double Urz = - 2.0 * cII / g1 * ( 1.0 - pow(( pI / pII ),gm) );
    double Uvk = - 2.0 * ( cII + cI ) / g1;
    double Udf = enI - enII;

    int il, ip;
    double p, r, te2, te3, en;

    if (Udf < Uvk)
    {
        il = -1;
        ip = -1;
    }
    else if  ( (Udf >= Uvk) && (Udf <= Urz) )
    {
        p = pI * pow(( ( Udf - Uvk ) / ( Urz - Uvk ) ),(1.0/gm));
        il = 0;
        ip = 0;
    }
    else if ( (Udf > Urz) && (Udf <= Uud))
    {
        devtwo(enI ,pI ,rI, enII,pII,rII, w ,p);
        il = 1;
        ip = 0;
    }
    else if (Udf > Uud)
    {
        newton(enI ,pI ,rI, enII,pII,rII, w ,p);
        il=1;
        ip=1;
    }

    //*********TWO SHOCKS**********************************************
        if ( (il == 1) && (ip == 1) )
        {
            //cout << "TWO SHOCKS" << endl;
            double aI = sqrt(rI *( g2/2.0*p+g1/2.0*pI  ));
            double aII = sqrt(rII*( g2/2.0*p+g1/2.0*pII ));

            double u =  (aI*enI+aII*enII+pI-pII)/(aI+aII);
            double dI = enI -aI /rI;
            double dII = enII+aII/rII;
            dsl=dI;
            dsp=dII;
            dsc = u;

            double UU = max( fabs(dsl), fabs(dsp));
            if (UU > eps8)
            {
                time = kurant * DIST/UU;
            }
            else
            {
                time = kurant * DIST/eps8;
            }
            
            
            if(w <= dI )
            {
                en = enI;
                p = pI;
                r = rI;
                te2 = teI2;
                te3 = teI3;
            }
            else if( (w > dI)&&(w <= u) )
            {
                en=u;
                p =p;
                r =rI*aI/(aI-rI*(enI-u));
                te2=teI2;
                te3=teI3;
            }
            else if( (w > u) && (w < dII) )
            {
                en=u;
                p =p;
                r =rII*aII/(aII+rII*(enII-u));
                te2=teII2;
                te3=teII3;
            }
            else if (w >= dII)
            {
                en = enII;
                p = pII;
                r = rII;
                te2 = teII2;
                te3 = teII3;
            }
        }


    //*********LEFT - SHOCK, RIGHT - EXPANSION FAN*******************
        if ( (il == 1) && (ip == 0) )
        {
           // cout << "LEFT - SHOCK, RIGHT - EXPANSION FAN" << endl;
            double aI = sqrt(rI*( g2/2.0*p+g1/2.0*pI ));
            double aII;
            if(fabs(p-pII) < eps)
            {
                aII = rII*cII;
            }
            else
            {
                aII = gm*rII*cII*(1.0-p/pII)/(1.0- pow((p/pII),gm));
            }

            double u = (aI*enI+aII*enII+pI-pII)/(aI+aII);
            double dI  = enI-aI/rI;
            double dII =enII +cII;
            double ddII =u + cII -g1*(enII-u)/2.0;
            dsl = dI;
            dsp=dII;
            dsc = u;

            double UU = max( fabs(dsl), fabs(dsp));
            UU = max(UU, fabs(ddII));
            if (UU > eps8)
            {
                time = kurant * DIST/UU;
            }
            else
            {
                time = kurant * DIST/eps8;
            }

                if(w <= dI )
                {
                    en=enI;
                    p =pI;
                    r =rI;
                    te2=teI2;
                    te3=teI3;
                }
                if( (w > dI) && (w <= u) )
                {
                    en=u;
                    p =p;
                    r =rI*aI/(aI-rI*(enI-u));
                    te2=teI2;
                    te3=teI3;
                }
                if( (w > u) && (w <= ddII))
                {
                    double ce = cII-g1/2.0*(enII-u);
                    en=u;
                    p =p;
                    r =ga*p/ce/ce;
                    te2=teII2;
                    te3=teII3;
                }
                if ( (w > ddII) && (w < dII))
                {
                    double ce=-g1/g2*(enII-w)+2.0/g2*cII;
                    en = w - ce;
                    p = pII* pow((ce/cII),(1.0/gm));
                    r =ga*p/ce/ce;
                    te2=teII2;
                    te3=teII3;
                }
                if(w >= dII)
                {
                    en=enII;
                    p =pII;
                    r =rII;
                    te2=teII2;
                    te3=teII3;
                }
        }
    //*********TWO EXPANSION FANS**************************************
        if ((il == 0) && (ip == 0))
        {
        //cout << "TWO EXPANSION FANS" << endl;
        //printf("p = %lf\n", p);
        double aI;
        if(fabs(p-pI) < eps)
        {
        aI =rI*cI;
        }
        else
        {
        aI = gm*rI*cI*(1.0-p/pI)/(1.0- pow((p/pI),gm));
        }

        //printf("aI = %lf\n", aI);
        double aII;
        if(fabs(p-pII) < eps)
        {
        aII =rII*cII;
        }
        else
        {
        aII=gm*rII*cII*(1.0-p/pII)/(1.0- pow((p/pII),gm));
        }
        
        //printf("aII = %lf\n", aI);
        double u=(aI*enI+aII*enII+pI-pII)/(aI+aII);
        double dI  =enI  -cI  ;
        double ddI =u   -cI  -g1*(enI -u)/2.0;
        double dII =enII +cII ;
        double ddII=u   +cII -g1*(enII-u)/2.0;
        dsl=dI;
        dsp=dII;
        dsc = u;
        // printf("enII = %lf\n", enII);
        // printf("cII = %lf\n", cII);
        // printf("u = %lf\n", u);
        // printf("dI = %lf\n", dI);
        // printf("dII = %lf\n", dII);
        // printf("ddI = %lf\n", ddI);
        // printf("ddII = %lf\n", ddII);

        double UU = max( fabs(dsl), fabs(dsp));
        UU = max(UU, fabs(ddII));
        UU = max(UU, fabs(ddI));
        if (UU > eps8)
        {
            time = kurant * DIST/UU;
        }
        else
        {
            time = kurant * DIST/eps8;
        }


            if(w <= dI )
            {
                en=enI;
                p =pI;
                r =rI;
                te2=teI2;
                te3=teI3;
            }
            if ((w > dI) && (w < ddI))
            {
                        double ce=g1/g2*(enI-w)+2.0/g2*cI;
                        en=w+ce;
                        p =pI*pow((ce/cI),(1.0/gm));
                        r =ga*p/ce/ce;
                        te2=teI2;
                        te3=teI3;
            }
            if( (w >= ddI) && (w <= u))
             {
                        double ce=cI+g1/2.0*(enI-u);
                        en=u;
                        p =p;
                        r =ga*p/ce/ce;
                        te2=teI2;
                        te3=teI3;
            }
            if( (w > u) && (w <= ddII))
             {
                        double ce=cII-g1/2.0*(enII-u);
                        en=u;
                        p =p;
                        r =ga*p/ce/ce;
                        te2=teII2;
                        te3=teII3;
            }
            if( (w > ddII) && (w < dII))
             {
                        double ce=-g1/g2*(enII-w)+2.0/g2*cII;
                        en=w-ce;
                        p =pII*pow((ce/cII),(1.0/gm));
                        r =ga*p/ce/ce;
                        te2=teII2;
                        te3=teII3;
            }
            if(w >= dII)
            {
                        en=enII;
                        p =pII;
                        r =rII;
                        te2=teII2;
                        te3=teII3;
            }
          }

    //*********VAKUUM ************************************************
         if ( (il == -1) && (ip == -1))
         {

            double dI  =enI  -cI ;
            double ddI =enI  +2.0/gg1*cI;
            double dII =enII +cII ;
            double ddII=enII -2.0/gg1*cII ;

            dsl=dI;
            dsp=dII;
            dsc = (dI+dII)/2.0;

            double UU = max( fabs(dsl), fabs(dsp));
            UU = max(UU, fabs(ddII));
            UU = max(UU, fabs(ddI));
            if (UU > eps8)
            {
                time = kurant * DIST/UU;
            }
            else
            {
                time = kurant * DIST/eps8;
            }


            if(w <= dI )
            {
                en=enI;
                p =pI;
                r =rI;
                te2=teI2;
                te3=teI3;
            }
            if ( (w > dI) && (w < ddI))
            {
                double ce = gg1/gg2*(enI-w)+2.0/gg2*cI;
                en = w+ce;
                p = pI* pow((ce/cI),(1.0/gm));
                r =gga*p/ce/ce;
                te2=teI2;
                te3=teI3;
            }
            if ( (w >= ddI) && (w <= ddII))
            {
                en=w;
                p =0.0;
                r =0.0;
                te2=0.0;
                te3=0.0;
            }
            if( (w > ddII) && (w < dII))
            {
                double ce=-gg1/gg2*(enII-w)+2.0/gg2*cII;
                en=w-ce;
                p =pII* pow((ce/cII),(1.0/gm));
                r =gga*p/ce/ce;
                te2=teII2;
                te3=teII3;
            }
            if(w >= dII)
            {
                en=enII;
                p =pII;
                r =rII;
                te2=teII2;
                te3=teII3;
            }
        }


        if(ipiz == 1)
        {
            en = - en;
            double dsl1 = dsl;
            double dsp1 = dsp;
            dsl =   -dsp1;
            dsp =   -dsl1;
            dsc = -dsc;
            w = -w;
        }

        double uo = al * en  + al2 * te2 + al3 * te3;
        double vo = be * en  + be2 * te2 + be3 * te3;
        double wo = ge * en  + ge2 * te2 + ge3 * te3;


        double eo = p / g1 + 0.5 * r *(uo*uo + vo*vo + wo*wo);
        en = al*uo + be*vo + ge*wo;

        qqq[0] = (r*(en-w));
        qqq[1] = (r*(en-w)*uo+al*p);
        qqq[2] = (r*(en-w)*vo+be*p);
        qqq[3] = (r*(en-w)*wo+ge*p);
        qqq[4] = (  (en-w)*eo+en*p);


         return time;

}

void Solvers::newton(const double& enI , const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
  const double& w , double& p)
{
    double fI, fIs, fII, fIIs;
    double cI = sqrt(ga*pI /rI );
    double cII = sqrt(ga*pII/rII);
    double pn = pI*rII*cII+pII*rI*cI+(enI-enII)*rI*cI*rII*cII;
    pn = pn/(rI*cI+rII*cII);

    double pee = pn;

    int kiter = 0;
    a1:
    p = pn;
    if(p <= 0.0)
    {
        ER_S;
        cout << "negative pressure, newton" << endl;
        exit(-1);
    }

    kiter = kiter + 1;

    fI=(p-pI)/(rI*cI*sqrt(gp*p/pI+gm));
    fIs=(ga+1.0)*p/pI+(3.0*ga-1.0);
    fIs = fIs/(4.0*ga*rI*cI* pow((gp*p/pI+gm),(3.0/2.0)) );

    fII=(p-pII)/(rII*cII*sqrt(gp*p/pII+gm));
    fIIs=(ga+1.0)*p/pII+(3.0*ga-1.0);
    fIIs = fIIs/(4.0*ga*rII*cII* pow((gp*p/pII+gm),(3.0/2.0)) );


    if ( kiter == 1100)
    {
        ER_S;
        cout << "zaciklilsya v raspade,i,j,k,KOBL,kdir:" << endl;
        watch(enI);
        watch(pI);
        watch(rI);
        watch(enII);
        watch(pII);
        watch(rII);
        exit(-1);
    }

    pn = p - (fI+fII-(enI-enII))/(fIs+fIIs);

    if(fabs(pn/pee - p/pee) >= eps)
    {
        goto a1;
    }

    p = pn;

    return;
}

void Solvers::devtwo(const double& enI , const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
 const double& w , double& p)
{
    const double epsil = 10e-10;
    double kl, kp, kc, ksi, ksir, um, ksit;
    int kpizd;

    kl =  pI;
    kp =  pII;


	this->lev(enI,pI,rI,enII,pII,rII,kl,ksi);
    this->lev( enI,pI,rI,enII,pII,rII,kp,ksir);

    if (fabs(ksi) <= epsil)
    {
        um = kl;
        goto a1;
    }

    if (fabs(ksir) <= epsil)
    {
        um = kp;
        goto a1;
    }

    kpizd = 0;

    a2:
    kpizd = kpizd + 1;

    if(kpizd == 1100)
    {
        ER_S;
        cout << "zaciklilsya, devtwo.f i,j,k,KOBL,kdir:" << endl;
        watch(enI);
        watch(pI);
        watch(rI);
        watch(enII);
        watch(pII);
        watch(rII);
        exit(-1);
    }


    kc = (kl+kp)/2.0;

    this->lev( enI,pI,rI,enII,pII,rII,kc,ksit);

    if (fabs(ksit) <= epsil)
    {
        goto a3;
    }

    if ((ksi*ksit) <= 0.0)
    {
        kp  = kc;
        ksir = ksit;
    }
    else
    {
        kl  = kc;
        ksi = ksit;
    }

	goto a2;

    a3:
    um = kc;
    a1:

    p = um;

    return;
}

void Solvers::lev (const double& enI, const double& pI, const double& rI, const double& enII,//
    const double& pII, const double& rII,double& uuu, double& fee)
{
    double cI = sqrt(ga*pI /rI );
    double cII = sqrt(ga*pII/rII);

    double fI = (uuu - pI)/(rI * cI * sqrt(gp * uuu / pI + gm));

    double fII = 2.0/g1 * cII * ( pow((uuu/pII),gm) - 1.0);

	double f1= fI+fII;
    double f2= enI-enII;
	fee = f1 - f2;
    return;
}


double Solvers::HLLC_Aleksashov(const double& ro_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
    const double& ro_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
    double* P, const double& n1, const double& n2, const double& n3, const double& rad)
{
    double n[3];
    n[0] = n1;
    n[1] = n2;
    n[2] = n3;
    //int id_bn = 1;
    //int n_state = 1;
    double FR[8], FL[8];
    double UL[8], UZ[8], UR[8];
    double UZL[8], UZR[8];

    double vL[3], vR[3], bL[3], bR[3];
    double vzL[3], vzR[3], bzL[3], bzR[3];
    double qv[3];
    double aco[3][3];

    double wv = 0.0;
    double r1 = ro_L;
    double u1 = v1_L;
    double v1 = v2_L;
    double w1 = v3_L;
    double p1 = p_L;
    double bx1 = 0.0;
    double by1 = 0.0;
    double bz1 = 0.0;


    double r2 = ro_R;
    double u2 = v1_R;
    double v2 = v2_R;
    double w2 = v3_R;
    double p2 = p_R;
    double bx2 = 0.0;
    double by2 = 0.0;
    double bz2 = 0.0;

    double ro = (r2 + r1) / 2.0;
    double ap = (p2 + p1) / 2.0;
    double abx = (bx2 + bx1) / 2.0;
    double aby = (by2 + by1) / 2.0;
    double abz = (bz2 + bz1) / 2.0;


    double bk = abx * n[0] + aby * n[1] + abz * n[2];
    double b2 = kv(abx) + kv(aby) + kv(abz);

    double d = b2 - kv(bk);
    aco[0][0] = n[0];
    aco[1][0] = n[1];
    aco[2][0] = n[2];
    if (d > eps)
    {
        d = sqrt(d);
        aco[0][1] = (abx - bk * n[0]) / d;
        aco[1][1] = (aby - bk * n[1]) / d;
        aco[2][1] = (abz - bk * n[2]) / d;
        aco[0][2] = (aby * n[2] - abz * n[1]) / d;
        aco[1][2] = (abz * n[0] - abx * n[2]) / d;
        aco[2][2] = (abx * n[1] - aby * n[0]) / d;
    }
    else
    {
        double aix, aiy, aiz;
        if ((fabs(n[0]) < fabs(n[1])) && (fabs(n[0]) < fabs(n[2])))
        {
            aix = 1.0;
            aiy = 0.0;
            aiz = 0.0;
        }
        else if (fabs(n[1]) < fabs(n[2]))
        {
            aix = 0.0;
            aiy = 1.0;
            aiz = 0.0;
        }
        else
        {
            aix = 0.0;
            aiy = 0.0;
            aiz = 1.0;
        }

        double aik = aix * n[0] + aiy * n[1] + aiz * n[2];
        d = sqrt(1.0 - kv(aik));
        aco[0][1] = (aix - aik * n[0]) / d;
        aco[1][1] = (aiy - aik * n[1]) / d;
        aco[2][1] = (aiz - aik * n[2]) / d;
        aco[0][2] = (aiy * n[2] - aiz * n[1]) / d;
        aco[1][2] = (aiz * n[0] - aix * n[2]) / d;
        aco[2][2] = (aix * n[1] - aiy * n[0]) / d;
    }

    for (int i = 0; i < 3; i++)
    {
        vL[i] = aco[0][i] * u1 + aco[1][i] * v1 + aco[2][i] * w1;
        vR[i] = aco[0][i] * u2 + aco[1][i] * v2 + aco[2][i] * w2;
        bL[i] = aco[0][i] * bx1 + aco[1][i] * by1 + aco[2][i] * bz1;
        bR[i] = aco[0][i] * bx2 + aco[1][i] * by2 + aco[2][i] * bz2;
    }

    double aaL = bL[0] / sqrt(r1);
    double b2L = kv(bL[0]) + kv(bL[1]) + kv(bL[2]);
    double b21 = b2L / r1;
    double cL = sqrt(ga * p1 / r1);
    double qp = sqrt(b21 + cL * (cL + 2.0 * aaL));
    double qm = sqrt(b21 + cL * (cL - 2.0 * aaL));
    double cfL = (qp + qm) / 2.0;
    double ptL = p1 + b2L / 2.0;

    double aaR = bR[0] / sqrt(r2);
    double b2R = kv(bR[0]) + kv(bR[1]) + kv(bR[2]);
    double b22 = b2R / r2;
    double cR = sqrt(ga * p2 / r2);
    qp = sqrt(b22 + cR * (cR + 2.0 * aaR));
    qm = sqrt(b22 + cR * (cR - 2.0 * aaR));
    double cfR = (qp + qm) / 2.0;
    double ptR = p2 + b2R / 2.0;

    double aC = (aaL + aaR) / 2.0;
    double b2o = (b22 + b21) / 2.0;
    double cC = sqrt(ga * ap / ro);
    qp = sqrt(b2o + cC * (cC + 2.0 * aC));
    qm = sqrt(b2o + cC * (cC - 2.0 * aC));
    double cfC = (qp + qm) / 2.0;
    double vC1 = (vL[0] + vR[0]) / 2.0;

    double SL = min((vL[0] - cfL), (vR[0] - cfR));
    double SR = max((vL[0] + cfL), (vR[0] + cfR));

    double suR = SR - vR[0];
    double suL = SL - vL[0];
    double SM = (suR * r2 * vR[0] - ptR + ptL - suL * r1 * vL[0]) / (suR * r2 - suL * r1);

    if (SR <= SL)
    {
        printf("231\n");
    }

    double SM00 = SM;
    double SR00 = SR;
    double SL00 = SL;
    double SM01, SR01, SL01;
    if ((SM00 >= SR00) || (SM00 <= SL00))
    {
        SL = min((vL[0] - cfL), (vR[0] - cfR));
        SR = max((vL[0] + cfL), (vR[0] + cfR));
        suR = SR - vR[0];
        suL = SL - vL[0];
        SM = (suR * r2 * vR[0] - ptR + ptL - suL * r1 * vL[0]) / (suR * r2 - suL * r1);
        SM01 = SM;
        SR01 = SR;
        SL01 = SL;
        if ((SM01 >= SR01) || (SM01 <= SL01))
        {
            printf("251\n");
        }
    }


    double UU = max(fabs(SL), fabs(SR));
    double time = kurant * rad / UU;

    double upt1 = (kv(u1) + kv(v1) + kv(w1)) / 2.0;
    double sbv1 = u1 * bx1 + v1 * by1 + w1 * bz1;

    double upt2 = (kv(u2) + kv(v2) + kv(w2)) / 2.0;
    double sbv2 = u2 * bx2 + v2 * by2 + w2 * bz2;

    double e1 = p1 / g1 + r1 * upt1 + b2L / 2.0;
    double e2 = p2 / g1 + r2 * upt2 + b2R / 2.0;

    FL[0] = r1 * vL[0];
    FL[1] = r1 * vL[0] * vL[0] + ptL - kv(bL[0]);
    FL[2] = r1 * vL[0] * vL[1] - bL[0] * bL[1];
    FL[3] = r1 * vL[0] * vL[2] - bL[0] * bL[2];
    FL[4] = (e1 + ptL) * vL[0] - bL[0] * sbv1;
    FL[5] = 0.0;
    FL[6] = vL[0] * bL[1] - vL[1] * bL[0];
    FL[7] = vL[0] * bL[2] - vL[2] * bL[0];

    FR[0] = r2 * vR[0];
    FR[1] = r2 * vR[0] * vR[0] + ptR - kv(bR[0]);
    FR[2] = r2 * vR[0] * vR[1] - bR[0] * bR[1];
    FR[3] = r2 * vR[0] * vR[2] - bR[0] * bR[2];
    FR[4] = (e2 + ptR) * vR[0] - bR[0] * sbv2;
    FR[5] = 0.0;
    FR[6] = vR[0] * bR[1] - vR[1] * bR[0];
    FR[7] = vR[0] * bR[2] - vR[2] * bR[0];

    UL[0] = r1;
    UL[4] = e1;
    UR[0] = r2;
    UR[4] = e2;


    for (int ik = 0; ik < 3; ik++)
    {
        UL[ik + 1] = r1 * vL[ik];
        UL[ik + 5] = bL[ik];
        UR[ik + 1] = r2 * vR[ik];
        UR[ik + 5] = bR[ik];
    }

    for (int ik = 0; ik < 8; ik++)
    {
        UZ[ik] = (SR * UR[ik] - SL * UL[ik] + FL[ik] - FR[ik]) / (SR - SL);
    }

    double suRm = suR / (SR - SM);
    double suLm = suL / (SL - SM);
    double rzR = r2 * suRm;
    double rzL = r1 * suLm;
    vzR[0] = SM;
    vzL[0] = SM;
    double ptzR = ptR + r2 * suR * (SM - vR[0]);
    double ptzL = ptL + r1 * suL * (SM - vL[0]);
    double ptz = (ptzR + ptzL) / 2.0;
    bzR[0] = UZ[5];
    bzL[0] = UZ[5];

    vzR[1] = UZ[2] / UZ[0];
    vzR[2] = UZ[3] / UZ[0];
    vzL[1] = vzR[1];
    vzL[2] = vzR[2];

    vzR[1] = vR[1] + UZ[5] * (bR[1] - UZ[6]) / suR / r2;
    vzR[2] = vR[2] + UZ[5] * (bR[2] - UZ[7]) / suR / r2;
    vzL[1] = vL[1] + UZ[5] * (bL[1] - UZ[6]) / suL / r1;
    vzL[2] = vL[2] + UZ[5] * (bL[2] - UZ[7]) / suL / r1;

    bzR[1] = UZ[6];
    bzR[2] = UZ[7];
    bzL[1] = bzR[1];
    bzL[2] = bzR[2];

    double sbvz = (UZ[5] * UZ[1] + UZ[6] * UZ[2] + UZ[7] * UZ[3]) / UZ[0];

    double ezR = e2 * suRm + (ptz * SM - ptR * vR[0] + UZ[5] * (sbv2 - sbvz)) / (SR - SM);
    double ezL = e1 * suLm + (ptz * SM - ptL * vL[0] + UZ[5] * (sbv1 - sbvz)) / (SL - SM);

    if (fabs(UZ[5]) < eps)
    {
        vzR[1] = vR[1];
        vzR[2] = vR[2];
        vzL[1] = vL[1];
        vzL[2] = vL[2];
        bzR[1] = bR[1] * suRm;
        bzR[2] = bR[2] * suRm;
        bzL[1] = bL[1] * suLm;
        bzL[2] = bL[2] * suLm;
    }
    UZL[0] = rzL;
    UZL[4] = ezL;
    UZR[0] = rzR;
    UZR[4] = ezR;

    for (int ik = 0; ik < 3; ik++)
    {
        UZL[ik + 1] = vzL[ik] * rzL;
        UZL[ik + 5] = bzL[ik];
        UZR[ik + 1] = vzR[ik] * rzR;
        UZR[ik + 5] = bzR[ik];
    }

    if (SL > wv)
    {
        P[0] = FL[0] - wv * UL[0];
        P[4] = FL[4] - wv * UL[4];
        for (int ik = 1; ik < 4; ik++)
        {
            qv[ik - 1] = FL[ik] - wv * UL[ik];
        }
    }
    else if ((SL <= wv) && (SM >= wv))
    {
        P[0] = FL[0] + SL * (rzL - r1) - wv * UZL[0];
        P[4] = FL[4] + SL * (ezL - e1) - wv * UZL[4];
        for (int ik = 1; ik < 4; ik++)
        {
            qv[ik - 1] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
        }
    }
    else if ((SM <= wv) && (SR >= wv))
    {
        P[0] = FR[0] + SR * (rzR - r2) - wv * UZR[0];
        P[4] = FR[4] + SR * (ezR - e2) - wv * UZR[4];
        for (int ik = 1; ik < 4; ik++)
        {
            qv[ik - 1] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
        }
    }
    else if (SR < wv)
    {
        P[0] = FR[0] - wv * UR[0];
        P[4] = FR[4] - wv * UR[4];
        for (int ik = 1; ik < 4; ik++)
        {
            qv[ik - 1] = FR[ik] + -wv * UR[ik];
        }
    }
    else
    {
        printf("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\n");
        cout << SL << " " << SR << " " << wv << endl;
    }


    P[1] = aco[0][0] * qv[0] + aco[0][1] * qv[1] + aco[0][2] * qv[2];
    P[2] = aco[1][0] * qv[0] + aco[1][1] * qv[1] + aco[1][2] * qv[2];
    P[3] = aco[2][0] * qv[0] + aco[2][1] * qv[1] + aco[2][2] * qv[2];

    return time;
}


double Solvers::HLLDQ_Korolkov(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
    const double& Bx_L, const double& By_L, const double& Bz_L, const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
    const double& Bx_R, const double& By_R, const double& Bz_R, double* P, double& PQ, const double& n1, const double& n2, const double& n3, const double& rad, int metod)
{// Не работает, если скорость грани не нулевая
 // Нормаль здесь единичная по осям координат

    double bx_L = Bx_L / spi4;
    double by_L = By_L / spi4;
    double bz_L = Bz_L / spi4;

    double bx_R = Bx_R / spi4;
    double by_R = By_R / spi4;
    double bz_R = Bz_R / spi4;

    double t1 = 0.0;
    double t2 = 0.0;
    double t3 = 0.0;

    double m1 = 0.0;
    double m2 = 0.0;
    double m3 = 0.0;

    if (n1 > 0.1)
    {
        t2 = 1.0;
        m3 = 1.0;
    }
    else if (n2 > 0.1)
    {
        t3 = 1.0;
        m1 = 1.0;
    }
    else if (n3 > 0.1)
    {
        t1 = 1.0;
        m2 = 1.0;
    }
    else if (n1 < -0.1)
    {
        t3 = -1.0;
        m2 = -1.0;
    }
    else if (n2 < -0.1)
    {
        t1 = -1.0;
        m3 = -1.0;
    }
    else if (n3 < -0.1)
    {
        t1 = -1.0;
        m2 = -1.0;
    }
    else
    {
        printf("EROROR 1421  normal_error\n");
    }


    double u1, v1, w1, u2, v2, w2;
    u1 = v1_L * n1 + v2_L * n2 + v3_L * n3;
    v1 = v1_L * t1 + v2_L * t2 + v3_L * t3;
    w1 = v1_L * m1 + v2_L * m2 + v3_L * m3;
    u2 = v1_R * n1 + v2_R * n2 + v3_R * n3;
    v2 = v1_R * t1 + v2_R * t2 + v3_R * t3;
    w2 = v1_R * m1 + v2_R * m2 + v3_R * m3;

    double bn1, bt1, bm1, bn2, bt2, bm2;
    bn1 = bx_L * n1 + by_L * n2 + bz_L * n3;
    bt1 = bx_L * t1 + by_L * t2 + bz_L * t3;
    bm1 = bx_L * m1 + by_L * m2 + bz_L * m3;
    bn2 = bx_R * n1 + by_R * n2 + bz_R * n3;
    bt2 = bx_R * t1 + by_R * t2 + bz_R * t3;
    bm2 = bx_R * m1 + by_R * m2 + bz_R * m3;

    //cout << " = " << bt2 * bt2 + bm2 * bm2 << endl;

    double sqrtroL = sqrt(ro_L);
    double sqrtroR = sqrt(ro_R);
    double ca_L = bn1 / sqrtroL;
    double ca_R = bn2 / sqrtroR;
    double cL = sqrt(ggg * p_L / ro_L);
    double cR = sqrt(ggg * p_R / ro_R);

    double bb_L = kv(bx_L) + kv(by_L) + kv(bz_L);
    double bb_R = kv(bx_R) + kv(by_R) + kv(bz_R);

    double aL = (kv(bx_L) + kv(by_L) + kv(bz_L)) / ro_L;
    double aR = (kv(bx_L) + kv(by_L) + kv(bz_L)) / ro_L;

    double uu_L = (kv(v1_L) + kv(v2_L) + kv(v3_L)) / 2.0;
    double uu_R = (kv(v1_R) + kv(v2_R) + kv(v3_R)) / 2.0;

    double cfL = sqrt((ggg * p_L + bb_L + //
        sqrt(kv(ggg * p_L + bb_L) - 4.0 * ggg * p_L * kv(bn1))) / (2.0 * ro_L));
    double cfR = sqrt((ggg * p_R + bb_R + //
        sqrt(kv(ggg * p_R + bb_R) - 4.0 * ggg * p_R * kv(bn2))) / (2.0 * ro_R));


    double SL = min(u1, u2) - max(cfL, cfR);
    double SR = max(u1, u2) + max(cfL, cfR);

    double pTL = p_L + bb_L / 2.0;
    double pTR = p_R + bb_R / 2.0;

    double suR = (SR - u2);
    double suL = (SL - u1);

    double SM = (suR * ro_R * u2 - suL * ro_L * u1 - pTR + pTL) //
        / (suR * ro_R - suL * ro_L);

    double PTT = (suR * ro_R * pTL - suL * ro_L * pTR + ro_L * ro_R * suR * suL * (u2 - u1))//
        / (suR * ro_R - suL * ro_L);

    double UU = max(fabs(SL), fabs(SR));
    double time = kurant * rad / UU;

    double FL[9], FR[9], UL[9], UR[9];

    double e1 = p_L / g1 + ro_L * uu_L + bb_L / 2.0;
    double e2 = p_R / g1 + ro_R * uu_R + bb_R / 2.0;


    FL[0] = ro_L * u1;
    FL[1] = ro_L * u1 * u1 + pTL - kv(bn1);
    FL[2] = ro_L * u1 * v1 - bn1 * bt1;
    FL[3] = ro_L * u1 * w1 - bn1 * bm1;
    FL[4] = (e1 + pTL) * u1 - bn1 * (u1 * bn1 + v1 * bt1 + w1 * bm1);
    //cout << uu_L << endl;
    FL[5] = 0.0;
    FL[6] = u1 * bt1 - v1 * bn1;
    FL[7] = u1 * bm1 - w1 * bn1;
    FL[8] = Q_L * u1;

    FR[0] = ro_R * u2;
    FR[1] = ro_R * u2 * u2 + pTR - kv(bn2);
    FR[2] = ro_R * u2 * v2 - bn2 * bt2;
    FR[3] = ro_R * u2 * w2 - bn2 * bm2;
    FR[4] = (e2 + pTR) * u2 - bn2 * (u2 * bn2 + v2 * bt2 + w2 * bm2);
    FR[5] = 0.0;
    FR[6] = u2 * bt2 - v2 * bn2;
    FR[7] = u2 * bm2 - w2 * bn2;
    FR[8] = Q_R * u2;

    UL[0] = ro_L;
    UL[1] = ro_L * u1;
    UL[2] = ro_L * v1;
    UL[3] = ro_L * w1;
    UL[4] = e1;
    UL[5] = bn1;
    UL[6] = bt1;
    UL[7] = bm1;
    UL[8] = Q_L;

    UR[0] = ro_R;
    UR[1] = ro_R * u2;
    UR[2] = ro_R * v2;
    UR[3] = ro_R * w2;
    UR[4] = e2;
    UR[5] = bn2;
    UR[6] = bt2;
    UR[7] = bm2;
    UR[8] = Q_R;

    double bn = (SR * UR[5] - SL * UL[5] + FL[5] - FR[5]) / (SR - SL);
    double bt = (SR * UR[6] - SL * UL[6] + FL[6] - FR[6]) / (SR - SL);
    double bm = (SR * UR[7] - SL * UL[7] + FL[7] - FR[7]) / (SR - SL);
    double bbn = bn * bn;

    double ro_LL = ro_L * (SL - u1) / (SL - SM);
    double ro_RR = ro_R * (SR - u2) / (SR - SM);
    double Q_LL = Q_L * (SL - u1) / (SL - SM);
    double Q_RR = Q_R * (SR - u2) / (SR - SM);

    if (metod == 2)   // HLLC  + mgd
    {
        double sbv1 = u1 * bn1 + v1 * bt1 + w1 * bm1;
        double sbv2 = u2 * bn2 + v2 * bt2 + w2 * bm2;

        double UZ0 = (SR * UR[0] - SL * UL[0] + FL[0] - FR[0]) / (SR - SL);
        double UZ1 = (SR * UR[1] - SL * UL[1] + FL[1] - FR[1]) / (SR - SL);
        double UZ2 = (SR * UR[2] - SL * UL[2] + FL[2] - FR[2]) / (SR - SL);
        double UZ3 = (SR * UR[3] - SL * UL[3] + FL[3] - FR[3]) / (SR - SL);
        double UZ4 = (SR * UR[4] - SL * UL[4] + FL[4] - FR[4]) / (SR - SL);
        double vzL, vzR, vLL, wLL, vRR, wRR, ppLR, btt1, bmm1, btt2, bmm2, ee1, ee2;


        double suRm = suR / (SR - SM);
        double suLm = suL / (SL - SM);
        double rzR = ro_R * suRm;
        double rzL = ro_L * suLm;

        double ptzR = pTR + ro_R * suR * (SM - u2);
        double ptzL = pTL + ro_L * suL * (SM - u1);
        double ptz = (ptzR + ptzL) / 2.0;


        vRR = UZ2 / UZ0;
        wRR = UZ3 / UZ0;
        vLL = vRR;
        wLL = wRR;

        /*vRR = v2 + bn * (bt2 - bt) / suR / ro_R;
        wRR = w2 + bn * (bm2 - bm) / suR / ro_R;
        vLL = v1 + bn * (bt1 - bt) / suL / ro_L;
        wLL = w1 + bn * (bm1 - bm) / suL / ro_L;*/

        btt2 = bt;
        bmm2 = bm;
        btt1 = btt2;
        bmm1 = bmm2;

        double sbvz = (bn * UZ1 + bt * UZ2 + bm * UZ3) / UZ0;

        ee2 = e2 * suRm + (ptz * SM - pTR * u2 + bn * (sbv2 - sbvz)) / (SR - SM);
        ee1 = e1 * suLm + (ptz * SM - pTL * u1 + bn * (sbv1 - sbvz)) / (SL - SM);

        /*if (fabs(bn) < 0.000001 )
        {
            vRR = v2;
            wRR = w2;
            vLL = v1;
            wLL = w1;
            btt2 = bt2 * suRm;
            bmm2 = bm2 * suRm;
            btt1 = bt1 * suLm;
            bmm1 = bm1 * suLm;
        }*/

        /*ppLR = (pTL + ro_L * (SL - u1) * (SM - u1) + pTR + ro_R * (SR - u2) * (SM - u2)) / 2.0;

        if (fabs(bn) < 0.000001)
        {
            vLL = v1;
            wLL = w1;
            vRR = v2;
            wRR = w2;

            btt1 = bt1 * (SL - u1) / (SL - SM);
            btt2 = bt2 * (SR - u2) / (SR - SM);

            bmm1 = bm1 * (SL - u1) / (SL - SM);
            bmm2 = bm2 * (SR - u2) / (SR - SM);

            ee1 = ((SL - u1) * e1 - pTL * u1 + ppLR * SM) / (SL - SM);
            ee2 = ((SR - u2) * e2 - pTL * u2 + ppLR * SM) / (SR - SM);
        }
        else
        {
            btt2 = btt1 = (SR * UR[6] - SL * UL[6] + FL[6] - FR[6]) / (SR - SL);
            bmm2 = bmm1 = (SR * UR[7] - SL * UL[7] + FL[7] - FR[7]) / (SR - SL);
            vLL = v1 + bn * (bt1 - btt1) / (ro_L * (SL - u1));
            vRR = v2 + bn * (bt2 - btt2) / (ro_R * (SR - u2));

            wLL = w1 + bn * (bm1 - bmm1) / (ro_L * (SL - u1));
            wRR = w2 + bn * (bm2 - bmm2) / (ro_R * (SR - u2));

            double sks1 = u1 * bn1 + v1 * bt1 + w1 * bm1 - SM * bn - vLL * btt1 - wLL * bmm1;
            double sks2 = u2 * bn2 + v2 * bt2 + w2 * bm2 - SM * bn - vRR * btt2 - wRR * bmm2;

            ee1 = ((SL - u1) * e1 - pTL * u1 + ppLR * SM + bn * sks1) / (SL - SM);
            ee2 = ((SR - u2) * e2 - pTR * u2 + ppLR * SM + bn * sks2) / (SR - SM);
        }*/


        double  ULL[9], URR[9], PO[9];
        ULL[0] = ro_LL;
        ULL[1] = ro_LL * SM;
        ULL[2] = ro_LL * vLL;
        ULL[3] = ro_LL * wLL;
        ULL[4] = ee1;
        ULL[5] = bn;
        ULL[6] = btt1;
        ULL[7] = bmm1;
        ULL[8] = Q_LL;

        URR[0] = ro_RR;
        URR[1] = ro_RR * SM;
        URR[2] = ro_RR * vRR;
        URR[3] = ro_RR * wRR;
        URR[4] = ee2;
        URR[5] = bn;
        URR[6] = btt2;
        URR[7] = bmm2;
        URR[8] = Q_RR;

        if (SL >= 0.0)
        {
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i];
            }
        }
        else if (SL < 0.0 && SM >= 0.0)
        {
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i] + SL * ULL[i] - SL * UL[i];
            }
        }
        else if (SR > 0.0 && SM < 0.0)
        {
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i] + SR * URR[i] - SR * UR[i];
            }
        }
        else if (SR <= 0.0)
        {
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i];
            }
        }



        double SN = max(fabs(SL), fabs(SR));

        PO[5] = -SN * (bn2 - bn1);

        P[1] = n1 * PO[1] + t1 * PO[2] + m1 * PO[3];
        P[2] = n2 * PO[1] + t2 * PO[2] + m2 * PO[3];
        P[3] = n3 * PO[1] + t3 * PO[2] + m3 * PO[3];
        P[5] = spi4 * (n1 * PO[5] + t1 * PO[6] + m1 * PO[7]);
        P[6] = spi4 * (n2 * PO[5] + t2 * PO[6] + m2 * PO[7]);
        P[7] = spi4 * (n3 * PO[5] + t3 * PO[6] + m3 * PO[7]);
        P[0] = PO[0];
        P[4] = PO[4];
        PQ = PO[8];

        double SWAP = P[4];
        P[4] = P[5];
        P[5] = P[6];
        P[6] = P[7];
        P[7] = SWAP;
        return time;

    }
    else if (metod == 3)  // HLLD
    {

        double ttL = ro_L * suL * (SL - SM) - bbn;
        double ttR = ro_R * suR * (SR - SM) - bbn;

        double vLL, wLL, vRR, wRR, btt1, bmm1, btt2, bmm2;

        if (fabs(ttL) >= 0.000001)
        {
            vLL = v1 - bn * bt1 * (SM - u1) / ttL;
            wLL = w1 - bn * bm1 * (SM - u1) / ttL;
            btt1 = bt1 * (ro_L * suL * suL - bbn) / ttL;
            bmm1 = bm1 * (ro_L * suL * suL - bbn) / ttL;
        }
        else
        {
            vLL = v1;
            wLL = w1;
            btt1 = 0.0;
            bmm1 = 0.0;
        }

        if (fabs(ttR) >= 0.000001)
        {
            vRR = v2 - bn * bt2 * (SM - u2) / ttR;
            wRR = w2 - bn * bm2 * (SM - u2) / ttR;
            btt2 = bt2 * (ro_R * suR * suR - bbn) / ttR;
            bmm2 = bm2 * (ro_R * suR * suR - bbn) / ttR;
            //cout << "tbr = " << (ro_R * suR * suR - bbn) / ttR << endl;
            //cout << "bt2 = " << bt2 << endl;
        }
        else
        {
            vRR = v2;
            wRR = w2;
            btt2 = 0.0;
            bmm2 = 0.0;
        }

        double eLL = (e1 * suL + PTT * SM - pTL * u1 + bn * //
            ((u1 * bn1 + v1 * bt1 + w1 * bm1) - (SM * bn + vLL * btt1 + wLL * bmm1))) //
            / (SL - SM);
        double eRR = (e2 * suR + PTT * SM - pTR * u2 + bn * //
            ((u2 * bn2 + v2 * bt2 + w2 * bm2) - (SM * bn + vRR * btt2 + wRR * bmm2))) //
            / (SR - SM);

        double sqrtroLL = sqrt(ro_LL);
        double sqrtroRR = sqrt(ro_RR);
        double SLL = SM - fabs(bn) / sqrtroLL;
        double SRR = SM + fabs(bn) / sqrtroRR;

        double idbn = 1.0;
        if (fabs(bn) > 0.001)
        {
            idbn = 1.0 * sign(bn);
        }
        else
        {
            idbn = 0.0;
            SLL = SM;
            SRR = SM;
        }

        double vLLL = (sqrtroLL * vLL + sqrtroRR * vRR + //
            idbn * (btt2 - btt1)) / (sqrtroLL + sqrtroRR);

        double wLLL = (sqrtroLL * wLL + sqrtroRR * wRR + //
            idbn * (bmm2 - bmm1)) / (sqrtroLL + sqrtroRR);

        double bttt = (sqrtroLL * btt2 + sqrtroRR * btt1 + //
            idbn * sqrtroLL * sqrtroRR * (vRR - vLL)) / (sqrtroLL + sqrtroRR);

        double bmmm = (sqrtroLL * bmm2 + sqrtroRR * bmm1 + //
            idbn * sqrtroLL * sqrtroRR * (wRR - wLL)) / (sqrtroLL + sqrtroRR);

        double eLLL = eLL - idbn * sqrtroLL * ((SM * bn + vLL * btt1 + wLL * bmm1) //
            - (SM * bn + vLLL * bttt + wLLL * bmmm));
        double eRRR = eRR + idbn * sqrtroRR * ((SM * bn + vRR * btt2 + wRR * bmm2) //
            - (SM * bn + vLLL * bttt + wLLL * bmmm));
        //cout << " = " << bn << " " << btt2 << " " << bmm2 << endl;
        //cout << "sbvr = " << (SM * bn + vRR * btt2 + wRR * bmm2) << endl;
        double  ULL[9], URR[9], ULLL[9], URRR[9];

        ULL[0] = ro_LL;
        ULL[1] = ro_LL * SM;
        ULL[2] = ro_LL * vLL;
        ULL[3] = ro_LL * wLL;
        ULL[4] = eLL;
        ULL[5] = bn;
        ULL[6] = btt1;
        ULL[7] = bmm1;
        ULL[8] = Q_LL;

        URR[0] = ro_RR;
        //cout << ro_RR << endl;
        URR[1] = ro_RR * SM;
        URR[2] = ro_RR * vRR;
        URR[3] = ro_RR * wRR;
        URR[4] = eRR;
        URR[5] = bn;
        URR[6] = btt2;
        URR[7] = bmm2;
        URR[8] = Q_RR;

        ULLL[0] = ro_LL;
        ULLL[1] = ro_LL * SM;
        ULLL[2] = ro_LL * vLLL;
        ULLL[3] = ro_LL * wLLL;
        ULLL[4] = eLLL;
        ULLL[5] = bn;
        ULLL[6] = bttt;
        ULLL[7] = bmmm;
        ULLL[8] = Q_LL;

        URRR[0] = ro_RR;
        URRR[1] = ro_RR * SM;
        URRR[2] = ro_RR * vLLL;
        URRR[3] = ro_RR * wLLL;
        URRR[4] = eRRR;
        URRR[5] = bn;
        URRR[6] = bttt;
        URRR[7] = bmmm;
        URRR[8] = Q_RR;

        double PO[9];

        if (SL >= 0.0)
        {
            //cout << "SL >= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i];
            }
        }
        else if (SL < 0.0 && SLL >= 0.0)
        {
            //cout << "SL < 0.0 && SLL >= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i] + SL * ULL[i] - SL * UL[i];
            }
            //cout << ULL[0] << endl;
        }
        else if (SLL <= 0.0 && SM >= 0.0)
        {
            //cout << "SLL <= 0.0 && SM >= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FL[i] + SLL * ULLL[i] - (SLL - SL) * ULL[i] - SL * UL[i];
            }
        }
        else if (SM < 0.0 && SRR > 0.0)
        {
            //cout << "SM < 0.0 && SRR > 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i] + SRR * URRR[i] - (SRR - SR) * URR[i] - SR * UR[i];
            }
            //cout << "P4 = " << URRR[4] << endl;
        }
        else if (SR > 0.0 && SRR <= 0.0)
        {
            //cout << "SR > 0.0 && SRR <= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i] + SR * URR[i] - SR * UR[i];
            }
            //cout << URR[0] << endl;
        }
        else if (SR <= 0.0)
        {
            //cout << "SR <= 0.0" << endl;
            for (int i = 0; i < 9; i++)
            {
                PO[i] = FR[i];
            }
        }



        double SN = max(fabs(SL), fabs(SR));

        PO[5] = -SN * (bn2 - bn1);

        P[1] = n1 * PO[1] + t1 * PO[2] + m1 * PO[3];
        P[2] = n2 * PO[1] + t2 * PO[2] + m2 * PO[3];
        P[3] = n3 * PO[1] + t3 * PO[2] + m3 * PO[3];
        P[5] = spi4 * (n1 * PO[5] + t1 * PO[6] + m1 * PO[7]);
        P[6] = spi4 * (n2 * PO[5] + t2 * PO[6] + m2 * PO[7]);
        P[7] = spi4 * (n3 * PO[5] + t3 * PO[6] + m3 * PO[7]);
        P[0] = PO[0];
        P[4] = PO[4];
        PQ = PO[8];

        double SWAP = P[4];
        P[4] = P[5];
        P[5] = P[6];
        P[6] = P[7];
        P[7] = SWAP;
        return time;
    }

    return time;

}




