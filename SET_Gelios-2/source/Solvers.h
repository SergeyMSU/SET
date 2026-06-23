#ifndef SOLVERS_H
#define SOLVERS_H
#include <vector>
#include <math.h>
#include "Help.h"

const double eps = 10e-10;
const double eps8 = 10e-8;
const double pi = 3.14159265358979323846;
const double cpi4 = 4*pi;
const double cpi8 = 8*pi;
const double spi4 = sqrt(cpi4);
const double epsb = 1e-6;
const double eps_p = 1e-6;
const double eps_d = 1e-3;

class Solvers
{
    public:
        Solvers();
        virtual ~Solvers();

        double Godunov_Solver_Alexashov(const std::vector<double>& qqq1, const std::vector<double>& qqq2,//
                                  const std::vector<double>& n, std::vector<double>& qqq,//
                double& dsl, double& dsp, double& dsc, double w = 0.0, const double& DIST= 0.0);
            /**
            * Распадник переписан с кода на языке Фортран, предоставленного Дмитрием Борисовичем Алексашовым
            * большинство обозначений оставленны такими-же, как в исходном коде
            * эта функция находит потоки, записывая их в qqq, при этом потом их необходимо умножать на шаг по времени
            * и площадь грани (она этого не делает)
            * ga - показатель адиабаты
            */
        double HLLD_Solver_Alexashov(const std::vector<double>& qqq1, const std::vector<double>& qqq2,//
                                  const std::vector<double>& n, std::vector<double>& qqq,//
                double& dsl, double& dsp, double& dsc, double w = 0.0, const double& DIST= 0.0);

        double HLLC_Aleksashov(const double& ro_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
            const double& ro_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
            double* P, const double& n1, const double& n2, const double& n3, const double& rad);

        double HLLDQ_Korolkov(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L, const double& v3_L,//
            const double& Bx_L, const double& By_L, const double& Bz_L, const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& v3_R,//
            const double& Bx_R, const double& By_R, const double& Bz_R, double* P, double& PQ, const double& n1, const double& n2, const double& n3, const double& rad, int metod);


    protected:

    private:

    void newton(const double& enI , const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
         const double& w , double& p);

    void devtwo(const double& enI , const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
 const double& w , double& p);

    void lev (const double& enI, const double& pI, const double& rI, const double& enII,//
    const double& pII, const double& rII,double& uuu, double& fee);

    void perpendicular(const std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, bool t = false);
        /** Функция, находящая два вектора b,c - перпендикулярные a
         *   при этом считается размерность векторов = 3
         *   bool t = true - если вектор a - единичной длины - тогда будет немного быстрее нахождение этих векторов
         */

};

#endif // SOLVERS_H
