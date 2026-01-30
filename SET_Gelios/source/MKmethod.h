#pragma once
#include "Help.h"

class MKmethod
{
public:
	// Переменные
	double R_[8];       // Радиусы начальных зон (и их количество)
	double alpha_[8];       // Коэффициенты зон при перезарядке(и их количество)
	double gam_[8];       // все гамма для запуска с границы
	double A0_;         // Главный интегралл сразу посчитанный 
	double A1_;         // Главный интегралл сразу посчитанный 
	double A2_;         // Главный интегралл сразу посчитанный 
	int num_area;      // Количество зон (именно количество дополнительных зон)
	double Int_[501];
	double Int_002[51][50];
	double Int_00625[51][50];
	double Int_02[51][50];
	double Int_055[51][50];

	MKmethod(void);


    // Здесь сами алгоритмы розыгрыша

	// розыгрыш начальных параметров запуска
	bool Init_Parametrs(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_);
	// Возвращает false, если не нужно запускать основной атом

	bool Init_Parametrs_(std::mt19937& gen, std::uniform_real_distribution<double>& dis, vector <double>& mu_,
		vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, vector <double>& X_);

	bool Init_Parametrs2(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_,
		vector <double>& Wr_, vector <double>& X_);

	int Init(Sensor* sens, vector <double>& mu_, vector <double>& Wt_, vector <double>& Wp_, vector <double>& Wr_, double& X_);

	// Розыгрыш скорости при перезарядке
	bool Change_Velosity(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Второй алгоритм Маламы
	bool Change_Velosity2(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Самый первый алгоритм Маламы
	bool Change_Velosity3(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Мой алгоритм (по второму Маламы) с табличным вычислением весов
	bool Change_Velosity4(Sensor* sens, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Мой алгоритм (по первому Маламы) с табличным вычислением весов
	bool Change_Velosity4(int& s1, int& s2, int& s3, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);  // Мой алгоритм (по первому Маламы) с табличным вычислением весов
	// Возвращает false, если не нужно запускать основной атом

	template<typename Random_type, typename Distribution_type>
	bool Change_Velosity5(Random_type& gen, Distribution_type& dis, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr_, vector <double>& Wthe_,//
		vector <double>& Wphi_, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex = 0.0,//
		const double& y_ex = 0.0, const double& z_ex = 0.0);

	double for_Wr_1(const double& Z, const double& gam, const double& ur);
	double H_Wr_1(const double& gam1, const double& gam2, const double& V, const double& ur, const double& p, const double& ksi);
	double for_Wr_2(const double& Z, const double& gam, const double& ur, const double& ut);
	double H_Wr_2(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& p, const double& ksi);


	void TEST(void);   // Функция - тестирующая разные вспомогательные функции, нужна была для отладки
	double play_mho(Sensor* sens, const double& c);
	double play_mho(Sensor& sens, const double& c);
	double play_mho(int& s1, int& s2, int& s3, const double& c);
	double play_mho2(Sensor* sens, const double& c, const double& p, const double& p1);
	double play_mho2(Sensor& sens, const double& c, const double& p, const double& p1);
	double play_mho2(int& s1, int& s2, int& s3, const double& c, const double& p, const double& p1);
	// Разыгрываем \mho 
	double norm_mho(const double& c);  // Ряд функции Бесселя
	double norm_mho2(const double& c); // Ряд функции Struvel
	double mho_H1(const double& c, const double& k, const double& ksi);
	double mho_H2(const double& c, const double& k, const double& ksi);
	double h_mho(const double& x, const double& c);

	template<typename Random_type, typename Distribution_type>
	low_type play_mho(Random_type& gen, Distribution_type& dis, const double& c);

	template<typename Random_type, typename Distribution_type>
	low_type play_mho2(Random_type& gen, Distribution_type& dis, const double& c, const double& p, const double& p1);


	// Интерполяция интеграллов

	inline double Int_cp_1(const double& x);
	inline double Int_cp_2(const double& x);

	// Табличное вычисление интеграллов

	double Get_Int(const double& uu);
	double Get_Int002(const double& Ur, const double& uu);
	double Get_Int00625(const double& Ur, const double& uu);
	double Get_Int02(const double& Ur, const double& uu);
	double Get_Int055(const double& Ur, const double& uu);

	// Табличное вычисление интеграллов (новое)
	inline double int_1(const double& x, const double& cp);
	inline double int_2(const double& x, const double& cp);
	inline double int_3(const double& x, const double& cp);
	inline double int_1_f1(const double& x);
	inline double int_1_f2(const double& x);
	inline double int_1_f3(const double& x);
	inline double int_2_f1(const double& x);
	inline double int_2_f2(const double& x);
	inline double int_2_f3(const double& x);
	inline double int_3_f1(const double& x);
	inline double int_3_f2(const double& x);
	inline double int_3_f3(const double& x);

	double Lin_Interpolate(const double& x1, const double& y1, const double& x2, const double& y2, const double& x);

	// Для розыгрыша на сфере
	double A0(const double& Y);  // Норма функции распределения (название совпадает с файлом)    ПРОВЕРЕНО!
	double G(const double& gam, const double& Y);  // Норма плотности     ПРОВЕРЕНО!
	double F(const double& X, const double& gam, const double& Y);  // Функция распределения для X   // ПРОВЕРЕНО!
	double F0(const double& X, const double& Y);  // Функция распределения для X   // ПРОВЕРЕНО!
	double FI(const double& Z, const double& X, const double& gam, const double& Y);   // ПРОВЕРЕНО!
	double R(const double& X, const double& Y);              // ПРОВЕРЕНО!

	double FF(const double& gam, const double& Yr);

	// Для розыгрыша перезарядки внутри области
	double f2(const double& V, const double& gam, const double& ur, const double& ut);
	double f2k(const double& V, const double& gam, const double& ur, const double& ut);  // учтено два члена ряда

	// Вспомогательные функции
	double Hx(const double& gam1, const double& gam2, const double& X, const double& Y, const double& ksi);
	double Hwr(const double& gam1, const double& gam2, const double& Z, const double& X, const double& Y, const double& ksi);

	double Hvr(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& ksi);
	double Hvrk(const double& gam1, const double& gam2, const double& V, const double& ur, const double& ut, const double& ksi);


private:
	

};


inline double MKmethod::Int_cp_1(const double& x)
{
	if (x < 0.0)
	{
		cout << "Error  1789 fuihferfrefdewed   " << x << endl;
		exit(-1);
	}
	else if (x <= 1.0)
	{
		return 1.0818657860101877 + 0.285218386210128 * kv(x) - //
			0.0013036140991678873 * pow(x, 3) - 0.027948426972744425 * pow(x, 4) - //
			0.02488624181111631 * pow(x, 5) + 0.05903086910650146 * pow(x, 6) - //
			0.07375065435261574 * pow(x, 7) + 0.061880698880672345 * pow(x, 8) - //
			0.02992298377445923 * pow(x, 9) + 0.006218848236075065 * pow(x, 10);
	}
	else if (x <= 10.0)
	{
		return 1.163413866981175 - 0.2948788733781214 * x + // 
			0.7222333466034818 * kv(x) - 0.33744822137483277 * pow(x, 3) + //
			0.1011120159809275 * pow(x, 4) - 0.020694315456416357 * pow(x, 5) + //
			0.0029240572198959364 * pow(x, 6) - 0.00028083329445766233 * pow(x, 7) + //
			0.000017508447325531657 * pow(x, 8) - 638.9792015974883 * pow((0.1 * x), 9) + //
			103.58193963558761 * pow((0.1 * x), 10);
	}
	else
	{
		cout << "Error  1805 fuihferfrefdewed   " << x << endl;
		exit(-1);
	}
}

inline double MKmethod::Int_cp_2(const double& x)
{
	if (x < 0.0)
	{
		cout << "Error  1771 dwcewrc   " << x << endl;
		exit(-1);
	}
	else if (x <= 1.0)
	{
		return 1.08187 + 0.0000148503 * x + 0.284859 * kv(x) + 0.00139356 * pow(x, 3) - //
			0.0392891 * pow(x, 4) + 0.00585294 * pow(x, 5) + 0.00170631 * pow(x, 6);
	}
	else if (x <= 10.0)
	{
		return 1.1427 - 0.227066 * x + 0.631423 * kv(x) - 0.271463 * pow(x, 3) + 0.0720498 * pow(x, 4) - //
			0.0125149 * pow(x, 5) + 0.00142189 * pow(x, 6) - 0.000101751 * pow(x, 7) + //
			415.925 * pow(0.1 * x, 8) - 74.0093 * pow(0.1 * x, 9);
	}
	else
	{
		cout << "Error  1893 vertrwx   " << x << endl;
		exit(-1);
	}
}

inline double MKmethod::int_1_f1(const double& x)
{
	if (x <= 1.0)
	{
		return 6.283185155644284 + 0.000024846677279866114 * x +//
			2.0934329078277405 * x * x + 0.008055998193903208 * x * x * x -//
			0.2355169235647438 * x * x * x * x + 0.03820480582423355 * x * x * x * x * x +//
			0.006992274370591744 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return 6.437524091973454 - 0.6331520099380095 * x +//
			3.1348881317268997 * x * x - 0.8454201478027856 * x * x * x +//
			0.1004702004260311 * x * x * x * x + 0.0009895488638964746 * x * x * x * x * x -//
			0.000920750276197054 * x * x * x * x * x * x;
	}
	else if (x <= 5)
	{
		return 4.4920780630505135 + 2.5133093267020654 * x +//
			1.1327223176567935 * x * x - 0.24648691152318875 * x * x * x +//
			0.031326738629523766 * x * x * x * x - 0.0021366031960331384 * x * x * x * x * x +//
			0.00005954097505746697 * x * x * x * x * x * x;
	}
	else if (x <= 7)
	{
		return 1.9138683588136232 + 5.350374732905213 * x -//
			0.16380205801427633 * x * x + 0.06765898334856263 * x * x * x -//
			0.011071118267864083 * x * x * x * x + 0.0008673476933852199 * x * x * x * x * x -//
			0.00002691859374483661 * x * x * x * x * x * x;
	}
	else if (x <= 50.0)
	{
		return 1.3138472469154294 + 5.336877156136497 * x +//
			0.020286308991329983 * x * x - 0.9780973533544137 * pow(x / 10.0, 3) +//
			0.26354051936651874 * pow(x / 10.0, 4) - 0.03711733070841712 * pow(x / 10.0, 5) +//
			0.002120935433043921 * pow(x / 10.0, 6);
	}
	else
	{
		cout << "Error  int_f1 > 7  =  " << x << endl;
		return 1.3138472469154294 + 5.336877156136497 * x +//
			0.020286308991329983 * x * x - 0.9780973533544137 * pow(x / 10.0, 3) +//
			0.26354051936651874 * pow(x / 10.0, 4) - 0.03711733070841712 * pow(x / 10.0, 5) +//
			0.002120935433043921 * pow(x / 10.0, 6);
	}
	return 0.0;
}

inline double MKmethod::int_1_f2(const double& x)
{
	if (x <= 1.0)
	{
		return 1.328216167939543 - 0.000004545681954848391 * x +//
			2.537368073155103 * x * x - 0.0020584991728545624 * x * x * x -//
			0.03742568018912792 * x * x * x * x - 0.010312136385277346 * x * x * x * x * x +//
			0.002767736179209713 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return 1.2959616295629246 + 0.1533684067037866 * x +//
			2.2354849981206106 * x * x + 0.3113395567715921 * x * x * x -//
			0.21656309882941488 * x * x * x * x + 0.041957500887605075 * x * x * x * x * x -//
			0.0029978773724628604 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return 1.903643456971281 - 1.4801836911099535 * x + 3.973958664572268 * x * x -//
			0.6482729779428982 * x * x * x + 0.07665007314658864 * x * x * x * x -//
			0.005369758193338703 * x * x * x * x * x + 0.00016605531871992049 * x * x * x * x * x * x;
	}
	else if (x <= 7.0)
	{
		return -4.484415105552316 + 5.3747429756690766 * x +//
			0.8892806582308143 * x * x + 0.09767316152573671 * x * x * x -//
			0.025704749778475783 * x * x * x * x + 0.0021937998296249206 * x * x * x * x * x -//
			0.00006928845984076111 * x * x * x * x * x * x;
	}
	return 0.0;
}

inline double MKmethod::int_1_f3(const double& x)
{
	if (x <= 1.0)
	{
		return 1.2938345594193854 - 0.000031719847351174835 * x +//
			1.3183710041280094 * x * x - 0.014150512069488197 * x * x * x +//
			0.4226114681928129 * x * x * x * x - 0.06985750969880078 * x * x * x * x * x -//
			0.015347864048406958 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return 0.9667460440956788 + 1.336271810704016 * x -//
			0.8687355257991665 * x * x + 1.7676868273627229 * x * x * x -//
			0.2731222764016417 * x * x * x * x + 0.004801770033831665 * x * x * x * x * x +//
			0.001780776080720323 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return 4.760566734123174 - 5.048204299463048 * x + 3.332342585228025 * x * x +//
			0.47584339615235993 * x * x * x - 0.12072786272726124 * x * x * x * x +//
			0.011870955604980658 * x * x * x * x * x - 0.0004580199652402304 * x * x * x * x * x * x;//
	}
	else if (x <= 7.0)
	{
		return 9.370493362469261 - 10.848615619383615 * x + 6.423326878282571 * x * x -//
			0.4148977656870439 * x * x * x + 0.025300923044176957 * x * x * x * x -//
			0.0010108688120876522 * x * x * x * x * x + 0.00001864423130429156 * x * x * x * x * x * x;
	}
	cout << " ERRRROROROROROROO MKMETHOD 2251  -   " << x << endl;
	return 0.0;
}

inline double MKmethod::int_1(const double& x, const double& cp)
{
	double b = 1.0 - a_2 * log(cp);
	return (cp / (sqrtpi_ * sqrtpi_ * sqrtpi_)) * (b * b * int_1_f1(x / cp) - 2.0 * a_2 * b * int_1_f2(x / cp) + kv(a_2) * int_1_f3(x / cp));
}

inline double MKmethod::int_2_f1(const double& x)
{
	if (x <= 1.0)
	{
		return 8.377571213788123 * x + 0.00047608508679086725 * x * x +//
			1.6710478320575737 * x * x * x + 0.016857530811432916 * x * x * x * x -//
			0.15132474125724812 * x * x * x * x * x + 0.030723378194358945 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return -0.11788367995598747 + 8.937936705157014 * x -//
			1.119886471634001 * x * x + 2.8831031948885917 * x * x * x -//
			0.735250146386514 * x * x * x * x + 0.10356311378423572 * x * x * x * x * x -//
			0.006231417172309398 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return 2.9044497739429858 + 2.757712415967557 * x + 4.239161941189675 * x * x +//
			0.36198838294786784 * x * x * x - 0.05737777787138304 * x * x * x * x +//
			0.004956250079677106 * x * x * x * x * x - 0.0001809238236975877 * x * x * x * x * x * x;
	}
	else if (x <= 7.0)
	{
		return 41.6323689028892 - 38.118317864344135 * x + 22.211189528076645 * x * x -//
			3.8547348524931246 * x * x * x + 0.5000517174807501 * x * x * x * x -//
			0.03446294709493891 * x * x * x * x * x + 0.0009860204070962582 * x * x * x * x * x * x;
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

inline double MKmethod::int_2_f2(const double& x)
{
	if (x <= 1.0)
	{
		return 3.8653461103376667 * x + 0.0001975300512691014 * x * x +//
			2.4468141895384012 * x * x * x + 0.005987984681429616 * x * x * x * x -//
			0.06453987836713967 * x * x * x * x * x + 0.0066920981111229004 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return -0.10983889480661446 + 4.321087890898017 * x -//
			0.7707850845797619 * x * x + 3.1237901158486583 * x * x * x -//
			0.31485222316123385 * x * x * x * x + 0.010270249760261791 * x * x * x * x * x +//
			0.0008259803934338584 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return 1.8468011847509729 + 1.1986396743254275 * x +//
			1.1421489029509448 * x * x + 2.606316149569781 * x * x * x -//
			0.2788783468089509 * x * x * x * x + 0.019815317035281846 * x * x * x * x * x -//
			0.0006379970557448899 * x * x * x * x * x * x;
	}
	else if (x <= 7.0)
	{
		return 9.480707804348185 - 8.022988228952784 * x + 5.823555900242488 * x * x +//
			1.3277220473440972 * x * x * x - 0.08074921612981413 * x * x * x * x +//
			0.0033058587723954185 * x * x * x * x * x - 0.00006041810279926061 * x * x * x * x * x * x;
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

inline double MKmethod::int_2_f3(const double& x)
{
	if (x <= 1.0)
	{
		return 2.6106039258326 * x - 0.0008357997793049243 * x * x +//
			2.0764571907368174 * x * x * x - 0.03182275644841273 * x * x * x * x +//
			0.26310521962808975 * x * x * x * x * x - 0.06034325471643871 * x * x * x * x * x * x;
	}
	else if (x <= 3.0)
	{
		return 0.20784760901369737 + 1.5920325291316857 * x +//
			2.0985329535259014 * x * x - 0.26286255221171206 * x * x * x +//
			1.4610329096120886 * x * x * x * x - 0.25626862029131897 * x * x * x * x * x +//
			0.01684969647300594 * x * x * x * x * x * x;
	}
	else if (x <= 5.0)
	{
		return -6.284115352064703 + 15.665162343948523 * x //
			- 10.766105772158252 * pow(x, 2) + 6.0821074614870465 * pow(x, 3) -//
			0.3181501403196319 * pow(x, 4) + 0.01232319194701587 * pow(x, 5)//
			- 0.00017890661550876597 * pow(x, 6);
	}
	else if (x <= 7.0)
	{
		return -4.355962170454177 + 13.835332665069274 * x //
			- 10.12766071646978 * pow(x, 2) + 5.999392227482686 * pow(x, 3) -//
			0.32171318647989955 * pow(x, 4) + 0.014181987261856027 * pow(x, 5)//
			- 0.00030579035551447497 * pow(x, 6);
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

inline double MKmethod::int_2(const double& x, const double& cp)
{
	double b = 1.0 - a_2 * log(cp);
	return -(cp * cp / (sqrtpi_ * sqrtpi_ * sqrtpi_)) * (b * b * int_2_f1(x / cp) - 2.0 * a_2 * b * int_2_f2(x / cp) + kv(a_2) * int_2_f3(x / cp));
}

inline double MKmethod::int_3(const double& x, const double& cp)
{
	double b = 1.0 - a_2 * log(cp);
	//cout << "1 = " << int_3_f1(x / cp) << endl;
	//cout << "2 = " << int_3_f2(x / cp) << endl;
	//cout << "3 = " << int_3_f3(x / cp) << endl;
	return (pow(cp, 3.0) / (sqrtpi_ * sqrtpi_ * sqrtpi_)) * (b * b * int_3_f1(x / cp) - 2.0 * a_2 * b * int_3_f2(x / cp) + kv(a_2) * int_3_f3(x / cp));
}

inline double MKmethod::int_3_f1(const double& x)
{
	if (x <= 1.0)
	{
		return 12.566370586001975 - 0.00001944816384202852 * x + 12.567558607381049 * pow(x, 2) -//
			0.010507444068349692 * pow(x, 3) + 1.2911398125420694 * pow(x, 4) - 0.05048482363937502 * pow(x, 5) -//
			0.029947937607835744 * pow(x, 6);
	}
	else if (x <= 3.0)
	{
		return 12.451555448799724 + 0.40252674353016715 * x + 12.081033298182223 * pow(x, 2) + 0.12939193776331415 * pow(x, 3) +//
			1.478876561367302 * pow(x, 4) - 0.2237491583356496 * pow(x, 5) + 0.014474521138620033 * pow(x, 6);
	}
	else if (x <= 5.0)
	{
		return 8.281962852844913 + 9.783032429527339 * x + 3.165848344887614 * pow(x, 2) + 4.711462666119614 * pow(x, 3) +//
			0.13739453130827392 * pow(x, 4) - 0.012096015315889195 * pow(x, 5) + 0.0004514225555943018 * pow(x, 6);
	}
	else if (x <= 7.0)
	{
		return 17.29966669035025 + 2.1457227895928668 * x + 5.572082327920818 * pow(x, 2) + 4.4083748449004645 * pow(x, 3) +//
			0.13640200155890422 * pow(x, 4) - 0.00854302147917508 * pow(x, 5) + 0.00022205921430504255 * pow(x, 6);
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

inline double MKmethod::int_3_f2(const double& x)
{
	if (x <= 1.0)
	{
		return 5.798024979296493 - 0.00001772671478406096 * x + 9.98769025405073 * pow(x, 2) - 0.0073593516014156535 * pow(x, 3) +//
			2.27820901023418 * pow(x, 4) - 0.03135086958655956 * pow(x, 5) - 0.030403716978821237 * pow(x, 6);
	}
	else if (x <= 3.0)
	{
		return 5.864728705834779 - 0.34799875480550213 * x + 10.760249318358127 * pow(x, 2) - 0.9493738978205943 * pow(x, 3) +//
			2.948810798239708 * pow(x, 4) - 0.29690823082625284 * pow(x, 5) + 0.015284639719532296 * pow(x, 6);
	}
	else if (x <= 5.0)
	{
		return -3.21810405152587 + 17.08054504204813 * x - 3.43427364926184 * pow(x, 2) + 5.342946243936558 * pow(x, 3) +//
			1.3459970589832742 * pow(x, 4) - 0.07448103623442687 * pow(x, 5) + 0.0021616888853670958 * pow(x, 6);
	}
	else if (x <= 7.0)
	{
		return -59.42860988375287 + 79.24709914283142 * x - 32.304295129821256 * pow(x, 2) + 12.558114389999929 * pow(x, 3) +//
			0.32138208180756495 * pow(x, 4) + 0.003976327853989966 * pow(x, 5) - 0.0003703205838879336 * pow(x, 6);
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

inline double MKmethod::int_3_f3(const double& x)
{
	if (x <= 1.0)
	{
		return 3.915885322797866 + 0.000044284982651632276 * x + 7.778946280540595 * pow(x, 2) + 0.01990833982643192 * pow(x, 3) +//
			2.710531013102172 * pow(x, 4) + 0.09480498076917274 * pow(x, 5) + 0.026454483470905545 * pow(x, 6);
	}
	else if (x <= 3.0)
	{
		return 4.371710139648201 - 1.787430730965795 * x + 10.530386754306065 * pow(x, 2) - 1.9945949141813966 * pow(x, 3) +//
			3.310710912254515 * pow(x, 4) + 0.1356460090430464 * pow(x, 5) - 0.019853464614841342 * pow(x, 6);
	}
	else if (x <= 5.0)
	{
		return 4.6940208896644435 - 6.423081982183987 * x + 18.137713177954524 * pow(x, 2) - 7.298291873534797 * pow(x, 3) +//
			5.202095367169497 * pow(x, 4) - 0.20620743501245353 * pow(x, 5) + 0.005094092519200813 * pow(x, 6);
	}
	else if (x <= 7.0)
	{
		return 194.2826492092263 - 195.77327124311523 * x + 95.83253215419927 * pow(x, 2) - 23.99281397751645 * pow(x, 3) +//
			7.170549820159753 * pow(x, 4) - 0.32568903981020303 * pow(x, 5) + 0.007955090180048264 * pow(x, 6);
	}
	cout << "Error  int_f1 > 7  =  " << x << endl;
	return 0.0;
}

