#include <filesystem>
#include <string>
#include <utility>
#include <cstdint>
#include <sstream>
#include <iostream>
#include <vector>
#include <list>
#include <thread>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <algorithm>
#include <fstream>
#include <map>
#include <math.h>
#include <cmath>
#include <limits>
#include <iterator>
#include <cstdlib>
#include <Eigen/Dense>
#include <optional>
#include <omp.h>
#include <chrono>
#include <random>

// Boost библиотека (надо подключать к компилятору отдельно)

#include <boost/parameter.hpp>
#include "boost/multi_array.hpp"


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/barycenter.h>

#define ga (5.0/3.0)          // Показатель адиабаты
#define ggg (5.0/3.0)
#define g1 (ga - 1.0)
#define kv(x) ( (x)*(x) )
#define kv3(x) ( (x)*(x)*(x) )
#define kv4(x) ( (x)*(x)*(x)*(x) )
#define kv5(x) ( (x)*(x)*(x)*(x)*(x) )
#define kv6(x) ( (x)*(x)*(x)*(x)*(x)*(x) )
#define pow3(x) ( (x)*(x)*(x) )
#define pow4(x) ( (x)*(x)*(x)*(x) )
#define pow5(x) ( (x)*(x)*(x)*(x)*(x) )
#define pow6(x) ( (x)*(x)*(x)*(x)*(x)*(x) )
#define pow7(x) ( (x)*(x)*(x)*(x)*(x)*(x)*(x) )
#define kvv(x,y,z)  (kv(x) + kv(y) + kv(z))
#define DOT_PRODUCT(x, y) ( ((x[0]) * (y[0])) + ((x[1]) * (y[1])) + ((x[2]) * (y[2])) )

#define pi_ 3.14159265358979323846
#define sqrtpi_ 1.77245385

#define n_p_LISM_ 3.0 //(3.0) 
#define a_2 0.1307345665  // 0.102578  // 0.10263
#define n_H_LISM_ (1.0)
#define Kn_  0.4326569 //0.4326569 // 0.4326569808 // 0.4326569808 // 6.0
#define sigma(x) (kv(1.0 - a_2 * log(x)))               // Дифференциальное сечение перезарядки
//#define sigma2(x, y) (kv(1.0 - (a_2/(1.0 - a_2 * log(y))) * log(x)))  // Для другого обезразмеривания скорости на cp
#define sigma2(x, y) (kv(1.0 - a_2 * log((x) * (y))))  // Для другого обезразмеривания скорости на cp

typedef CGAL::Exact_predicates_inexact_constructions_kernel KKexact;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, KKexact> Vb;
typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, KKexact> Vb2;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Triangulation_data_structure_2<Vb2> Tds2;
typedef CGAL::Delaunay_triangulation_3<KKexact, Tds> Delaunay;
typedef CGAL::Delaunay_triangulation_2<KKexact, Tds2> Delaunay2;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay2::Vertex_handle Vertex_handle2;
typedef Delaunay2::Face_handle Face_handle;
typedef KKexact::FT FT;

typedef KKexact::Point_3 Point;
typedef KKexact::Vector_3 Vector;
typedef KKexact::Tetrahedron_3 Tetrahedron;
typedef KKexact::Point_2 Point2;

#define sqrtpi_ 1.77245385
const double const_pi = 3.14159265358979323846;
const double sqrt_pi = sqrt(const_pi);
const double cpi4 = 4.0 * const_pi;
const double cpi8 = 8.0 * const_pi;
const double spi4 = sqrt(cpi4);
const double eps = 1E-12;
const double epsb = 1E-4;
const double eps_p = 1E-6;
const double eps_d = 1E-3;
const double MF_meDmp = (1.0 / 1836.15);  // Отношения массы электрона к массе протона

using namespace std;


// Глобальные переменные (аналогично common-блоку Фортрана)
double host_time_all;
std::int32_t host_N_cell;   // integer(4) в Фортране
std::int32_t host_N_gran;   // читаем, но не используем
const int host_num_param = 16;

// Матрицы Eigen (хранятся по столбцам, как в Фортране)
Eigen::Matrix<double, host_num_param, Eigen::Dynamic> host_Cell_par;  // параметры ячеек
Eigen::Matrix<double, 3, Eigen::Dynamic> host_Cell_center;           // центры ячеек

Delaunay2* Delone2; // указатель на триангуляцию (инициализируется в BuildDelaunay2D)


// Считываем сначала всю сетку (центры ячеек и начальные переменные)
void Set_Storage(string name)
{
    // Открываем бинарный файл
    std::ifstream file(name, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open file FCMHD_4.bin" << std::endl;
        return;
    }

    // Читаем три скаляра: время, число ячеек, число граней
    file.read(reinterpret_cast<char*>(&host_time_all), sizeof(double));
    file.read(reinterpret_cast<char*>(&host_N_cell), sizeof(std::int32_t));
    file.read(reinterpret_cast<char*>(&host_N_gran), sizeof(std::int32_t));

    if (!file) {
        std::cerr << "ERROR: Failed to read header (time, N_cell, N_gran)" << std::endl;
        return;
    }

    std::cout << "AA = " << host_time_all << " " << host_N_cell << " " << host_N_gran << std::endl;

    // Изменяем размер матриц под реальное число ячеек
    host_Cell_par.resize(host_num_param, host_N_cell);
    host_Cell_center.resize(3, host_N_cell);

    // Читаем host_Cell_par (идут сразу после заголовка)
    file.read(reinterpret_cast<char*>(host_Cell_par.data()),
        host_N_cell * host_num_param * sizeof(double));
    if (!file) {
        std::cerr << "ERROR: Failed to read host_Cell_par" << std::endl;
        return;
    }

    // Читаем host_Cell_center (следует непосредственно за host_Cell_par)
    file.read(reinterpret_cast<char*>(host_Cell_center.data()),
        host_N_cell * 3 * sizeof(double));
    if (!file) {
        std::cerr << "ERROR: Failed to read host_Cell_center" << std::endl;
        return;
    }

    file.close();

    // Пример проверки, аналогичной фортрановской (можно добавить по желанию)
    // for (int i = 0; i < host_N_cell; ++i) {
    //     if (host_Cell_par(8, i) < 0.0) { // 9-й параметр (индекс 8)
    //         std::cerr << "ERROR: host_Cell_par 9 < 0 at cell " << i+1 << std::endl;
    //         // ...
    //     }
    // }

    std::cout << "Data successfully read: host_Cell_par and host_Cell_center" << std::endl;
}

void WriteCentersToTecplot(const std::string& filename) 
{
    if (host_N_cell == 0 || host_Cell_center.size() == 0) 
    {
        std::cerr << "No cell centers to write." << std::endl;
        return;
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "ERROR: Cannot open output file " << filename << std::endl;
        return;
    }

    // Заголовок Tecplot
    out << "TITLE = \"HP\"  VARIABLES = x, y, rhoH3, rhoH4\n";

    // Данные: каждая строка – x y (координаты центра ячейки)
    for (int i = 0; i < host_N_cell; ++i) {
        double x = host_Cell_center(0, i); // строка 0 (x)
        double y = host_Cell_center(1, i); // строка 1 (y)
        // z (строка 2) – нулевая, не выводим
        out << x << " " << y << " " << host_Cell_par(8, i) << " " << host_Cell_par(12, i) << "\n";
    }

    out.close();
    std::cout << "Cell centers written to " << filename << std::endl;
}

void ReadUpdatedStorage(string name)
{
    std::ifstream file(name, std::ios::binary);
    if (!file.is_open()) 
    {
        std::cerr << "ERROR: Cannot open file FCMHD_4.7_out.bin" << std::endl;
        return;
    }

    // Читаем время
    file.read(reinterpret_cast<char*>(&host_time_all), sizeof(double));
    if (!file) {
        std::cerr << "ERROR: Failed to read host_time_all from update file" << std::endl;
        return;
    }

    // Проверяем, что размер host_Cell_par соответствует ожидаемому
    if (host_N_cell == 0) {
        std::cerr << "ERROR: host_N_cell is not set. Call Set_Storage first." << std::endl;
        return;
    }

    // Читаем массив параметров ячеек прямо в данные Eigen-матрицы
    file.read(reinterpret_cast<char*>(host_Cell_par.data()),
        host_N_cell * host_num_param * sizeof(double));
    if (!file) {
        std::cerr << "ERROR: Failed to read host_Cell_par from update file" << std::endl;
        return;
    }

    // Читаем проверочную константу (можно проигнорировать)
    double cf;
    file.read(reinterpret_cast<char*>(&cf), sizeof(double));
    if (!file) 
    {
        std::cerr << "ERROR: Failed to read cf from update file" << std::endl;
        return;
    }

    std::cout << "Read cf = " << cf << " (expected 321.0)" << std::endl;
    file.close();
    std::cout << "Updated data successfully read." << std::endl;
}

void BuildDelaunay2D() 
{
    // Собираем точки с информацией (индекс ячейки)
    std::vector<std::pair<Point2, size_t>> points;
    points.reserve(host_N_cell);
    for (int i = 0; i < host_N_cell; ++i) 
    {
        double x = host_Cell_center(0, i);
        double y = host_Cell_center(1, i);
        points.emplace_back(Point2(x, y), i);
    }

    // Если триангуляция уже существовала, удаляем старую
    if (Delone2) 
    {
        delete Delone2;
    }
    // Строим новую триангуляцию
    Delone2 = new Delaunay2(points.begin(), points.end());
}

std::array<double, 3> barycentric_coordinates_2d(const Point2& p, const Point2& p0, const Point2& p1, const Point2& p2) 
{
    // Вычисляем площади треугольников
    double area = CGAL::area(p0, p1, p2); // площадь опорного треугольника
    // Предполагаем area != 0 (треугольник не вырожден)
    double area0 = CGAL::area(p, p1, p2);
    double area1 = CGAL::area(p0, p, p2);
    double area2 = CGAL::area(p0, p1, p);

    return { area0 / area, area1 / area, area2 / area };
}

// 2. Функция получения значений параметров по координате (x, y)
// Возвращает true, если точка найдена внутри выпуклой оболочки, и записывает интерполированные параметры в result
// hint — предыдущая грань для ускорения поиска (опционально)
bool InterpolateAtPoint(double x, double y, std::vector<double>& result,
    Face_handle hint, Face_handle& next_face) 
{
    if (!Delone2) 
    {
        std::cerr << "ERROR: Delaunay triangulation not built. Call BuildDelaunay2D first." << std::endl;
        return false;
    }
    if (host_N_cell == 0) return false;

    Point2 query(x, y);
    Face_handle face = Delone2->locate(query, hint);
    if (face == nullptr) 
    {
        return false; // не должно происходить, но на всякий случай
    }

	if (Delone2->is_infinite(face)) 
	{
		return false; // точка вне области триангуляции
	}

    // Проверяем, лежит ли точка внутри треугольника (или на границе)
    // Получаем вершины треугольника
    Point2 p0 = face->vertex(0)->point();
    Point2 p1 = face->vertex(1)->point();
    Point2 p2 = face->vertex(2)->point();

    // Вычисляем барицентрические координаты
    auto coords = barycentric_coordinates_2d(query, p0, p1, p2);

    // Проверяем, что все координаты неотрицательны (с учётом погрешности)
    const double eps = 1e-12;
    if (coords[0] < -eps || coords[1] < -eps || coords[2] < -eps) 
    {
        return false; // точка снаружи треугольника (возможно, снаружи выпуклой оболочки)
    }

    // Получаем индексы ячеек, соответствующих вершинам
    size_t i0 = face->vertex(0)->info();
    size_t i1 = face->vertex(1)->info();
    size_t i2 = face->vertex(2)->info();

	if (i0 >= host_N_cell || i1 >= host_N_cell || i2 >= host_N_cell || i0 < 0 || i1 < 0 || i2 < 0)
	{
		cout << "Error 0w394t7y37gfoe4t4t3ef   " << i0 << " " << i1 << " " << i2 << " " << host_N_cell << endl;
		exit(-1);
	}

    // Интерполяция параметров (линейная)
    result.resize(host_num_param);
    for (int p = 0; p < host_num_param; ++p) 
    {
        double val0 = host_Cell_par(p, i0); // предполагается, что host_Cell_par — глобальная матрица параметров
        double val1 = host_Cell_par(p, i1);
        double val2 = host_Cell_par(p, i2);
        result[p] = coords[0] * val0 + coords[1] * val1 + coords[2] * val2;
    }

    next_face = face; // сохраняем найденный треугольник для следующего вызова
    return true;
}

double Velosity_1(const double& u, const double& cp)
{
    if (u < 0.00001)
    {
        return 2.0 * cp / sqrtpi_ + 2.0 * u * u / (3.0 * cp * sqrtpi_) - u * u * u * u / (15.0 * cp * cp * cp * sqrtpi_);
    }
    else
    {
        return  exp(-u * u / kv(cp)) * cp / sqrtpi_ + (u + kv(cp) / (2.0 * u)) * erf(u / cp);
    }
}

double maxwell_distribution(double rho, double cp,
    double u, double v, double w,
    double vx, double vy, double vz)
{
    // Квадрат относительной скорости
    double dvx = vx - u;
    double dvy = vy - v;
    double dvz = vz - w;
    double v_sq = dvx * dvx + dvy * dvy + dvz * dvz;

    // Нормировочная константа
    double norm = rho / (const_pi * sqrt(const_pi) * cp * cp * cp);

    // Экспоненциальная часть
    double exponent = exp(-v_sq / (cp * cp));

    return norm * exponent;
}

// Функция линейной интерполяции
double linear_interpolate(double x, const std::vector<double>& x_arr,
    const std::vector<double>& y_arr) 
{
	int N1 = x_arr.size();
	int N2 = y_arr.size();

	if (N1 != N2)
	{
		cout << "error ieguhyuifgew347yr " << endl;
		exit(-1);
	}

    // Если x за пределами массива, возвращаем крайние значения
    if (x <= x_arr[0]) return y_arr[0];
    if (x >= x_arr.back()) return y_arr.back();

    // Находим индекс i такой, что x_arr[i] <= x < x_arr[i+1]
    auto it = std::upper_bound(x_arr.begin(), x_arr.end(), x);
    int i = std::distance(x_arr.begin(), it) - 1;

	if (i >= N1 - 1 || i < 0)
	{
		cout << "error gertgerfgerferf     " << i << "   " << N1 << endl;
		exit(-1);
	}

    // Линейная интерполяция
    double x0 = x_arr[i];
    double x1 = x_arr[i + 1];
    double y0 = y_arr[i];
    double y1 = y_arr[i + 1];

    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

double velocity_to_energy(double v) // безразмерная скорость в эВ
{
    double kms = 0.0963358;    // 1 км/с
    return kv(v / kms) * 5.226 * 0.001;
}

// Функция для создания массива энергий из массива скоростей
std::vector<double> create_energy_from_velocity(const std::vector<double>& v_mas, int N)
{
    std::vector<double> energy_mas(N);
    for (int i = 0; i < N; i++)
    {
        energy_mas[i] = velocity_to_energy(v_mas[i]);
    }
    return energy_mas;
}


inline double int_1_f1(const double& x)
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

inline double int_1_f2(const double& x)
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

inline double int_1_f3(const double& x)
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

inline double int_1(const double& x, const double& cp)
{
	double b = 1.0 - a_2 * log(cp);
	return (cp / (sqrtpi_ * sqrtpi_ * sqrtpi_)) * (b * b * int_1_f1(x / cp) - 2.0 * a_2 * b * int_1_f2(x / cp) + kv(a_2) * int_1_f3(x / cp));
}

void Culc_ENA()
{
	cout << "Start: Culc_ENA" << endl;
	std::ofstream outFile2("jENA_angle_test.txt");

	for (double the = 0.008; the < const_pi - 0.008; the = the + const_pi / 180.0)
	{
		cout << "the = " << the * 180 / const_pi << "  degree " << endl;
		//double the = pi / 2.0 + pi/1000.0;
		double LOS_x = cos(the);   // Направление интегрирования
		double LOS_y = sin(the);

		// Задаём значения
		double ae = 0.00313166;    // 1 ае
		double kms = 0.0963358;    // 1 км/с

		double x0 = 0.0;
		double y0 = 0.0;

		double a = 0.0;   // Вспомогательная переменная
		double S = 0.0;   // Поглощение ENA
		double jENA = 0.0;
		int b;
		double rho_p, p_p, u_p, v_p, cp_p;
		double rho_H3, p_H3, u_H3, v_H3, cp_H3;
		double rho_H4, p_H4, u_H4, v_H4, cp_H4;
		double u, uz, nu_ex, nu_H, nu_H3, nu_H4;

		// На всякий случай нормируем LOS
		a = sqrt(kvv(0.0, LOS_x, LOS_y));
		LOS_x /= a;
		LOS_y /= a;

		// Задаём начальное положение (на 2 АЕ)
		x0 = x0 + 5.0 * LOS_x * ae;
		y0 = y0 + 5.0 * LOS_y * ae;

		int N = 1000;
		double vL = 1 * kms;   // Скорость, для которой считаем поток
		double vR = 1000 * kms;   // Скорость, для которой считаем поток
		std::vector<double> v_mas(N);
		std::vector<double> jENA_mas(N);
		std::vector<double> S_mas(N);
		//double v = 100 * kms;   // Скорость, для которой считаем поток
		double ds = ae / 30.0;   // Шаг интегрирования
		//S = 1.0;
		//jENA = 0.0;
		x0 = x0 + LOS_x * ds / 2.0;
		y0 = y0 + LOS_y * ds / 2.0;

		double step = (vR - vL) / N;

		// Заполняем массив координатами центров ячеек
		for (int i = 0; i < N; ++i)
		{
			v_mas[i] = vL + (i + 0.5) * step;
			S_mas[i] = 1.0;
			jENA_mas[i] = 0.0;

		}

		// Цикл интегрирования

        Face_handle current_face; // инициализируется nullptr (пустой)
        std::vector<double> result;

		//cout << "B1" << endl;
		while (true)
		{
			//cout << "A1" << endl;
            bool bb = InterpolateAtPoint(x0, y0, result, current_face, current_face);
			//cout << "A2" << endl;
            if (bb == false) break;

			if (true)
			{
				// Получаем параметры плазмы и водорода в ячейке
                rho_p = result[0]; //A->par[0].ro;
				p_p = result[4]; //A->par[0].p;
				u_p = result[1]; //A->par[0].u;
				v_p = result[2]; //A->par[0].v;
				cp_p = sqrt(p_p / rho_p);

				rho_H3 = result[8] * 3.0; //A->par[0].H_n[2] * 3.0;
				p_H3 = 0.5 * rho_H3 * result[11]; //0.5 * rho_H3 * A->par[0].H_T[2];
				u_H3 = result[9]; //A->par[0].H_u[2];
				v_H3 = result[10]; //A->par[0].H_v[2];
				cp_H3 = sqrt(p_H3 / rho_H3);

				rho_H4 = result[12] * 3.0; // A->par[0].H_n[3] * 3.0;
				p_H4 = 0.5 * rho_H4 * result[15]; // 0.5 * rho_H4 * A->par[0].H_T[3];
				u_H4 = result[13]; // A->par[0].H_u[3];
				v_H4 = result[14]; // A->par[0].H_v[3];
				cp_H4 = sqrt(p_H4 / rho_H4);

				// Считаем частоты

				for (int i = 0; i < N; ++i)
				{
					double v = v_mas[i];

					u = sqrt(kvv(-v * LOS_x - u_p, -v * LOS_y - v_p, 0.0));
					if (u / cp_p > 7.0)
					{
						uz = Velosity_1(u, cp_p);
						nu_ex = (rho_p * uz * sigma(uz)) / Kn_;
					}
					else
					{
						nu_ex = (rho_p * int_1(u, cp_p)) / Kn_;        // Пробуем вычислять интеграллы численно
					}

					u = sqrt(kvv(-v * LOS_x - u_H3, -v * LOS_y - v_H3, 0.0));
					if (u / cp_H3 > 7.0)
					{
						uz = Velosity_1(u, cp_H3);
						nu_H3 = (rho_H3 * uz * sigma(uz)) / Kn_;
					}
					else
					{
						nu_H3 = (rho_H3 * int_1(u, cp_H3)) / Kn_;        // Пробуем вычислять интеграллы численно
					}

					u = sqrt(kvv(-v * LOS_x - u_H4, -v * LOS_y - v_H4, 0.0));
					if (u / cp_H4 > 7.0)
					{
						uz = Velosity_1(u, cp_H4);
						nu_H4 = (rho_H4 * uz * sigma(uz)) / Kn_;
					}
					else
					{
						nu_H4 = (rho_H4 * int_1(u, cp_H4)) / Kn_;        // Пробуем вычислять интеграллы численно
					}

					nu_H = nu_H3 + nu_H4;


					double r = sqrt(x0 * x0 + y0 * y0);
					double nu_ph = 770.291 * kv(ae / r);

					S_mas[i] = S_mas[i] * exp(-(nu_ex + nu_ph) * ds / v);

					jENA_mas[i] = jENA_mas[i] + nu_H * maxwell_distribution(rho_p, cp_p, u_p, v_p, 0.0, -v * LOS_x, -v * LOS_y, 0.0) * v * S_mas[i] * ds;
				}
			}

			/*cout << rho_p << " " << rho_H3 << " " << rho_H4 << "  " << nu_H << " " << nu_ex << endl;
			cout << maxwell_distribution(rho_p, cp_p, u_p, v_p, 0.0, -v * LOS_x, -v * LOS_y, 0.0) << " " << S << endl;
			cout << "------------------" << endl;*/

			// Продвигаемся
			x0 = x0 + LOS_x * ds;
			y0 = y0 + LOS_y * ds;
		}


		double flux = 6.75186 * 1E7;

		//if (the * 180 / const_pi > 9.0 && the * 180 / const_pi < 11.0)
		if (false)
		{
			std::ofstream outFile("jENA_(10)_1.6.txt");

			if (!outFile.is_open()) 
			{
				std::cerr << "Error iueghuioeyhgiegrergr" << std::endl;
				exit(-1);
			}

			for (int i = 0; i < N; ++i)
			{
				outFile << kv(v_mas[i] / kms) * 5.226 * 0.000001 << " " << v_mas[i] / kms << " " << jENA_mas[i] * flux << std::endl;
			}

			outFile.close();
		}

		//cout << "B2" << endl;

		//cout << "END ENA   " << x0 << " " << y0 << endl;
		//cout << v << " " << jENA << endl;

		// Создаем массив энергий для jENA_mas
		std::vector<double> energy_jENA = create_energy_from_velocity(v_mas, N);

		double E1 = 15.0;
		double j1 = linear_interpolate(E1, energy_jENA, jENA_mas);

		double E2 = 29.0;
		double j2 = linear_interpolate(E2, energy_jENA, jENA_mas);

		double E3 = 55.0;
		double j3 = linear_interpolate(E3, energy_jENA, jENA_mas);

		double E4 = 110.0;
		double j4 = linear_interpolate(E4, energy_jENA, jENA_mas);

		double E5 = 209.0;
		double j5 = linear_interpolate(E5, energy_jENA, jENA_mas);

		double E6 = 439.0;
		double j6 = linear_interpolate(E6, energy_jENA, jENA_mas);

		double E7 = 872.0;
		double j7 = linear_interpolate(E7, energy_jENA, jENA_mas);

		double E8 = 1821.0;
		double j8 = linear_interpolate(E8, energy_jENA, jENA_mas);

		// ---------------------

		double E9 = 450.0;
		double j9 = linear_interpolate(E8, energy_jENA, jENA_mas);

		double E10 = 710.0;
		double j10 = linear_interpolate(E8, energy_jENA, jENA_mas);

		double E11 = 1110.0;
		double j11 = linear_interpolate(E8, energy_jENA, jENA_mas);

		double E12 = 1740.0;
		double j12 = linear_interpolate(E8, energy_jENA, jENA_mas);

		double E13 = 2730.0;
		double j13 = linear_interpolate(E8, energy_jENA, jENA_mas);

		double E14 = 4290.0;
		double j14 = linear_interpolate(E8, energy_jENA, jENA_mas);

		outFile2 << the << " " << j1 * flux << " " << j2 * flux << " " << j3 * flux << " " << 
			j4 * flux << " " << j5 * flux << " " << j6 * flux << " " << j7 * flux << " " << 
			j8 * flux << " " << j9 * flux << " " << j10 * flux << " " << j11 * flux << " " <<
			j12 * flux << " " << j13 * flux << " " << j14 * flux << std::endl;
	}

	outFile2.close();
}

void Culc_ENA_the(double the, string name)
{
	cout << "Start: Culc_ENA" << endl;
	cout << "the = " << the * 180 / const_pi << "  degree " << endl;
	//double the = pi / 2.0 + pi/1000.0;
	double LOS_x = cos(the);   // Направление интегрирования
	double LOS_y = sin(the);

	// Задаём значения
	double ae = 0.00313166;    // 1 ае
	double kms = 0.0963358;    // 1 км/с

	double x0 = 0.0;
	double y0 = 0.0;

	double a = 0.0;   // Вспомогательная переменная
	double S = 0.0;   // Поглощение ENA
	double jENA = 0.0;
	int b;
	double rho_p, p_p, u_p, v_p, cp_p;
	double rho_H3, p_H3, u_H3, v_H3, cp_H3;
	double rho_H4, p_H4, u_H4, v_H4, cp_H4;
	double u, uz, nu_ex, nu_H, nu_H3, nu_H4;

	// На всякий случай нормируем LOS
	a = sqrt(kvv(0.0, LOS_x, LOS_y));
	LOS_x /= a;
	LOS_y /= a;

	// Задаём начальное положение (на 2 АЕ)
	x0 = x0 + 5.0 * LOS_x * ae;
	y0 = y0 + 5.0 * LOS_y * ae;

	int N = 1000;
	double vL = 1 * kms;   // Скорость, для которой считаем поток
	double vR = 1000 * kms;   // Скорость, для которой считаем поток
	std::vector<double> v_mas(N);
	std::vector<double> jENA_mas(N);
	std::vector<double> S_mas(N);
	//double v = 100 * kms;   // Скорость, для которой считаем поток
	double ds = ae / 30.0;   // Шаг интегрирования
	//S = 1.0;
	//jENA = 0.0;
	x0 = x0 + LOS_x * ds / 2.0;
	y0 = y0 + LOS_y * ds / 2.0;

	double step = (vR - vL) / N;

	// Заполняем массив координатами центров ячеек
	for (int i = 0; i < N; ++i)
	{
		v_mas[i] = vL + (i + 0.5) * step;
		S_mas[i] = 1.0;
		jENA_mas[i] = 0.0;

	}

	// Цикл интегрирования

	Face_handle current_face; // инициализируется nullptr (пустой)
	std::vector<double> result;

	//cout << "B1" << endl;
	while (true)
	{
		//cout << "A1" << endl;
		bool bb = InterpolateAtPoint(x0, y0, result, current_face, current_face);
		//cout << "A2" << endl;
		if (bb == false) break;

		if (true)
		{
			// Получаем параметры плазмы и водорода в ячейке
			rho_p = result[0]; //A->par[0].ro;
			p_p = result[4]; //A->par[0].p;
			u_p = result[1]; //A->par[0].u;
			v_p = result[2]; //A->par[0].v;
			cp_p = sqrt(p_p / rho_p);

			rho_H3 = result[8] * 3.0; //A->par[0].H_n[2] * 3.0;
			p_H3 = 0.5 * rho_H3 * result[11]; //0.5 * rho_H3 * A->par[0].H_T[2];
			u_H3 = result[9]; //A->par[0].H_u[2];
			v_H3 = result[10]; //A->par[0].H_v[2];
			cp_H3 = sqrt(p_H3 / rho_H3);

			rho_H4 = result[12] * 3.0; // A->par[0].H_n[3] * 3.0;
			p_H4 = 0.5 * rho_H4 * result[15]; // 0.5 * rho_H4 * A->par[0].H_T[3];
			u_H4 = result[13]; // A->par[0].H_u[3];
			v_H4 = result[14]; // A->par[0].H_v[3];
			cp_H4 = sqrt(p_H4 / rho_H4);

			// Считаем частоты

			for (int i = 0; i < N; ++i)
			{
				double v = v_mas[i];

				u = sqrt(kvv(-v * LOS_x - u_p, -v * LOS_y - v_p, 0.0));
				if (u / cp_p > 7.0)
				{
					uz = Velosity_1(u, cp_p);
					nu_ex = (rho_p * uz * sigma(uz)) / Kn_;
				}
				else
				{
					nu_ex = (rho_p * int_1(u, cp_p)) / Kn_;        // Пробуем вычислять интеграллы численно
				}

				u = sqrt(kvv(-v * LOS_x - u_H3, -v * LOS_y - v_H3, 0.0));
				if (u / cp_H3 > 7.0)
				{
					uz = Velosity_1(u, cp_H3);
					nu_H3 = (rho_H3 * uz * sigma(uz)) / Kn_;
				}
				else
				{
					nu_H3 = (rho_H3 * int_1(u, cp_H3)) / Kn_;        // Пробуем вычислять интеграллы численно
				}

				u = sqrt(kvv(-v * LOS_x - u_H4, -v * LOS_y - v_H4, 0.0));
				if (u / cp_H4 > 7.0)
				{
					uz = Velosity_1(u, cp_H4);
					nu_H4 = (rho_H4 * uz * sigma(uz)) / Kn_;
				}
				else
				{
					nu_H4 = (rho_H4 * int_1(u, cp_H4)) / Kn_;        // Пробуем вычислять интеграллы численно
				}

				nu_H = nu_H3 + nu_H4;


				double r = sqrt(x0 * x0 + y0 * y0);
				double nu_ph = 770.291 * kv(ae / r);

				S_mas[i] = S_mas[i] * exp(-(nu_ex + nu_ph) * ds / v);

				jENA_mas[i] = jENA_mas[i] + nu_H * maxwell_distribution(rho_p, cp_p, u_p, v_p, 0.0, -v * LOS_x, -v * LOS_y, 0.0) * v * S_mas[i] * ds;
			}
		}

		/*cout << rho_p << " " << rho_H3 << " " << rho_H4 << "  " << nu_H << " " << nu_ex << endl;
		cout << maxwell_distribution(rho_p, cp_p, u_p, v_p, 0.0, -v * LOS_x, -v * LOS_y, 0.0) << " " << S << endl;
		cout << "------------------" << endl;*/

		// Продвигаемся
		x0 = x0 + LOS_x * ds;
		y0 = y0 + LOS_y * ds;
	}


	double flux = 6.75186 * 1E7;

	if (true)
	{
		std::ofstream outFile("jENA_" + name + ".txt");

		if (!outFile.is_open())
		{
			std::cerr << "Error iueghuioeyhgiegrergr" << std::endl;
			exit(-1);
		}

		for (int i = 0; i < N; ++i)
		{
			outFile << kv(v_mas[i] / kms) * 5.226 * 0.000001 << " " << v_mas[i] / kms << " " << jENA_mas[i] * flux << std::endl;
		}

		outFile.close();
	}

	

}

int main()
{
	Set_Storage("FCMHD_1.bin"); // "FCMHD_1.bin"
    ReadUpdatedStorage("FCMHD_1.8_out.bin"); // "FCMHD_1.6_out.bin"

	//WriteCentersToTecplot("centers.txt");

    BuildDelaunay2D();

    // std::vector<double> result;
	// Face_handle current_face;
    // InterpolateAtPoint(0.2, 0.2, result, current_face, current_face);

	Culc_ENA_the(3 * const_pi / 180.0, "1.8(3)");
	//Culc_ENA_the(52.62 * const_pi / 180.0, "5.5(52.62)");
	//Culc_ENA_the(103.48 * const_pi / 180.0, "5.5(103.48)");
	//Culc_ENA();
}
