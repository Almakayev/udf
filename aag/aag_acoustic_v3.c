#include "udf.h"
#include "mem.h"
#include "stdio.h"
#include <math.h>

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)
#define PRESSURE_DL(p) (p / GAMMA/ P0) // возвращает б/р давление
#define TEMPER_DL(temp) (temp / T0 - 1.0) // возвращает б/р относительную температуру
//Количество элементов в сетке
const int Nx = 1002,  Ny = 32;
//физические величины
const double OMEGA_FACTOR = 0.86;
const double ampX0 = 0.01;// Амплитуда хождения поршня
double U0, OMEGA, PERIOD, SOUND_SPEED, DENSITY_0, U_R;
const double P0 = 109511.7; //начальное давление
const double LENGTH = 0.005;// длина резонатора
const double HEIGHT = 0.0001; //высота резонатора
const double T0 = 300.0; // Начальная температура
const double R = 286.8; //Удельная газовая постоянная
const double GAMMA = 1.4;// Адиабата Пуассона
const double VISCOSITY = 1.935e-05; // Динамическая вязкость газа
int phaseCounterQuarter = 1;
int phaseCounterOther = 1;
int ELEMENTS_COUNT, PISTON_INDEX;
double *x_array, *y_array, *u, *v, *pressure, *temp, *u_acc, *v_acc, *density, *temp_acc;
double dx, dy, pMaxPiston, pMinPiston, velocityMax, tempAvgMax;

/**
Используется для инициализации параметров счета, массивов и сохр-ия мета-данных в файл
Необходимо вызывать вручную в GUI Fluent и перед заданием ГУ (DEFINE_PROFILE)
*/
DEFINE_ON_DEMAND(init_on_demand)
{
	double f, acousticBoundaryLayer;
	
	FILE *file = fopen("meta data.txt", "w");
	fprintf(file, "Length: %f, Height: %f\n", LENGTH, HEIGHT);
	SOUND_SPEED = sqrt(GAMMA * R * T0);
	fprintf(file, "Sound speed: %f\n", SOUND_SPEED);
	OMEGA = OMEGA_FACTOR * SOUND_SPEED / LENGTH;// see AA Gubaidullin, Pyatkova
	fprintf(file, "Cyclic frequency OMEGA: %f, omega factor: %f\n", OMEGA, OMEGA_FACTOR);
	PERIOD = 2.0 * M_PI / OMEGA;// период колебаний поршня
	fprintf(file, "Period Tw: %f\n", PERIOD);
	fprintf(file, "Piston amplitude ampX0: %f\n", ampX0);
	U0 = ampX0 * OMEGA;
	fprintf(file, "Piston velocity amplitude U0: %f\n", U0);
	fprintf(file, "Initial pressure: %f\n", P0);
	fprintf(file, "Initial temperature: %f\n", T0);
	DENSITY_0 = P0 / (R * T0);
	fprintf(file, "Initial density: %f\n", DENSITY_0);
	acousticBoundaryLayer = sqrt(2 * VISCOSITY / DENSITY_0 / OMEGA);
	fprintf(file, "Acoustic boundary layer thickness : %f\n", acousticBoundaryLayer);
	dx = LENGTH / Nx;
	dy = HEIGHT / Ny;
	fprintf(file, "AcousBoundLayer / dy: %f\n", acousticBoundaryLayer / dy);
	ELEMENTS_COUNT = Nx * Ny; //Общее количество элементов
	fprintf(file, "Nx: %d\t  Ny: %d\t Total elements: %d\n", Nx, Ny, ELEMENTS_COUNT);
	fprintf(file, "ampX0 / dx: %f\n", ampX0 / dx);
	
	allocateArrayMemory();
	initCoordinates();
	PISTON_INDEX = Ny - 1;//Индекс для середины поршня
	fprintf(file, "Piston index = %d\t x_piston: %f\t  y_piston: %f\t\n", PISTON_INDEX, x_array[PISTON_INDEX], y_array[PISTON_INDEX]);
	initMinAndMaxCharacteristics();
	fclose(file);
	printf("My UDF for initialization was executed...\n");
	fflush(stdout);
}

void findMaxCharacteristics() {
	int i;
	for(i = 0; i < ELEMENTS_COUNT; i++) {
		double velocity = sqrt(u[i] * u[i] + v[i] * v[i]);
		if(velocity > velocityMax) {
			velocityMax = velocity;
			U_R = 3.0 / 16.0 * SQR(velocityMax) / SOUND_SPEED;//скорость Рэлея
		}
		
		if(temp_acc[i] > tempAvgMax) {
			tempAvgMax = temp_acc[i];
		}
	}
}

/*
Данная фукнция вызывается в конце каждого шага по времени. 
Здесь происходит вызов обработичиков в промежутках/частях периода (см. признак isChanged)
*/
DEFINE_EXECUTE_AT_END(execute_at_end) 
{	
	int isChanged, i;
	double tmp, periodPart, countPeriods;
	double currentTime = RP_Get_Real("flow-time");
	Domain *domain;
	Thread *c_thread;
	cell_t cell;
	domain=Get_Domain(1);
	
	/*Пробегаем по элементам сетки и получаем значения в текущий момент времени*/
	i = 0;
	thread_loop_c(c_thread, domain) {
		begin_c_loop(cell, c_thread){
			//мгновенная скорость
			u[i] = C_U(cell, c_thread) / SOUND_SPEED;
			v[i] = C_V(cell, c_thread) / SOUND_SPEED;
			//акустическая скорость
			u_acc[i] += u[i] * CURRENT_TIMESTEP / PERIOD / U_R;
			v_acc[i] += v[i] * CURRENT_TIMESTEP / PERIOD / U_R;
			
			pressure[i] = PRESSURE_DL(C_P(cell, c_thread));
			density[i] = C_R(cell, c_thread) / DENSITY_0;
			
			temp[i] = TEMPER_DL(C_T(cell, c_thread));
			temp_acc[i] += temp[i] * CURRENT_TIMESTEP / PERIOD;
			
			i++;
		}
		end_c_loop(cell, c_thread);
	}
	
	countPeriods =  currentTime / PERIOD;
	pMaxPiston = MAX(pMaxPiston, pressure[PISTON_INDEX]);
	pMinPiston = MIN(pMinPiston, pressure[PISTON_INDEX]);
	printf("Tw: %f\t, pMaxPiston: %f\t, pMinPiston: %f\t, velocityMax: %f\t, tempAvgMax: %f\n", countPeriods, pMaxPiston, pMinPiston, velocityMax, tempAvgMax);
	fflush(stdout);
	
	//После каждой 1/30 части периода
	isChanged = isPeriodPartChanged(currentTime, PERIOD / 30.0);
	if(isChanged) {
		findMaxCharacteristics();
		handlePistonCharacteristics(countPeriods);
	}
	
	//После каждой 1/4 части периода
	isChanged = isPeriodPartChanged(currentTime, 0.25 * PERIOD);
	if(isChanged) {
		int phase = phaseCounterQuarter % 4;
		handleQuarterOfPeriod(currentTime, phase, countPeriods);
		phaseCounterQuarter++;
	}
}

int isPeriodPartChanged(double currentTime, double periodPart) {
	double tmp = floor(currentTime / periodPart) * periodPart;
	return currentTime < (tmp + periodPart) && (currentTime + CURRENT_TIMESTEP) >= (tmp + periodPart) ? 1 : 0;
}

void allocateArrayMemory() {
	x_array = malloc(sizeof(double) * ELEMENTS_COUNT);
	y_array = malloc(sizeof(double) * ELEMENTS_COUNT);
	u = malloc(sizeof(double) * ELEMENTS_COUNT);
	v = malloc(sizeof(double) * ELEMENTS_COUNT);
	pressure = malloc(sizeof(double) * ELEMENTS_COUNT);
	temp = malloc(sizeof(double) * ELEMENTS_COUNT);
	temp_acc = malloc(sizeof(double) * ELEMENTS_COUNT);
	u_acc = malloc(sizeof(double) * ELEMENTS_COUNT);
	v_acc = malloc(sizeof(double) * ELEMENTS_COUNT);
	density = malloc(sizeof(double) * ELEMENTS_COUNT);
	initArrays();
}

void initArrays() {
	int i;
	for(i = 0; i < ELEMENTS_COUNT; i++) {
		u[i] = 0.0;
		u_acc[i] = 0.0;
		v[i] = 0.0;
		v_acc[i] = 0.0;
		temp[i] = 0.0;
		temp_acc[i] = 0.0;
		pressure[i] = 0.0;
		density[i] = 0.0;
		
	}
}

void initMinAndMaxCharacteristics() {
	pMaxPiston = -1000000000.0;
	pMinPiston = 1000000000.0;
	tempAvgMax = -1000000000.0;
	velocityMax = -1000000000.0;
}

void initCoordinates() {
	int i = 0;
	Domain *domain;
	Thread *c_thread;
	cell_t cell;
	domain=Get_Domain(1);
	thread_loop_c(c_thread, domain) {
		begin_c_loop(cell, c_thread){
			double coordinates[ND_ND];
			C_CENTROID(coordinates, cell, c_thread);
			x_array[i] = coordinates[0];
			y_array[i] = coordinates[1];
			i++;
		}
		end_c_loop(cell, c_thread);
	}
}

/*
Обработка четвертинок периода (1/4Pi) колебания поршня.
Здесь сохраняются значения мгоновенной скорости газа, давления
и переход к обработке целого периода, если на него попали.
*/
void handleQuarterOfPeriod(double currentTime, int phase, double countPeriods) {
	char filename[20];
	sprintf(filename, "%d_4T.txt", phase);
	saveHorizontalToFile(filename);
	saveAtHalfLength(countPeriods, phase);
	if(phase == 0) { // т.е. целый период
		handlePeriod(currentTime, countPeriods);
	}
}
/*
Обработка полного периода (2Pi) колебания поршня. 
Здесь вычисляются и записываются значения скорости акустического течения, 
среднее и амлитуда размаха колебания давления газа на поршне и т.п. 
После обработки характеристик течения, в конце периода скидываем значения
*/
void handlePeriod(double currentTime, double countPeriods) {
	saveAcouticAll();
	saveAcousticAt34Length(currentTime, countPeriods);
	savePressureAmplitude(countPeriods);
	//скидываем
	initMinAndMaxCharacteristics();
	clearAccumulatedValues();
}
/*
Сохранение значений акустической скорости по всему резонатору
*/
void saveAcouticAll() {
	FILE *x_file = fopen("x.txt", "w");
	FILE *y_file = fopen("y.txt", "w");
	FILE *u_file = fopen("u_acoustic.txt", "w");
	FILE *v_file = fopen("v_acoustic.txt", "w");
	FILE *temp_acc_file = fopen("temp_avg_period.txt", "w");

	int i;
	for(i = 0; i < ELEMENTS_COUNT; i++) {
		//Осреднение и обезразмерование для компонент скорости акустического течения
		double u_acoustic = u_acc[i];
		double v_acoustic = v_acc[i];
		fprintf(x_file, "%1.10f\n", x_array[i] / LENGTH);
		fprintf(y_file, "%1.10f\n", y_array[i] / HEIGHT);
		fprintf(u_file, "%1.10f\n", u_acoustic);
		fprintf(v_file, "%1.10f\n", v_acoustic);		
		fprintf(temp_acc_file, "%f\n", temp_acc[i]);		
	}
	fclose(x_file);		fclose(y_file);		fclose(u_file);		fclose(v_file);		fclose(temp_acc_file);
}
/**
Cохранение значений продольной акустической скорости на 3/4 длины резонатора, т.е. u_acoustic(3/4L, y)
*/
void saveAcousticAt34Length(double currentTime, double countPeriods) {
	FILE *vertical = fopen("u_acoustic(y).txt", "w");		
	int i, from, to;
	//вертикаль на 3/4 длины резонатора	
	from = 0.75 * Nx * Ny;
	to = Ny * (0.75 * Nx + 1);
	fprintf(vertical, "\n Flow time: %f, periods: %f\n", currentTime, countPeriods);
	for(i = from; i < to; i++) {
		double y = y_array[i] / HEIGHT;
		fprintf(vertical, "%f\t %f\n", y, u_acc[i]);
	}
	fclose(vertical);
}
/*
скидываем аккумулированные значения скорости
*/
void clearAccumulatedValues() {
	int i;
	for(i = 0; i < ELEMENTS_COUNT; i++) {
		u_acc[i] = 0.0;
		v_acc[i] = 0.0;
		temp_acc[i] = 0.0;
	}
}
/*
Сохраняем размах колебания - амплитуду давления газа на поршне в зав-ти от количество периодов
*/
void savePressureAmplitude(double countPeriods) {
	double pAmpl = (pMaxPiston - pMinPiston) / P0;// размах колебаний давления газа
	FILE *file = fopen("pressure_amplitude.txt", "a");
	fprintf(file, "%f\t %f\n", countPeriods, pAmpl);
	fclose(file);
}

/*
Запись характеристик в середине резонатора (т.е. x = L/2, y = [R, 0])
*/
void saveAtHalfLength(double countPeriods, int phase) {
	char filename[20];
	FILE *file;
	double y, u_vel, p;
	int i, from, to;
	
	sprintf(filename, "05L_%d_4T.txt", phase);
	file = fopen(filename, "w");
	fprintf(file, "Count period: %f\t, phase: %d\n", countPeriods, phase);
	from = 0.5 * Nx * Ny;
	to = Ny * (0.5 * Nx + 1);
	for(i = from; i < to; i++) {
		u_vel = u[i];
		y = y_array[i] / HEIGHT;
		p = (pressure[i + Ny] - pressure[i - Ny]) / P0 / (2.0 * dx);// dp/dx
		fprintf(file, "%f\t %f\t %f\t %f\t %f\n", y, u_vel, p, temp[i], temp_acc[i]);
	}
	fclose(file);
}

void saveHorizontalToFile(char *filename) {
	FILE *file = fopen(filename, "w");
	double p, x, u_vel, u_acoustic;
	int i, k;
	for (i = 0; i < ELEMENTS_COUNT; i++) {
		k = (i + 1) % Ny;
		if (k == 0) {
			x = x_array[i] / LENGTH;
			u_vel = u[i];
			u_acoustic = u_acc[i];
			fprintf(file, "%f\t %f\t %f\t %f\t %f\t %f\n", x, u_vel, u_acoustic, pressure[i], temp[i], temp_acc[i]);
		}
	}
	fclose(file);
}
/*
Обработка значений характеристик на поршне и средней точке трубы
*/
void handlePistonCharacteristics(double countPeriods) {
	FILE *p_piston = fopen("on_piston.txt","a");
	FILE *u_file = fopen("middle_point.txt","a");
	double u_vel, p_left, p_right;
	int i, k;
	//средняя точка на поршне
	i = PISTON_INDEX;//левый конец
	p_left = pressure[i];
	k = Nx * Ny - 1;//правый конец
	p_right = pressure[k];
	fprintf(p_piston, "%f\t %f\t %f\t %f\t %f\t %f\n", countPeriods, u[i], p_right, p_left, temp[i], temp_acc[i]);
	fclose(p_piston);
	
	//точка в середине резонатора
	i = Ny * Nx / 2 + Ny - 1;//axisymmetric
	fprintf(u_file, "%f\t %f\t %f\t %f\t %f\n", countPeriods, u[i], pressure[i], temp[i], temp_acc[i]);
	fclose(u_file);
}
/**
Дополнительный инерционный член (A*w^2*cos(w*t)*rho) в уравнении сохр. импульса по Z
*/
DEFINE_SOURCE(x_mom_source, c, t, dS, eqn)
 {
	double time = RP_Get_Real("flow-time");
	double source = ampX0 * SQR(OMEGA) * cos(OMEGA * time) * C_R(c, t);
    dS[eqn] = 0.0;
	return source;
 } 