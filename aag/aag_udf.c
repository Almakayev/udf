#include "udf.h"
#include "mem.h"
#include "stdio.h"
#include <math.h>

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)
#define PRESSURE_DL(p) (p / GAMMA/ P0)
#define TEMPER_DL(temp) (temp / T0 - 1.0)

const int Nx = 1002,  Ny = 22;
const double OMEGA_FACTOR = 0.14;
const double ampX0 = 0.01;// piston amplitude = 2 * L
double U0, OMEGA, PERIOD, SOUND_SPEED, DENSITY_0, U_R;
const double P0 = 109511.7; // atmospheric pressure
const double LENGTH = 0.005;// resonator length
const double HEIGHT = 0.0001; //height
const double T0 = 300.0; // initial temperature
const double R = 286.8; 
const double GAMMA = 1.4;
const double VISCOSITY = 1.935e-05; // dynamic viscosity of the gas
int phaseCounterQuarter, ELEMENTS_COUNT, PISTON_INDEX;
double *x_array, *y_array, *u, *v, *pressure, *temp, *u_acc, *v_acc, *density, *temp_acc;
double dx, dy, pMaxPiston, pMinPiston, velocityMax, tempAvgMax;

/**
Initializes the variables for computation and save some meta data to the file. 
The procedure should be called from GUI manually.
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
	PERIOD = 2.0 * M_PI / OMEGA;
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
	ELEMENTS_COUNT = Nx * Ny; //total elements count
	fprintf(file, "Nx: %d\t  Ny: %d\t Total elements: %d\n", Nx, Ny, ELEMENTS_COUNT);
	fprintf(file, "ampX0 / dx: %f\n", ampX0 / dx);
	
	allocateArrayMemory();
	initCoordinates();
	PISTON_INDEX = Ny - 1;//hardcoded index of the middle of the piston
	fprintf(file, "Piston index = %d\t x_piston: %f\t  y_piston: %f\t\n", PISTON_INDEX, x_array[PISTON_INDEX], y_array[PISTON_INDEX]);
	initMinAndMaxCharacteristics();
	fclose(file);
	printf("My UDF for initialization was executed...\n");
	fflush(stdout);
	
	phaseCounterQuarter = 1;
	pMaxPiston = -1e+09;
	pMinPiston = 1e+09;
	velocityMax = -1e+09;
	tempAvgMax = -1e+09;
}

/*
Handles all computation characteristics and consumes them after each "piece of the period" (see variable isChanged).   
*/
DEFINE_EXECUTE_AT_END(execute_at_end) 
{	
	int isChanged, i;
	double tmp, periodPart, tau, countPeriods;
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
	tau = currentTime * SOUND_SPEED / LENGTH;
	
	pMaxPiston = MAX(pMaxPiston, pressure[PISTON_INDEX]);
	pMinPiston = MIN(pMinPiston, pressure[PISTON_INDEX]);
	printf("tau: %f\t, Tw: %f\t, pMaxPiston: %f\t, pMinPiston: %f\t, velocityMax: %f\t, tempAvgMax: %f\n", tau, countPeriods, pMaxPiston, pMinPiston, velocityMax, tempAvgMax);
	fflush(stdout);
	
	//After each 1/30 of the period
	isChanged = isPeriodPartChanged(currentTime, PERIOD / 30.0);
	if(isChanged) {
		findMaxCharacteristics();
		handleCharacteristicsAtPoints(tau);
	}
	
	//After each 1/4 of the period
	isChanged = isPeriodPartChanged(currentTime, 0.25 * PERIOD);
	if(isChanged) {
		int phase = phaseCounterQuarter % 4;
		handleQuarterOfPeriod(currentTime, phase, tau);
		phaseCounterQuarter++;
	}
}

int isPeriodPartChanged(double currentTime, double periodPart) {
	double tmp = floor(currentTime / periodPart) * periodPart;
	return currentTime < (tmp + periodPart) && (currentTime + CURRENT_TIMESTEP) >= (tmp + periodPart) ? 1 : 0;
}

/*
Handles each 1/4 of the period. 
*/
void handleQuarterOfPeriod(double currentTime, int phase, double tau) {
	char filename[20];
	sprintf(filename, "%d_4T.txt", phase);
	saveHorizontalToFile(filename);
	saveAtHalfLength(tau, phase);
	if(phase == 0) { // т.е. целый период
		handlePeriod(currentTime, tau);
	}
}
/*
Handles after the whole period 2*PI. 
After that, we need to clear some arrays, counters etc.
*/
void handlePeriod(double currentTime, double tau) {
	saveAcousticAll();
	saveAcousticAt34Length(currentTime, tau);
	savePressureAmplitude(tau);
	//скидываем
	initMinAndMaxCharacteristics();
	clearAccumulatedValues();
}
/*
Saves acoustic velocity for whole resonator.
*/
void saveAcousticAll() {
	FILE *x_file = fopen("x.txt", "w");
	FILE *y_file = fopen("y.txt", "w");
	FILE *u_file = fopen("u_acoustic.txt", "w");
	FILE *v_file = fopen("v_acoustic.txt", "w");
	FILE *temp_acc_file = fopen("temp_avg_period.txt", "w");

	int i;
	for(i = 0; i < ELEMENTS_COUNT; i++) {
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
Saves the values for 3/4 of the length of the resonator
*/
void saveAcousticAt34Length(double currentTime, double tau) {
	FILE *vertical = fopen("u_acoustic(y).txt", "w");		
	int i, from, to;	
	from = 0.75 * Nx* Ny;
	to = from + Ny;
	fprintf(vertical, "\n Flow time: %f, periods: %f\n", currentTime, tau);
	for(i = from; i < to; i++) {
		double y = y_array[i] / HEIGHT;
		fprintf(vertical, "%f\t %f\n", y, u_acc[i]);
	}
	fclose(vertical);
}

/*
Saves the current time and the amplitude of the total pressure
*/
void savePressureAmplitude(double tau) {
	double pAmpl = (pMaxPiston - pMinPiston) / P0;// размах колебаний давления газа
	FILE *file = fopen("pressure_amplitude.txt", "a");
	fprintf(file, "%f\t %f\n", tau, pAmpl);
	fclose(file);
}

/*
Saves the values for cross-section x = L/2, y = [R, 0] 
*/
void saveAtHalfLength(double tau, int phase) {
	char filename[20];
	FILE *file;
	double y, u_vel, p;
	int i, from, to;
	
	sprintf(filename, "05L_%d_4T.txt", phase);
	file = fopen(filename, "w");
	fprintf(file, "Count period: %f\t, phase: %d\n", tau, phase);
	from = 0.5 * Nx * Ny;
	to = from + Ny;
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
	int i;
	for (i = Ny - 1; i < ELEMENTS_COUNT; i+= Ny) {
			x = x_array[i] / LENGTH;
			u_vel = u[i];
			u_acoustic = u_acc[i];
			fprintf(file, "%f\t %f\t %f\t %f\t %f\t %f\t %f\n", x, u_vel, u_acoustic, pressure[i], density[i], temp[i], temp_acc[i]);
	}
	fclose(file);
}
/*
Handles the values at the 3 different points.
*/
void handleCharacteristicsAtPoints(double tau) {
	//middle point of the left boundary
	int i = PISTON_INDEX;
	saveValues("0L.txt", i, tau);
	
	//point 0.25*L
	i = Ny * Nx / 4 + Ny - 1;
	saveValues("025L.txt", i, tau);
	
	//point 0.5*L
	i = Ny * Nx / 2 + Ny - 1;//axisymmetric
	saveValues("05L.txt", i, tau);
}

void saveValues(char *filename, int i, double time) {
	FILE *file = fopen(filename, "a");
	fprintf(file, "%f\t %f\t %f\t %f\t %f\t %f\t %f\n", x_array[i], y_array[i], time, u[i], pressure[i], temp[i], temp_acc[i]);
	fclose(file);
}

/**
The additinal source term for the x-momentum equation.
The unit should be N/m^3 or kg/(m^2*s^2)
*/
DEFINE_SOURCE(x_mom_source, c, t, dS, eqn)
{
	double time = RP_Get_Real("flow-time");
	double source = ampX0 * SQR(OMEGA) * cos(OMEGA * time) * C_R(c, t);
    dS[eqn] = 0.0;
	return source;
}

/**
Computes and returns the acceleration of the force for the discrete phase.
The unit should be m/s^2
 */ 
DEFINE_DPM_BODY_FORCE(particle_body_force, tp, cartesian_comp_index)
{
    if(cartesian_comp_index == 0) {//for x-momentum 
		double time = RP_Get_Real("flow-time");
		double acceleration = ampX0 * SQR(OMEGA) * cos(OMEGA * time);
		return acceleration;
    }
    else return 0.0;
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
	pMaxPiston = -1e+09;
	pMinPiston = 1e+09;
	tempAvgMax = -1e+09;
	velocityMax = -1e+09;
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

void clearAccumulatedValues() {
	int i;
	for(i = 0; i < ELEMENTS_COUNT; i++) {
		u_acc[i] = 0.0;
		v_acc[i] = 0.0;
		temp_acc[i] = 0.0;
	}
}