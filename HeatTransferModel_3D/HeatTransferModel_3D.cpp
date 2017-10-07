#define Section 12
#define CoolSection 8
#define MoldSection 4
#include <ctime>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

struct Meantemperature 
{
	float *Mean_TSurfacecore;
	float *Mean_TSurfacemiddle;
	float *Mean_TCentral;
};

struct physicalpara
{
	float lamd;
	float Ce;
	float pho;
};

physicalpara Physicial_Parameters(float T);
float Boundary_Condition(int j, int ny, float Ly, float *, float *);
Meantemperature Calculation_MeanTemperature(int nx, int ny, int nz, float dy, float *ccml, float ***T);
void differencecalculation(float ***T_New, float ***T_Last, int nx, int ny, int nz, float Ly, float dx, float dy, float dz, float tao, float *ccml, float T_Cast, float *H_Init, float Vcast, int pstep, bool & disout);
void initcondition(float***T_Last, int nx, int ny, int nz, float T_Cast, bool & disout);

int main()
{
	float ***T_New, ***T_Last;
	float H_Init[Section] = { 1380.0f, 1170.0f, 980.0f, 800.0f, 1223.16f, 735.05f, 424.32f, 392.83f, 328.94f, 281.64f, 246.16f, 160.96f };
	float ccml[Section + 1] = { 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0925f, 2.27f, 4.29f, 5.831f, 9.6065f, 13.6090f, 19.87014f, 28.599f };
	const int nx = 31, ny = 3000, nz = 31, tnpts = 10001;
	float Lx = 0.25f, Ly = 28.599f, Lz = 0.25f, dx, dy, dz, t_final = 2000.0f, tao, T_Cast = 1558.0f, Vcast = -0.02f;
	int inter_num = 100, count = 0;
	bool disout = true;
	clock_t begin, duration;

	begin = clock();
	dx = Lx / (nx - 1);
	dy = Ly / (ny - 1);
	dz = Lz / (nz - 1);
	tao = t_final / (tnpts - 1);

	T_New = new float**[nx];
	for (int i = 0; i < nx; i++)
	{
		T_New[i] = new float*[ny];
		for (int j = 0; j < ny; j++)
			T_New[i][j] = new float[nz];
	}
	T_Last = new float**[nx];
	for (int i = 0; i < nx; i++)
	{
		T_Last[i] = new float*[ny];
		for (int j = 0; j < ny; j++)
			T_Last[i][j] = new float[nz];
	}


	cout << "Casting Temperature " << T_Cast << endl;
	cout << "The length of steel billets(m) " << Ly << endl;
	cout << "The width of steel billets(m) " << Lz << endl;
	cout << "The thick of steel billets(m) " << Lx << endl;
	cout << "dx(m) " << dx << ", ";
	cout << "dy(m) " << dy << ", ";
	cout << "dz(m) " << dz << ", ";
	cout << "tao(s) " << tao << ", ";
	cout << "simulation time(s) " << t_final << endl;
	
	int pstep = 1;
	initcondition(T_Last, nx, ny, nz, T_Cast, disout);
	for (int k = 0; k < tnpts - 1; k++)
	{
		differencecalculation(T_New, T_Last, nx, ny, nz, Ly, dx, dy, dz, tao, ccml, T_Cast, H_Init, Vcast, pstep, disout);
		if (k % inter_num == 0) 
		{
			count++;
			cout << endl;
			cout << "Output time step = " << k <<", ";
			cout << "Simulation time = " << k * tao << endl;
			Meantemperature meantemperature = Calculation_MeanTemperature(nx, ny, nz, dy, ccml, T_New);
			for (int temp = 0; temp < CoolSection; temp++)
				cout << meantemperature.Mean_TSurfacecore[temp + MoldSection]<<", ";
		}
	}

	duration = fabs(begin - clock());
	cout << endl;
	cout << "running time = " << duration / CLK_TCK << endl;
	cout << "count = " << count << endl;

	delete T_New;
	delete T_Last;
	for (int i = 0; i < nx; i++)
	{
		delete T_New[i];
		delete T_Last[i];
		for (int j = 0; j < ny; j++)
		{
			delete T_New[i][j];
			delete T_Last[i][j];
		}
	}

}


physicalpara Physicial_Parameters(float T)
{
	physicalpara steel;
	float Ts = 1462.0, Tl = 1518.0, lamds = 30, lamdl = 50, phos = 7000, phol = 7500, ce = 540.0, L = 265600.0, fs = 0.0;
	if (T<Ts)
	{
		fs = 0;
		steel.pho = phos;
		steel.lamd = lamds;
		steel.Ce = ce;
	}

	if (T >= Ts && T <= Tl)
	{
		fs = (T - Ts) / (Tl - Ts);
		steel.pho = fs*phos + (1 - fs)*phol;
		steel.lamd = fs*lamds + (1 - fs)*lamdl;
		steel.Ce = ce + L / (Tl - Ts);
	}

	if (T>Tl)
	{
		fs = 1;
		steel.pho = phol;
		steel.lamd = lamdl;
		steel.Ce = ce;
	}
	return steel;
}

float Boundary_Condition(int j, int ny, float Ly, float *ccml_zone, float *H_Init)
{
	float YLabel, h;
	YLabel = (j*Ly) / float(ny - 1);

	for (int i = 0; i < Section; i++)
	{
		if (YLabel >= *(ccml_zone + i) && YLabel <= *(ccml_zone + i + 1))
		{
			h = *(H_Init + i);
		}
	}
	return h;
}

Meantemperature Calculation_MeanTemperature(int nx, int ny, int nz, float dy, float *ccml, float ***T)
{
	Meantemperature meantemperature;
	meantemperature.Mean_TSurfacecore = new float[Section];
	meantemperature.Mean_TSurfacemiddle = new float[Section];
	meantemperature.Mean_TCentral = new float[Section];
	float y;
	int count = 0;
	for (int i = 0; i < Section; i++)
	{
		meantemperature.Mean_TSurfacecore[i] = 0.0;
		for (int j = 0; j < ny; j++)
		{
			y = j * dy;
			if (y > *(ccml + i) && y <= *(ccml + i + 1))
			{
				meantemperature.Mean_TSurfacecore[i] = meantemperature.Mean_TSurfacecore[i] + T[nx - 1][j][0];
				meantemperature.Mean_TSurfacemiddle[i] = meantemperature.Mean_TSurfacemiddle[i] + T[nx - 1][j][int((nz - 1) / 2)];
				meantemperature.Mean_TCentral[i] = meantemperature.Mean_TCentral[i] + T[int((nx - 1) / 2)][j][int((nz - 1) / 2)];
				count++;
			}
		}
		meantemperature.Mean_TSurfacecore[i] = meantemperature.Mean_TSurfacecore[i] / float(count);
		meantemperature.Mean_TSurfacemiddle[i] = meantemperature.Mean_TSurfacemiddle[i] / float(count);
		meantemperature.Mean_TCentral[i] = meantemperature.Mean_TCentral[i] / float(count);
		count = 0;
	}

	return meantemperature;
}

void differencecalculation(float ***T_New, float ***T_Last, int nx, int ny, int nz, float Ly, float dx, float dy, float dz, float tao, float *ccml, float T_Cast, float *H_Init, float Vcast, int pstep, bool & disout)
{
	float T_Up, T_Down, T_Right, T_Left, T_Forw, T_Back, a, Tw = 30.0f, h;
	if(disout)
	{
		for (int j = 0; j < ny; j++)
		{
			h = Boundary_Condition(j, ny, Ly, ccml, H_Init);
			for (int i = 0; i < nx; i++)
				for (int m = 0; m < nz; m++)
				{
					physicalpara steel = Physicial_Parameters(T_Last[i][j][m]);
					a = steel.lamd / (steel.pho * steel.Ce);
					if (j == 0 && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //1
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == 0 && i == 0 && m != 0 && m != (nz - 1)) //2
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == 0 && i == (nx - 1) && m != 0 && m != (nz - 1))//3
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == 0 && i != 0 && i != (nx - 1) && m == 0) //4
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == 0 && i != 0 && i != (nx - 1) && m == (nz - 1)) //5
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == 0 && i == 0 && m == 0)  //6
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == 0 && i == 0 && m == (nz - 1))  //7
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == 0 && i == (nx - 1) && m == 0)  //8
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == 0 && i == (nx - 1) && m == (nz - 1)) //9
					{
						T_New[i][j][m] = T_Cast;
					}

					else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //10
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == 0 && m != 0 && m != (nz - 1)) //11
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i + 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //12
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i - 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //13
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m + 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //14
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m - 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == 0 && m == 0)  //15
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i + 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m + 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == 0 && m == (nz - 1))  //16
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i + 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m - 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == (nx - 1) && m == 0)  //17
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i - 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m + 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == (nx - 1) && m == (nz - 1))  //18
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i - 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j - 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m - 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //19
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m + 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //20
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m - 1] - 2 * dz*h*(T_Last[i][j][m] - Tw) / steel.lamd;
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == 0 && m == 0) //21
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i + 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m + 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == 0)  //22
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i - 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m + 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == 0 && m == (nz - 1)) //23
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i + 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m - 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == (nz - 1)) //24
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i - 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m - 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == 0 && m != 0 && m != (nz - 1))  //25
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i + 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //26
					{
						//T[i][j][m][k + 1] = T_Cast;
						T_Up = T_Last[i - 1][j][m] - 2 * dx * h * (T_Last[i][j][m] - Tw) / steel.lamd;
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else  //27
					{
						T_Up = T_Last[i + 1][j][m];
						T_Down = T_Last[i - 1][j][m];
						T_Right = T_Last[i][j + 1][m];
						T_Left = T_Last[i][j - 1][m];
						T_Forw = T_Last[i][j][m + 1];
						T_Back = T_Last[i][j][m - 1];
						T_New[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_Last[i][j][m]
							+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}
				}
		}
	}

	else
		{
			for (int j = 0; j < ny; j++)
			{
				h = Boundary_Condition(j, ny, Ly, ccml, H_Init);
				for (int i = 0; i < nx; i++)
					for (int m = 0; m < nz; m++)
					{
						physicalpara steel = Physicial_Parameters(T_New[i][j][m]);
						a = steel.lamd / (steel.pho * steel.Ce);
						if (j == 0 && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //1
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == 0 && i == 0 && m != 0 && m != (nz - 1)) //2
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == 0 && i == (nx - 1) && m != 0 && m != (nz - 1))//3
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == 0 && i != 0 && i != (nx - 1) && m == 0) //4
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == 0 && i != 0 && i != (nx - 1) && m == (nz - 1)) //5
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == 0 && i == 0 && m == 0)  //6
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == 0 && i == 0 && m == (nz - 1))  //7
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == 0 && i == (nx - 1) && m == 0)  //8
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == 0 && i == (nx - 1) && m == (nz - 1)) //9
						{
							T_Last[i][j][m] = T_Cast;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //10
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m != 0 && m != (nz - 1)) //11
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i + 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //12
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i - 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //13
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m + 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //14
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m - 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m == 0)  //15
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i + 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m + 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m == (nz - 1))  //16
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i + 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m - 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m == 0)  //17
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i - 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m + 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m == (nz - 1))  //18
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i - 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j - 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m - 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //19
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m + 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //20
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m - 1] - 2 * dz*h*(T_New[i][j][m] - Tw) / steel.lamd;
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m == 0) //21
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i + 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m + 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == 0)  //22
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i - 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m + 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m == (nz - 1)) //23
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i + 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m - 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == (nz - 1)) //24
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i - 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m - 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m != 0 && m != (nz - 1))  //25
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i + 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //26
						{
							//T[i][j][m][k + 1] = T_Cast;
							T_Up = T_New[i - 1][j][m] - 2 * dx * h * (T_New[i][j][m] - Tw) / steel.lamd;
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else  //27
						{
							T_Up = T_New[i + 1][j][m];
							T_Down = T_New[i - 1][j][m];
							T_Right = T_New[i][j + 1][m];
							T_Left = T_New[i][j - 1][m];
							T_Forw = T_New[i][j][m + 1];
							T_Back = T_New[i][j][m - 1];
							T_Last[i][j][m] = a*(tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*Vcast / dy))*T_New[i][j][m]
								+ a*(tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*Vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}
					}
			}
		}
		disout = !disout;
}

void initcondition(float***T_Last, int nx, int ny, int nz, float T_Cast, bool & disout)
{
	disout = true;
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			for (int m = 0; m < nz; m++)
				T_Last[i][j][m] = T_Cast;
}