#define Section 12
#define CoolSection 8
#define MoldSection 4
#include "time.h"
#include "stdio.h"
#include "math.h"
void Physicial_Parameters(float T);
float Boundary_Condition(int j, int ny, float Ly, float *, float *);
void OutputTemperature(int nx, int ny, int nz, int tnpts, int inter_num);
float T[21][3001][21][200] = { 0.0 };
float T_Last[21][3001][21];
float T_New[21][3001][21];
float Vcast = -0.02, h = 100.0, lamd = 50.0, Ce = 540.0, pho = 7000.0, a;

int main()
{
	FILE *fp = NULL;
	float H_Init[Section] = { 1380.0, 1170.0, 980.0, 800.0, 1223.16, 735.05, 424.32, 392.83, 328.94, 281.64, 246.16, 160.96 };
	float ccml[Section + 1] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0925, 2.27, 4.29, 5.831, 9.6065, 13.6090, 19.87014, 28.599 };
	float T_SumZone[CoolSection][3001] = {};
	float y[3001];
	int nx = 21, ny = 3001, nz = 21, tnpts = 10001;
	float Lx = 0.25, Ly = 28.599, Lz = 0.25, dx, dy, dz, t_final = 2000.0, tao, T_Cast = 1558, Tw = 30.0, T_Up, T_Down, T_Right, T_Left, T_Forw, T_Back;
	int inter_num = 500, count = 0;
	clock_t begin, duration;

	begin = clock();
	dx = Lx / (nx - 1);
	dy = Ly / (ny - 1);
	dz = Lz / (nz - 1);
	tao = t_final / (tnpts - 1);

	for (int j = 0; j < ny; j++)
		y[j] = j*dy;

	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			for (int m = 0; m < nz; m++)
				    T_Last[i][j][m] = T_Cast;

	printf("Casting Temperature = %f ", T_Cast);
	printf("\n");
	printf("The length of steel billets(m) = %f ", Ly);
	printf("\n");
	printf("The width of steel billets(m) = %f ", Lz);
	printf("\n");
	printf("The thick of steel billets(m) = %f ", Lx);
	printf("\n");
	printf("dx(m) = %f ", dx);
	printf("dy(m) = %f ", dy);
	printf("dz(m) = %f ", dz);
	printf("tao(s) = %f ", tao);
	printf("\n");
	printf("simulation time(s) = %f ", t_final);

	for (int k = 0; k < tnpts-1; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			h = Boundary_Condition(j, ny, Ly, ccml, H_Init);
			for (int i = 0; i < nx; i++)
				for (int m = 0; m < nz; m++)
				{
					Physicial_Parameters(T_Last[i][j][m]);
					a = lamd / (pho*Ce);
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
						T_New[i][j][m]= T_Cast;
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
						T_Forw = T_Last[i][j][m - 1] - 2 * dz*h*(T_Last[i][j][m] - Tw) / lamd;
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
						T_Up = T_Last[i - 1][j][m] - 2 * dx * h * (T_Last[i][j][m] - Tw) / lamd;
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

		if (k % inter_num == 0)
		{
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					for (int m = 0; m < nz; m++)
						T[i][j][m][count] = T_Last[i][j][m];
			count++;
			printf("\n");
			printf("\n");
			printf(" Output time step = %d", k);
			printf(" Simulation time = %f", k*tao);

			int Y_Label[Section + 1] = {0};
			float temp_XZ[Section] = {0.0}, temp_XY[Section] = { 0.0 }, temp_central[Section] = { 0.0 };
			Y_Label[Section] = ny-1;
			for (int j = 0; j < ny - 1; j++)
			{
				for (int i = 1; i < Section; i++)
				{
					if (y[j] <= ccml[i] && y[j + 1] >= ccml[i])
					{
						Y_Label[i] = j + 1;
					}
				}
			}	
			//for (int i = 0; i < Section+1; i++)
				//printf("Label %d = %d  ", i + 1, Y_Label[i]);

			for (int j = 0; j < ny; j++)
				for (int i = 0; i < Section; i++)
					if (j > Y_Label[i] && j <= Y_Label[i + 1])
					{
						temp_XZ[i] = temp_XZ[i] + T_Last[nx - 1][j][0];
						temp_XY[i] = temp_XY[i] + T_Last[5][j][nz-1];
						temp_central[i] = temp_central[i] + T_Last[0][j][0];
					}

			printf("\n Surface temperature\n");
			float Mean_Temperature_XZ[Section] = {}, Mean_Temperature_XY[Section] = {}, Mean_Temperature_central[Section] = {};
			for (int i = 0; i < Section; i++)
			{
				Mean_Temperature_XZ[i] = temp_XZ[i] / (Y_Label[i + 1] - Y_Label[i]);
				//Mean_Temperature_central[i] = temp_central[i] / (Y_Label[i + 1] - Y_Label[i]);
				//Mean_Temperature_XY[i] = temp_XY[i] / (Y_Label[i + 1] - Y_Label[i]);
				printf("zone %d = %f  ", i + 1, Mean_Temperature_XZ[i]);
				//printf("zone %d = %f  ", i + 1, Mean_Temperature_XY[i]);
				//printf("Label %d = %d  ", i + 1, Y_Label[i]);
			}

			printf("\n Central temperature\n");
			for (int i = 0; i < Section; i++)
			{
				Mean_Temperature_central[i] = temp_central[i] / (Y_Label[i + 1] - Y_Label[i]);
				printf("zone %d = %f  ", i + 1, Mean_Temperature_central[i]);
			}

		}

		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				for (int m = 0; m < nz; m++)
					T_Last[i][j][m] = T_New[i][j][m];
	}

	duration = fabs(begin - clock());
	printf("\n");
	printf("running time = %ld\n", duration / CLK_TCK);
	printf("\n");
	printf("count = %d", count);
	printf("\n");
	printf("Please wait");

	OutputTemperature(nx, ny, nz, tnpts, inter_num);
}


void Physicial_Parameters(float T)
{
	float Ts = 1462.0, Tl = 1518.0, lamds = 30, lamdl = 50, phos = 7000, phol = 7500, ce = 540.0, L = 265600.0, fs = 0.0;
	if (T<Ts)
	{
		fs = 0;
		pho = phos;
		lamd = lamds;
		Ce = ce;
	}

	if (T >= Ts && T <= Tl)
	{
		fs = (T - Ts) / (Tl - Ts);
		pho = fs*phos + (1 - fs)*phol;
		lamd = fs*lamds + (1 - fs)*lamdl;
		Ce = ce + L / (Tl - Ts);
	}

	if (T>Tl)
	{
		fs = 1;
		pho = phol;
		lamd = lamdl;
		Ce = ce;
	}
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

void OutputTemperature(int nx, int ny, int nz, int tnpts, int inter_num)
{
	FILE*fp = NULL;
	int i, j, k, m;

	fp = fopen("C:\\Temperature3D_Static.txt", "w");
	k = (tnpts - 1)/inter_num;
	printf("\n");
	printf("Output time step is %d",k);
	for (j = 0; j < ny; j++)
	{
		for (i = 0; i < nx; i++)
		{
			for (m = 0; m < nz; m++)
			{
				fprintf(fp, " %f", T[i][j][m][k]);
				fprintf(fp, " %d, %d, %d,", i, j, m);
			} 
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	/*fp = fopen("C:\\Temperature3Dxyn.txt", "w");
	m = nz - 1;
	for (k = 0; k < tnpts - 1; k++)
	{
		for (j = 0; j < ny; j++)
		{
			for (i = 0; i < nx; i++)
			{
				fprintf(fp, " %f", T[i][j][m][k]);
				fprintf(fp, " %d, %d, %d,", i, j, m);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("C:\\Temperature3Dxy1.txt", "w");
	m = 0;
	for (k = 0; k < tnpts - 1; k++)
	{
		for (j = 0; j < ny; j++)
		{
			for (i = 0; i < nx; i++)
			{
				fprintf(fp, " %f", T[i][j][m][k]);
				fprintf(fp, " %d, %d, %d,", i, j, m);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("C:\\Temperature3Dyzn.txt", "w");
	i = nx - 1;
	for (k = 0; k < tnpts - 1; k++)
	{
		for (j = 0; j < ny; j++)
		{
			for (m = 0; m < nz; m++)
			{
				fprintf(fp, " %f", T[i][j][m][k]);
				fprintf(fp, " %d, %d, %d,", i, j, m);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("C:\\Temperature3Dyz1.txt", "w");
	i = 0;
	for (k = 0; k < tnpts - 1; k++)
	{
		for (j = 0; j < ny; j++)
		{
			for (m = 0; m < nz; m++)
			{
				fprintf(fp, " %f", T[i][j][m][k]);
				fprintf(fp, " %d, %d, %d,", i, j, m);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("C:\\Temperature3Dxzn.txt", "w");
	j = ny - 1;
	for (k = 0; k < tnpts - 1; k++)
	{
		for (i = 0; i < nx; i++)
		{
			for (m = 0; m < nz; m++)
			{
				fprintf(fp, " %f", T[i][j][m][k]);
				fprintf(fp, " %d, %d, %d,", i, j, m);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("C:\\Temperature3Dxz1.txt", "w");
	j = 0;
	for (k = 0; k < tnpts - 1; k++)
	{
		for (i = 0; i < nx; i++)
		{
			for (m = 0; m < nz; m++)
			{
				fprintf(fp, " %f", T[i][j][m][k]);
				fprintf(fp, " %d, %d, %d,", i, j, m);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);*/
}