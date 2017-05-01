# include "Header.h"  
using namespace std;

int main()
{

	
	double  EndTime = 0.1; 
		double dt = 0.00001; 
		double Lx = 1.0; 
	    double Ly = 1.0; 
		double miu = 1.5; 
		double PrdRate = 100; 
		const int	nx = 51;
		const int	ny = 51;
		const int   nz = 2;
		double growRate = 4;
		double WingCntrX = 0.5;
		double WingCntrY = 0.5;
		double SrcWidthRatio = 0.2;
		

	
	double cold[nx+1][ny+1];
	double c[nx+1][ny+1];
	double cSd[nx + 1][ny + 1];
	double Prd[nx + 1][ny + 1];
	double x[nx+1]; 
	double y[ny+1]; 
	double dx, dy; 
	double RWing ;
	double distX; 
	
	

	

	dx = Lx / (nx - 1); 
	dy = Ly / (ny - 1);

	for (int i = 0; i <= nx+1; i++) {
		x[i] = 0.5*dx + (i - 1)*dx; 
	}

	for (int j = 0;  j <= ny; j++) {
		y[j] = 0.5*dy + (j - 1)*dy; 
	}

	int CountEnd= EndTime / dt; 

	//initialize
	for (int i = 0; i <= nx; i++) {
		for (int j = 0; j <= ny; j++) {

			cold[i][j] = 0;
			c[i][j] = 0;
			Prd[i][j] = 0;

		}

	}


	//initialize source at the center
	int DPPcntr = floor((nx - 1) / 2);
	for (int j = 0; j <= ny; j++) {

	//	cold[DPPcntr][j] = 1.0;
	//	c[DPPcntr][j] = cold[DPPcntr][j];
	}


	RWing = 0.05;

	int plotCounter; 
		plotCounter = 0; 
	//main loop
	for (int t = 0; t <= CountEnd; t++) {


		RWing = RWing + growRate*dt;

		

		plotCounter = plotCounter + 1;

		 // update the values
		for (int i = 0; i <= nx; i++) {
			for (int j = 0; j <= ny; j++) {

				cold[i][j] = c[i][j];
			}

		}

		cout << "Time is" << t *dt << endl;

	
		

		double dist;
		for (int i = 0; i <= nx; i++) {
			for (int j = 0; j <= ny; j++) {

				dist = sqrt((x[i] - WingCntrX)*(x[i] - WingCntrX) + (y[j] - WingCntrY)*(y[j] - WingCntrY));

				if (dist < RWing) {

					cSd[i][j] = 1;
				}
				else {
					cSd[i][j] = 0;
				}

			}

		}
		
		for (int i = 0; i <= nx; i++) {
			for (int j = 0; j <= ny; j++) {

				if (cSd[i][j] == 1) {

					distX = sqrt((x[i] - WingCntrX)*(x[i] - WingCntrX));
						//	cold[DPPcntr][j] = 1.0;
						if (distX / RWing < SrcWidthRatio) {

							Prd[i][j] = PrdRate;
						}
						else {
							Prd[i][j] = 0;
						}
				}


					else {
						Prd[i][j] = 0;
					}
				}
			}
		


		for (int i = 1; i <= nx-1; i++) {
			for (int j = 1; j <= ny-1; j++) {

				c[i][j] = cold[i][j] + dt*Prd[i][j]+ dt*miu*(  (cold[i + 1][j] - 2 * cold[i][j] + cold[i - 1][j]) / (dx*dx) +
					                                           (cold[i][j + 1] - 2 * cold[i][j] + cold[i][j - 1]) / (dy*dy) +

					                      );

			}


		}

		// treating boundary cells 
		for (int i = 1; i <= nx - 1; i++) {
			c[i][0] = c[i][1]; 
			c[i][ny] = c[i][ny-1]; 
		}

		for (int j = 0; j <= ny; j++) {
			c[0][j] = c[1][j];
			c[nx][j] = c[nx-1][j];
		}

		//z is just for output purpose
		double z[2]; 
		z[1] = 1;

		

		if (plotCounter == 50) {




			std::string vtkFileName = "DPP_" + std::to_string(t) + ".vtk";
			ofstream SignalOut;
			SignalOut.open(vtkFileName.c_str());
			SignalOut << "# vtk DataFile Version 2.0" << endl;
			SignalOut << "Result for paraview 2d code" << endl;
			SignalOut << "ASCII" << endl;
			SignalOut << "DATASET RECTILINEAR_GRID" << endl;
			SignalOut << "DIMENSIONS" << " " << nx - 1 << " " << " " << ny - 1 << " " << nz - 1 << endl;








			SignalOut << "X_COORDINATES " << nx - 1 << " float" << endl;
			//write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
			for (int i = 1; i <= nx - 1; i++) {
				SignalOut << x[i] << endl;
			}

			SignalOut << "Y_COORDINATES " << ny - 1 << " float" << endl;
			//write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
			for (int j = 1; j <= ny - 1; j++) {
				SignalOut << y[j] << endl;
			}

			SignalOut << "Z_COORDINATES " << nz - 1 << " float" << endl;
			//write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
			for (int k = 1; k <= nz - 1; k++) {
				SignalOut << z[k] << endl;
			}

			SignalOut << "POINT_DATA " << (nx - 1)*(ny - 1)*(nz - 1) << endl;
			SignalOut << "SCALARS DPP float 1" << endl;
			SignalOut << "LOOKUP_TABLE default" << endl;

			for (int k = 1; k <= nz - 1; k++) {
				for (int j = 1; j <= ny - 1; j++) {
					for (int i = 1; i <= nx - 1; i++) {
						SignalOut << c[i][j]* cSd[i][j] << endl;

					}
				}
			}

			SignalOut << "SCALARS Wing float 1" << endl;
			SignalOut << "LOOKUP_TABLE default" << endl;

			for (int k = 1; k <= nz - 1; k++) {
				for (int j = 1; j <= ny - 1; j++) {
					for (int i = 1; i <= nx - 1; i++) {
						SignalOut << cSd[i][j] << endl;

					}
				}
			}

			plotCounter = 0; 
		}

		}
    

	

	int b;

	cin >> b;

return 0; 
}