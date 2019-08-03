#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#define G 100  // gravitational constant
#define N 6  // number of bodies
#include <string.h>
#include "suzaku.h" // Suzaku routines



int main(int argc, char *argv[] ) {
	double dt;
	int i, j, a, T, t;
	double fx, fy, vx, vy, x, y, dx, dy;
	double x_diff, y_diff, r, F;
	double fxsum[N], fysum[N];
	int end, flag[N];
	double ct[N][N];  // used to distroy body
	
	int P, rank; // number of processes

	double time1, time2;  // for measuring time
	double time;

	

	
	

	//Hardcoding the initial data
	double A[6][5]= {
		{25.0, 400.0, 400.0, 0.0, 0.0},
		{20.0, 200.0, 400.0, 3.0, 4.0},
		{30.0, 50.0, 600.0, 1.0, 0.0},
		{50.0, 400.0, 200.0, 1.0, 0.0},
		{40.0, 700.0, 700.0, -1.0, 0.0},
		{70.0, 200.0, 100.0, -1.0, 0.0}
	};

	SZ_Init(P); // initialize MPI message-passing environment
	//average of masses    // used to destroy body
	for (i=0;i<N;i++) {
		for (j=0;j<N;j++) {
			ct[i][j] = (A[i][0] + A[j][0])*0.5;
		}
		
	}

	// taking input from user
	printf("Enter value of timestep:");
	fflush(stdout);
	scanf("%lf",&dt);
	
	printf("Enter total number of iterations:");
	fflush(stdout);
	scanf("%d",&T);

	//printing initial table
	printf("Initial data table:\n");
	fflush(stdout);
	printf("Body\tMass\tx pos\ty pos\tx vel\ty vel\n");
	fflush(stdout);
	for (i=0;i<N;i++) {
		printf("%2d\t",i);
		for (j=0;j<5;j++) {
			printf("%6.2f\t", A[i][j]);
		}
		printf("\n");
	}


	
	time1 = SZ_Wtime();          //time start
	
	SZ_Parallel_begin   // Parallel start

	a = SZ_Get_process_num();  // geting process number
	SZ_Broadcast(&T); 
	SZ_Broadcast(&dt);
	SZ_Broadcast(A);  // broadcasting initial array and other values
	SZ_Broadcast(ct);
	
	//printf("rank %d\n",a);

	for (t=0; t<T; t++ ) {    //Time loop

		//for (a=0; a<N; a++) {   // This loop is removed 
			fxsum[a] = 0;
			fysum[a] = 0;
			for (i=0; i<N; i++) {
				if ((a != i) && (A[i][0] != 0) && (A[a][0] != 0)) {   // if body is destroyed no calculations done
				x_diff = A[i][1]-A[a][1];
				
				y_diff = A[i][2]-A[a][2];
				
				r = sqrt((x_diff*x_diff)+(y_diff*y_diff));	//Calculate distance
				
				F = (G*A[a][0]*A[i][0])/(r*r);           //Calculate force
				
				fx = F*x_diff/r;						//Resolve force
				
				fy = F*y_diff/r;
				
				fxsum[a] = fxsum[a] + fx;
				
				fysum[a] = fysum[a] + fy;
				
				if (r < ct[a][i]) {

					flag[0] = a;
					flag[1] = i;
					A[a][0] = 0.0;
					A[i][0] = 0.0;
				}

				}
			}
			
			 //Updating the values

			if (A[a][0] == 0) {

				//Body is destroyed so not updating anything

			}
			else {
				vx = A[a][3] + (fxsum[a]*dt/A[a][0]);

				vy = A[a][4] + (fysum[a]*dt/A[a][0]);
			
				x = A[a][1] + vx*dt;
			
				y = A[a][2] + vy*dt;
			

				A[a][1] = x;
				A[a][2] = y;
				A[a][3] = vx;
				A[a][4] = vy;
			}

		SZ_AllBroadcast(A);
		
		

	}
	SZ_Parallel_end; //parallel end

	time2 = SZ_Wtime();      //time end
	time = time2 - time1;
	
	for (i=0;i<N;i++) {
		if (A[i][0] == 0.0) {
			printf("Body %d is destroyed.\n", i);
		}
	}

	printf("Final data table:\n");
	printf("Body\tMass\tx pos\ty pos\tx vel\ty vel\n");
	for (i=0;i<N;i++) {
		printf("%2d\t",i);
		for (j=0;j<5;j++) {
			printf("%6.2f\t", A[i][j]);
		}
		printf("\n");
		
	}

	printf("Time elapsed: %f s\n",time);


	




	SZ_Finalize(); 
	return 0;
}
