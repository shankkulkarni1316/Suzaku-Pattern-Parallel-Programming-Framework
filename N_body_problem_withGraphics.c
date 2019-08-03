#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#define G 100  // gravitational constant
#define N 6  // number of bodies   // number of tasks
#include <string.h>
#include "suzaku.h" // Suzaku routines

// required Suzaku constants
#define D 6 // number of data items in each task, doubles only
#define R 6 // number of data items in result of each task, doubles only

//global parameters
double dt;
double ct[N][N];  // used to distroy body
double A[N][N][2];
int a, i;
double fxsum, fysum;
double fx, fy, vx, vy, x, y, dx, dy;
double x_diff, y_diff, r, F;

void init(int *tasks, int *data_items, int *result_items) {
	*tasks = N;
	*data_items = D;
	*result_items = R;
}

void diffuse(int taskID,double output[D]) { // taskID not used in computation
	a = taskID;	

	output[0] = A[a][0][0];
	output[1] = A[a][1][0];
	output[2] = A[a][2][0];
	output[3] = A[a][3][0];
	output[4] = A[a][4][0];
	output[5] = A[a][5][0]; 
}

void compute(int taskID, double input[D], double output[R]) {
	double r2;
	
	//printf("check 2\n" );
	a = taskID;
	
	//printf("task id %d\n", a);
	A[a][0][0]=input[0];
	A[a][1][0]=input[1];
	A[a][2][0]=input[2];
	A[a][3][0]=input[3];
	A[a][4][0]=input[4];
	A[a][5][0]=input[5];
	
	
	fxsum = 0;
	fysum = 0;
	for (i=0; i<N; i++) {
		if ((a != i) && (A[i][0][0] != 0) && (A[a][0][0] != 0)) {   // if body is destroyed no calculations done
			x_diff = A[i][1][0]-A[a][1][0];
		
			y_diff = A[i][2][0]-A[a][2][0];
			r2 = (x_diff*x_diff)+(y_diff*y_diff);	
			r = sqrt(r2);	//Calculate distance
				
			F = (G*A[a][0][0]*A[i][0][0])/(r*r);           //Calculate force
				
			fx = F*x_diff/r;						//Resolve force
				
			fy = F*y_diff/r;
				
			fxsum = fxsum + fx;
				
			fysum = fysum + fy;
				
			if (r < ct[a][i]) {

				
				A[a][0][0] = 0.0;
				A[i][0][0] = 0.0;
				A[a][0][1] = 0.0;
				A[i][0][1] = 0.0;
			}

		}
	}
			
			 //Updating the values

	if (A[a][0][0] == 0) {

		//Body is destroyed so not updating anything

	}
	else {
		vx = A[a][3][0] + (fxsum*dt/A[a][0][0]);

		vy = A[a][4][0] + (fysum*dt/A[a][0][0]);
			
		x = A[a][1][0] + vx*dt;
			
		y = A[a][2][0] + vy*dt;
			

		A[a][1][1] = x;
		A[a][2][1] = y;
		A[a][3][1] = vx;
		A[a][4][1] = vy;
	}
	output[0] = A[a][0][1];
	output[1] = A[a][1][1];
	output[2] = A[a][2][1];
	output[3] = A[a][3][1];
	output[4] = A[a][4][1];
	output[5] = A[a][5][1];
}

void gather(int taskID, double input[R]) {
	
	a = taskID;
	A[a][0][1] = input[0];
	A[a][1][1] = input[1];
	A[a][2][1] = input[2];
	A[a][3][1] = input[3];
	A[a][4][1] = input[4];
	A[a][5][1] = input[5];
}



int main(int argc, char *argv[] ) {
	//double dt;
	int  j, T, t;

	
	int P, rank; // number of processes

	double time1, time2;  // for measuring time
	double time;


	

	
	

	//Hardcoding the initial data
	double A1[6][5]= {
		{25.0, 400.0, 400.0, 0.0, 0.0},
		{20.0, 200.0, 400.0, 3.0, 4.0},
		{30.0, 50.0, 600.0, 1.0, 0.0},
		{50.0, 400.0, 200.0, 1.0, 0.0},
		{40.0, 700.0, 700.0, -1.0, 0.0},
		{70.0, 200.0, 100.0, -1.0, 0.0}
	};
	
	//Initializing array A
	for (i=0;i<N;i++) {
		for (j=0;j<N;j++) {
			A[i][j][0] = A1[i][j];
			A[i][j][1] = A1[i][j];
			
		}
		
	}	

	SZ_Init(P); // initialize MPI message-passing environment
	//average of masses    // used to destroy body
	for (i=0;i<N;i++) {
		for (j=0;j<N;j++) {
			ct[i][j] = (A[i][0][0] + A[j][0][0])*0.5;
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
			printf("%6.2f\t", A[i][j][0]);
		}
		printf("\n");
	}


	
	time1 = SZ_Wtime();          //time start
	
	printf("check 1\n");
	SZ_Parallel_begin   // Parallel start

	a = SZ_Get_process_num();  // geting process number
	SZ_Broadcast(&T); 
	SZ_Broadcast(&dt);
	SZ_Broadcast(A);  // broadcasting initial array and other values
	SZ_Broadcast(ct);
	
	

	for (t=0; t<T; t++ ) {    //Time loop

	
		SZ_Workpool(init,diffuse,compute,gather);   // Workpool here
		
		
		for (i=0;i<N;i++) {
			for (j=0;j<N;j++) {
				A[i][j][0]=A[i][j][1];
			}
		}
		SZ_Broadcast(A);   // broadcasting updated array again
		

	}
	SZ_Parallel_end; //parallel end

	time2 = SZ_Wtime();      //time end
	time = time2 - time1;
	
	for (i=0;i<N;i++) {
		if (A[i][0][0] == 0.0) {
			printf("Body %d is destroyed.\n", i);
		}
	}

	printf("Final data table:\n");
	printf("Body\tMass\tx pos\ty pos\tx vel\ty vel\n");
	for (i=0;i<N;i++) {
		printf("%2d\t",i);
		for (j=0;j<5;j++) {
			printf("%6.2f\t", A[i][j][0]);
		}
		printf("\n");
		
	}

	printf("Time elapsed: %f s\n",time);


	




	SZ_Finalize(); 
	return 0;
}
