#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#define G 100  // gravitational constant
#define N 6  // number of bodies

#include "X11Macros.h" // X11 macros
#define X_RESN 1000 // x resolution
#define Y_RESN 1000 // y resolution


int main(int argc, char **argv ) {
	double dt;
	int i, j, a, T, t;
	double fx, fy, vx, vy, x, y, dx, dy;
	double x_diff, y_diff, r, F;
	double fxsum[N], fysum[N];
	int end, flag[N];
	double ct[N][N];  // used to distroy body
	int x1, y1; // for graphical output

	clock_t time1, time2;
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
	
	//average of masses    // used to destroy body
	for (i=0;i<N;i++) {
		for (j=0;j<N;j++) {
			ct[i][j] = (A[i][0] + A[j][0])*0.5;
		}
		
	}

	// taking input from user
	printf("Enter value of timestep:");
	scanf("%lf",&dt);
	
	printf("Enter total number of iterations:");
	scanf("%d",&T);

	//printing initial table
	printf("Initial data table:\n");
	printf("Body\tMass\tx pos\ty pos\tx vel\ty vel\n");
	for (i=0;i<N;i++) {
		printf("%2d\t",i);
		for (j=0;j<5;j++) {
			printf("%6.2f\t", A[i][j]);
		}
		printf("\n");
	}

	/* --------------------------- X11 graphics setup ------------------------------ */
	initX11(X_RESN,Y_RESN); // includes the X11 initialization code

	time1 = clock();          //time start
	for (t=0; t<T; t++ ) {    //Time loop
		for (a=0; a<N; a++) {
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

					//flag[0] = a;
					//flag[1] = i;
					A[a][0] = 0.0;
					A[i][0] = 0.0;
				}

				}
			}
			
		}
		//printf("check 1");
		for (a=0; a<N; a++) {  //Updating the values

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


		}
	XClearWindow(display, win);  				// clear window for next drawing
	XSetForeground(display,gc,(long) 0xDC143C);  		// color of foreground (applies to object to be drawn)
	

	XFillArc (display,win,gc,A[0][1]-12.5,A[0][2]-12.5,A[0][0],A[0][0],0,23040);

	XFillArc (display,win,gc,A[1][1]-10,A[1][2]-10,A[1][0],A[1][0],0,23040);

	XFillArc (display,win,gc,A[2][1]-15,A[2][2]-15,A[2][0],A[2][0],0,23040);

	XFillArc (display,win,gc,A[3][1]-25,A[3][2]-25,A[3][0],A[3][0],0,23040);

	XFillArc (display,win,gc,A[4][1]-20,A[4][2]-20,A[4][0],A[4][0],0,23040);

	XFillArc (display,win,gc,A[5][1]-35,A[5][2]-35,A[5][0],A[5][0],0,23040);
	


 
	
	XFlush(display);					// necessary to write to display
	
	usleep(10000);						// provide a delay beween each drawing


	}
	time2 = clock();      //time end
	time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	
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

	





	return 0;
}
