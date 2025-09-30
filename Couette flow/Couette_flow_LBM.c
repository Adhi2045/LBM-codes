#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){

    int NX = 201, NY = 21;
    float **f[9], **ft[9];  //f is dist. func. and ft is post collision dist. func.
    float feq[9];
    int time = 20000;
    int i, j, a, a1, ts, ia, ja;
    float visc, H;
    float tau = 0.8f;
    float U_wall = 0.1f;
    float rho0 = 1.0f;
    float ux, uy, rho, wt[9];
    float u2, term1, term2;
    int ex[9], ey[9], kb[9];
    float delta = 0.5f;
    float ux_exact[NY];
    float y, y2;
    float c1, c2;
    float l2;
    FILE *out, *data;

    visc = (tau-0.5)/3.0;
    H = NY-1-2.0*delta;

    ex[0] = 0 ; ey[0] = 0;
    ex[1] = 1 ; ey[1] = 0;
    ex[2] = 0 ; ey[2] = 1;
    ex[3] = -1 ; ey[3] = 0;
    ex[4] = 0 ; ey[4] = -1;
    ex[5] = 1 ; ey[5] = 1;
    ex[6] = -1 ; ey[6] = 1;
    ex[7] = -1 ; ey[7] = -1;
    ex[8] = 1 ; ey[8] = -1;

    for (a=0; a<9; a++){
        if (a==0){wt[a]=4.0f/9.0;}
        if (a>=1 & a<=4){wt[a]=1.0f/9.0;}
        if (a>=5 & a<=8){wt[a]=1.0f/36.0;}
    }

    for (a=0; a<9; a++){
        f[a] = (float**)malloc(NX*sizeof(float*));
        ft[a] = (float**)malloc(NX*sizeof(float*));

        for (i=0; i<NX; i++){
            f[a][i] = (float*)malloc(NY*sizeof(float));
            ft[a][i] = (float*)malloc(NY*sizeof(float));
        }
    }

    for (a=0; a<9; a++){
        for (a1=a; a1<9; a1++){
            if (ex[a]+ex[a1]==0 & ey[a]+ey[a1]==0){
                kb[a] = a1;
                kb[a1] = a;
            }
        }
    }

    //Initialize the f = feq(rho0,u=0)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            for (a=0; a<9; a++){
                f[a][i][j] = wt[a]*rho0;
                ft[a][i][j] = wt[a]*rho0;
            }
        }
    }

    //Start the time-loop
    for(ts=1; ts<=time; ts++){
        for (i=0; i<NX; i++){
            for (j=1; j<NY-1; j++){
                ux = 0.0; uy = 0.0; rho = 0.0;
                for (a=0; a<9; a++){
                    rho += f[a][i][j];
                    ux += f[a][i][j] * ex[a];
                    uy += f[a][i][j] * ey[a];
                }
                if (rho > 1e-12) {
                    ux /= rho;
                    uy /= rho;
                } 
                else {
                    ux = uy = 0.0;
                }
                u2 = ux*ux + uy*uy;

                //Collision step
                for (a=0; a<9; a++){
                    term1 = ux*ex[a] + uy*ey[a];
                    term2 = term1*term1;
                    feq[a] = wt[a]*rho*(1.0f+3.0f*term1+4.5f*term2-1.5f*u2);
                    ft[a][i][j] = f[a][i][j] - (f[a][i][j]-feq[a])/tau;
                }
            }
        }

        if (fmod(ts,1000)==0){
            printf("ts = %d Completed\n", ts);
        }

        //Boundary Conditions
        for (i=0; i<NX; i++){

            j = 0; //Bottom wall stationary
            ft[2][i][j] = ft[4][i][j];
            ft[5][i][j] = ft[7][i][j];
            ft[6][i][j] = ft[8][i][j];

            j = NY-1; //Top wall moving
            rho = 0.0;
            for (a=0; a<9; a++){rho += f[a][i][j];}
            ft[4][i][j] = ft[2][i][j];
            ft[7][i][j] = ft[5][i][j] - (1.0f/12.0f) * rho * U_wall; 
            ft[8][i][j] = ft[6][i][j] + (1.0f/12.0f) * rho * U_wall;
        }

        //Streaming
        for (i=0; i<NX; i++){
            for (j=0; j<NY; j++){
                for (a=0; a<9; a++){
                    ia = i-ex[a];
                    ja = j-ey[a];
                    if (ia<0){ia = NX-1;}
                    else if (ia>NX-1){ia = 0;}

                    if (ja >= 0 && ja < NY){
                        f[a][i][j] = ft[a][ia][ja];
                    }
                }
            }
        }
    }

    for (j=0; j<NY; j++){
        y = (j-delta)/H;
        y2 = y*y;
        ux_exact[j] = U_wall*y;
        if (j==0){ux_exact[j]=0;}
        if (j==NY-1){ux_exact[j]=0.1;}
    }
    out = fopen("ux_prof_couette.dat","w");
    i = NX-1;
    l2 = 0.0; c1 = 0.0; c2 = 0.0;
    for (j=0; j<NY; j++){
        rho = 0.0; ux = 0.0; uy = 0.0;
        if (j!=0 && j!=NY-1){
            for (a=0; a<9; a++){
                rho += f[a][i][j];
                ux += f[a][i][j]*ex[a];
                uy += f[a][i][j]*ey[a];
            }
            ux /= rho; uy /= rho;
        }
        if (j==NY-1){ux = U_wall;}

        c1 += ux_exact[j]*ux_exact[j];
        c2 += (ux_exact[j]-ux)*(ux_exact[j]-ux);
        fprintf(out, "%d %d %12.8f %12.8f %12.8f\n", NX-1, j, ux, ux_exact[j], uy, rho);
    }

    printf("c1 = %12.8e c2 = %12.8e\n", c1, c2);
    l2 = pow((c2/c1), 0.5);
    printf("l2 = %12.8e\n", l2);
    fclose(out);
}