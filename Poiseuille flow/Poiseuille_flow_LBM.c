#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    int nx = 201, ny = 21;
    int i, j, a, a1, ts, ia, ja;
    float **f[9], **ft[9];
    int isn[nx][ny];
    float wt[9], ux, uy, rho, term1, term2, u2;
    int ex[9], ey[9], kb[9];
    float source[9];
    int time = 20000;
    float dpdx = 1.0e-5;
    float tau = 0.8;
    float rho0 = 1.0;
    float visc, H;
    float feq[9];
    float rho_av;
    float ux_exact[ny];
    int kpor;
    float delta = 0.5;
    float y, y2;
    float c1, c2;
    float l2;
    FILE *out, *data;

    visc = (tau-0.5)/3.0;
    H = ny-1-2.0*delta;
    
    ex[0] = 0; ey[0] = 0;
    ex[1] = 1; ey[1] = 0;
    ex[2] = 0; ey[2] = 1;
    ex[3] = -1; ey[3] = 0;
    ex[4] = 0; ey[4] = -1;
    ex[5] = 1; ey[5] = 1;
    ex[6] = -1; ey[6] = 1;
    ex[7] = -1; ey[7] = -1;
    ex[8] = 1; ey[8] = -1;

    for (a=0; a<9; a++){
        f[a] = (float**)malloc(nx*sizeof(float*));
        ft[a] = (float**)malloc(nx*sizeof(float*));

        for (i=0; i<nx; i++){
            f[a][i] = (float*)malloc(ny*sizeof(float));
            ft[a][i] = (float*)malloc(ny*sizeof(float));
        }
    }

    for (a=0; a<9; a++){
        if (a==0){wt[a] = 4.0/9.0;}
        if (a>=1 & a<=4){wt[a] = 1.0/9.0;}
        if (a>=5 & a<=8){wt[a] = 1.0/36.0;}
    }

    for (a=0; a<9; a++){
        for (a1=a; a1<9; a1++){
            if ((ex[a]+ex[a1])==0 & (ey[a]+ey[a1])==0){
                kb[a] = a1;
                kb[a1] = a;
            }
        }
    }

    for (a=0; a<9; a++){
        printf("a = %d, kb = %d\n", a, kb[a]);
    }

    for (i=0; i<nx; i++){
        for (j=0; j<ny; j++){
            isn[i][j] = 0;
            if(j==0|j==ny-1){isn[i][j]=1;}
            for (a=0; a<9; a++){
                f[a][i][j] = wt[a]*rho0;
            }
        }
    }

    for (ts=1; ts<=time; ts++){         //this is the time loop
        rho_av = 0.0; kpor = 0;
        for (i=0; i<nx; i++){
            for(j=1; j<ny-1; j++){
                ux = 0.0; uy = 0.0; rho = 0.0;
                for (a=0; a<9; a++){
                    rho += f[a][i][j];
                    ux += f[a][i][j]*ex[a];
                    uy += f[a][i][j]*ey[a];
                }

                rho_av += rho;
                kpor += 1;
                ux += dpdx/2.0;
                ux /= rho; uy /= rho;
                u2 = ux*ux+uy*uy;

                for (a=0; a<9; a++){
                    term1 = ux*ex[a] + uy*ey[a];
                    term2 = term1*term1;
                    source[a] = (1.0-0.5/tau)*wt[a]*(3.0*(ex[a]-ux)+9.0*(ex[a]*ux+ey[a]*uy)*ex[a])*dpdx;
                    feq[a] = wt[a]*rho*(1.0+3.0*term1+4.5*term2-1.5*u2);
                    ft[a][i][j] = f[a][i][j] - (f[a][i][j]-feq[a])/tau + source[a];
                }
            }    
        }

        if (fmod(ts,100)==0){
            printf("ts = %d\t rho_av = %12.8f\n", ts, rho_av/kpor);
        }

        for (i=0; i<nx; i++){
            for (j=1; j<ny-1; j++){
                for (a=0; a<9; a++){
                    ia = i+ex[a];
                    ja = j+ey[a];
                    if (ia<0){ia = nx-1;}
                    if (ia>nx-1){ia = 0;}
                    f[a][ia][ja] = ft[a][i][j];
                }
            }
        }

        for (i=0; i<nx; i++){
            for (j=1; j<ny-1; j++){
                if (isn[i][j]==0){
                    for (a=0; a<9; a++){
                        ia = i-ex[a];
                        ja = j-ey[a];
                        if (ia<0){ia = nx-1;}
                        if (ia>nx-1){ia = 0;}

                        if (isn[ia][ja]==1){
                            f[a][i][j] = f[kb[a]][ia][ja];
                        }
                    } 
                }
            }
        }
    }

    for (j=0; j<ny; j++){
        y = (j-delta)/H;
        y2 = y*y;
        ux_exact[j] = 0.5*dpdx*H*H*(y-y2)/visc;
        if (j==0|j==ny-1){ux_exact[j]=0;}
    }

    out = fopen("ux_prof.dat","w");
    i = nx-1;
    l2 = 0.0; c1 = 0.0; c2 = 0.0;
    for (j=0; j<ny; j++){
        rho = 0.0; ux = 0.0; uy = 0.0;
        if (isn[i][j]==0){
            for (a=0; a<9; a++){
                rho += f[a][i][j];
                ux += f[a][i][j]*ex[a];
                uy += f[a][i][j]*ey[a];
            }
            ux /= rho; uy /= rho;
        }
        c1 += ux_exact[j]*ux_exact[j];
        c2 += (ux_exact[j]-ux)*(ux_exact[j]-ux);
        fprintf(out, "%d %d %12.8f %12.8f %12.8f\n", nx-1, j, ux, ux_exact[j], uy, rho);
    }

    printf("c1 = %12.8e c2 = %12.8e\n", c1, c2);
    l2 = pow((c2/c1), 0.5);
    printf("l2 = %12.8e\n", l2);
    fclose(out);
}