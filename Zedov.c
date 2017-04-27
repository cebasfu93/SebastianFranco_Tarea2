#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "struct.h"

//------------------Constantes----------------
#define FLOAT double

#define GAM 1.4
#define NDIM 3
#define Lx 256.0
#define Ly 256.0
#define Lz 256.0
#define delx 2.0
#define dely 2.0
#define delz 2.0
#define T 300.0
#define E 10
#define P_atm 1.0
#define core 63

#define deltat 0.01
#define deltax 0.1
#define deltay 0.1
#define deltaz 0.1

//------------------Variables globales----------------

int Nx = (int) Lx/delx;
int Ny = (int) Ly/dely;
int Nz = (int) Lz/delz;
int N_tot= 2097152; //(int) Nx*Ny*Nz;
FLOAT time_step;

FLOAT ener;
FLOAT pres;
FLOAT entalpia;
FLOAT v_sound;
FLOAT v_max;
FLOAT v_temp;

FLOAT rad;
FILE *rad_dat, *dens_dat;

int i, j, k, m;

//------------------Declaraciones Funciones----------------
void init_to_zero(FLOAT *p, int n_points);
physics_grid * create_physics_grid(void);
U_grid * create_U_grid(void);
F_grid * create_F_grid(void);
void init_problem(physics_grid *P, U_grid *U, F_grid *Fx, F_grid *Fy, F_grid *Fz);
void init_zedov(physics_grid *P, U_grid *U, F_grid *Fx, F_grid *Fy, F_grid *Fz);
int ndx(int i, int j, int k, FLOAT length, FLOAT height, int cube);
void update(physics_grid *P, U_grid *U, F_grid *Fx, F_grid *Fy, F_grid *Fz, FLOAT delt);
FLOAT fromU2Fx(FLOAT u0, FLOAT u1, FLOAT u2, FLOAT u3, FLOAT u4, int m);
FLOAT fromU2Fy(FLOAT u0, FLOAT u1, FLOAT u2, FLOAT u3, FLOAT u4, int m);
FLOAT fromU2Fz(FLOAT u0, FLOAT u1, FLOAT u2, FLOAT u3, FLOAT u4, int m);
FLOAT presion(FLOAT u0, FLOAT u1, FLOAT u2, FLOAT u3, FLOAT u4);
void fromUtoRho(U_grid *U);

//------------------Main----------------
int main(){
  rad_dat=fopen("rad_dat.txt", "w");
  dens_dat=fopen("dens_dat.txt", "w");
  physics_grid * P_state;
  U_grid * U_state;
  F_grid * Fx_state;
  F_grid * Fy_state;
  F_grid * Fz_state;

  P_state = create_physics_grid();
  U_state = create_U_grid();
  Fx_state = create_F_grid();
  Fy_state = create_F_grid();
  Fz_state = create_F_grid();

  init_problem(P_state, U_state, Fx_state, Fy_state, Fz_state);
  init_zedov(P_state, U_state, Fx_state, Fy_state, Fz_state);
  fromUtoRho(U_state);
  /*for(i=0;i<Ny;i++){
    for(j=0;j<Nz;j++){
      printf("%f ", P_state->P[ndx(i,j,core, Ny, Nz,3)]);
    }
    printf("\n");
  }*/

  return 0;
}

//------------------Funciones----------------
void init_to_zero(FLOAT *p, int n_points){
  for(i=0;i<n_points;i++){
    p[i]=0.0;
  }
}
physics_grid * create_physics_grid(void){
  physics_grid *G;
  if(!(G = malloc(sizeof(physics_grid)))){
   fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
   exit(0);
  }
  G->L_x=0.0;
  G->L_y=0.0;
  G->L_z=0.0;
  G->delta_x=0.0;
  G->delta_y=0.0;
  G->delta_z=0.0;
  G->N_x=0;
  G->N_y=0;
  G->N_z=0;
  G->N_cells=0;
  G->P=NULL;
  return G;
}
U_grid * create_U_grid(void){
  U_grid *G;
  if(!(G = malloc(sizeof(U_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  }
  G->N_x=0.0;
  G->N_y=0.0;
  G->N_z=0.0;
  G->N_cells=0.0;
  G->U=NULL;
  return G;
}
F_grid * create_F_grid(void){
  F_grid *G;
  if(!(G = malloc(sizeof(F_grid)))){
    fprintf(stderr, "Problem with data allocation\n");fflush(stdout);
    exit(0);
  }
  G->N_x=0.0;
  G->N_y=0.0;
  G->N_z=0.0;
  G->N_cells=0.0;
  G->F=NULL;
  return G;
}
void init_problem(physics_grid *P, U_grid *U, F_grid *Fx, F_grid *Fy, F_grid *Fz){

  P->L_x = Lx;
  P->L_y = Ly;
  P->L_z = Lz;
  P->delta_x = delx;
  P->delta_y = dely;
  P->delta_z = delz;
  P->N_x = Nx;
  P->N_y = Ny;
  P->N_z = Nz;
  P->N_cells = N_tot;

  U->N_x = Nx;
  U->N_y = Ny;
  U->N_z = Nz;
  U->N_cells = N_tot;

  Fx->N_x = Nx;
  Fx->N_y = Ny;
  Fx->N_z = Nz;
  Fx->N_cells = N_tot;

  Fy->N_x = Nx;
  Fy->N_y = Ny;
  Fy->N_z = Nz;
  Fy->N_cells = N_tot;

  Fz->N_x = Nx;
  Fz->N_y = Ny;
  Fz->N_z = Nz;
  Fz->N_cells = N_tot;

  if(!(P->P = malloc(P->N_cells * (NDIM +2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with pressure allocation");
    exit(1);
  }
  init_to_zero(P->P, P->N_cells * (NDIM +2));

  if(!(U->U = malloc(U->N_cells * (NDIM +2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with U_1 allocation");
    exit(1);
  }
  init_to_zero(U->U, U->N_cells * (NDIM +2));

  if(!(Fx->F = malloc(Fx->N_cells * (NDIM) * (NDIM + 2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F allocation");
    exit(1);
  }
  init_to_zero(Fx->F, Fx->N_cells * NDIM * (NDIM +2));

  if(!(Fy->F = malloc(Fy->N_cells * (NDIM) * (NDIM + 2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F allocation");
    exit(1);
  }
  init_to_zero(Fy->F, Fy->N_cells * NDIM * (NDIM +2));

  if(!(Fz->F = malloc(Fz->N_cells * (NDIM) * (NDIM + 2) * sizeof(FLOAT)))){
    fprintf(stderr, "Problem with F allocation");
    exit(1);
  }
  init_to_zero(Fz->F, Fz->N_cells * NDIM * (NDIM +2));
}
void init_zedov(physics_grid *P, U_grid *U, F_grid *Fx, F_grid *Fy, F_grid *Fz){

  FLOAT dens=1.0/E*(GAM-1);
  /*P->P[ndx(core,core,core, P->N_x, P->N_y, 3)]=P_atm;
  P->P[ndx(core,core,core, P->N_x, P->N_y, 4)]=dens;

  U->U[ndx(core,core,core, U->N_x, U->N_y, 0)]=dens;
  U->U[ndx(core,core,core, U->N_x, U->N_y, 4)]=dens*E;

  Fx->F[ndx(core,core,core, Fx->N_x, Fx->N_y, 1)]=P_atm;

  Fy->F[ndx(core,core,core, Fy->N_x, Fy->N_y, 2)]=P_atm;

  Fz->F[ndx(core,core,core, Fz->N_x, Fz->N_y, 3)]=P_atm;*/

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
        U->U[ndx(i,j,k,Nx,Ny,0)]= (FLOAT) pow(pow(i-core,2)+pow(j-core,2)+pow(j-core,2),0.5)*delx/2;
      }
    }
  }

}
int ndx(int i, int j, int k, FLOAT length, FLOAT height, int cube){
  int res;
  res = ((i*length*height)+j*height+k)+cube*N_tot;
  return res;
}
void update(physics_grid *P, U_grid *U, F_grid *Fx, F_grid *Fy, F_grid *Fz, FLOAT delt){

  FLOAT * U_halfs;
  FLOAT * F_halfs;

  U_halfs=malloc(sizeof(FLOAT)*2*NDIM*(NDIM+2));
  F_halfs=malloc(sizeof(FLOAT)*2*NDIM*(NDIM+2));

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
        for(m=0;m<(NDIM+2);m++){
          if (i==0 || j==0 || k==0 || i==(Nx-1) || j==(Ny-1) || k==(Nz-1)){
            U->U[ndx(i, j, k, Nx, Ny, m)]=U->U[ndx(i, j, k, Nx, Ny, m)];
          }
          else{
            U_halfs[0+m*(2*NDIM)]=(U->U[ndx(i-1,j,k,Nx,Ny,m)]+U->U[ndx(i,j,k,Nx,Ny,m)])*0.5;
            U_halfs[1+m*(2*NDIM)]=(U->U[ndx(i,j,k,Nx,Ny,m)]+U->U[ndx(i+1,j,k,Nx,Ny,m)])*0.5;
            U_halfs[2+m*(2*NDIM)]=(U->U[ndx(i,j-1,k,Nx,Ny,m)]+U->U[ndx(i,j,k,Nx,Ny,m)])*0.5;
            U_halfs[3+m*(2*NDIM)]=(U->U[ndx(i,j,k,Nx,Ny,m)]+U->U[ndx(i,j+1,k,Nx,Ny,m)])*0.5;
            U_halfs[4+m*(2*NDIM)]=(U->U[ndx(i,j,k-1,Nx,Ny,m)]+U->U[ndx(i,j,k,Nx,Ny,m)])*0.5;
            U_halfs[5+m*(2*NDIM)]=(U->U[ndx(i,j,k,Nx,Ny,m)]+U->U[ndx(i,j,k+1,Nx,Ny,m)])*0.5;
          }
        }
        for(m=0;m<(NDIM+2);m++){
          F_halfs[0+m*(2*NDIM)]=fromU2Fx(U_halfs[0], U_halfs[6], U_halfs[12], U_halfs[18], U_halfs[24], m);
          F_halfs[1+m*(2*NDIM)]=fromU2Fx(U_halfs[1], U_halfs[7], U_halfs[13], U_halfs[19], U_halfs[25], m);
          F_halfs[2+m*(2*NDIM)]=fromU2Fy(U_halfs[2], U_halfs[8], U_halfs[14], U_halfs[20], U_halfs[26], m);
          F_halfs[3+m*(2*NDIM)]=fromU2Fy(U_halfs[3], U_halfs[9], U_halfs[15], U_halfs[21], U_halfs[27], m);
          F_halfs[4+m*(2*NDIM)]=fromU2Fz(U_halfs[4], U_halfs[10], U_halfs[16], U_halfs[22], U_halfs[28], m);
          F_halfs[5+m*(2*NDIM)]=fromU2Fz(U_halfs[5], U_halfs[11], U_halfs[17], U_halfs[23], U_halfs[29], m);
        }
        for(m=0;m<(NDIM+2);m++){
          U->U[ndx(i,j,k,Nx,Ny,m)]=U->U[ndx(i,j,k,Nx,Ny,m)]+delt/deltax*(F_halfs[0+m*(2*NDIM)]-F_halfs[1+m*(2*NDIM)])+delt/deltay*(F_halfs[2+m*(2*NDIM)]-F_halfs[3+m*(2*NDIM)])+delt/deltaz*(F_halfs[4+m*(2*NDIM)]-F_halfs[5+m*(2*NDIM)]);
        }
      }
    }
  }


}
FLOAT fromU2Fx(FLOAT u0, FLOAT u1, FLOAT u2, FLOAT u3, FLOAT u4, int m){
  if(m==0){
    return u1;
  }
  else if(m==1){
    FLOAT P=presion(u0, u1, u2, u3, u4);
    return pow(u1,2)/u0+P;
  }
  else if(m==2){
    return u1*u2/u0;
  }
  else if(m==3){
    return u1*u3/u0;
  }
  else{
    FLOAT P=presion(u0, u1, u2, u3, u4);
    return u4+P*u1/u0;
  }
}
FLOAT fromU2Fy(FLOAT u0, FLOAT u1, FLOAT u2, FLOAT u3, FLOAT u4, int m){
  if(m==0){
    return u2;
  }
  else if(m==1){
    return u1*u2/u0;
  }
  else if(m==2){
    FLOAT P=presion(u0, u1, u2, u3, u4);
    return pow(u2,2)*(u0)+P;
  }
  else if(m==3){
    return u2*u3/u0;
  }
  else{
    FLOAT P=presion(u0, u1, u2, u3, u4);
    return u4+P*u2/u0;
  }
}
FLOAT fromU2Fz(FLOAT u0, FLOAT u1, FLOAT u2, FLOAT u3, FLOAT u4, int m){
  if(m==0){
    return u3;
  }
  else if(m==1){
    return u1*u3/u0;
  }
  else if(m==2){
    return u2*u3/u0;
  }
  else if(m==3){
    FLOAT P=presion(u0, u1, u2, u3, u4);
    return pow(u3,2)/u0+P;
  }
  else{
    FLOAT P=presion(u0, u1, u2, u3, u4);
    return u4+P*u3/u0;
  }
}
FLOAT presion(FLOAT u0, FLOAT u1, FLOAT u2, FLOAT u3, FLOAT u4){
  return (GAM-1)*(u4-0.5*(pow(u1,2)+pow(u2,2)+pow(u3,2)));
}
FLOAT t_step(physics_grid *P, U_grid *U){
  v_max=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
        ener=U->U[ndx(i,j,k,Nx,Ny,4)]/U->U[ndx(i,j,k,Nx,Ny,0)];
        pres=presion(U->U[ndx(i,j,k,Nx,Ny,0)], U->U[ndx(i,j,k,Nx,Ny,1)], U->U[ndx(i,j,k,Nx,Ny,2)], U->U[ndx(i,j,k,Nx,Ny,3)], U->U[ndx(i,j,k,Nx,Ny,4)]);
        entalpia=ener+pres/U->U[ndx(i,j,k,Nx,Ny,0)];
        v_sound=pow((GAM-1)*entalpia,0.5);
        v_temp=fabs(P->P[ndx(i,j,k,Nx,Ny,0)])+v_sound;
        if(v_temp > v_max){
          v_max=v_temp;
        }
        v_temp=fabs(P->P[ndx(i,j,k,Nx,Ny,1)])+v_sound;
        if(v_temp > v_max){
          v_max=v_temp;
        }
        v_temp=fabs(P->P[ndx(i,j,k,Nx,Ny,2)])+v_sound;
        if(v_temp > v_max){
          v_max=v_temp;
        }
      }
    }
  }
  time_step=2.0*delx/v_max;
  return time_step;
}
void fromUtoRho(U_grid *U){
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
        rad=(FLOAT) pow(pow(i-core,2)+pow(j-core,2)+pow(j-core,2),0.5)*delx;
        fprintf(rad_dat, "%f ", rad);
        fprintf(dens_dat, "%f ", U->U[ndx(i,j,k,Nx,Ny,0)]);
      }
      fprintf(rad_dat, "\n");
      fprintf(dens_dat, "\n");
    }
  }
}
