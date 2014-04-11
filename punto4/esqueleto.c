#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1

FLOAT * get_memory(int n_points);
void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points, FLOAT radius);
void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT vel, FLOAT radius);
void initialize_mass(FLOAT *mass, int n_points, FLOAT unit_mass);
void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points);
void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points);
void get_kineticenergy(FLOAT Ktotal, FLOAT *K, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *mass, int n_points);
void get_potentialenergy(FLOAT Utotal, FLOAT *U,FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points);
void initialize_kineticenergy(FLOAT *K, int n_points);
void initialize_potentialenergy(FLOAT *mass, FLOAT UTotal,FLOAT *U, FLOAT *x,Float*y, FLOAT *z, int n_points);
int main(int argc, char **argv){
 
  /*positions of all particles*/
  FLOAT *x;
  FLOAT *y;
  FLOAT *z;
 
  /*velocities of all particles*/
  FLOAT *v_x;
  FLOAT *v_y;
  FLOAT *v_z;

  /*accelerations of all particles*/
  FLOAT *a_x;
  FLOAT *a_y;
  FLOAT *a_z;

  /*Energia potencial*/
  FLOAT *U;
  FLOAT Utotal;

  /*Energia cinetica*/
  FLOAT *K;
  FLOAT Ktotal;

  /*masses*/
  FLOAT *mass;

  /*timestep variables*/
  FLOAT delta_t= 0.001;
  int n_steps = (int)(1000000000/delta_t);
  int n_points = 10000;
  FLOAT unit_mass = 1.0;
  
  int i,j;
 
  /*memory allocation*/
  x = get_memory(n_points);
  y = get_memory(n_points);
  z = get_memory(n_points);
  v_x = get_memory(n_points);
  v_y = get_memory(n_points);
  v_z = get_memory(n_points);
  a_x = get_memory(n_points);
  a_y = get_memory(n_points);
  a_z = get_memory(n_points);
  mass = get_memory(n_points);
  U = get_memory(n_points);
  K = get_memory(n_points);

  initialize_pos(x,y,z, n_points);
  initialize_vel(v_x,v_y,v_z, n_points);
  initialize_mass(mass, n_points, unit_mass);
  initialize_potentialenergy(mass,U,Utotal,x,y,z, n_points);
  initialize_kineticenergy(K, n_points);

  /*implementation of a simple Euler integration*/
  for(i=0;i<n_steps;i++){
    get_acceleration(a_x, a_y, a_z, x, y, z, mass, n_points);
    for(j=0;j<n_points;j++){   
      
      x[j] = x[j] + delta_t * v_x[j];
      y[j] = y[j] + delta_t * v_y[j];
      z[j] = z[j] + delta_t * v_z[j];

      v_x[j] = v_x[j] + delta_t * a_x[j];
      v_y[j] = v_y[j] + delta_t * a_y[j];
      v_z[j] = v_z[j] + delta_t * a_z[j];
        
      
   }
    print_status(x,y,z,v_x,v_y,v_z, a_x, a_y, a_z, n_points);
  }

  /*Creo una funcion que grafique el cociente entre el cambio de la energia y la energia inicial y otra que grafique 2K+U */
}

void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points){
  int i,j;
  FLOAT r_ij;
  for(i=0;i<n_points;i++){
    ax[i]=0.0;
    ay[i]=0.0;
    az[i]=0.0;
   
    for(j=0;j<n_points;j++){
      if(j!=i){
    r_ij = (pow((x[i] - x[j]),2.0) +
        pow((y[i] - y[j]),2.0) +
        pow((z[i] - z[j]),2.0));
    r_ij = sqrt(r_ij);
    /*modifico la ley de la gravedad*/
    ax[i] += -G_GRAV *mass[j]/ pow(r_ij+0.01,1.5) * (x[i] - x[j]);
    ay[i] += -G_GRAV *mass[j]/ pow(r_ij+0.01,1.5) * (y[i] - y[j]);
    az[i] += -G_GRAV *mass[j] / pow(r_ij+0.01,1.5) * (z[i] - z[j]);
      }
    }   
  } 
}
void get_potentialenergy(FLOAT Utotal, FLOAT *U,FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points){
  int i,j;
  FLOAT r_ij;
  

  for(i=0;i<n_points;i++){
    
   
    for(j=0;j<n_points;j++){
      if(j!=i){
    r_ij = (pow((x[i] - x[j]),2.0) +
        pow((y[i] - y[j]),2.0) +
        pow((z[i] - z[j]),2.0));
    r_ij = sqrt(r_ij);
    /*modifico la ley de la gravedad*/
    U[j] += -G_GRAV *mass[j]*mass[i]/ (r_ij+0.01);
    
    
      }
    }
    Utotal+= U[i]/2.0;   
  }
}
void get_kineticenergy(FLOAT Ktotal, FLOAT *K, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *mass, int n_points){
  int i;
  
  for(i=0;i<n_points;i++){
    K[i]=mass[i]*(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0))*0.5;
    Ktotal+= K[i];
   
    
      }
    }   
  } 
}



void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points, FLOAT radius){
  int i,j;
  FLOAT theta;
  FLOAT phi;
  FLOAT rad;

 
  /*genero los puntos aleatoriamente en una esfera de radio 20 parsecs*/
  for(i=0;i<n_points;i++){
    rad=pow(drand48(),0.33333)*20.0;
    phi=2.0*PI*drand48();
    theta=arccos((2.0*drand48())-1.0);
    x[i] = sin(theta) * cos(phi) * rad;
    y[i] = sin(theta) * sin(phi) * rad;
    z[i] = cos(theta) * rad;
    for(j=0;j<n_points;j++){
     if 
   } 
 }
}

void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT vel, FLOAT radius){
  int i;
  /*establezco las velocidades iniciales iguales a cero*/
 
  for(i=0;i<n_points;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
  } 

}

void initialize_mass(FLOAT *mass, int n_points, FLOAT unit_mass){
  int i;
  for (i=0;i<n_points;i++){
    mass[i] = 1.0;
  }
}

void initialize_kineticenergy(FLOAT *K, int n_points){
  int i;
  for (i=0;i<n_points;i++){
    K[i]=0.0;
 }
}
   


FLOAT * get_memory(int n_points){
  FLOAT * x;
  if(!(x = malloc(sizeof(FLOAT) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;
}

void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points){
  int i;
  for(i=0;i<n_points;i++){
    printf("%f %f %f ",
       x[i], y[i], z[i]);
  }
  printf("\n");
} 
