#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1

/* El siguiente codigo esta basado en el codigo de Jaime Forero. Las diferencias estan en el cambio de velocidades iniciales para los tres cuerpos y en que se implementa el metodo de Runge kutta de cuarto orden en lugar del metodo de Euler para poder enocntrar las posiciones y velocidades de los tres cuerpos en el tiempo
 */

FLOAT * get_memory(int n_points);
void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points, FLOAT radius);
void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT vel, FLOAT radius);
void initialize_mass(FLOAT *mass, int n_points, FLOAT unit_mass);
void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points);
void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points);

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
    
  /*terminos usados en la implementacion del metodo de Runge Kutta*/
  
  FLOAT *a2_x;
  FLOAT *a2_y;
  FLOAT *a2_z;
  FLOAT *a3_x;
  FLOAT *a3_y;
  FLOAT *a3_z;
  FLOAT *a4_x;
  FLOAT *a4_y;
  FLOAT *a4_z;

  FLOAT *dx1;
  FLOAT *dy1;
  FLOAT *dz1;
  FLOAT *dx2;
  FLOAT *dy2;
  FLOAT *dz2;
  FLOAT *dx3;
  FLOAT *dy3;
  FLOAT *dz3;
  FLOAT *dx4;
  FLOAT *dy4;
  FLOAT *dz4;

  FLOAT *dv1_x;
  FLOAT *dv1_y;
  FLOAT *dv1_z;
  FLOAT *dv2_x;
  FLOAT *dv2_y;
  FLOAT *dv2_z;
  FLOAT *dv3_x;
  FLOAT *dv3_y;
  FLOAT *dv3_z;
  FLOAT *dv4_x;
  FLOAT *dv4_y;
  FLOAT *dv4_z;

  FLOAT *x2;
  FLOAT *y2;
  FLOAT *z2;
  FLOAT *x3;
  FLOAT *y3;
  FLOAT *z3;
  FLOAT *x4;
  FLOAT *y4;
  FLOAT *z4;  
    
  /*masses*/
  FLOAT *mass;

  /*timestep variables*/
  FLOAT h= 0.001;
  int n_steps = (int)(100.0/h);
  int n_points = 3;
  FLOAT radius = 100.0;
  FLOAT unit_mass = 1.0; 
  FLOAT vel_initial = sqrt((11.0/3.0) * G_GRAV * unit_mass / (sqrt(3.0)*radius));
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
  a2_x = get_memory(n_points);
  a2_y = get_memory(n_points);
  a2_z = get_memory(n_points);
  a3_x = get_memory(n_points);
  a3_y = get_memory(n_points);
  a3_z = get_memory(n_points);
  a4_x = get_memory(n_points);
  a4_y = get_memory(n_points);
  a4_z = get_memory(n_points);
  dx1 = get_memory(n_points);
  dy1 = get_memory(n_points);
  dz1 = get_memory(n_points);
  dx2 = get_memory(n_points);
  dy2 = get_memory(n_points);
  dz2 = get_memory(n_points);
  dx3 = get_memory(n_points);
  dy3 = get_memory(n_points);
  dz3 = get_memory(n_points);
  dx4 = get_memory(n_points);
  dy4 = get_memory(n_points);
  dz4 = get_memory(n_points);
  dv1_x = get_memory(n_points);
  dv1_y = get_memory(n_points);
  dv1_z = get_memory(n_points);
  dv2_x = get_memory(n_points);
  dv2_y = get_memory(n_points);
  dv2_z = get_memory(n_points);
  dv3_x = get_memory(n_points);
  dv3_y = get_memory(n_points);
  dv3_z = get_memory(n_points);
  dv4_x = get_memory(n_points);
  dv4_y = get_memory(n_points);
  dv4_z = get_memory(n_points);
  x2 = get_memory(n_points);
  y2 = get_memory(n_points);
  z2 = get_memory(n_points);
  x3 = get_memory(n_points);
  y3 = get_memory(n_points);
  z3 = get_memory(n_points);
  x4 = get_memory(n_points);
  y4 = get_memory(n_points);
  z4 = get_memory(n_points);

    
  initialize_pos(x,y,z, n_points, radius);
  initialize_vel(v_x,v_y,v_z, n_points, vel_initial, radius);
  initialize_mass(mass, n_points, unit_mass);

  /*Aca se hace una implementacion del metodo de Runge Kutta de cuarto orden*/
    
    FILE *in;
    in = fopen("3cuerpos_posiciones_velocidades.dat","w");
    
    for(i=0;i<n_steps;i++){
      
    get_acceleration(a_x, a_y, a_z, x, y, z, mass, n_points);
    for(j=0;j<n_points;j++){
       /*primer paso del Runge Kutta*/
        dx1[j] = h*v_x[j];
        dy1[j] = h*v_y[j];
        dz1[j] = h*v_z[j];
        
        dv1_x[j] = h*a_x[j];
        dv1_y[j] = h*a_y[j];
        dv1_z[j] = h*a_z[j];
        
        x2[j] = x[j]+(dx1[j]/2.0);
        y2[j] = y[j]+(dy1[j]/2.0);
        z2[j] = z[j]+(dz1[j]/2.0);
        
    }
      
      get_acceleration(a2_x,a2_y,a2_z,x2,y2,z2, mass, n_points);
      
      for(j=0;j<n_points;j++){
      /*segundo paso del runge Kutta*/  
        dx2[j] = h*(v_x[j]+(dv1_x[j]/2.0));
        dy2[j] = h*(v_y[j]+(dv1_y[j]/2.0));
        dz2[j] = h*(v_z[j]+(dv1_z[j]/2.0));
        
        dv2_x[j] = h*a2_x[j];
        dv2_y[j] = h*a2_y[j];
        dv2_z[j] = h*a2_z[j];
        
        x3[j] = x[j]+(dx2[j]/2.0);
        y3[j] = y[j]+(dy2[j]/2.0);
        z3[j] = z[j]+(dz2[j]/2.0);
          
        }
      
      get_acceleration(a3_x,a3_y,a3_z,x3,y3,z3, mass, n_points);
      
      for(j=0;j<n_points;j++){
      /*tercer paso del Runge Kutta*/    
        dx3[j] = h*(v_x[j]+(dv2_x[j]/2.0));
        dy3[j] = h*(v_y[j]+(dv2_y[j]/2.0));
        dz3[j] = h*(v_z[j]+(dv2_z[j]/2.0));
        
        dv3_x[j] = h*a3_x[j];
        dv3_y[j] = h*a3_y[j];
        dv3_z[j] = h*a3_z[j];
        
        x4[j] = x[j]+dx3[j];
        y4[j] = y[j]+dy3[j];
        z4[j] = z[j]+dz3[j];
          
      }
      
      get_acceleration(a4_x,a4_y,a4_z,x4,y4,z4, mass, n_points);
      
      for(j=0;j<n_points;j++){
      /*cuarto paso del Runge Kutta*/    
        dx4[j] = h*(v_x[j]+dv3_x[j]);
        dy4[j] = h*(v_y[j]+dv3_y[j]);
        dz4[j] = h*(v_z[j]+dv3_z[j]);
        
        dv4_x[j] = h*a4_x[j];
        dv4_y[j] = h*a4_y[j];
        dv4_z[j] = h*a4_z[j];
        
        
      }
      
      for(j=0;j<n_points;j++){
 /*Donde las aproximaciones para las velocidades y posiciones en cada instante de tiempo vendran dadas por: */  
          
	x[j] = x[j] + (1.0/6.0)*(dx1[j]+(2*dx2[j])+(2*dx3[j])+dx4[j]);
        y[j] = y[j] + (1.0/6.0)*(dy1[j]+(2*dy2[j])+(2*dy3[j])+dy4[j]);
        z[j] = z[j] + (1.0/6.0)*(dz1[j]+(2*dz2[j])+(2*dz3[j])+dz4[j]);
        
        v_x[j] = v_x[j] + (1.0/6.0)*(dv1_x[j]+(2*dv2_x[j])+(2*dv3_x[j])+dv4_x[j]);
        v_y[j] = v_y[j] + (1.0/6.0)*(dv1_y[j]+(2*dv2_y[j])+(2*dv3_y[j])+dv4_y[j]);
        v_z[j] = v_z[j] + (1.0/6.0)*(dv1_z[j]+(2*dv2_z[j])+(2*dv3_z[j])+dv4_z[j]);
	fprintf(in," %f %f %f %f %f %f\n ", x[j], y[j], z[j], v_x[j], v_y[j], v_z[j]);
      }
      
      
  }
    fclose(in);
    
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
          r_ij = (pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0));
          r_ij = sqrt(r_ij);
          ax[i] += -G_GRAV *mass[j]/ pow(r_ij,1.5) * (x[i] - x[j]);
          ay[i] += -G_GRAV *mass[j]/ pow(r_ij,1.5) * (y[i] - y[j]);
          az[i] += -G_GRAV *mass[j] / pow(r_ij,1.5) * (z[i] - z[j]);
          
      }
    }    
  }  
}
   

void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points, FLOAT radius){
  int i; 
  FLOAT delta_theta;
  delta_theta = 2.0*PI/n_points;
  
  for(i=0;i<n_points;i++){
    x[i] = cos(delta_theta * i) * radius;
    y[i] = sin(delta_theta * i) * radius;
    z[i] = 0.0;
  }
}
/*Las velocidades iniciales son perpendiculares a el plano que definen los tres cuerpos */
void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT vel, FLOAT radius){
  int i; 
  FLOAT delta_theta;
  delta_theta = 2.0*PI/n_points;
  
  for(i=0;i<3;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = vel;
  }  

    
}

void initialize_mass(FLOAT *mass, int n_points, FLOAT unit_mass){
  int i;
  for (i=0;i<n_points;i++){
    mass[i] = (i+1) * unit_mass;
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
