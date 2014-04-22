//
//  volterra_lotka.c
//
//
//  Created by Juliana Ayala and David Aleman
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*separo espacio en el computador para todo */

float * get_memory(int n_steps);

/*Remplazando los valores dados para A=20, B=1, C=30 y D=1 encontramos los valores de x(0)=30 y y(0)=20 para que exista el equilibrio*/

/*----MAIN-----*/
int main(){
    /* Variables a usar*/
    
  float *x;
  float *y;
  float *xprima;
  float *yprima;
  float delta_t=0.001;
  int n_steps= (int)(1.0/delta_t);
  int i,j;
 
  x = get_memory(n_steps);
  y = get_memory(n_steps);
  xprima = get_memory(n_steps);
  yprima = get_memory(n_steps);

    /* creo el archivo donde voy a escribir los datos*/
    FILE *in;
    in=fopen("datosxy.dat","w");

    
  /*Implementacion del metodo de Euler para resolver el problema*/
    for(i=30; i>0; i--){
    x[0]=i;
    y[0]=20.0;
    fprintf(in,"%f %f\n", x[0],y[0]);
    for(j=1; j<n_steps; j++){
      x[j]= x[j-1]+delta_t*(20*x[j-1]-x[j-1]*y[j-1]);
      y[j]= y[j-1]+delta_t*(-30*y[j-1]+x[j-1]*y[j-1]);
      fprintf(in,"%f %f\n",x[j],y[j]); 
  }
  
}
    return 0;
  
}

float * get_memory(int n_steps){
  float * x;
  if(!(x=malloc(sizeof(float)*n_steps))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;

}
