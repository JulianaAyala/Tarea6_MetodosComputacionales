//
//  IC.c
//  
//
//  Created by Juliana Ayala and David Aleman
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//Defino las constantes que se van a usar

float const PI= 3.14159265;
float const G= 4.296E-6; //Constante gravitacional en Km2 * kpc/(mSolares * s2)


//Defino la funcion que me va a devolver la magnitud de la velocidades
float magn_v(float m,float g, float r);

//-------------Main---------------
int main(int argc, char **argv){
    
    if(argc!=10){
        printf("Se requieren 9 datos para introducir para el centro de masa; \n-posicion inicial en eje x:\n-posicion inicial en eje y:\n-posicion inicial en eje z:\n-velocidad inicial en el eje x:\n-velocidad inicial en el eje y:\n-velocidad inicial en el eje z:\n-masa: \n-radio de la orbita externa: \n-numero de particulas:\n");
        exit(1);
    }
    
    //Datos de entrada iniciales en kpc y km/s
    float xo=atof(argv[1]);
    float yo=atof(argv[2]);
    float zo=atof(argv[3]);
    float V_xo=atof(argv[4]);
    float V_yo=atof(argv[5]);
    float V_zo=atof(argv[6]);
    
    //Datos de entrada de la masa, radio externo y numero de particulas
    
    float M=atof(argv[7]);
    float R=atof(argv[8]);
    int N=atoi(argv[9]);
    
    /*------------------------------------------------------------------------
     Checkpoint; Los datos de posicion y velocidades fueron introducidos
     --------------------------------------------------------------------------*/
    
    /*------------------------------------------------------------------------
     Variables
     Output
     Distancia de las orbitas para encontrar las orbitas de acuerdo a R introducido
     --------------------------------------------------------------------------*/
    FILE *output;
    float dist_orbi=10.0*3.0*pow(10.0,16.0); //en km
    char filename[100]="IC_output.dat";
    output=fopen(filename,"a");
    /*---------------------------------------------------------
     llama el archivo producido, convierte y coloca en distintos arrays cada columna
     convierte a km
     -----------------------------------------------------------*/
    float *x_axis;
    float *y_axis;
    float *z_axis;
    float *Vx_axis;
    float *Vy_axis;
    float *Vz_axis;
    
  
    x_axis= malloc(N*sizeof(float));
    y_axis= malloc(N*sizeof(float));
    z_axis= malloc(N*sizeof(float));

    Vx_axis= malloc(N*sizeof(float));
    Vy_axis= malloc(N*sizeof(float));
    Vz_axis= malloc(N*sizeof(float));

    
    /*---------------------------------------------------------
     se tiene la magnitud de la velocidad dependiente de radio, masa y la constante G
     -----------------------------------------------------------*/

    int i=0;
    float *r=malloc(N*sizeof(float));
    float *v=malloc(N*sizeof(float));
    
    //inicio las orbitas multiplicando 2PI por numeros random dentro del radio externo que di, teniendo el angulo, se le da una ubicacion en coordenadas cartesianas considerando el radio y se les suma la posicion de la masa central
    
    for (i=0;i<N;i++){
        
        float orb=2*PI*drand48();
        
        x_axis[i]=cos(orb)+xo;
        y_axis[i]=sin(orb)+yo;
        z_axis[i]=zo;
        
        r[i]=sqrt(pow(x_axis[i]-xo,2)+pow(y_axis[i]-yo,2)+pow(z_axis[i]-zo,2));
    }

    /*---------------------------------------------------------
     luego se calcula la magnitud de la velocidad de cada particula y se separa la componente vertical de la horizontal y se guarda en los array de velocidades, luego se suma la velocidad de la particula central.
     -----------------------------------------------------------*/
    for(i=0;i<N;i++){
        
        v[i]= magn_v(M,G,r[i]);
        
        Vx_axis[i]=V_xo - v[i]*((y_axis[i]-yo)/r[i]);
        Vy_axis[i]=V_yo + v[i]*((x_axis[i]-xo)/r[i]);
        Vz_axis[i]=V_zo;
    }
    
    /*---------------------------------------------------------
     se guarda la primera fila de datos introducidos en posicion -1 y se guarda tambien los arrays en orden de posiciones y velocidades de 0 hasta el numero de particulas
     -----------------------------------------------------------*/
    fprintf(output,"-1 %f %f %f %f %f %f\n",xo,yo,zo,V_xo,V_yo,V_zo);
    
    float data1,data2,data3,data4,data5,data6;
    for(i=0;i<N;i++){
        
        data1=x_axis[i];
        data2=y_axis[i];
        data3=z_axis[i];
        data4=Vx_axis[i];
        data5=Vy_axis[i];
        data6=Vz_axis[i];
        
        fprintf(output,"%d %f %f %f %f %f %f\n",i,data1,data2,data3,data4,data5,data6);
    }
    
    /*---------------------------------------------------------
     FIN
     -----------------------------------------------------------*/
    fclose(output);
    return 0;
}



/*---------------------------------------------------------
 Funciones
 -----------------------------------------------------------*/
//entran valores de masa, constante gravitacional y radio, devuelve la magnitud de la velocidad.
//Se requiere salida en terminos de km/s
float magn_v(float m,float g, float r){
    float dato=m*g/r;
    float dato2=sqrt(dato);
    return dato2;
}
