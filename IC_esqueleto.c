//
//  IC_esqueleto.c
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
float const G= 1.3267297*pow(10,11); //Constante gravitacional en Km³/(mSolares * s²)


//Defino la funcion que me va a devolver la magnitud de la velocidades
float magn_v(float m,float g, float r);

//-------------Main---------------
int main(int argc, char **argv){
    
    //Datos de entrada iniciales en kpc y km/s
    float xo=atof(argv[1]);
    float yo=atof(argv[2]);
    float zo=atof(argv[3]);
    float V_xo=atof(argv[4]);
    float V_yo=atof(argv[5]);
    float V_zo=atof(argv[6]);
    
    //Datos de entrada de la masa, radio y numero de particulas
    
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
    
    /*---------------------------------------------------------
     se les da un valor espacial a las particulas en cada orbita
     ej:
     primera: 12
     segunda: 18
     tercera: 24
     cuarta: 30
     quinta: 36
     -----------------------------------------------------------*/
    
    //se hace todo lo siguiente de acuerdo al radio externo y las orbitas creadas equiespaciadas con dist_orbi
    
    int (par_1)=12;
    
    //inicio las orbitas dividiendo 2PI entre el numero de particulas en cada orbita, teniendo el angulo, se le da una ubicacion en coordenadas cartesianas considerando el radio y se les suma la posicion de la masa central
    
    float firstorb=2*PI/(par_1);
    
    int alpha=0;
    int i=0;
    int k=0;
    
    //radio para cada orbita
    float radio1=10.0*3.0*pow(10.0,16.0);//en km
    
    for (i=0;i<par_1;i++){
        x_axis[i]=radio1*(cosf(firstorb*i))+xo;
        y_axis[i]=radio1*(sinf(firstorb*i))+yo;
    }
    k+=par_1;
    
    
    /*---------------------------------------------------------
     luego se calcula la magnitud de la velocidad de cada particula y se separa la componente vertical de la horizontal y se guarda en los array de velocidades, luego se suma la velocidad de la particula central. Lo mismo para cada orbita
     -----------------------------------------------------------*/
    float V_orb_1;
    V_orb_1= magn_v(M, G, radio1);
    
    /*---------------------------------------------------------
     se guarda la primera fila de datos introducidos en posicion -1 y se guarda tambien los arrays en orden de posiciones y velocidades de 0 hasta el numero de particulas
     -----------------------------------------------------------*/
    fprintf(output,"-1 %f %f %f %f %f %f\n",xo/(3.0*pow(10.0,16.0)),yo/(3.0*pow(10.0,16.0)),zo/(3.0*pow(10.0,16.0)),V_xo/(3.1536*pow(10,13)),V_yo/(3.1536*pow(10,13)),V_zo/(3.1536*pow(10,13)));
    
    float data1,data2,data3,data4,data5,data6;
    for(i=0;i<N;i++){
        
        data1=x_axis[i]/(3.0*pow(10.0,16.0));
        data2=y_axis[i]/(3.0*pow(10.0,16.0));
        data3=z_axis[i]/(3.0*pow(10.0,16.0));
        data4=Vx_axis[i]/(3.1536*pow(10,13));
        data5=Vy_axis[i]/(3.1536*pow(10,13));
        data6=Vz_axis[i]/(3.1536*pow(10,13));
        
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
