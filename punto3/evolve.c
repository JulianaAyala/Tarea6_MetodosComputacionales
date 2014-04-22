//
//  evolve.c
//
//
//  Created by Juliana Ayala and David Aleman
//
//
//

/*---------------------------------------------------------
 codigo que evoluciona un archivo con condiciones iniciales
 -----------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define FLOAT float

//Definimos las constantes iguales a las de las condiciones iniciales
const float G=4.296E-6; //Constante gravitacional en Kpc3/(mSolares *kyear2)

void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT M, int N);

/*----------Main-----------*/

int main(int argc, char **argv){
    
    //defino valores de entrada como masa, tiempo, momentos y el numero de particulas antes usado
    
    if(argc!=6){
        printf("Se requieren 5 datos a introducir para inicializar; \n-Nombre del archivo:\n-Masa:\n-Tiempo:\n-Momentos:\n-Numero de datos(igual al usado antiguamente):\n");
        exit(1);
    }
    
    
    float M=atof(argv[2]);
    float T=atof(argv[3]);
    float m=atof(argv[4]);
    int N=atoi(argv[5]);
    int num_lineas = N+1;
    
    //llamo el archivo con las condiciones iniciales para evolucionar
    FILE *input;
    input=fopen(argv[1],"r");

    /*---------------------------------------------------------
     llama el archivo producido, convierte y coloca en distintos arrays cada columna
     -----------------------------------------------------------*/
    float *x;
    float *y;
    float *z;
    float *v_x;
    float *v_y;
    float *v_z;
    int *ID;
    
    ID=malloc(num_lineas*sizeof(int));
    x= malloc(num_lineas*sizeof(float));
    y= malloc(num_lineas*sizeof(float));
    z= malloc(num_lineas*sizeof(float));
    v_x= malloc(num_lineas*sizeof(float));
    v_y= malloc(num_lineas*sizeof(float));
    v_z= malloc(num_lineas*sizeof(float));

    
    /*---------------------------------------------------------
     Registro de Datos en los arrays
     -----------------------------------------------------------*/
    rewind(input);
    int j;
    for(j=0;j<num_lineas;j++){
        fscanf(input, "%d %f %f %f %f %f %f\n",&ID[j], &x[j],&y[j],&z[j],&v_x[j],&v_y[j], &v_z[j]);//en kpc, km/s
        
    /*-----------------------------------------------------------------------------------------
        Crea grupos de archivos de salida
    --------------------------------------------------------------------------------------*/
        
    int maximo=100000;
    int p;
    FILE *output1, *out;
        
    for (p=1;p<10;p++){
            
        char filename1[20];
        sprintf(filename1,"Datos%d.txt",p);
        
        out=fopen(filename1,"r");
        if (out){
            fclose(out);
        }
        else{
            output1=fopen(filename1,"a");
            break;
            
        }
    }


        
    /*---------------------------------------------------------
        Usando RUNGEKUTTA:
        se aplica a cada particula alrededor de la particula central el cambio de posicion y
        velocidades dado por la interaccion con cada particula central presente por cada cambio de tiempo
    -----------------------------------------------------------*/
    int part,centros;
    
    //Conteo de todas las particulas centrales
    centros=0;
    for (j=0;j<num_lineas;j++){
        if( ID[j]==-1){
            centros++;
        }
    }
        
        /*-----------------------------------------------------------------------------------------
         Crea una lista de la ubicacion y velocidades de los nuevos centros de masa
         --------------------------------------------------------------------------------------*/
    int *CID;
    float *Cx;
    float *Cy;
    float *Cz;
    float *CVx;
    float *CVy;
    float *CVz;
    int i=0;
        
    CID=malloc(num_lineas*sizeof(int));
    Cx=malloc(num_lineas*sizeof(float));
    Cy=malloc(num_lineas*sizeof(float));
    Cz=malloc(num_lineas*sizeof(float));
    CVx=malloc(num_lineas*sizeof(float));
    CVy=malloc(num_lineas*sizeof(float));
    CVz=malloc(num_lineas*sizeof(float));
        
    for (j=-1;j<num_lineas;j++){
            CID[i]=j;
            Cx[i]=x[j];
            Cy[i]=y[j];
            Cz[i]=z[j];
            CVx[i]=v_x[j];
            CVy[i]=v_y[j];
            CVz[i]=v_z[j];
            i++;
        fprintf(output1,"%f %f %f %f %f %f\n ", Cx[j], Cy[j], Cz[j], CVx[j], CVy[j], CVz[j]);
    }
        
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
    
    /*timestep variables*/
    FLOAT h= 0.001;
    int n_steps = (int) T/(100*(m-1));
    
    /*memory allocation*/
    a2_x = malloc(num_lineas*sizeof(float));
    a2_y = malloc(num_lineas*sizeof(float));
    a2_z = malloc(num_lineas*sizeof(float));
    a3_x = malloc(num_lineas*sizeof(float));
    a3_y = malloc(num_lineas*sizeof(float));
    a3_z = malloc(num_lineas*sizeof(float));
    a4_x = malloc(num_lineas*sizeof(float));
    a4_y = malloc(num_lineas*sizeof(float));
    a4_z = malloc(num_lineas*sizeof(float));
    dx1 = malloc(num_lineas*sizeof(float));
    dy1 = malloc(num_lineas*sizeof(float));
    dz1 = malloc(num_lineas*sizeof(float));
    dx2 = malloc(num_lineas*sizeof(float));
    dy2 = malloc(num_lineas*sizeof(float));
    dz2 = malloc(num_lineas*sizeof(float));
    dx3 = malloc(num_lineas*sizeof(float));
    dy3 = malloc(num_lineas*sizeof(float));
    dz3 = malloc(num_lineas*sizeof(float));
    dx4 = malloc(num_lineas*sizeof(float));
    dy4 = malloc(num_lineas*sizeof(float));
    dz4 = malloc(num_lineas*sizeof(float));
    dv1_x = malloc(num_lineas*sizeof(float));
    dv1_y = malloc(num_lineas*sizeof(float));
    dv1_z = malloc(num_lineas*sizeof(float));
    dv2_x = malloc(num_lineas*sizeof(float));
    dv2_y = malloc(num_lineas*sizeof(float));
    dv2_z = malloc(num_lineas*sizeof(float));
    dv3_x = malloc(num_lineas*sizeof(float));
    dv3_y = malloc(num_lineas*sizeof(float));
    dv3_z = malloc(num_lineas*sizeof(float));
    dv4_x = malloc(num_lineas*sizeof(float));
    dv4_y = malloc(num_lineas*sizeof(float));
    dv4_z = malloc(num_lineas*sizeof(float));
    x2 = malloc(num_lineas*sizeof(float));
    y2 = malloc(num_lineas*sizeof(float));
    z2 = malloc(num_lineas*sizeof(float));
    
    /*Aca se hace una implementacion del metodo de Runge Kutta de cuarto orden*/
    
    for(i=0;i<n_steps;i++){
        
        get_acceleration(CVx, CVy, CVz, Cx, Cy, Cz, M, N);
        
        for(j=0;j<N;j++){
            /*primer paso del Runge Kutta*/
            dx1[j] = h*v_x[j];
            dy1[j] = h*v_y[j];
            dz1[j] = h*v_z[j];
            
            dv1_x[j] = h*CVx[j];
            dv1_y[j] = h*CVy[j];
            dv1_z[j] = h*CVz[j];
            
            x2[j] = x[j]+(dx1[j]/2.0);
            y2[j] = y[j]+(dy1[j]/2.0);
            z2[j] = z[j]+(dz1[j]/2.0);
            
        }
        
        get_acceleration(a2_x,a2_y,a2_z,x2,y2,z2, M, N);
        
        for(j=0;j<N;j++){
            /*segundo paso del runge Kutta*/
            dx2[j] = h*(v_x[j]+(dv1_x[j]/2.0));
            dy2[j] = h*(v_y[j]+(dv1_y[j]/2.0));
            dz2[j] = h*(v_z[j]+(dv1_z[j]/2.0));
            
            dv2_x[j] = h*a2_x[j];
            dv2_y[j] = h*a2_y[j];
            dv2_z[j] = h*a2_z[j];
            
            x2[j] = x[j]+(dx2[j]/2.0);
            y2[j] = y[j]+(dy2[j]/2.0);
            z2[j] = z[j]+(dz2[j]/2.0);
            
        }
        
        get_acceleration(a3_x,a3_y,a3_z,x2,y2,z2, M, N);
        
        for(j=0;j<N;j++){
            /*tercer paso del Runge Kutta*/
            dx3[j] = h*(v_x[j]+(dv2_x[j]/2.0));
            dy3[j] = h*(v_y[j]+(dv2_y[j]/2.0));
            dz3[j] = h*(v_z[j]+(dv2_z[j]/2.0));
            
            dv3_x[j] = h*a3_x[j];
            dv3_y[j] = h*a3_y[j];
            dv3_z[j] = h*a3_z[j];
            
            x2[j] = x[j]+dx3[j];
            y2[j] = y[j]+dy3[j];
            z2[j] = z[j]+dz3[j];
            
        }
        
        get_acceleration(a4_x,a4_y,a4_z,x2,y2,z2, M, N);
        
        for(j=0;j<N;j++){
            /*cuarto paso del Runge Kutta*/
            dx4[j] = h*(v_x[j]+dv3_x[j]);
            dy4[j] = h*(v_y[j]+dv3_y[j]);
            dz4[j] = h*(v_z[j]+dv3_z[j]);
            
            dv4_x[j] = h*a4_x[j];
            dv4_y[j] = h*a4_y[j];
            dv4_z[j] = h*a4_z[j];
            
            
        }
        
        for(j=0;j<N;j++){
            /*Donde las aproximaciones para las velocidades y posiciones en cada instante de tiempo vendran dadas por: */
            
            x[j] = x[j] + (1.0/6.0)*(dx1[j]+(2*dx2[j])+(2*dx3[j])+dx4[j]);
            y[j] = y[j] + (1.0/6.0)*(dy1[j]+(2*dy2[j])+(2*dy3[j])+dy4[j]);
            z[j] = z[j] + (1.0/6.0)*(dz1[j]+(2*dz2[j])+(2*dz3[j])+dz4[j]);
            
            v_x[j] = v_x[j] + (1.0/6.0)*(dv1_x[j]+(2*dv2_x[j])+(2*dv3_x[j])+dv4_x[j]);
            v_y[j] = v_y[j] + (1.0/6.0)*(dv1_y[j]+(2*dv2_y[j])+(2*dv3_y[j])+dv4_y[j]);
            v_z[j] = v_z[j] + (1.0/6.0)*(dv1_z[j]+(2*dv2_z[j])+(2*dv3_z[j])+dv4_z[j]);
            fprintf(output1," %f %f %f %f %f %f\n ", x[j], y[j], z[j], v_x[j], v_y[j], v_z[j]);
        }
        
        
    }
    fclose(output1);

}
}

void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT M, int N){
    int i,j;
    FLOAT r_ij;
    for(i=0;i<N;i++){
        ax[i]=0.0;
        ay[i]=0.0;
        az[i]=0.0;
        
        for(j=0;j<N;j++){
            if(j!=i){
                r_ij = (pow((x[i] - x[j]),2.0) + pow((y[i] - y[j]),2.0) + pow((z[i] - z[j]),2.0));
                r_ij = sqrt(r_ij);
                ax[i] += -G *M/ pow(r_ij,1.5) * (x[i] - x[j]);
                ay[i] += -G *M/ pow(r_ij,1.5) * (y[i] - y[j]);
                az[i] += -G *M/ pow(r_ij,1.5) * (z[i] - z[j]);
                
            }
        }    
    }  
}

    