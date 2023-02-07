#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 9  

double a(double x_j, double x_k, double y_j, double y_k, double m){
	return (m*(x_k-x_j))/(pow(pow(x_k-x_j,2)+pow(y_k-y_j,2),1.5));
}

double V(double m_j, double m_k, double x_j, double x_k, double y_j, double y_k){
	return -(m_j*m_k)/(pow(pow(x_k-x_j,2)+pow(y_k-y_j,2),0.5));
}

double Temp_1(double R_s, double T_s, double d_s){
	return pow((0.9*pow(R_s, 2)*pow(T_s, 4))/(4*d_s*d_s), 0.25);
}

double Temp_2(double R_s, double T_s, double x_t, double x_s, double y_t, double y_s){
	return pow((0.9*pow(R_s, 2)*pow(T_s, 4))/(4*(pow(x_t-x_s, 2)+pow(y_t-y_s, 2))), 0.25);
}


int main(){
	
int i, j, k, N, l;
double T_nep, dt, dr, a_x, a_y, v, K, v_i, K_i, s_E_tot, H_s, H_a, H_k, H_p;

dt=0.000016615384; //el pas de temps més petit que podem agafar sense que el metode divergeixi, equival a 1 dia terrestre
T_nep=1; // periode orbital de Neptú normalitzat
N = T_nep/dt;
dr=0.000153846; 
	
double k1_x[M], k2_x[M], k3_x[M], k4_x[M], l1_x[M], l2_x[M], l3_x[M], l4_x[M], k1_y[M], k2_y[M], k3_y[M], k4_y[M], l1_y[M], l2_y[M], l3_y[M], l4_y[M];
double x_1[M], y_1[M], v_x_1[M], v_y_1[M], x_2[M], y_2[M], v_x_2[M], v_y_2[M],x_3[M], y_3[M], v_x_3[M], v_y_3[M]; 
double En[M], En_i[M], E_tot[M];
double T_est[5]={1, 2.1666667, 1.058333333, 0.595, 1.0768333};
double R_est[5]={0.000155556, 0.00000155556, 0.00032666667, 0.00004557777, 0.0000124444};
double m[M]={39.21390206, 0.00000649604, 0.00009565815, 0.0001184339, 0.0000126752, 85*0.03762950198, 0.01184339062, 0.00171709359, 0.0020795251}; //jupiter 2 = 0.081 Ms =85*M_jup
double T[M]={0.0,0.001461538, 0.0037307692, 0.006076923, 0.0114230769, 0.071923076, 0.178846153, 0.51153846, 1.0};
double T_2[M], T_3[M], T_4[M], T_s[N], T_e[N], T_tot[N];
double x[M]={0.0, 0.012867, 0.024, 0.03311, 0.05067, 0.17289, 0.3178, 0.6378, 1};
double y[M]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double v_x[M]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double v_y[M]={0.0,2*M_PI*x[1]/T[1],2*M_PI*x[2]/T[2],2*M_PI*x[3]/T[3],2*M_PI*x[4]/T[4],2*M_PI*x[5]/T[5],2*M_PI*x[6]/T[6],2*M_PI*x[7]/T[7],2*M_PI*x[8]};

for (j=0; j<M; j++){
	T_2[j]=T[j]*pow(1/0.63, 0.5); //AE Aquarii A
	T_3[j]=T[j]*pow(1/1.42, 0.5); //Procyon A
	T_4[j]=T[j]*pow(1/0.274, 0.5); //Estrella de Kapteyn 
	
}

double v_y2[M]={0.0,2*M_PI*x[1]/T_2[1],2*M_PI*x[2]/T_2[2],2*M_PI*x[3]/T_2[3],2*M_PI*x[4]/T_2[4],2*M_PI*x[5]/T_2[5],2*M_PI*x[6]/T_2[6],2*M_PI*x[7]/T_2[7],2*M_PI*x[8]/T_2[8]};
double v_y3[M]={0.0,2*M_PI*x[1]/T_3[1],2*M_PI*x[2]/T_3[2],2*M_PI*x[3]/T_3[3],2*M_PI*x[4]/T_3[4],2*M_PI*x[5]/T_3[5],2*M_PI*x[6]/T_3[6],2*M_PI*x[7]/T_3[7],2*M_PI*x[8]/T_3[8]};
double v_y4[M]={0.0,2*M_PI*x[1]/T_4[1],2*M_PI*x[2]/T_4[2],2*M_PI*x[3]/T_4[3],2*M_PI*x[4]/T_4[4],2*M_PI*x[5]/T_4[5],2*M_PI*x[6]/T_4[6],2*M_PI*x[7]/T_4[7],2*M_PI*x[8]/T_4[8]};

	//FILE* fichero;
	//fichero= fopen("Hsol", "w");
for(l=1; l<N; l++){
	H_s=Temp_1(R_est[0], T_est[0], l*dr);
	H_a=Temp_1(R_est[1], T_est[1], l*dr);
	H_p=Temp_1(R_est[2], T_est[2], l*dr);
	H_k=Temp_1(R_est[3], T_est[3], l*dr);
	//printf("%lf \t %lf \t %lf \t %lf\n", H_s*6000, H_a*6000, H_p*6000, H_k*6000);
	//fprintf(fichero,"%lf \t %lf \t %lf \t %lf \t %lf \n", H_s*6000, H_a*6000, H_p*6000, H_k*6000, l*dr);
}
    //fclose(fichero);
for(j=0; j<M; j++){
	
	K_i=0.5*m[j]*(pow(v_x[j], 2)+pow(v_y[j], 2));
	v_i=0;

	
	for(k=0; k<M; k++){
		if(k==j){
		}
		else{
			v_i+=V(m[j], m[k], x[j], x[k], y[j], y[k]);
		}
	}
	
	
	En_i[j]=K_i+v_i;
	//printf("%lf\n", En_i[j]);
}
	FILE* fichero;
	fichero= fopen("T real d", "w");
//Runge-Kutta 

for(i=1; i<N; i++){
    for(j=0; j<M; j++){
    	
		k1_x[j]=v_x[j];
        k1_y[j]=v_y[j];
		a_x=0;
        a_y=0;
		
		for(k=0; k<M; k++){
            if(k==j){
            }
            else{
                a_x+=a(x[j],x[k],y[j],y[k],m[k]);
				a_y+=a(y[j],y[k],x[j],x[k],m[k]);
            }
        }		

       	l1_x[j]=a_x;
       	l1_y[j]=a_y;
		x_1[j]=x[j]+0.5*dt*k1_x[j];
        y_1[j]=y[j]+0.5*dt*k1_y[j];
		v_x_1[j]=v_x[j]+0.5*dt*l1_x[j];
		v_y_1[j]=v_y[j]+0.5*dt*l1_y[j];
		
    }
    
    for(j=0; j<M; j++){

        k2_x[j]= v_x_1[j];
        k2_y[j] = v_y_1[j];
        a_x=0;
        a_y=0;

        for(k=0; k<M; k++){
            if(k==j){
            }
            else{
            	a_x+=a(x_1[j],x_1[k],y_1[j],y_1[k],m[k]);
				a_y+=a(y_1[j],y_1[k],x_1[j],x_1[k],m[k]);
            }
        }

        l2_x[j]=a_x;
        l2_y[j]=a_y;
		x_2[j]=x[j]+0.5*dt*k2_x[j];
        y_2[j]=y[j]+0.5*dt*k2_y[j];
        v_x_2[j]=v_x[j]+0.5*dt*l2_x[j];
       	v_y_2[j]=v_y[j]+0.5*dt*l2_y[j];
    }
    
	for(j=0; j<M; j++){

        k3_x[j]=v_x_2[j];
        k3_y[j]=v_y_2[j];
        a_x=0;
        a_y=0;

        for(k=0; k<M; k++){
            if(k==j){
        	}
            else{
                a_x+=a(x_2[j],x_2[k],y_2[j],y_2[k],m[k]);
				a_y+=a(y_2[j],y_2[k],x_2[j],x_2[k],m[k]);
        	}
        }

        l3_x[j]=a_x;
        l3_y[j]=a_y;
        x_3[j]=x[j]+dt*k3_x[j];
        y_3[j]=y[j]+dt*k3_y[j];
        v_x_3[j]=v_x[j]+dt*l3_x[j];
        v_y_3[j]=v_y[j]+dt*l3_y[j];
    }
    
    for(j=0; j<M; j++){

        k4_x[j]=v_x_3[j];
        k4_y[j]=v_y_3[j];
        a_x=0;
        a_y=0;

        for(k=0; k<M; k++){
            if(k==j){
            }
            else{
                a_x+=a(x_3[j],x_3[k],y_3[j],y_3[k],m[k]);
				a_y+=a(y_3[j],y_3[k],x_3[j],x_3[k],m[k]);
            }
    }

    l4_x[j]=a_x;
    l4_y[j]=a_y;
    x[j]+=(dt/6)*(k1_x[j]+2.0*k2_x[j]+2.0*k3_x[j]+k4_x[j]);
    y[j]+=(dt/6)*(k1_y[j]+2.0*k2_y[j]+2.0*k3_y[j]+k4_y[j]);
    v_x[j]+=(dt/6)*(l1_x[j]+2.0*l2_x[j]+2.0*l3_x[j]+l4_x[j]);
    v_y[j]+=(dt/6)*(l1_y[j]+2.0*l2_y[j]+2.0*l3_y[j]+l4_y[j]);
    		}
			//printf("%lf \t %lf \t %lf \t %lf \n", x[8], y[8], v_x[8], v_y[8]);


	for(j=0; j<M; j++){
		
		K=0.5*m[j]*(pow(v_x[j], 2)+pow(v_y[j], 2));
		v=0;
		
		for(k=0; k<M; k++){
			
			if(k==j){
			}
			else{
				v+=V(m[j], m[k], x[j], x[k], y[j], y[k]);
			}
			
		}
		
	En[j]=K+v;
	E_tot[j]=En[j]-En_i[j];
	s_E_tot+=E_tot[j];

	}
	

	//printf("%lf \t %lf\n", s_E_tot, v);
	s_E_tot=0;

	T_s[i]=Temp_2(R_est[0], T_est[0], x[3], x[0], y[3], y[0]);
	T_e[i]=Temp_2(R_est[4], T_est[4], x[3], x[5], y[3], y[5]);
	T_tot[i]=T_s[i]+T_e[i];
	
//	printf("%lf \t %lf \t %lf \t %lf \t %lf \n", x[5], y[5], v_x[5], v_y[5], i*dt);
//	fprintf(fichero,"%lf \t %lf \t %lf \t %lf \t %lf \n", x[5], y[5], v_x[5], v_y[5], i*dt);
printf("%lf \t %lf \n", ((T_tot[i]*6000)-273.15), i*dt);
fprintf(fichero,"%lf \t %lf \n", ((T_tot[i]*6000)-273.15), i*dt*60000);
	}
	
	fclose(fichero);
	
	return 0;
	
}
