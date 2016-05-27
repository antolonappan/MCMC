#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <thread>
#include <unistd.h> 
# define intr 100

using namespace std;
unsigned int microseconds=1000000;
float sum[2][3]={};
float intial_value[2][3]={},df=1.03;
float chi_sq[2][1000000] ={},A[2][1000000] ={},B[2][1000000] ={},C[2][1000000] ={},x[10]={},y[10]={},dy[10]={},sa,sb,sc;
int ter = 0,st=1,en=intr,sto[2],ki[2],h=0;
float N=intr,M=2;


int chain(int itr)
	{
        int n,t,k=0;
        float a,b,m1,p,m2,yf,mean,mean1,mean2,sd,sd1,sd2,c,c1,c2,l,x2,f,avg,avg1;
        float s,SUMa=0,SUMb=0,SUMc=0,ye[10]={},q[1000000]={};
        float pbO,pbN,pr,pb;
        
        const char* xdata = "x.txt";
        const char* ydata = "y.txt";
        const char* dydata = "dy.txt";
  	ifstream inX(xdata);
        ifstream inY(ydata);
        ifstream inDY(dydata);
        int i=0;
     	while(!inX.eof() && i<10 )
  	                           { 
  	  	                   inX >>setprecision(10)>> x[i];
                                   i++;
                                   } 
        i=0;
        while(!inY.eof() && i<10 )
  	                           { 
  	  	                   inY >>setprecision(10)>> y[i];
                                   i++;
                                   } 
         i=0;
         while(!inDY.eof() && i<10 )
  	                           { 
  	  	                   inDY >>setprecision(10)>> dy[i];
                                   i++;
                                   } 
     

                             
         unsigned int total = n; // no of iteration 
         n = 1000000;
        // usleep(microseconds);
         default_random_engine generator; // normal distribution
         
         random_device rd;//uniform random generator
         mt19937 re(rd());//uniform random twister
         
         uniform_real_distribution<double> ud(0, 1); 
         uniform_real_distribution<double> od(-100, 100);    
         
         mean =od(re);         intial_value[itr][0]=mean;
         mean1 =od(re);        intial_value[itr][1]=mean1;
         mean2 =od(re);        intial_value[itr][2]=mean2; 

         
         sd=.01;
         
         q[0]=0;
         for (i=0;i<10;i++) // setting intial value of chi.sq with entered parameter
                        {
                       ye[i]=(mean * pow(x[i],2)) + (mean1 * x[i]) + mean2;
                       //   ye[i]=(mean*x[i])+mean2;
                        f = (ye[i]-y[i])/dy[i];
                        x2 = pow(f,2);
                        q[0]=q[0]+x2;
                        }


         for (unsigned int j(0); j < n;++j) 
                                             {  
      					     if(ter ==1)
 							{
							terminate();
							n=sto[itr];
							}
                                             normal_distribution<double> da(mean ,.01); a=da(generator);
                                             normal_distribution<double> db(mean1,.02); b=db(generator);
                                             normal_distribution<double> dc(mean2,.03); c=dc(generator);
                                             
                                             s=0;//intializing the sum of chi.sq
               
                                             for (i=0;i<10;i++)
                                                            { 
                                                            ye[i]=(a * pow(x[i],2)) + (b * x[i]) + c;
                                                            //ye[i]=(a*x[i])+b;
                                                            f = (ye[i]-y[i])/dy[i];
                                                            x2 = pow(f,2);
                                                            s=s+x2;//summing chi.sq
                                                            }   
                               		     
 					                    q[k+1]=s;
    					
                                        	            if (q[k+1]<q[k])
                                                       			    { 
                                                                       
                                                             		    SUMa = SUMa + a;
                                                            		    SUMb = SUMb + b;
                                                             		    SUMc = SUMc + c;
                                                            
							   		    sum[itr][0]=SUMa;
							    		    sum[itr][1]=SUMb;
 							    		    sum[itr][2]=SUMc;
       		 
                		                                            chi_sq[itr][k+1] = q[k+1];
						      			    A[itr][k+1] = a;
							   		    B[itr][k+1] = b;
							   		    C[itr][k+1] = c;
                                                              
                                                           		    mean = a ; // using MH algorithm
                                                           		    mean1 = b;
                                                         		    mean2 =c;  // "     "     "
                                                         		    ki[itr]=k+1;
							  		    k++;
                                                                            } 
                                             
					             	    else if (q[k+1]>q[k])
                                                                                {
                                                                		pr =ud(re);
                                                                		pbO = exp (-(pow(q[k],2 )/2));
                                                                		pbN = exp (-(pow(q[k+1],2 )/2));
                                                                		pb = pbN/pbO;
                                                                		if(pr<pb)
                                                                       			{
 									 
 									       	         SUMa = SUMa + a;
 									                 SUMb = SUMb + b;
  										         SUMc = SUMc + c;
											 sum[itr][0]=SUMa;
											 sum[itr][1]=SUMb;
											 sum[itr][2]=SUMc;
 
                                                                       			 chi_sq[itr][k+1] = q[k+1];
							              			 A[itr][k+1] = a;
							              			 B[itr][k+1] = b;
							               			 C[itr][k+1] = c;
                                                                      
								     			 mean = a ; // using MH algorithm
                                                                       			 mean1 = b;  // "     "     "
                                                                       		         mean2 =c;
                                                                   			 ki[itr]=k+1;
								       			 k++;

                                                                       			}
                                                              else 
                                                                   j--;
                               
                                                              }
                        
                   
                   
                                              }
        return 0;
        }

//***************************************************************************************************************************************
int R()
       {
	float theta_bar_1_a,theta_bar_1_b,theta_bar_1_c,theta_bar_2_a,theta_bar_2_b,theta_bar_2_c,theta_bar_a,theta_bar_b,theta_bar_c,B_a,B_b,B_c;
	float Sa,Sb,Sc,S_1_a=0,S_1_b=0,S_1_c=0,S_2_a=0,S_2_b=0,S_2_c=0,S_1_af,S_1_bf,S_1_cf,S_2_af,S_2_bf,S_2_cf,Wa,Wb,Wc,Nf,Mf,Ra,Rb,Rc;
ofstream printR;
printR.open("R.csv");
printR<<"#Ra,#Rb,#Rc\n";
	while(1) 
		{

    		if(A[0][en]==0 || A[1][en]==0){
					      }
		else
		   {
		   sa=0;sb=0;sc=0;
                   for(int z=(en/2);z<en;z++)
                                        {
                                         sa=sa+ A[0][z];
                                         }
                   for(int z=(en/2);z<en;z++)
                                        {
                                         sb=sb+ B[0][z];
                                         }
                   for(int z=(en/2);z<en;z++)
                                        {
                                         sc=sc+ C[0][z];
                                         }

	           theta_bar_1_a = sa/(N/2);
		   theta_bar_1_b = sb/(N/2);
		   theta_bar_1_c = sc/(N/2);

                   sa=0;sb=0;sc=0;
                   for(int z=(en/2);z<en;z++)
                                        {
                                         sa=sa+ A[1][z];
                                         }
                   for(int z=(en/2);z<en;z++)
                                        {
                                         sb=sb+ B[1][z];
                                         }
                   for(int z=(en/2);z<en;z++)
                                        {
                                         sc=sc+ C[1][z];
                                         }
	
	 	   theta_bar_2_a = sa/(N/2);
  	           theta_bar_2_b = sb/(N/2);
	           theta_bar_2_c = sc/(N/2);
	
	           theta_bar_a = (theta_bar_1_a + theta_bar_2_a) / M;	
	           theta_bar_b = (theta_bar_1_b + theta_bar_2_b) / M;
	           theta_bar_c = (theta_bar_1_c + theta_bar_2_c) / M;

	           B_a = pow((theta_bar_1_a - theta_bar_a),2) + pow((theta_bar_2_a - theta_bar_a),2);
	           B_b = pow((theta_bar_1_b - theta_bar_b),2) + pow((theta_bar_2_b - theta_bar_b),2);
	           B_c = pow((theta_bar_1_c - theta_bar_c),2) + pow((theta_bar_2_c - theta_bar_c),2);

	           for(int u = (en/2); u<=en; u++)
			{
			 Sa =pow( (A[0][u] - theta_bar_1_a),2);
			 S_1_a = S_1_a + Sa;
                
			 Sb =pow( (B[0][u] - theta_bar_1_b),2);
			 S_1_b = S_1_b + Sb;
          
                         Sc =pow( (C[0][u] - theta_bar_1_c),2);
		         S_1_c = S_1_c + Sc; 
                         }

			 S_1_af = S_1_a/((N/2)-1);
			 S_1_bf = S_1_b/((N/2)-1);
			 S_1_cf = S_1_c/((N/2)-1);


		   for(int l = (en/2); l<=en; l++)
			{
			 Sa =pow( (A[1][l] - theta_bar_2_a),2);
			 S_2_a = S_2_a + Sa;
                
			 Sb =pow( (B[1][l] - theta_bar_2_b),2);
			 S_2_b = S_2_b + Sb;
          
                         Sc =pow( (C[1][l] - theta_bar_2_c),2);
			 S_2_c = S_2_c + Sc; 
                         }

		  S_2_af = S_2_a/((N/2)-1);
	          S_2_bf = S_2_b/((N/2)-1);
	          S_2_cf = S_2_c/((N/2)-1);


  	   	  Wa = (S_1_af + S_2_af)/2;
	          Wb = (S_1_bf + S_2_bf)/2;
	          Wc = (S_1_cf + S_2_cf)/2;
	
	
	
	
	          Nf = ((N/2)-1)/(N/2);
	          Mf = 1+(1/M);
	

	          Ra= ( (Nf*Wa) + (B_a*Mf))/Wa;
	          Rb= ( (Nf*Wb) + (B_b*Mf))/Wb;
	          Rc= ( (Nf*Wc) + (B_c*Mf))/Wc;
         
    

        if(Ra<df && Rb<df && Rc<df)
               {
                cout<<"\n"<<Ra<<"\t\t"<<Rb<<"\t\t"<<Rc/*<<"\t"<<N/2<<"\t"<<en*/;
                printR<<Ra<<","<<Rb<<","<<Rc<<endl;

                 float var,var1,var2;
                 var = Nf*Wa + B_a/(N);
                 var1= Nf*Wb + B_b/(N);
                 var2= Nf*Wc + B_c/(N);
                

                  float v,v1,v2;
                  v  = sqrt(var/(N*N));
                  v1 = sqrt(var1/(N*N));
                  v2 = sqrt(var2/(N*N));


ofstream printv;
printv.open("variance.csv");

                  cout<<"***********mean variance of a,b,c**********"<<endl;
                  cout<<"a\t\tb\t\tc\n";                
                  cout<<v<<"\t\t"<<v1<<"\t"<<v2<<endl;
 printv<<"#a,#b,#c\n"; 
 printv<<v<<","<<v1<<","<<v2<<endl;
printv.close();
		ofstream print;
		
	
		print.open("chain1.csv");
		print<<"#chisq,#a,#b,#c\n";

	        for(int i=1; i<100000; i++)
			 {
			  print<<chi_sq[0][i]<<","<<A[0][i]<<","<<B[0][i]<<","<<C[0][i]<<"\n";
     			 if (i>3 &&  chi_sq[0][i]==0 && A[0][i]==0 /*&& B[0][i]==0 && C[0][i]==0*/) {break;}
			      			 
			 }
		print.close();

		print.open("chain2.csv");
		print<<"#chisq,#a,#b,#c\n";
 		for(int i=1;i<100000;i++)
    			{ 
			print<<chi_sq[1][i]<<","<<A[1][i]<<","<<B[1][i]<<","<<C[1][i]<<"\n";
     			if ( i>2 && chi_sq[1][i]==0 && A[1][i]==0 /*&& B[1][i]==0 && C[1][i]==0*/){break;}
                             			
			}

		print.close();
	        ter=1;
		break;
		}
    
		cout<<"\n"<<Ra<<"\t\t"<<Rb<<"\t\t"<<Rc/*<<"\t"<<N/2<<"\t"<<en*/;  
printR<<Ra<<","<<Rb<<","<<Rc<<endl;  
           en=en+intr;
		N=N+intr;                             
		 }


	}
}


    




//****************************************************************************************************************************************
          

  int main()
            {
            cout<<"\nChains started to execute concurrently\n";


	    thread t1(chain,0);
	    thread t2(chain,1);
	    thread t3(R);

	    t1.join();
	    t2.join();
	    t3.join();

	    cout<<"\nTerminated : See profile.txt\n\n";
	    return 0;

	    }
