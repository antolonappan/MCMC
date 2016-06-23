#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <thread>
#include <unistd.h> 
#include <stdio.h> 
# define intr 100

using namespace std;
float sum[2][3]={};
float intial_value[2][3]={},df=1;
float chi_sq[2][10000000] ={},A[2][50000000] ={},B[2][50000000] ={},C[2][50000000] ={},D[2][50000000] ={},E[2][50000000] ={},x[25]={},y[25]={};
float dy[25]={},sa,sb,sc,sd,se;
int ter = 0,st=1,en=intr,sto[2],ki[2],h=0;
float N=intr,M=2;


int chain(int itr)
	{
        int n,t,k=0;
        float a,b,d,e,m1,p,m2,yf,mean,mean1,mean2,mean3,mean4,sd,sd1,sd2,sd3,sd4,c,c1,c2,l,x2,f,avg,avg1;
        float s,SUMa=0,SUMb=0,SUMc=0,SUMd=0,SUMe=0,ye[10]={},q[1000000]={};
        float pbO,pbN,pr,pb;
        
        const char* xdata = "z.txt";
        const char* ydata = "H.txt";
        const char* dydata = "dH.txt";
  	ifstream inX(xdata);
        ifstream inY(ydata);
        ifstream inDY(dydata);
        int i=0;
     	while(!inX.eof() && i<25 )
  	                           { 
  	  	                   inX >>setprecision(10)>> x[i];
                                   i++;
                                   } 
        i=0;
        while(!inY.eof() && i<25 )
  	                           { 
  	  	                   inY >>setprecision(10)>> y[i];
                                   i++;
                                   } 
         i=0;
         while(!inDY.eof() && i<25 )
  	                           { 
  	  	                   inDY >>setprecision(10)>> dy[i];
                                   i++;
                                   } 
     
         
                             
       
         n = 10000000;
        
         default_random_engine generator;
         
         random_device rd;
         mt19937 re(rd());
         
         uniform_real_distribution<double> ud(0, 1); 
         uniform_real_distribution<double> od(0, 3);    
         
         mean =od(re);         intial_value[itr][0]=mean;
         mean1 =od(re);        intial_value[itr][1]=mean1;
         mean2 =od(re);        intial_value[itr][2]=mean2; 
         //mean3 =od(re);        intial_value[itr][3]=mean3;
         //mean4 =od(re);        intial_value[itr][4]=mean4; 

         //float z,H0,ohmM0,ohmL
         
         q[0]=0;
         for (i=0;i<25;i++) 
                        {\*****************************************lambda-CDM equation*************************************\
                            
                         ye[i] = sqrt(pow(mean,2)*( (mean1*pow((1+x[i]),3)) + mean2 + ((1-mean1-mean2)*pow((1+x[i]),2)) ) );
                                              
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
                                             normal_distribution<double> da(mean ,1); a=da(generator);
                                             normal_distribution<double> db(mean1,1); b=db(generator);
                                             normal_distribution<double> dc(mean2,1); c=dc(generator);
                                            // normal_distribution<double> dd(mean3,.01); d=dd(generator);
                                            // normal_distribution<double> de(mean4,.01); e=de(generator);
                                             
                                             s=0;
               
                                             for (i=0;i<25;i++)
                                                            { 
                         /*Polynomial to be fitted*/         ye[i] = sqrt(pow(a,2)*( (b*pow((1+x[i]),3)) + c + ((1-b-c)*pow((1+x[i]),2)) ) );
                                                         
                                                            f = pow((ye[i]-y[i]),2)/pow(dy[i],2);
                                                            x2 = f;
                                                            s=s+x2;
                                                            }   
                               		  
 					                    q[k+1]=s;
    					
                                        	            if (q[k+1]<q[k])
                                                       			    { 
                                                                       
                                                             		    SUMa = SUMa + a;
                                                            		    SUMb = SUMb + b;
                                                             		    SUMc = SUMc + c;
                                                                           // SUMd = SUMd + d;
                                                                           // SUMe = SUMe + e;

							   		    sum[itr][0]=SUMa;
							    		    sum[itr][1]=SUMb;
              						    		    sum[itr][2]=SUMc;
       		                                                           // sum[itr][3]=SUMd;
       		                                                           // sum[itr][4]=SUMe;
       		 
                		                                            chi_sq[itr][k+1] = q[k+1];
						      			    A[itr][k+1] = a;
							   		    B[itr][k+1] = b;
							   		    C[itr][k+1] = c;
                                                                           // D[itr][k+1] = d;
                                                           	           // E[itr][k+1] = e;
	                                                                    mean = a ; 
                                                           		    mean1 = b;
                                                         		    mean2 =c;
                                                                           // mean3 =d;
                                                                           // mean4 =e;
  
                                                         		    ki[itr]=k+1;
					   		  		    k++;// cout<<k<<endl;

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
                                                                                       //  SUMd = SUMd + d;
                                                                                       //  SUMe = SUMe + e;


											 sum[itr][0]=SUMa;
											 sum[itr][1]=SUMb;
											 sum[itr][2]=SUMc;
                                                                                        // sum[itr][3]=SUMd;
                                                                      		        // sum[itr][4]=SUMe;	 
                                                                                        
                                                                                         chi_sq[itr][k+1] = q[k+1];
							              			 A[itr][k+1] = a;
							              			 B[itr][k+1] = b;
							               			 C[itr][k+1] = c;
                                                                                        // D[itr][k+1] = d;

                                                                                       //  E[itr][k+1] = e;
								     			 mean = a ;
                                                                       			 mean1 = b;
                                                                       		         mean2 =c;
                                                                                         //   mean3 =d;
                                                                                           //  mean4 =e;


                                                                   			 ki[itr]=k+1;
								       		    	 k++;
                                                                                  //       cout<<k<<endl;

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
	float theta_bar_1_a,theta_bar_1_b,theta_bar_1_c,theta_bar_1_d,theta_bar_1_e,theta_bar_2_a,theta_bar_2_b,theta_bar_2_c,theta_bar_2_d;
        float theta_bar_2_e, theta_bar_a,theta_bar_b,theta_bar_c,theta_bar_d,theta_bar_e,B_a,B_b,B_c,B_d,B_e;
	float Sa,Sb,Sc,Sd,Se,S_1_a=0,S_1_b=0,S_1_c=0,S_1_d=0,S_1_e=0;
        float S_2_a=0,S_2_b=0,S_2_c=0,S_2_d=0,S_2_e=0,S_1_af,S_1_bf,S_1_cf,S_1_df,S_1_ef,S_2_af,S_2_bf,S_2_cf,S_2_df,S_2_ef,Wa,Wb,Wc,Wd,We,Nf,Mf;
        float Ra,Rb,Rc,Rd,Re;
ofstream printR;
printR.open("R.csv");
printR<<"#Ra,#Rb,#Rc\n";
	while(1) 
		{

    		if(A[0][en]==0 || A[1][en]==0)
                     {
					      }
		else
		   {
		   sa=0;sb=0;sc=0;sd=0;se=0;
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
                  /* for(int z=(en/2);z<en;z++)
                                        {
                                         sd=sd+ D[0][z];
                                         }
                   for(int z=(en/2);z<en;z++)
                                        {
                                         se=se+ E[0][z];
                                         }*/
	           theta_bar_1_a = sa/(N/2);
		   theta_bar_1_b = sb/(N/2);
		   theta_bar_1_c = sc/(N/2);
                  // theta_bar_1_d = sd/(N/2);
                  // theta_bar_1_e = se/(N/2);

                   sa=0;sb=0;sc=0;sd=0;se=0;
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
	          /* for(int z=(en/2);z<en;z++)
                                        {
                                         sd=sd+ D[1][z];
                                         }
                   for(int z=(en/2);z<en;z++)
                                        {
                                         se=se+ E[1][z];
                                         }*/

	 	   theta_bar_2_a = sa/(N/2);
  	           theta_bar_2_b = sb/(N/2);
	           theta_bar_2_c = sc/(N/2);
	          // theta_bar_2_d = sd/(N/2);
                  // theta_bar_2_e = se/(N/2);
	
	
	           theta_bar_a = (theta_bar_1_a + theta_bar_2_a) / M;	
	           theta_bar_b = (theta_bar_1_b + theta_bar_2_b) / M;
	           theta_bar_c = (theta_bar_1_c + theta_bar_2_c) / M;
                   //theta_bar_d = (theta_bar_1_d + theta_bar_2_d) / M;
                   //theta_bar_e = (theta_bar_1_e + theta_bar_2_e) / M;
	           
                   B_a = pow((theta_bar_1_a - theta_bar_a),2) + pow((theta_bar_2_a - theta_bar_a),2);
	           B_b = pow((theta_bar_1_b - theta_bar_b),2) + pow((theta_bar_2_b - theta_bar_b),2);
	           B_c = pow((theta_bar_1_c - theta_bar_c),2) + pow((theta_bar_2_c - theta_bar_c),2);
                   //B_d = pow((theta_bar_1_d - theta_bar_d),2) + pow((theta_bar_2_d - theta_bar_d),2);
                   //B_e = pow((theta_bar_1_e - theta_bar_e),2) + pow((theta_bar_2_e - theta_bar_e),2);

	           for(int u = (en/2); u<=en; u++)
			{
			 Sa =pow( (A[0][u] - theta_bar_1_a),2);
			 S_1_a = S_1_a + Sa;
                
			 Sb =pow( (B[0][u] - theta_bar_1_b),2);
			 S_1_b = S_1_b + Sb;
          
                         Sc =pow( (C[0][u] - theta_bar_1_c),2);
		         S_1_c = S_1_c + Sc; 
                         
                     /*    Sd =pow( (D[0][u] - theta_bar_1_d),2);
 		         S_1_d = S_1_d + Sd; 

                         Se =pow( (E[0][u] - theta_bar_1_e),2);
		         S_1_e = S_1_e + Se;*/ 
                         }

			 S_1_af = S_1_a/((N/2)-1);
			 S_1_bf = S_1_b/((N/2)-1);
			 S_1_cf = S_1_c/((N/2)-1);
                        // S_1_df = S_1_d/((N/2)-1);
                        // S_1_ef = S_1_e/((N/2)-1);


		   for(int l = (en/2); l<=en; l++)
			{
			 Sa =pow( (A[1][l] - theta_bar_2_a),2);
			 S_2_a = S_2_a + Sa;
                
			 Sb =pow( (B[1][l] - theta_bar_2_b),2);
			 S_2_b = S_2_b + Sb;
          
                         Sc =pow( (C[1][l] - theta_bar_2_c),2);
			 S_2_c = S_2_c + Sc; 
                        /* Sd =pow( (D[1][l] - theta_bar_2_d),2);
			 S_2_d = S_2_d + Sd; 

                         Se =pow( (E[1][l] - theta_bar_2_e),2);
			 S_2_e = S_2_e + Se;   */                       

                        }

		  S_2_af = S_2_a/((N/2)-1);
	          S_2_bf = S_2_b/((N/2)-1);
	          S_2_cf = S_2_c/((N/2)-1);
                 // S_2_df = S_2_d/((N/2)-1);
                 // S_2_ef = S_2_e/((N/2)-1);


  	   	  Wa = (S_1_af + S_2_af)/2;
	          Wb = (S_1_bf + S_2_bf)/2;
	          Wc = (S_1_cf + S_2_cf)/2;
                 // Wd = (S_1_df + S_2_df)/2;
	         // We = (S_1_ef + S_2_ef)/2;
	
	          Nf = ((N/2)-1)/(N/2);
	          Mf = 1+(1/M);
	

	          Ra= ( (Nf*Wa) + (B_a*Mf))/Wa;
	          Rb= ( (Nf*Wb) + (B_b*Mf))/Wb;
	          Rc= ( (Nf*Wc) + (B_c*Mf))/Wc;
                 // Rd= ( (Nf*Wd) + (B_d*Mf))/Wd;
                 // Re= ( (Nf*We) + (B_e*Mf))/We;
          
              

        if(Ra<df && Rb<df && Rc<df)
               {cout<<"chains terminated"<<endl;
                cout<<"\n"<<Ra<<"\t"<<Rb<<"\t"<<Rc<<endl;
                
                 printR<<Ra<<","<<Rb<<","<<Rc<<endl;

                 float var,var1,var2,var3,var4;
                 var = Nf*Wa + B_a/(N/2);
                 var1= Nf*Wb + B_b/(N/2);
                 var2= Nf*Wc + B_c/(N/2);
                 //var3= Nf*Wd + B_d/(N/2);
                 //var4= Nf*We + B_e/(N/2);

                  float v,v1,v2,v3,v4;
                  v  = sqrt(var/(N/2));
                  v1 = sqrt(var1/(N/2));
                  v2 = sqrt(var2/(N/2));
                 // v3 = sqrt(var3/(N/2));
                 // v4 = sqrt(var4/(N/2));




ofstream printv;
printv.open("variance.csv");


                  cout<<"***********mean variance of a,b,c **********"<<endl;
                 // cout<<"a\t\tb\t\tc\n";                
                  cout<<v<<"\t\t"<<v1<<"\t\t"<<v2<<"\t\t"<<v3<<"\t"<<v4<<endl;
 printv<<"#a,#b,#c,#d,#e\n"; 
 printv<<v<<","<<v1<<","<<v2<<","<<v3<<","<<v4<<endl;
printv.close();
		ofstream print;
		
	
		print.open("chain1.csv");
		print<<"#chisq,#a,#b,#c\n";

	        for(int i=1; i<100000; i++)
			 {
			  print<<chi_sq[0][i]<<","<<A[0][i]<<","<<B[0][i]<<","<<C[0][i]<<"\n";
     			 if (i>3 &&  chi_sq[0][i]==0 && A[0][i]==0 ) 
                         {
                        break;
                         }
			      			 
			 }
		print.close();

	 print.open("chain2.csv");

		print<<"#chisq,#a,#b,#c\n";
 		for(int i=1;i<100000;i++){
               print<<chi_sq[1][i]<<","<<A[1][i]<<","<<B[1][i]<<","<<C[1][i]<<"\n";



     			if ( i>2 && chi_sq[1][i]==0 && A[1][i]==0)
                      {
                     break;
                       }
                           
  			
			}

		print.close();

	        ter=1;
		break;
		}
    
		cout<<"\n"<<Ra<<"\t"<<Rb<<"\t"<<Rc<<endl; 
                printR<<Ra<<","<<Rb<<","<<Rc<<endl;   
           en=en+intr;
		N=N+intr;                             
		 }


	}
}


    




//****************************************************************************************************************************************
          

  int main()
            {
            


	    thread t1(chain,0);
	    thread t2(chain,1);
	    thread t3(R);

	    t1.join();
	    t2.join();
	    t3.join();

	   
	    return 0;

	    }
